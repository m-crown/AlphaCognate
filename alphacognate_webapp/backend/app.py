from fastapi import FastAPI, Depends, HTTPException, Query
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse, StreamingResponse, Response
from typing import Annotated, List, Optional
from sqlmodel import Session, select, func
import pandas as pd

from models import CognateLigand, CognateLigandResponse, Structure, StructuresListResponse, Transplant, TransplantFull, Ligand, CognateLigandMapping, TransplantResponse, TransplantsListResponse
from database import create_db_engine

from dotenv import load_dotenv
import os
from fastapi.middleware.cors import CORSMiddleware
from fastapi import Query

from abc import ABC, abstractmethod
from typing import BinaryIO
from pathlib import Path
import boto3
import gzip 
import io 
import gemmi

load_dotenv()  # Load environment variables from .env file
DB_PASSWORD = os.getenv("DB_PASSWORD")

if DB_PASSWORD is None:
    raise ValueError("DB_PASSWORD environment variable is not set!")

engine = create_db_engine(username = "admin", password = DB_PASSWORD, host = "postgresdb", database = "postgres")

def get_session():
    with Session(engine) as session:
        yield session


SessionDep = Annotated[Session, Depends(get_session)]

app = FastAPI()

origins = [
    "http://localhost",
    "http://localhost:8080",
    "http://localhost:5173",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
async def root():
    return {"message": "Hello World"}

@app.get("/structures/", response_model=StructuresListResponse)
def read_structures(
        session: Session = Depends(get_session),
        offset: int = 0,
        limit: Annotated[int, Query(le=100)] = 100,
) -> list[Structure]:
    total_count = session.exec(select(func.count()).select_from(Structure)).one()  # Get total row count
    structures = session.exec(select(Structure).offset(offset).limit(limit)).all()
    
    return StructuresListResponse(
        data=structures,
        meta={"totalRowCount": total_count}
    )

@app.get("/transplants/", response_model=TransplantsListResponse)
def read_transplants(
    session: Session = Depends(get_session),
    structure_name: str = Query(..., description="Filter transplants by structure name"),
    best: bool = Query(False, description="Return only the best matching cognate ligand per transplant")
) -> TransplantsListResponse:
    
    structure = session.exec(select(Structure).where(Structure.name == structure_name)).first()
    
    transplants = session.exec(
        select(Transplant).where(Transplant.structure_name == structure_name)
    ).all()

    transplant_full_list = []

    for transplant in transplants:
        ligand = session.get(Ligand, transplant.ligand_id)
        transplant_response = TransplantResponse(
            transplant_id=transplant.id,
            structure_name=transplant.structure_name,
            ligand_id=transplant.ligand_id,
            ligand_chain=transplant.ligand_chain,
            ligand_residues=transplant.ligand_residues,
            ec_list=transplant.ec_list,
            tcs=transplant.tcs,
            foldseek_rmsd=transplant.foldseek_rmsd,
            global_rmsd=transplant.global_rmsd,
            local_rmsd=transplant.local_rmsd,
            hetcode=ligand.het_code,
            name=ligand.name,
            smiles=ligand.smiles
        )

        # Get all mappings
        cognate_ligands = session.exec(
            select(CognateLigandMapping).where(
                CognateLigandMapping.ligand_instance_id == transplant.id
            )
        ).all()

        # If best=True, keep only the mapping with the highest similarity
        if best and cognate_ligands:
            max_sim = max(cl.similarity for cl in cognate_ligands)
            cognate_ligands = [cl for cl in cognate_ligands if cl.similarity == max_sim]

        cognate_ligand_response_list = []
        for cognate_ligand in cognate_ligands:
            cognate_ligand_detail = session.get(CognateLigand, cognate_ligand.cognate_ligand_id)
            cognate_ligand_response = CognateLigandResponse(
                id=cognate_ligand_detail.id,
                name=cognate_ligand_detail.name,
                smiles=cognate_ligand_detail.smiles,
                xref=cognate_ligand_detail.xref,
                similarity=cognate_ligand.similarity
            )
            cognate_ligand_response_list.append(cognate_ligand_response)

        transplant_full_list.append(
            TransplantFull(
                transplant=transplant_response,
                cognate_ligands=cognate_ligand_response_list
            )
        )

    response_data = TransplantsListResponse(data=transplant_full_list, structure_data=structure)
    return response_data

@app.get("/search/", response_model=List[str])
def search_structures(session: Session = Depends(get_session), query: str = Query("")):
    results = session.exec(select(Structure.name).where(Structure.name.ilike(f"%{query}%")).limit(5)).all()
    return results

class StructureStorage(ABC):
    @abstractmethod
    def get(self, structure_id: str) -> BinaryIO:
        pass

class LocalStructureStorage(StructureStorage):
    def __init__(self, base_path: Path):
        self.base_path = base_path

    def get(self, structure_id: str) -> BinaryIO:
        path = self.base_path / f"{structure_id}_transplants.cif.gz"
        print(f"Loading structure from TEST {path}", flush=True)
        return path.open("rb")

class S3StructureStorage(StructureStorage):
    def __init__(self, bucket: str):
        self.bucket = bucket
        self.s3 = boto3.client("s3")

    def get(self, structure_id: str) -> BinaryIO:
        import io
        obj = self.s3.get_object(Bucket=self.bucket, Key=f"{structure_id}_transplants.cif.gz")
        return io.BytesIO(obj["Body"].read())

USE_S3 = os.getenv("USE_S3", "false").lower() == "true"

def get_storage() -> StructureStorage:
    if USE_S3:
        return S3StructureStorage(bucket="your-bucket-name")
    return LocalStructureStorage(base_path=Path("/app/cif-files/"))

@app.get("/structure/{structure_id}")
def get_structure(
    structure_id: str,
    chains: Optional[list[str]] = Query(default=None),
    storage: StructureStorage = Depends(get_storage),
):
    try:
        file_obj = storage.get(structure_id)
        print(f"Loading structure from {getattr(storage, 'base_path', 'unknown')}", flush=True)
    except FileNotFoundError:
        raise HTTPException(status_code=404, detail=f"Structure not found {getattr(storage, 'base_path', 'unknown')}")

    if chains:
        return StreamingResponse(
            filter_structure_chains(file_obj, chains),
            media_type="chemical/x-cif"
        )
    print(f"Loading structure from {getattr(storage, 'base_path', 'unknown')}", flush=True)
    return StreamingResponse(file_obj, media_type="chemical/x-cif")

def filter_structure_chains(file_obj: BinaryIO, allowed_chains: list[str]) -> BinaryIO:
    #allowed_chains = ["A"]
    with gzip.open(file_obj, 'rt') as f:
        doc = gemmi.cif.read_string(f.read())

    block = doc.sole_block()
    structure = gemmi.make_structure_from_block(block)
    model = structure[0]
    struct_conf = block.get_mmcif_category('_struct_conf.')
    struct_conf_type = block.get_mmcif_category('_struct_conf_type.')
    # Filter to select only requested chains
    to_be_deleted = []

    for idx, chain in enumerate(model):
        if chain.name not in allowed_chains:
            to_be_deleted.append(idx)

    for index in reversed(to_be_deleted):
        del model[index]

    # Convert back to mmCIF
    structure.update_mmcif_block(block)
    block.set_mmcif_category("_struct_conf.", struct_conf)
    block.set_mmcif_category("_struct_conf_type.", struct_conf_type)
    #remove problematic chem comp category
    block.find_mmcif_category('_chem_comp.').erase()
    # Compress and return as BinaryIO
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, mode="wb") as gz:
        gz.write(block.as_string().encode("utf-8"))
    output.seek(0)
    return output