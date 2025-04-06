from fastapi import FastAPI, Depends, HTTPException, Query
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse
from typing import Annotated, List
from sqlmodel import Session, select, func
from .database import create_db_engine
from .models import CognateLigand, CognateLigandResponse, Structure, StructuresListResponse, Transplant, TransplantFull, Ligand, CognateLigandMapping, TransplantResponse, TransplantsListResponse
from dotenv import load_dotenv
import os
from fastapi.middleware.cors import CORSMiddleware


load_dotenv()  # Load environment variables from .env file
DB_PASSWORD = os.getenv("DB_PASSWORD")

if DB_PASSWORD is None:
    raise ValueError("DB_PASSWORD environment variable is not set!")

engine = create_db_engine(username = "admin", password = DB_PASSWORD, host = "localhost", database = "postgres")

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
    structure_name: str = Query(..., description="Filter transplants by structure name")
) -> TransplantsListResponse:
    
    structure = session.exec(select(Structure).where(Structure.name == structure_name)).first()
    
    transplants = session.exec(
        select(Transplant).where(Transplant.structure_name == structure_name)
    ).all()

    transplant_full_list = []
    for transplant in transplants:
        ligand = session.get(Ligand, transplant.ligand_id)
        transplant_response = TransplantResponse(
            transplant_id = transplant.id,
            structure_name = transplant.structure_name,
            ligand_id = transplant.ligand_id,
            ligand_chain = transplant.ligand_chain,
            ligand_residues = transplant.ligand_residues,
            ec_list = transplant.ec_list,
            tcs = transplant.tcs,
            foldseek_rmsd = transplant.foldseek_rmsd,
            global_rmsd = transplant.global_rmsd,
            local_rmsd = transplant.local_rmsd,
            hetcode = ligand.het_code,
            name = ligand.name,
            smiles = ligand.smiles
        )
        cognate_ligands = session.exec(
            select(CognateLigandMapping).where(CognateLigandMapping.ligand_instance_id == transplant.id)
        ).all()
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