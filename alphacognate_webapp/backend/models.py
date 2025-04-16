from typing import Optional, List
from pydantic import BaseModel
from sqlmodel import Field, SQLModel, Relationship


#TODO : DRAW OUT THE MODEL FOR THIS INTRERACTIONS TO UNDERSTAND IT BETTER

from typing import Optional, List
from sqlmodel import SQLModel, Field, Relationship


class Structure(SQLModel, table=True):
    name: str = Field(primary_key=True)
    url: str
    runtime: float
    num_transplants: int

    transplants: List["Transplant"] = Relationship(back_populates="structure")


class Ligand(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    het_code: str
    smiles: str

    instances: List["Transplant"] = Relationship(back_populates="ligand")


class Transplant(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    structure_name: str = Field(foreign_key="structure.name")
    ligand_id: int = Field(foreign_key="ligand.id")
    ligand_chain: str
    ligand_residues: str
    ec_list: str

    tcs: float
    transplant_error: str
    foldseek_rmsd: float
    global_rmsd: float
    local_rmsd: float

    structure: Optional[Structure] = Relationship(back_populates="transplants")
    ligand: Optional[Ligand] = Relationship(back_populates="instances")
    mappings: List["CognateLigandMapping"] = Relationship(back_populates="ligand_instance")


class CognateLigand(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    smiles: str
    xref: str

    mappings: List["CognateLigandMapping"] = Relationship(back_populates="cognate_ligand")


class CognateLigandMapping(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    ligand_instance_id: int = Field(foreign_key="transplant.id")
    cognate_ligand_id: int = Field(foreign_key="cognateligand.id")
    similarity: float

    ligand_instance: Optional[Transplant] = Relationship(back_populates="mappings")
    cognate_ligand: Optional[CognateLigand] = Relationship(back_populates="mappings")

# Define a Pydantic model for the meta information
class MetaResponse(BaseModel):
    totalRowCount: int

class StructuresListResponse(BaseModel):
    data: List[Structure]
    meta: MetaResponse


# Define a Pydantic model for the cognate ligand response - combines cognate ligand and cognate ligand mapping.
class CognateLigandResponse(BaseModel):
    id: int
    name: str
    smiles: str
    xref: str
    similarity: float

class TransplantResponse(BaseModel):
    transplant_id: int
    structure_name: str
    ligand_id: int
    ligand_chain: str
    ligand_residues: str
    ec_list: str
    tcs: float
    foldseek_rmsd: float
    global_rmsd: float
    local_rmsd: float
    hetcode: str
    name: str
    smiles: str
    
class TransplantFull(BaseModel):
    transplant: TransplantResponse
    cognate_ligands: List[CognateLigandResponse]

class TransplantsListResponse(BaseModel):
    data: List[TransplantFull]
    structure_data: Structure