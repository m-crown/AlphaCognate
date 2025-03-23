from typing import Optional, List
from pydantic import BaseModel
from sqlmodel import Field, SQLModel, Relationship

class Transplant(SQLModel, table=True):
    __tablename__ = "transplants"

    id: Optional[int] = Field(default=None, primary_key=True)
    structure_name: str = Field(foreign_key="structures.name")  # Foreign key to Structure
    date: str  # Example field (you might want to use datetime)
    success: bool  # Example field to track transplant success

    # Relationship: Each transplant is linked to one structure
    structure: Optional["Structure"] = Relationship(back_populates="transplants")

class Structure(SQLModel, table=True):
    __tablename__ = "structures"

    name: str = Field(primary_key=True)
    url: str
    transplanted: bool
    #eventually we could have new metadata on structure added here?

    # Relationship: A structure can have multiple transplants
    transplants: List["Transplant"] = Relationship(back_populates="structure")

    def __repr__(self) -> str: #repr is a special instance that defines how an instance of a class is represented when you print it
        return f"<Structure(name={self.name}, url={self.url}, transplanted={self.transplanted})>"
    

# Define a Pydantic model for the meta information
class MetaResponse(BaseModel):
    totalRowCount: int

# Define a Pydantic model for the entire API response
class StructuresListResponse(BaseModel):
    data: List[Structure]
    meta: MetaResponse

# Define a Pydantic model for the entire API response
class TransplantsListResponse(BaseModel):
    data: List[Transplant]