from typing import Optional, List
from pydantic import BaseModel
from sqlmodel import Field, SQLModel

class Structure(SQLModel, table=True):
    __tablename__ = "structures"

    id: Optional[int] = Field(default = None, primary_key=True)
    name: str
    url: str
    transplanted: bool
    #eventually we could have new metadata on structure added here?

    def __repr__(self) -> str: #repr is a special instance that defines how an instance of a class is represented when you print it
        return f"<Structure(id={self.id}, name={self.name}, url={self.url}, transplanted={self.transplanted})>"
    

# Define a Pydantic model for the meta information
class MetaResponse(BaseModel):
    totalRowCount: int

# Define a Pydantic model for the entire API response
class StructuresListResponse(BaseModel):
    data: List[Structure]
    meta: MetaResponse