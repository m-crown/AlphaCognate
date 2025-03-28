from fastapi import FastAPI, Depends, HTTPException, Query
from typing import Annotated, List
from sqlmodel import Session, select, func
from .database import create_db_engine
from .models import Structure, StructuresListResponse, Transplant, TransplantsListResponse
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
        session: SessionDep,
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
    session: SessionDep,
    structure_name: str = Query(..., description="Filter transplants by structure name")
) -> TransplantsListResponse:
    transplants = session.exec(
        select(Transplant).where(Transplant.structure_name == structure_name)
    ).all()
    structure = session.exec(select(Structure).where(Structure.name == structure_name)).first()
    
    return TransplantsListResponse(transplant_data=transplants,
                                   structure_data=structure)

@app.get("/search/", response_model=List[str])
def search_structures(session: Session = Depends(get_session), query: str = Query("")):
    results = session.exec(select(Structure.name).where(Structure.name.ilike(f"%{query}%")).limit(5)).all()
    return results