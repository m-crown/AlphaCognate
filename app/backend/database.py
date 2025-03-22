from typing import Optional
from sqlmodel import Field, Session, SQLModel, create_engine, select
from sqlalchemy import URL
from .models import Structure
from dotenv import load_dotenv
import os

def create_db_engine(username, password, host, database, dbtype = "postgresql"):
    
    url_object = URL.create(
        dbtype,
        username=username,
        password=password,
        host=host,
        database=database,
    )

    engine = create_engine(url_object)
    return engine

def create_db_and_tables(engine): #we pass an engine instance to the function
    SQLModel.metadata.create_all(engine)

def main() -> None:
    """
    database.py provides functions for creating a db connection and building a database. When run directly, creates the database and tables.
    """

    load_dotenv()  # Load environment variables from .env file
    DB_PASSWORD = os.getenv("DB_PASSWORD")

    if DB_PASSWORD is None:
        raise ValueError("DB_PASSWORD environment variable is not set!")
    
    engine = create_db_engine(username = "admin", password = DB_PASSWORD, host = "localhost", database = "postgres")
    create_db_and_tables(engine)

    structures = [Structure(name = f"test{x}", url=f"Matt{x}", transplanted=True) for x in range(0,100)]

    with Session(engine) as session:
        for structure in structures:
            session.add(structure)
        session.commit()
        statement = select(Structure)
        results = session.exec(statement)
        for structure in results:
            print(structure)



if __name__ == "__main__":
    main()





