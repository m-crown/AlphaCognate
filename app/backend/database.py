from typing import Optional
from sqlmodel import Session, SQLModel, create_engine, select, text
from sqlalchemy import URL

#fastapi needs relative import, running directly needs direct.
try:
    from models import Structure, Transplant  # For standalone script execution
except ImportError:
    from .models import Structure, Transplant  # For FastAPI when run inside a package

from dotenv import load_dotenv
import os
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
import pandas as pd
from pathlib import Path

def create_db_engine(username, password, host, database, dbtype="postgresql"):
    url_object = URL.create(
        dbtype,
        username=username,
        password=password,
        host=host,
        database=database,
    )
    engine = create_engine(url_object)
    return engine

def drop_all_tables(dbname, user, password, host="localhost"):
    """Drops all tables in the database without dropping the database itself."""
    conn = psycopg2.connect(dbname=dbname, user=user, password=password, host=host)
    conn.autocommit = True
    cur = conn.cursor()

    # Execute the SQL block to drop all tables
    cur.execute("""
        DO $$ 
        DECLARE 
            r RECORD;
        BEGIN 
            FOR r IN (SELECT tablename FROM pg_tables WHERE schemaname = 'public') 
            LOOP 
                EXECUTE 'DROP TABLE ' || quote_ident(r.tablename) || ' CASCADE'; 
            END LOOP; 
        END $$;
    """)

    cur.close()
    conn.close()

def create_db_and_tables(engine):
    SQLModel.metadata.create_all(engine)  # Ensures both Structure and Transplant tables are created

def clear_tables(session):
    """Truncate tables to remove all existing records and reset primary keys."""
    session.exec(text("TRUNCATE TABLE transplants RESTART IDENTITY CASCADE;"))
    session.exec(text("TRUNCATE TABLE structures RESTART IDENTITY CASCADE;"))
    session.commit()

def main() -> None:
    """
    database.py provides functions for creating a db connection and building a database.
    When run directly, creates the database and tables.
    """
    
    load_dotenv()  # Load environment variables from .env file
    DB_USER = "admin"
    DB_HOST = "localhost"
    DB_NAME = "postgres"
    TRANSPLANTS_FILE = Path("/Users/matthewcrown/GitHub/AlphaCognate/app/cif-files/combined_transplants.tsv.gz")
    STRUCTURE_SUMMARY = Path("/Users/matthewcrown/GitHub/AlphaCognate/app/cif-files/combined_structure_summaries.tsv.gz")
    DB_PASSWORD = os.getenv("DB_PASSWORD")

    if DB_PASSWORD is None:
        raise ValueError("DB_PASSWORD environment variable is not set!")
    
    drop_all_tables(DB_NAME, DB_USER, DB_PASSWORD, host = DB_HOST)

    engine = create_db_engine(username=DB_USER, password=DB_PASSWORD, host=DB_HOST, database=DB_NAME)

    create_db_and_tables(engine)


    #here we need to use the alphacognate tsv files from the gzipped files to build the database. Find and process these from a specified directory.
    structures = []
    transplants = []
    transplants_df = pd.read_csv(TRANSPLANTS_FILE, sep = "\t")
    #once we work out why structure info is in transplants can remove this - error should nto be present in trasnplants file.
    if "error" in transplants_df.columns:
        transplants_df = transplants_df.loc[transplants_df.error.isna()]
    for index, row in transplants_df.iterrows():
        if "tcs" in row.index:
            if pd.notna(row.tcs):
                tcs = row.tcs
            else:
                tcs = 0
        else:
            tcs = 0
        transplant = Transplant(
            structure_name=row.accession,
            ligand = row.ligand,
            tcs = tcs,
            struct_asym_id = row.ligand_chain
        )
        transplants.append(transplant)

    print("done transplants")

    structures_df = pd.read_csv(STRUCTURE_SUMMARY, sep = "\t")
    for index, row in structures_df.iterrows():
        structures.append(Structure(
            name = row.accession,
            url = row.accession + "_transplants.cif.gz",
            runtime = row.runtime,
            num_transplants = int(row.num_transplants)
        ))

    with Session(engine) as session:
        clear_tables(session)  # Clear existing records before adding new data
        structure_names = []
        for structure in structures:
            structure_names.append(structure.name)
            if structure.name == -1:
                print(structure)
            session.add(structure)
        session.commit()

        # Fetch structures to create transplants for them
        statement = select(Structure)
        results = session.exec(statement).all()

        # Create some transplants for each structure
        for transplant in transplants:
            if transplant.structure_name not in structure_names:
                print(transplant)
            if not transplant.structure_name:
                print(transplant)
            session.add(transplant)

        session.commit()

        # Print all structures and their transplants
        # for structure in results:
        #     print(structure.name)
        #     transplants = session.exec(select(Transplant).where(Transplant.structure_name == structure.name)).all()
        #     for transplant in transplants:
        #         print("  ->", transplant)

if __name__ == "__main__":
    main()
