from typing import Optional
from sqlmodel import Session, SQLModel, create_engine, select, text
from sqlalchemy import URL
from .models import Structure, Transplant  # Import both models
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
    TRANSPLANT_FILES_DIR = Path("/Users/matthewcrown/GitHub/AlphaCognate/app/cif-files")
    DB_PASSWORD = os.getenv("DB_PASSWORD")

    if DB_PASSWORD is None:
        raise ValueError("DB_PASSWORD environment variable is not set!")
    
    drop_all_tables(DB_NAME, DB_USER, DB_PASSWORD, host = DB_HOST)

    engine = create_db_engine(username=DB_USER, password=DB_PASSWORD, host=DB_HOST, database=DB_NAME)

    create_db_and_tables(engine)


    #here we need to use the alphacognate tsv files from the gzipped files to build the database. Find and process these from a specified directory.
    structures = []
    transplants = []
    for file in TRANSPLANT_FILES_DIR.glob("*.tsv.gz"):
        transplants_df = pd.read_csv(TRANSPLANT_FILES_DIR / file, sep = "\t")
        print(transplants_df.num_transplants.values[0].astype(int))
        print(transplants_df.columns)
        if "error" in transplants_df.columns:
            transplants_df = transplants_df.loc[transplants_df.error.isna()]
        if len(transplants_df) == 0:
            continue
        structures.append(Structure(
            name = transplants_df.accession.values[0],
            url = transplants_df.accession.values[0] + "_transplants.cif.gz",
            transplanted = True,
            num_transplants = int(transplants_df.num_transplants.fillna(0).values[0])
        ))

        for index, row in transplants_df.iterrows():
            print(row.transplanted_ligand_chain)
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
                struct_asym_id = row.transplanted_ligand_chain
            )
            transplants.append(transplant)

    print("done")

    with Session(engine) as session:
        clear_tables(session)  # Clear existing records before adding new data

        for structure in structures:
            session.add(structure)
        session.commit()

        # Fetch structures to create transplants for them
        statement = select(Structure)
        results = session.exec(statement).all()

        # Create some transplants for each structure
        for transplant in transplants:
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
