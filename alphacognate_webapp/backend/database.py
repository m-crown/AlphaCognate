from typing import Optional
from sqlmodel import Session, SQLModel, create_engine, select, text
from sqlalchemy import URL

#fastapi needs relative import, running directly needs direct.
try:
    from models import Ligand, Structure, Transplant, CognateLigand, CognateLigandMapping  # For standalone script execution
except ImportError:
    from .models import Ligand, Structure, Transplant, CognateLigand, CognateLigandMapping   # For FastAPI when run inside a package

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

# Function to drop all tables in the database
def drop_all_tables(engine):
    with engine.connect() as conn:
        conn.execute(text("DROP SCHEMA public CASCADE;"))
        conn.execute(text("CREATE SCHEMA public;"))
        conn.commit()
        conn.close()

def create_db_and_tables(engine):
    SQLModel.metadata.create_all(engine)

def clear_tables(session):
    """Truncate tables to remove all existing records and reset primary keys."""
    session.exec(text("TRUNCATE TABLE transplants RESTART IDENTITY CASCADE;"))
    session.exec(text("TRUNCATE TABLE structures RESTART IDENTITY CASCADE;"))
    session.exec(text("TRUNCATE TABLE cognate_ligands RESTART IDENTITY CASCADE;"))
    session.commit()

def load_data(session, df_transplants, df_structures):
    #load structures into db
    for _, row in df_structures.iterrows():
        structure = Structure(name=row['accession'], 
                              url = row['accession'] + "_transplants.cif.gz",
                              runtime=row['runtime'], 
                              num_transplants=row['num_transplants'])
        session.add(structure)

    session.commit()  # Commit to save the structures in the database
    
    #group dataframe into each translpant and its associated coglig mappings

    grouped_transplants = df_transplants.groupby(['accession', 'ligand'])
    for _, group in grouped_transplants:
        ligand_data = group.iloc[0]  # Take the first row to get ligand transplant and ligand details
        ligand = session.exec(select(Ligand).where(Ligand.smiles == ligand_data['ligand_smiles'])).first()
        
        if not ligand:
            # If the ligand doesn't exist, insert it
            ligand = Ligand(
                het_code=ligand_data['ligand_het_code'],
                name = ligand_data['ligand_name'],
                smiles=ligand_data['cognate_mapping_smiles'],
            )
            session.add(ligand)
            session.commit()  # Commit to generate ligand id before referencing it in transplant

        #insert transplant data.
        transplant = Transplant(
            structure_name=ligand_data['accession'], #ref to structure
            ligand_id=ligand.id,  # Reference to the Ligand
            ligand_chain=ligand_data['ligand_chain'],
            ligand_residues=ligand_data['ligand_residues'],
            ec_list=ligand_data['cognate_mapping_ec_list'],
            tcs=float(ligand_data['tcs']),
            transplant_error=float(ligand_data['transplant_error']),
            foldseek_rmsd=float(ligand_data['foldseek_rmsd']),
            global_rmsd=float(ligand_data['global_rmsd']),
            local_rmsd=float(ligand_data['local_rmsd'])
        )
        session.add(transplant)

        for _, row in group.iterrows():
            #check if cognate ligand exists
            coglig = session.exec(select(CognateLigand).where(CognateLigand.smiles == row['cognate_mapping_smiles'])).first()
            if not coglig:
                # If the cognate ligand doesn't exist, insert it
                coglig = CognateLigand(
                    id = row['cognate_mapping_id'], #we already have ids from procoggraph so can specify these for consistency.
                    name=row['cognate_mapping_name'],
                    smiles=row['cognate_mapping_smiles'],
                    xref=row['cognate_mapping_xref']
                )
                session.add(coglig)
                session.commit() # Commit to generate cognate ligand id before referencing it in mapping
            #add the cognateligandmapping
            cognate_mapping = CognateLigandMapping(
                ligand_instance_id=transplant.id,  # Reference to the Transplant
                cognate_ligand_id=coglig.id,  # Reference to the Cognate Ligand
                similarity=row['cognate_mapping_similarity'],
            )
            session.add(cognate_mapping)

    session.commit()  # Commit to save all Transplants and Cognate Ligands

def main() -> None:
    """
    database.py provides functions for creating a db connection and building a database.
    When run directly, creates the database and tables.
    """
    
    load_dotenv()  # Load environment variables from .env file
    DB_USER = "admin"
    DB_HOST = "localhost"
    DB_NAME = "postgres"
    DB_PASSWORD = os.getenv("DB_PASSWORD")

    TRANSPLANTS_FILE = Path("/Users/matthewcrown/GitHub/AlphaCognate/app/cif-files/combined_transplants.tsv.gz")
    STRUCTURE_SUMMARY = Path("/Users/matthewcrown/GitHub/AlphaCognate/app/cif-files/combined_structure_summaries.tsv.gz")

    if DB_PASSWORD is None:
        raise ValueError("DB_PASSWORD environment variable is not set!")

    engine = create_db_engine(username=DB_USER, password=DB_PASSWORD, host=DB_HOST, database=DB_NAME)
    
    drop_all_tables(engine)
    create_db_and_tables(engine)

    # Create session to interact with DB
    with Session(engine) as session:
        # Load data from TSV files into pandas DataFrames
        df_transplants = pd.read_csv(TRANSPLANTS_FILE, sep="\t", compression='gzip', dtype={'ligand_residues': str})
        df_structures = pd.read_csv(STRUCTURE_SUMMARY, sep="\t", compression='gzip')

        # # Optionally clear the tables if you want to load fresh data
        # clear_tables(session)

        # Load data into the database
        load_data(session, df_transplants, df_structures)
    
    print("Data loaded successfully!")

if __name__ == "__main__":
    main()