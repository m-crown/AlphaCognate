import pandas as pd
from typing import Optional
from pydantic import BaseModel
import numpy as np
import pandas as pd

class TransplantResult(BaseModel):
    accession: Optional[str] = None
    query_structure: Optional[str] = None
    transplant_structure: Optional[str] = None
    foldseek_rmsd: Optional[float] = None
    global_rmsd: Optional[float] = None
    ligand: Optional[str] = None
    interaction: Optional[str] = None
    error: Optional[str] = None
    center_of_mass: Optional[float] = None

transplant_results = []  # Store all instances

for i in range(5):  # Example loop
    transplant = TransplantResult()  # Create a new instance

    # Simulate processing logic (sometimes fields are missing)
    if i % 2 == 0:
        transplant.accession = f"P12345_{i}"
        transplant.query_structure = f"query_{i}"
        transplant.transplant_structure = f"target_{i}"
        transplant.foldseek_rmsd = np.random.uniform(0.5, 2.0)
        transplant.ligand = f"LIG{i}"
        transplant.interaction = "H-bond"

    # Sometimes an error occurs (no data)
    else:
        transplant.error = "Processing failed"

    transplant_results.append(transplant)

#based on 

transplant_df = pd.DataFrame([s.__dict__ for s in transplant_results])
transplant_df = pd.DataFrame([s.__dict__ for s in [TransplantResult(error = "Processing failed")]])
print(transplant_df)