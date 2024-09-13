
import pandas as pd

cath_af = pd.read_csv("/raid/MattC/repos/ProCogGraphData/analysis/notebooks/cath-v4_3_0.alphafold-v2.2022-11-22 (1).tsv", sep = "\t")
cath_af["af_accession"] = cath_af["domain_ID"].str.extract("^af_([A-Za-z0-9]+)_")
cath_af[["domain_start", "domain_end"]] = cath_af["domain_ID"].str.extract("^af_[A-Za-z0-9]+_(\d+)_(\d+)", expand = True)
cath_af[["domain_start", "domain_end"]] = cath_af[["domain_start", "domain_end"]].astype("int")
cath_af["domain_profile"] = cath_af["sfam_id"] + ":" + cath_af["domain_start"].astype("str") + "-" + cath_af["domain_end"].astype("str")
cath_af = cath_af.sort_values(["af_accession", "domain_start"])
cath_af_grouped = cath_af.groupby("af_accession").agg({"domain_profile": list}).reset_index()
cath_af_grouped["domain_profile"] = cath_af_grouped["domain_profile"].str.join(";")
cath_af_grouped.to_csv("cath_alphafold_domains.tsv.gz", sep = "\t", compression = "gzip", index = False)#.head(20)

