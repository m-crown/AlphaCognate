# AlphaCognate

AlphaCognate is a tool to identify the cognate ligand binding potential of an AlphaFold structure.

It adapts the AlphaFill methodology first described by Hekkelman et al. in AlphaFill (2021) - updated to use FoldSeek and an expanded set of ligands for transplanting, based upon those with identified cognate ligand matches in the ProCogGraph database.

AlphaCognate processes single chain AlphaFold predictions, and as such utilises a subset of the ProCogGraph database in which domain interactions are detected within a single chain.

AlphaCognate also makes use of the CATH domain annotations of AlphaFold structures from (XYZ) and will be updated to utilise information from The Encyleopedia of Domains (TED) when this becomes available.