#!/usr/bin/env bash

find /raid/MattC/repos/ProCogGraphData/structures/ -name '*_bio-h.cif.gz' -print0 | xargs -0 -I {} ln -s {} .
foldseek createdb --chain-name-mode 1 pcg_assemblies/ pcgDB
foldseek createindex pcgDB pcgDB_index
