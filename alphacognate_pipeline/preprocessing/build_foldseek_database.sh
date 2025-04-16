#!/usr/bin/env bash
mkdir -p "$2"/pcgDB
foldseek createdb --chain-name-mode 1 "$1" "$2"/pcgDB/pcgDB
foldseek createindex "$2"/pcgDB/pcgDB "$2"/pcgDB/pcgDB_index
