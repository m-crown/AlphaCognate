#!/usr/bin/env bash
mkdir -p "$2"/pcgDB
foldseek createdb --chain-name-mode 1 "$1" "$3"
foldseek createindex "$3" "$4"
