#!/usr/bin/env bash

#first check if venv exists - if it does activate it if not create it and install requirements
if [ -d ".venv" ]; then
    echo "Activating existing virtual environment"
    source .venv/bin/activate
else
    echo "Creating new virtual environment"
    python3 -m venv .venv --prompt alphacognateBackend
    echo "Activating alphacognateApp virtual environment"
    source .venv/bin/activate
    echo "Upgrading pip and installing depedencies"
    pip install --upgrade pip
    pip install -r dependencies.txt
fi
