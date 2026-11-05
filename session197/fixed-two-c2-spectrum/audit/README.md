# Replaying the independent audit

The three checker sources and all receipts here are unchanged from the independent audit bound in ../REVIEW.md. The mathematical source snapshot is the parent directory of this same immutable repository checkpoint. To provide the checkers' expected local source/ directory without duplicating that snapshot in the repository, run from this audit directory:

    mkdir -p source
    cp ../twoport_core.py ../CONTROL.json ../spectrum.py ../families.py ../core.py source/
    python check_language.py
    python -O check_language.py
    c++ -std=c++17 -O2 check_graphs.cpp -o check_graphs
    ./check_graphs inverse-replay.csv
    cmp inverse-replay.csv inverse-normal.csv
    python check_inverse.py
    python -O check_inverse.py

The C++ program's stdout is compared to graphs-normal.json. The two Python programs' stdout is compared to language-normal.json and inverse-check-normal.json respectively. Debug/optimized and normal/optimized reference receipts are included. No binary or unrelated working-tree data is required.
