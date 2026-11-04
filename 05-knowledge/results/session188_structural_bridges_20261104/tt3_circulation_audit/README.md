# Independent TT3 circulation audit

Read REVIEW.md for the exact reviewed scope and source pins. Run `python audit.py` and `python -O audit.py` for independent physical reconstruction/counting; run `python formal.py` and `python -O formal.py` for the universal free-symbol Schur contraction. The programs require Python 3 and g++, install nothing, and use explicit checks.

The physical audit reads the frozen sibling `tt3_circulation/normal.json` and pins that packet's mathematical sources and paired receipts. The formal audit is standalone. Neither imports producer algorithms.

Exact frozen allowlist: REVIEW.md, README.md, audit.py, physical.cpp, normal.json, optimized.json, formal.py, formal.normal.json, formal.optimized.json, SHA256SUMS.
