# Reproduce the independent endpoint-obstruction audit

Run `python3 audit.py > normal.json` and `python3 -O audit.py > optimized.json`, then compare the outputs. The script uses the sibling `even_endpoint_obstruction` packet as frozen data, but imports none of its code. The installed C++17 compiler is used in a temporary directory.

The audit independently checks all endpoint vectors, proper modules, complete small-even census histograms, all 2,818 actual response values and every retained cubic/quartic coefficient. See `REVIEW.md` for scope.

Frozen allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.
