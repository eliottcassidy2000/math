# Independent audit

[Aggregate verdict and boundaries](REPORT.md) ·
[Local-closure checkpoint audit](LOCAL_CLOSURE.md).

Run `python3 -B run_checks.py` from this directory or any other directory.
It replays every independent script normally and optimized and compares the
frozen receipts. No producer modules, repository data, or external packages
are imported. Only `kernel.py` is shared among these independent scripts.

`REVIEWED_INPUTS.json` pins producer manuscripts and associated code at final
audit. It is not a claim that all producer scripts were independently rerun;
our separate computation universes are listed in the report. `PUBLIC_FILES.txt`
is the curated deliverable list; bytecode and exploratory filesystem debris
are excluded.
