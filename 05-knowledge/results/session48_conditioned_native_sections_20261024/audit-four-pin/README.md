# Four-cheap-freeze independent audit

Run with Python3 and a C++17 compiler (`g++`), using only the standard library:

    python -B audit.py
    python -O -B audit.py
    python -B verify_manifest.py

`audit.py` reads only the copied pinned inputs within this directory and compiles
`enumerate.cpp` in a temporary directory. It independently reconstructs the full
native menus, electrical rows, exact census, cost envelopes, and endpoint routes,
then compares against the saved exact results. `--write` regenerates those results.
No producer modules or result tables are imported. Checks are explicit exceptions,
so Python optimization does not disable them.

`AUDIT.md` gives the proof, quantified scope, and hostile controls.
`receipt.json` contains native, budget, witness, and endpoint summaries.
`all_four_freeze_cases.json` gives all4845 passing-source sets and failure witnesses.
`producer_comparison.json` records a post-computation comparison with the producer
snapshot, not a dependency needed to reproduce this audit.
`input_provenance.json` and `copied_archives.txt` identify all copied inputs.
`MANIFEST.json` authenticates every packet file except itself. Nested manifests,
if present, are included as ordinary files; none are globally ignored.

Normal, optimized, and copied-directory runs must all succeed. Original harmonic
saturation is an explicitly inherited mathematical theorem, as detailed in the
audit; finite original Laplacian rows are recomputed here.
