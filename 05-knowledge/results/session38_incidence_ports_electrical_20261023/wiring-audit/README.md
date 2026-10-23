# Reproduce the independent Session38 audit

Read AUDIT.md for exact scope and results. Run from this directory:

    g++ -O3 -std=c++17 independent_wiring.cpp -o independent_wiring
    ./independent_wiring .
    python3 -B compare_ledgers.py
    python3 -B audit_graph_inventory.py
    python3 -B audit_class_metadata.py
    python3 -B independent_signatures.py
    python3 -B audit_scalar_repair.py
    python3 -B audit_local_transport.py
    python3 -B audit_circulants.py
    python3 -B audit_electrical.py
    python3 -B audit_deeper_lifts.py
    python3 -B audit_small_switch.py
    python3 -B audit_bridge_core.py
    python3 -B scan_converse_menus.py
    python3 -B audit_reflected_counterfeits.py

For independent optimized checks, repeat each Python command with `-O -B`.
Every check uses explicit exceptions, not removable assertions. Python scripts
write only their own output receipts. The C++ audit regenerates its own complete
class map and uncompressed ledgers. Those large intermediate files are omitted
from the publication whitelist and should not be committed after a replay.

The comparison scripts look for adjacent `wiring/` and `electrical/` packets;
for the original working directory they also accept adjacent `session38-wiring/`
and `session38-electrical/`. The input representative TSV is included here and
is independently certified, not trusted as a complete list without checking.
No producer Python module is imported. Primary ledger and electrical witness
files are data inputs only, compared against independent reconstruction or
integer linear algebra as explained in AUDIT.md.

The reflected-counterfeit verifier consumes explicit witness data from adjacent
`theory/residual_counterfeits.json` (working-directory fallback `session38-root/`).
It does not import the producer signature checker. The converse-menu scan
independently recomputes both orientations' signatures on all representatives.
The nine-symbol result is specific to the full mirrored transportation model;
the twelve-symbol result concerns ordinary class separation.
