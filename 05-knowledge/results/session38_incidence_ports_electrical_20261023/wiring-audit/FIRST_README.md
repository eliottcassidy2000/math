# Independent full-wiring checkpoint

Requirements: C++17 compiler and Python3 standard library. From this directory:

    g++ -O3 -std=c++17 independent_wiring.cpp -o independent_wiring
    ./independent_wiring .
    python3 -B compare_ledgers.py
    python3 -O -B compare_ledgers.py
    python3 -B audit_graph_inventory.py
    python3 -O -B audit_graph_inventory.py

The C++ audit generates all its own uncompressed class assignments and ledgers.
The comparator reads the adjacent wiring/class_of_mask.u16le.gz with gzip
decompression and compares every one of its4,194,304 bytes. It also consumes
the adjacent compressed classes, ports, edges and endpoint-capacity TSVs and
the producer receipt. The local-workspace fallback directory is session38-wiring.
No raw producer map is needed. No producer Python module is imported.

Generated raw audit files and the executable are deliberately excluded from
publication. They are reproducible intermediates, pinned in the freeze receipt.
Do not commit them after a replay. Every Python check uses an explicit exception
and survives -O. See FIRST_WIRING_AUDIT.md for exact coverage and boundaries.
