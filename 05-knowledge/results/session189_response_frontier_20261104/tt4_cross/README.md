# Crossed TT4 exact quartic packet

Read PROOF.md for the general two-kernel replacement identity, four constant symmetric certificates, positive fixed-ambient transfer, and exact source-relative scope. The full all-transitive-child integration remains conditional on separately reviewed incoming common-end separation.

Rebuild with standard installed tools:

    python3 prepare.py
    g++ -O3 -std=c++17 bank.cpp -o /tmp/s189_tt4_bank
    /tmp/s189_tt4_bank bank.json < graphs.txt
    python3 prepare_counts.py
    g++ -O3 -std=c++17 predecessor.cpp -o /tmp/s189_tt4_predecessor
    /tmp/s189_tt4_predecessor < independent_input.txt > independent_values.txt
    /tmp/s189_tt4_predecessor < physical_input.txt > physical_values.txt
    python3 prepare_certificates.py
    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json
    cmp normal.json receipt.json

inner_kernels.json preserves every L/R coefficient for all12 marked triples in each of the five physical TT4 orders. bank.json preserves all36 complete named response vectors. individual_quartics.json retains all768 physical local-replacement differences; capacity_controls.json retains all four complete capacity matrices. physical_queries.json records all544 distinct actual grouped queries and their full augmented rows and physical masks. Source, quotient target, and local replacement are never identified without the explicit restoration test.

The generated large independent_input.txt is excluded from the publication allowlist and is recreated by prepare_counts.py. Binaries, caches, and logs are excluded. The included signed64-bit counting engines are safe here because the largest augmented tournament has12 names, so every Hamilton count and finite-difference sum is far below the signed64-bit bound. No software installation is needed.
