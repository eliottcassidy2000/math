# Crossed TT4 independent audit

Read REVIEW.md for the all-order argument and the approved source-relative synthesis. All full named coefficient banks, including zero entries, are retained.

From this directory, with the frozen producer in ../session189-tt4-cross:

    python3 audit.py > normal.stdout
    cp receipt.json receipt.normal.json
    python3 -O audit.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp receipt.json receipt.normal.json
    python3 boundary.py > boundary.normal.stdout
    python3 -O boundary.py > boundary.optimized.stdout
    cmp boundary.normal.stdout boundary.optimized.stdout
    sha256sum -c SHA256SUMS

A relocated producer path can be given as audit.py's sole argument. Python standard library and installed g++ suffice; no downloads or installations are needed. The verifier replays the producer in a temporary directory and does not modify it. The exact publication allowlist excludes executables, caches, and temporary files.

literal_banks.json reproduces all36 producer physical graphs. arbitrary_replacement_banks.json and controls retain twelve nontransitive replacement graphs and378 separate general-identity checks. q2_boundary_banks.json covers the minimal two-child quotient. capacity_transfer.json retains full matrices and the nontrivial cyclic ambient's exact excursion weights. physical_queries.in and physical_query_values.txt retain every independently counted actual grouped query.
