# Reproduce the independent odd-transitive and conversed-child review

REVIEW.md approves the two frozen producer proofs and states the exact scope. Source pins bind every producer publication file. Complete independently recounted physical banks and all additional family data are retained, including zeros and every named cubic/quartic coordinate.

Use standard Python and the existing C++17 compiler. Install nothing. The two source arguments below are the actual frozen producer directories; the example publication names can be adjusted without editing an audit source.

    python3 bank_audit.py ../odd_transitive ../collision_next > bank.normal.stdout
    python3 -O bank_audit.py ../odd_transitive ../collision_next > bank.optimized.stdout
    cmp bank.normal.stdout bank.optimized.stdout
    python3 kernel_audit.py > kernel.normal.stdout
    python3 -O kernel_audit.py > kernel.optimized.stdout
    cmp kernel.normal.stdout kernel.optimized.stdout
    python3 family_audit.py > family.normal.stdout
    python3 -O family_audit.py > family.optimized.stdout
    cmp family.normal.stdout family.optimized.stdout
    python3 queries_audit.py ../odd_transitive ../collision_next > queries.normal.stdout
    python3 -O queries_audit.py ../odd_transitive ../collision_next > queries.optimized.stdout
    cmp queries.normal.stdout queries.optimized.stdout
    python3 replay_audit.py ../odd_transitive ../collision_next > replay.normal.stdout
    python3 -O replay_audit.py ../odd_transitive ../collision_next > replay.optimized.stdout
    cmp replay.normal.stdout replay.optimized.stdout

The first program independently recounts all 101 producer graphs, not just selected coordinates. The kernel program uses subset-path polynomials rather than the producer's permutation/cut enumeration. The family program verifies all 4,096 root-containing original-module families on C3[TT5,TT2,TT2] in two ambient contexts. Its relabelings transport full named graphs and queries together and fix the exterior. The query program separately checks the actual grouped and anchor masks/counts. The final program replays both producer modes in temporary copies and verifies every regenerated allowlisted hash.

PUBLICATION_FILES.json is the exact audit allowlist. SHA256SUMS binds those files. Binaries, temporary replay copies, benchmarks, caches, and source-reading copies under input/ are excluded. INHERITED.json identifies the immutable inherited sources reviewed. No later producer changes or exploratory files are implicitly approved.
