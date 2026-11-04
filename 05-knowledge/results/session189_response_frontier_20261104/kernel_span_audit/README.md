# Independent kernel-span audit

REVIEW.md explains the all-order argument, exact novelty relative to incoming E4, and the physical local-converse scope. The audit retains full source/target response polynomials, not only their contrasts.

The standard-library scripts require installed gcc/g++ but no downloads or package installation. From this directory, with the two frozen producers at the default sibling paths:

    python3 replay_audit.py > replay.normal.stdout
    python3 -O replay_audit.py > replay.optimized.stdout
    cmp replay.normal.stdout replay.optimized.stdout
    python3 audit.py > normal.stdout
    cp receipt.json receipt.normal.json
    python3 -O audit.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp receipt.json receipt.normal.json
    python3 query_audit.py > query.normal.stdout
    python3 -O query_audit.py > query.optimized.stdout
    cmp query.normal.stdout query.optimized.stdout
    python3 boundaries.py > boundary.normal.stdout
    python3 -O boundaries.py > boundary.optimized.stdout
    cmp boundary.normal.stdout boundary.optimized.stdout
    sha256sum -c SHA256SUMS

For relocated producer packets, pass the kernel-span directory as the first argument to audit.py and replay_audit.py; pass the odd-transitive directory as replay_audit.py's second argument. The latter replays a byte-identical copy in a temporary directory with the dependency's original relative name, without modifying either producer.

polynomial_certificates.json retains exact rational span inverses and finite-jet certificates. moment_bank.json retains all independently calculated raw jets and maximum matchings. physical_banks.json and boundary_banks.json retain complete named coefficient banks, including zeros. grouped_queries.txt, grouped_query_values.txt, and grouped_certificates.json retain every actual queried graph, physical mask, count, and finite difference.

PUBLICATION_FILES.json is the exact release allowlist. Executables and temporary directories are excluded. The E4 reference is a prior theorem, not a newly established claim of this packet. Generic detectability does not automatically imply this restricted quartic battery, nor an arbitrary-ambient root-family result.
