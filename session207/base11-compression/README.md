# Base-eleven witness compression

For every attained Hamiltonian path count h>1, a tournament witness of order m
exists with 3^5 * 11^m <= 11^3 * h^5. PROOF.md gives the all-order argument,
explicit forward-spectrum stabilization window, cutoff algorithm and limits.
This is scalar replacement, not compression of an unchanged SCC core.

Dependencies: the cited Busch minimum and 13 explicit small witnesses.
Witness encoding: lexicographically ordered pairs i<j; bit one means i->j.
REVIEW.md records the independent audit. All arithmetic is integral.

From this directory:

    python3 verify.py > /tmp/session207-base11-result.json
    cmp RESULT.json /tmp/session207-base11-result.json
    python3 independent_audit.py > /tmp/session207-base11-independent.out
    cmp independent_audit.out /tmp/session207-base11-independent.out
    sha256sum -c SHA256SUMS

No external data or packages are required for replay. Both verification paths
check physical Hamiltonian paths directly by permutation enumeration.
