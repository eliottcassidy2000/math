# Session190: a nonovershooting local Hamilton-count patch

For every k>=1, the packet constructs strong tournaments with Hamilton
counts 29+157*2^k and 31+157*2^k. Four specified old arcs change, the physical
root and exterior remain fixed, and every intermediate count is below the
starting count. The two inner endpoints both have one starting Hamilton
path and the same named E2 response157, while their other responses differ.

`PROOF.md` gives the all-order construction, exact scope and two fixed-parent
obstructions. `verify.py` independently checks the literal patch, every
shortest rooted relabeling/flip path, and the larger twin controls.
`run_census.py` and `certify.cpp` reproduce the finite size boundary and the
all-order parent-profile exclusions using the inherited audited order-eight
bank. The bank is included as `bank8.json`; its original receipt and review
are included and cited.

This is not an arbitrary-target constructor, an indefinitely repeatable +2
operation, or a closed ternary return grammar. The global spectrum is open.
