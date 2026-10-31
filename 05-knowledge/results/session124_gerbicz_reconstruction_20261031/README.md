# Gerbicz square-sum construction: primary-certificate replication

**Status: independently verified classical result, not a new solution.**
Robert Gerbicz's2018 construction supplies a Hamiltonian cycle in the graph
on[1,N], with distinct vertices adjacent when their sum is square, for every
N>=32. It also supplies Hamiltonian paths for every N>=25. This packet
independently checks its finite construction certificates and explains the
universal induction. It does not infer small-order nonexistence from empty
table rows. Other ordinary/seeded proof chains already in this repository
remain distinct; this packet specifically replicates the primary data.

Read [the proof](GERBICZ_RECONSTRUCTION.md) and [independent audit](AUDIT.md).
The original method and data are Gerbicz's; the parsers, checkers, proof
exposition and integration controls in this packet were independently written.

## Required external input — not offline self-contained

The original archive/C source and extracted certificate arrays are **not
included**: no redistribution license was found. Obtain the source separately
from the public source pinned in [SOURCE.json](SOURCE.json). The provenance
chain is Berlekamp–Guy, *Fibonacci Plays Billiards*,
[arXiv:2002.03705, appendix/reference2](https://arxiv.org/html/2002.03705),
which links Gerbicz's archive. His author comment at
[OEIS A090461](https://oeis.org/A090461) records the2018 theorem.

Inspect the ZIP and extract its single file `squares.c` as data. Do not compile
or run it for this verification. Both checkers require its exact SHA256 before
parsing only the three integer initializer arrays. They never import, compile,
or execute external C. Neither checker trusts the original program's runtime
assertions or weaker aggregate-permutation check.

## Portable replay

Python3 standard library only. From any directory, use explicit paths:

```
python /path/to/packet/verify_gerbicz_data.py /path/to/squares.c
python -O /path/to/packet/verify_gerbicz_data.py /path/to/squares.c
python /path/to/packet/check_independent.py /path/to/squares.c
python -O /path/to/packet/check_independent.py /path/to/squares.c
python /path/to/packet/construct_primary.py /path/to/squares.c
python -O /path/to/packet/construct_primary.py /path/to/squares.c
```

The producer parses literal integer arrays through a strict character whitelist
and JSON. The independent auditor instead uses a recursive brace/comma/integer
grammar and imports no producer code. Both verify1992 base pairs at41..2032,
196 glue words and their98 paired parity conditions. These finite obligations,
not extrapolation, prove the all-order induction starting2033. Literal cycles
at32..40 close the lower finite boundary. The auditor also tests206 expanded
paths and five deliberately corrupted transition tables. Our constructor
checks four expanded pairs including one at1,000,000. All checks use explicit
exceptions and remain active with optimization.

Six JSON replay receipts record normal/-O results. The portable replay was
run from `/tmp`, not from the packet directory. `release-check.json` records
file-whitelist and source-rejection checks; `SHA256SUMS` pins packet artifacts.
No external archive, C source, or extracted table is part of this packet.
