# Session 140 — Hadwiger–Nelson recovery and exact representation bridges

## Main deliverables

- [`PROOF.md`](PROOF.md): a functional graph is a three-colorable pseudoforest; the standard Collatz graph is exactly three-colorable; every countable locally finite simple pseudoforest has an induced Euclidean unit-distance realization (crossings and arbitrarily close vertices allowed); continuum many noncongruent labeled Collatz placements; and a necessary-and-sufficient cycle-displacement certificate for lifting modular edge labels to exact positions.
- [`HISTORICAL_RECOVERY.md`](HISTORICAL_RECOVERY.md): seven actual historical HN methods, not seven new discoveries, with the corrected lattice cap, Q(ζ55) prime classification, square-oracle/descent boundary, polygon-moment status split, exact hexagonal metric lift, and current external-source status.
- [`verify.py`](verify.py): standard-library-only exact controls; normal and optimized outputs are identical.
- [`SOURCE_PINS.json`](SOURCE_PINS.json): thirteen targeted small Git sources, their original commit and blob identities, sizes and SHA-256 hashes. Copies are in `sources/`.

The strongest lattice-growth bound belongs to the parallel Collatz lane:
`|C^{-5r}(m)|≥2^r` for `3∤m`. Our slower odd-preimage construction is retained
only as an independently derived control, not a duplicate primary advance.

## Verified scope

`verification.json` and `verification.optimized.json` record **50,685 checks**:
all 3,413 functional maps on at most five vertices; all six odd-inverse
residue classes; inverse levels through depth12; the triangular ball formula
through radius12; all30 rooted K7 quotient triangles (12 exactly closed,
18 closed only modulo the ideal); the forty classes modulo55; exact small-
field torus/Kloosterman trace distributions and open-chain quartic counts
at5,7,11,13; all36 pair distances of the recovered Γ9 placement; and all
source hashes. No optimization-removable asserts are used.

The infinite results have written proofs; finite checks are controls.
Inherited Selmer computations, higher modularity, and the external HN six
lower bound were not independently rerun or formalized here. No result
here proves Collatz convergence or determines the exact plane chromatic
number. Generic realization is an existence proof, not an effective
procedure discovering unknown components/cycles.

## Replay

```sh
python3 verify.py > verification.replayed.json
python3 -O verify.py > verification.replayed.optimized.json
cmp verification.replayed.json verification.replayed.optimized.json
```

Produced in `/tmp/session140_research/hadwiger_nelson`, using repository reads
from `/tmp/session139-source`; no Git mutation or remote write was performed
by this lane.
