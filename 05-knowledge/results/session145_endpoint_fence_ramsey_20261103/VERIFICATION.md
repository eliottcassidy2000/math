# Verification and scope

## Checkpoint 1: pooled interior-parity carrier

The self-contained discovery checker enumerates all old edge subsets and expands all
Hamiltonian words through order six; it checks all distinct signature pairs
through six and every exterior subset with two new vertices through old order
four. The independent audit re-enumerates all 32,768 six-vertex edge subsets
and verifies every one-birth pullback through old order five.

The finite lower-rank certificate has 571 exterior contexts, each using at most
12 new vertices. Its independent checker regenerates its own 571 columns and
reconstructs all matrix entries. It obtains rank 571 modulo 101 and 103, with
determinants 40 and 84. Deleting the final row or replacing it by the first gives
rank 570. Invalid exterior edges, duplicated edges, and loops are rejected.
The coefficient formula is independently checked by direct Hamiltonian-word
expansion through order five.

An independent C++ subset-DP counter evaluates 589,824 actual Boolean graphs:
18 contexts on all 32,768 six-vertex sources. Exact Walsh transforms validate
all signed quotient relations and forbidden coefficients. Separate exhaustive
three-vertex controls exhibit the failures of conditioning and adaptive choices.

Run from `carrier/`:

```sh
python3 -B verify_frozen_packet.py
python3 -B -O verify_frozen_packet.py
```

The verifier checks every frozen byte, then runs all pinned exact replays. The
actual-graph replay builds its executable in a temporary directory and requires
`g++` with C++17 support. The modular matrix verifier requires NumPy; its matrix
arithmetic is bounded integer modular arithmetic, not floating-point rank.
No script uses Python `assert` as a validity gate. Normal and optimized stdout
must match their pinned hashes exactly. The manifest excludes mutable overviews
and navigation; it does not include or claim to hash itself.

## Repository checks

Publication preparation uses an explicit file map and tree-only Git objects,
without an index, reset, or full repository hydration. Exact owned blobs,
predecessors, changed-path set, protected old packets, relative Markdown links,
and unchanged targeted startup blobs are checked. S111-plus is strict append.

The mandatory startup surface remains 179,983/180,000 bytes. This is a bounded
startup-budget and owned-document check, **not** a claim that the complete
repository-wide documentation checker ran. The verification scope is limited
to the explicit source pins, owned files, links, and startup budgets.

Published commit and tree identities are external receipts. They will be
verified against the actual remote by the publisher, not inferred from a local
preparation result.
