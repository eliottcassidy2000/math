# The refined exact-six cycle core: one k2 three-cycle and an acyclic k3 graph

## Status

**FINITE-EXACT UNNUMBERED SIDECAR.**  This note freezes a structural fact about
the current refined six-clock mutation relations.  It is not a theorem-ID
dependency, an LRC(14) terminal, a physical cover, or an iteration theorem.

## Inheritance and reduction

The raw exact-six companion freezes the full mutation relation on all `3003`
six-clock bodies.  In both `k=2` and `k=3` it has exactly two nontrivial SCCs,
each a three-cycle:

```text
A = {(1,3,7,8,9,10), (2,4,5,7,9,12), (4,5,7,8,9,10)},
B = {(2,3,5,7,9,12), (3,7,8,9,10,12), (4,5,6,7,9,10)}.
```

The current refined relation uses the same strict open-cell target and the
same six-clock coverage predicate, but retains only body/divisor rows having
positive current occurrence weight.  Every refined edge is therefore a raw
edge: refinement deletes rows and cannot create or merge SCCs.  Consequently
every nontrivial refined SCC must lie inside `A` or `B`.  It is enough to
reconstruct every raw nonself edge internal to these two components and test
the six corresponding current row keys.

This reduction is exact and avoids replaying either full refined ledger.

## Targeted exact rows and weights

The sidecar independently reconstructs the six raw internal edges from the
strict endpoints and all open atoms.  It also verifies that their targets have
no cover by at most five pool clocks, so they genuinely belong to the
exact-six relation.  Each internal edge occurs at a unique divisor row:

| raw edge `F -> C` | `L` | `D` | `|S_D|` | refined k2 occurrence weight | refined k3 occurrence weight |
|---|---:|---:|---:|---:|---:|
| `(1,3,7,8,9,10) -> (4,5,7,8,9,10)` | 35280 | 17640 | 9984 | 58269 | 0 |
| `(2,3,5,7,9,12) -> (3,7,8,9,10,12)` | 17640 | 4410 | 3280 | 0 | 0 |
| `(2,4,5,7,9,12) -> (1,3,7,8,9,10)` | 17640 | 4410 | 3000 | 3199 | 0 |
| `(3,7,8,9,10,12) -> (4,5,6,7,9,10)` | 35280 | 17640 | 9984 | 58269 | 0 |
| `(4,5,6,7,9,10) -> (2,3,5,7,9,12)` | 17640 | 8820 | 4848 | 145541 | 0 |
| `(4,5,7,8,9,10) -> (2,4,5,7,9,12)` | 35280 | 17640 | 10116 | 522861 | 0 |

For `k=2`, all three edges of `A` survive.  In `B`, the first source row has
weight zero, leaving only the directed path

```text
(3,7,8,9,10,12) -> (4,5,6,7,9,10)
                    -> (2,3,5,7,9,12).
```

Thus the exact list of nontrivial refined `k=2` SCCs is `(A,)`.  For `k=3`,
all six cycle-source rows have weight zero, so neither raw SCC retains an
internal edge; the refined `k=3` relation is acyclic.

The zero weights are current-key absence, not a claim that the raw geometric
edge disappeared.  Conversely, every positive row above already survived the
at-most-five complement-clock subtraction because its target is exact-six.

## Scope and lost information

The SCC conclusion is now reproducible exact evidence, but it does not repair
the type loss of mutation.  An edge keeps the present body, divisor row, and
pointwise coverage witness while forgetting the next sector and divisor,
quotient-tail labels and residues, endpoint ownership transport, collision
data, and physical time.  Even the surviving `k=2` cycle therefore supplies
no lawful physical recursion.  The six complements plus seven inherited
clocks still reach the first open `7+6=13` boundary, so no LRC(14) conclusion
is claimed.

## Reproduction and security

Run

```bash
python3 -B 04-computation/lrc14_refined_exact_six_scc_sidecar_20260815.py
python3 -B -O 04-computation/lrc14_refined_exact_six_scc_sidecar_20260815.py
```

Both runs must reproduce the stored transcript byte for byte.  The sidecar
pins the raw relation script/output, fixed arrangement and support programs,
both refined current-ledger compositions, their direct row generators, and
the targeted k2 C++ engine.  It builds only the six body/divisor inputs.  The
k2 engine receives only their finite queries and must agree under `-O2` and
`-O3`; k3 evaluates only divisors `4410`, `8820`, and `17640`.  There is no
floating point, assertion-dependent truth gate, user input, or network call.
MISTAKE-401 repaired the original platform-dependent path serialization by
printing dependency paths in POSIX form; the graph data and semantic digest
did not change.

The LF-normalized hashes are

```text
script  cb0c02afe0d9bcf44f9790ea51c2499de41ce7f3b2df4571f56809153bade2f2
output  c7e711150aded2449064ac633c3316b7141081e232fb2553af55dd2fef9e28b1
```

The targeted engine transcript hash is
`859d75a7b7870c777438d6abee483c3458f761aca1e1e9a54179574385a399ef`,
and the SCC semantic hash is
`d3be3507d47946f2fa688f70faa572aad1312b8e2ecee30330af9e684bbbed31`.
