# Exact universes, independent paths and reproduction

The scripts use Python's standard library; the inherited metagraph replay
also uses a C++17 compiler. Assertions are explicit RuntimeError checks,
not removable Python `assert` statements. There is no randomness.

Run from the repository root:

```
python -B 05-knowledge/results/frobenius_lca_metagraph_20261028/verify.py --replay --replay-metagraph
python -B -O 05-knowledge/results/frobenius_lca_metagraph_20261028/verify.py
```

Each producer also runs directly and prints its JSON receipt, e.g.
`python -B 05-knowledge/results/frobenius_lca_metagraph_20261028/fields.py`.
The `--output PATH` option saves it. `counting.py` defaults to orders 0..50.
The replay does not fetch the network or alter the stored receipts.

## Scope and controls

| computation | complete universe and inherited filters | positive / hostile control | independent path |
|---|---|---|---|
| field forms | all Frobenius-invariant alternating forms in F16,F64,F256; all subfield lines; nondegeneracy explicitly checked | trace form / base line does not select it; F64 period-two ambiguity | polynomial-basis nullspace versus normal-basis circulant polynomial gcd |
| field symmetry | every Frobenius power on all nonzero points of these fields | F64 order-three powers / all nonidentity F256 powers | odd-cycle formula versus direct oriented-pair orbit traversal |
| phase queries | every phase and every individual basis-port entry; every normal generator: 8,24,128 | explicit parity decoders / raw-query lower bounds | exact decision-tree recursion plus elementary weight bounds in F2 |
| different orbit pairs | all4,16,256 ordered normal-orbit pairs in F16,F64,F256 | cross-basis parity decoders / full-period nonnormal collisions | replay incoming nonadaptive census; add adaptive optimization |
| mixed periods | all 64 ordered triples of degrees dividing6, and all64 dividing8; all8 pair graphs; every phase tuple | maximum-period anchor / degree-two anchor loses a lift | complete common-shift orbit enumeration versus prime-layer connectivity |
| affine cells | every nonempty affine subset of the fixed-path spaces at orders1..4; every permutation guard; all229 joins with child orders<=4 and total<=6 | symbolic affine rules / ab,bc lose global parity; blue-loop H=11/13 | direct sets compared with affine spans |
| 81/27 pair | all8 C3 decorations by A=C3->1 or its converse; actual groups for000/001 | constant / nonconstant decorations; equal local stabilizers but27/0 transporters | Hamilton subset DP; adjacency-checked automorphism search versus cyclic stabilizer formula; every element of G minus H checked fixed-point-free |
| reset count | every fixed-path presentation through n6 and every path permutation | classwise constant rows / trace(R tau) differs from class count | reset trace versus independent odd-partition Burnside sum |
| scalar count | every n=0..50 | exact rational integrality; strong-component recurrence | frozen OEIS numeric table plus metagraph overlap |
| actual metagraph | all presentations through n8; all legal right relocations; every literal complement pair | ordinary/converse lift and blue/black totals | rebuild existing Session67 C++; compare its original census and independent Burnside counts |
| incoming order-nine bank | every saved representative's relation and exact encoding; counts through9; entire fibre histograms through8 | no diagonal, exactly one arc per pair, no duplicate encodings | incoming root-switch construction versus E2 counts and native metagraph histograms; no independent n9 isomorphism census |

All five Python producers have equal normal and optimized receipts. The
separate field audit uses the same elementary field multiplication routine,
but a different classification and coordinate system; it is not an independent
field-arithmetic library. The phase-query lower bounds are proved directly
in the report. This is a same-session audit, not a claim of review by another
agent or mathematician.

`SOURCES.json` pins the inherited repository texts and enumerator by Git blob
and LF-normalized SHA256. `SHA256SUMS` seals every file in this packet except
itself, also normalizing CRLF to LF. The external numeric table records its
source and retrieval date. External counts are a validation oracle, never
an input to `counting.py` or the C++ metagraph enumerator.

## What these checks do not establish

- No literal enumeration at n>=9 was performed in this session. The incoming
  order-nine construction is cited and its saved bank inspected, not regenerated.
- No all-degree single-line reconstruction theorem is inferred from three fields.
- No identification of arbitrary tournament port structures with field forms.
- No all-order classification by strong-atom response histograms.
- No polynomial-size closure of affine mixtures under arbitrary reset,
  complement and variable prime quotient substitution.
- No weakening of the shared-alignment or run-multiplicity requirements.

The all-order mathematical statements are the proofs in [THEORY.md](THEORY.md),
with inherited inputs identified there. Finite receipts have their exact
universes above and do not replace those proofs.
