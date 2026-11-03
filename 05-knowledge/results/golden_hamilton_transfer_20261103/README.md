# Golden complements: a three-coordinate Hamilton recurrence

**PROVED:** the single-bridge composition formula, its coarsest future-H
observer, and the seven golden transfer matrices. **FINITE-EXACT:** complete
bounded-H reachability at all vertex orders, with small closed certificates
for the first family obstructions. Same-agent proof audit with independent
literal-graph algorithms; no claim of independent human review or priority.
The unrestricted spectrum conjecture, positive odd integers minus {7,21},
remains **OPEN**.

The requested recurrence closes in three nonnegative integers `(A,B,C)`.
Here H=A+B, B counts paths ending at the retained boundary vertex, and C
counts paths starting there plus spanning two-path covers with a path
starting there. A singleton changes the state by

    (A,B,C) -> (C,A,A+B+C).

There is a different exact matrix for each of Q4's five marked presentations.
See [RECURRENCE](RECURRENCE.md) for the matrices and all-order proof. The
three observations H now, after one singleton, and after two singletons
recover the state. Thus this is the coarsest future-H quotient for the
specified joins and has minimum rational linear dimension three. It is
not a three-element state space or an arbitrary-attachment carrier.

The principal discoveries are a usable recurrence and a precise obstruction:

| Source alphabet before actual path complement | Exact all-order bounded-H result |
|---|---|
| Singleton, C3, all five Q4 charts | First missing prime 37; first missing 7p is 133=7*19. Also misses 63 and343 at every order. |
| Above plus all50 strong five-vertex charts | Through10001 the only missing odd values are7,21,47,83,103,491,709,833,3227,4619. |

For the second alphabet the missing prime lane is47,83,103,491,709; its
only missing7p through10001 is3227=7*461. Absence here is relative to this
restricted construction, not the unrestricted tournament spectrum.
The [reachability proof and receipts](REACHABILITY.md) explain why these
are all-order conclusions despite a finite computation.

There are exactly two marked golden H55 sources. Neither has any
one-vertex attachment with H111. The inherited successful55->111 arrow
uses a different actual source. [CONNECTIONS](CONNECTIONS.md) relates this
to the earlier56-class/55-leakage/111-memory calculation and the square-sum
graph's111 edges at56 vertices, without identifying the different objects.

The aggregate Hamilton-path sum over marked golden sources has a rational
generating function with a degree12 denominator. This supplies a second,
scalar recurrence; it counts with presentation multiplicity, not isotypes.

The hostile12-vertex pair still has identical67 labelled cyclic triples
and H462915, automorphism orders81/27, and different named attachment
responses2390013/1754919. It is replayed, not discarded by the new quotient.

## Reproduction and inheritance

Run from the repository root, with Python and SymPy available:

```text
python -X utf8 05-knowledge/results/golden_hamilton_transfer_20261103/transfer.py
python -X utf8 05-knowledge/results/golden_hamilton_transfer_20261103/reach.py --limit 100000
python -X utf8 05-knowledge/results/golden_hamilton_transfer_20261103/extension.py
python -X utf8 05-knowledge/results/golden_hamilton_transfer_20261103/verify.py
python -O -X utf8 05-knowledge/results/golden_hamilton_transfer_20261103/verify.py --output verification_optimized.json
```

The proof and independent checks are detailed in [AUDIT](AUDIT.md).
The closest inherited sources, concept board, connection maps, and next
decisive test are in [HANDOFF](HANDOFF.md). No scarce theorem ID is reserved.
