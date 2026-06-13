---
id: HYP-1774
status: CONFIRMED_SMALL_N
source: even_graph_endpoint_transfer_s95.py
session: codex-2026-05-29
---

# HYP-1774: Even-Graph Endpoint Rank Defect

## Statement

The endpoint-transfer row/column theorem is quotient-agnostic, but full
`GF(2)` row rank is not.

Formalized in `01-canon/theorems/THM-266-endpoint-transfer-boundary.md`.

For the even-graph quotient of the fixed-path tiling cube:

```text
Σ_D A_n(E,D) = 2^(n-1) fiber_n(E)
Σ_E A_n(E,D) = fiber_{n+1}(D)
```

still holds exactly. However, unlike the tournament quotient, the parity
transfer matrix loses row rank once `n >= 3`.

## Evidence

Exact endpoint-transfer enumeration:

```text
transitions n->n+1: 2->3, 3->4, 4->5, 5->6, 6->7
even parent counts: [1, 2, 3, 7, 16]
GF(2) ranks:        [1, 1, 2, 6, 8]
rank defects:       [0, 1, 1, 1, 8]
```

The child boundary is not "all child classes" as in the tournament quotient.
It is the subset of even-graph classes with odd labeled orbit size:

```text
odd-orbit child boundary: [2, 2, 4, 4, 16]
```

## Interpretation

Tournament fixed-path fibers are all odd because:

```text
F(C) = H(C)/|Aut(C)|
```

and both `H(C)` and `|Aut(C)|` are odd.

Even-graph fibers are orbit sizes:

```text
n! / |Aut(G)|.
```

Their parity is controlled by whether `Aut(G)` contains the full 2-adic
valuation of `n!`. Thus even-graph transfer parity sees highly symmetric graph
orbits, not every child class.

## Consequence

The endpoint-transfer recursion separates the two lenses:

- tournament quotient: parity-injective through tested levels;
- even-graph quotient: coarser cycle-space lens forgets parity information.

So recursive parity-injectivity is a tournament-side structure, not a generic
property of tiling-cube quotients.
