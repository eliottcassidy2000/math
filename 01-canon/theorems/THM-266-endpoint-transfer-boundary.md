# THM-266: Endpoint Transfer Boundary for Fixed-Path Quotients

**Type:** Theorem
**Certainty:** 5 — PROVED for row/column and parity-boundary laws
**Status:** PROVED; even-graph rank defect verified computationally
**Filed by:** codex-2026-05-29-S1
**Tags:** #endpoint-transfer #tiling-cube #quotient #parity #even-graphs

---

## Statement

Let `X_n = Q_{m(n)}` be the fixed-path tiling cube for tournaments on
vertices `0,...,n-1`, where `m(n)=C(n-1,2)`. Endpoint insertion sends

```text
X_n x Q_{n-1} -> X_{n+1}
```

by keeping the old tiling, fixing the new path edge `n-1 -> n`, and choosing
the `n-1` new endpoint arcs between `n` and `0,...,n-2`.

For any finite quotient maps `q_n : X_n -> C_n`, define the fiber size

```text
F_n(C) = |q_n^{-1}(C)|
```

and let `A_n(C,D)` count endpoint-extension children of parent class `C`
landing in child class `D`.

Then:

```text
sum_D A_n(C,D) = 2^(n-1) F_n(C)
sum_C A_n(C,D) = F_{n+1}(D)
```

Consequently over `GF(2)`, for `n >= 2`,

```text
XOR_C A_n(C, -) = (D |-> F_{n+1}(D) mod 2).
```

That is, the parity boundary of the endpoint-transfer rows is exactly the
odd-fiber indicator at the child level.

## Proof

For the row sum, each tiling in parent fiber `C` has exactly `2^(n-1)`
endpoint signatures: the orientations of the new endpoint arcs to
`0,...,n-2`. Summing over all parent tilings gives
`2^(n-1) F_n(C)`.

For the column sum, every child tiling in class `D` has a unique endpoint
deletion parent: delete vertex `n` and forget the `n-1` endpoint choices.
Partitioning the child fiber by the quotient class of this unique parent gives
`sum_C A_n(C,D) = F_{n+1}(D)`.

Reducing modulo 2, each row sum is even for `n >= 2`. Therefore the XOR of all
rows has coordinate `D` equal to the column sum modulo 2, namely
`F_{n+1}(D) mod 2`.

## Even-Graph Corollary

For the even-graph quotient, the child fiber of an unlabeled even graph `G` is
the labeled orbit size:

```text
F_{n+1}(G) = (n+1)! / |Aut(G)|.
```

Therefore the endpoint parity boundary is not all child even-graph classes.
It is the subset whose automorphism group absorbs the full 2-adic valuation of
`(n+1)!`.

Exact enumeration through transition `6 -> 7` gives:

```text
transitions n->n+1: 2->3, 3->4, 4->5, 5->6, 6->7
even parent counts: [1, 2, 3, 7, 16]
GF(2) ranks:        [1, 1, 2, 6, 8]
rank defects:       [0, 1, 1, 1, 8]
```

So full row rank fails for the even-graph quotient in the tested range, unlike
the tournament quotient and complement-merged tournament quotient from
HYP-1773.

## Verification Record

Exact row/column sums, parity boundaries, and even-graph rank data are stored
in:

- `04-computation/even_graph_endpoint_transfer_s95.py`
- `05-knowledge/results/even_graph_endpoint_transfer_s95.out`
- `05-knowledge/hypotheses/HYP-1774-even-graph-endpoint-rank-defect.md`
