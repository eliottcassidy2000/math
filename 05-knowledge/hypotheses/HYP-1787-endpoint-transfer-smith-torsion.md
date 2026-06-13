---
id: HYP-1787
status: OPEN
source: endpoint_transfer_smith_s95.py
session: codex-2026-05-29
---

# HYP-1787: Endpoint Transfer Smith Torsion Split

## Statement

The tournament endpoint-transfer matrices may have no 2-primary Smith torsion:
all nonzero Smith invariant factors are odd.

Equivalently, the full GF(2) row rank phenomenon may lift to an integral
statement:

```text
coker(A_n^tournament) has no 2-primary obstruction on the row-rank part.
```

The even-graph endpoint-transfer matrices do have 2-primary Smith torsion, and
the number of odd Smith factors is exactly their GF(2) rank.

## Evidence

Exact Smith normal forms:

```text
tournament invariant factors:
2->3: [1]
3->4: [1, 1]
4->5: [1, 1, 1, 1]
5->6: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 15]

merged invariant factors:
2->3: [1]
3->4: [1, 1]
4->5: [1, 1, 1]
5->6: [1, 1, 1, 1, 1, 1, 1, 1, 3, 15]

even-graph invariant factors:
2->3: [1]
3->4: [1, 2]
4->5: [1, 1, 2]
5->6: [1, 1, 1, 1, 1, 1, 8]
6->7: [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 8]
```

Odd Smith factor counts:

```text
tournament: [1, 2, 4, 12]
merged:     [1, 2, 3, 10]
even graph: [1, 1, 2, 6, 8]
```

## Interpretation

The parity-rank separation is not merely a mod-2 artifact. In the even-graph
quotient it is visible as explicit 2-primary Smith factors. The tournament
quotient appears to route endpoint transfer through odd invariant factors only.

This matches the automorphism story:

- tournament automorphism groups have odd order;
- even-graph automorphism groups commonly carry powers of two;
- endpoint-transfer parity is preserved on the tournament side and torsionized
  on the graph side.

## Next Step

Compute tournament and merged Smith forms for `6->7` if a faster integer-SNF
route is available. The matrices are `56x456` and `34x272`; sympy may be too
slow, but modular/minor methods might detect whether any even invariant factor
appears.
