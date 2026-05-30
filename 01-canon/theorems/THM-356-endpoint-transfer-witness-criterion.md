---
id: THM-356
name: endpoint-transfer-witness-criterion
status: PROVED
date: 2026-05-29
session: codex-2026-05-29-s95
depends_on:
  - THM-266
  - THM-347
---

# THM-356: Endpoint Transfer Witness Criterion

## Terminology

Let `A` be an endpoint-transfer matrix over `F_2`, with rows indexed by parent
quotient classes and columns indexed by child quotient classes.

An odd child column `d` is **private** for row `c` if:

```text
A(c,d) = 1
A(c',d) = 0 for every c' != c.
```

The odd-support graph of `A` is the bipartite graph with edge `c-d` when
`A(c,d)=1`.

## Statement

1. If every parent row has a distinct private child column, then `A` has full
   row rank over `F_2`.

2. More generally, full row rank implies the odd-support graph has a matching
   covering all rows only after replacing support by a suitable nonzero minor.
   A full support matching alone does not imply full row rank over `F_2`.

3. Therefore endpoint parity-injectivity has three increasingly strong
   possible proof mechanisms:

   ```text
   private child witnesses  =>  triangular minor  =>  full row rank.
   support matching alone   does not suffice.
   ```

## Proof

For (1), order the rows and their distinct private columns together. The
submatrix on these rows and columns is the identity matrix, so the rows are
linearly independent.

For (2), support matching only proves the existence of a permutation of nonzero
entries. Over `F_2`, the determinant of the corresponding square submatrix may
still vanish because an even number of nonzero permutation terms can cancel.

For example,

```text
1 1
1 1
```

has a perfect support matching but rank `1`.

Thus a matching is a useful obstruction test but not a rank proof.

## Boundary Corollary

If `A` is an endpoint-transfer matrix satisfying the boundary law of THM-266
and `d` is a private odd child column, then `d` lies in the child fiber
boundary.

Indeed, the column sum of `d` is odd because exactly one row has an odd entry
in that column. By THM-266 this column-sum parity is the child fiber parity.

Consequently, for the complement-merged tournament quotient, every private odd
child witness is necessarily self-complementary by THM-347/HYP-1772. For the
unmerged tournament quotient this condition is vacuous, since every tournament
class has odd fixed-path tiling fiber.

## Computed Context

`endpoint_transfer_witnesses_s95.py` found:

```text
unmerged tournament private rows: [1, 2, 4, 12, 56]
unmerged tournament row counts:   [1, 2, 4, 12, 56]

merged private rows:              [1, 2, 3, 9, 28]
merged row counts:                [1, 2, 3, 10, 34]
merged ranks:                     [1, 2, 3, 10, 34]

even support matchings:           [1, 1, 2, 7, 8]
even ranks:                       [1, 1, 2, 6, 8]
even row counts:                  [1, 2, 3, 7, 16]
```

So the unmerged tournament quotient currently has the strongest visible proof
mechanism: every tested parent has a private odd child. The merged tournament
quotient needs a subtler rank proof by `n=5->6`; private witnesses no longer
cover all rows, but full row rank persists. The even-graph quotient shows both
zero-row losses and a genuine matching-versus-rank cancellation at `n=5->6`.
