# Endpoint Transfer Smith S95

The GF(2) transfer rank was only the first shadow. Smith normal form shows the
same split integrally.

For tournaments through `5->6`, all nonzero Smith factors are odd:

```text
tournament 5->6:
[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 15]

merged 5->6:
[1, 1, 1, 1, 1, 1, 1, 1, 3, 15]
```

For even graphs, 2-primary factors appear immediately:

```text
even 3->4: [1, 2]
even 4->5: [1, 1, 2]
even 5->6: [1, 1, 1, 1, 1, 1, 8]
even 6->7: [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4, 8]
```

This is a better explanation than "the mod-2 rank drops." The even-graph
quotient has real 2-adic torsion in the endpoint-transfer cokernel. The
tournament quotient, at least in tested range, does not.

## Mathematical Picture

Endpoint transfer has three layers:

```text
integer matrix A
Smith factors of A
mod-2 row rank of A
```

The number of odd Smith factors equals the GF(2) rank. The even graph rank
defect is exactly the count of even Smith factors.

This suggests a refined conjecture:

```text
tournament endpoint transfer has only odd Smith factors.
```

The reason should be the same reason every tournament fiber is odd: the
tournament quotient is protected from 2-torsion by Redei oddness and odd
automorphism groups. The even-graph quotient is not protected because graph
automorphism groups have Sylow-2 structure.

## New Proof Target

The private-child witness conjecture may prove full GF(2) rank, but it will not
by itself prove odd Smith factors beyond mod 2. To prove the Smith statement,
we need either:

1. an integer triangular minor with odd determinant;
2. a signed incidence interpretation where all row-rank minors have odd gcd;
3. a 2-local chain homotopy showing endpoint transfer is split-injective over
   `Z_(2)`.

The third option is the most structural: tournament endpoint transfer might be
an injection in a 2-local quotient category, while even graphs introduce
2-torsion in the quotient map.
