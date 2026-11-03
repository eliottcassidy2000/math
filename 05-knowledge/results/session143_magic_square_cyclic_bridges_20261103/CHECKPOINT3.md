# Checkpoint 3: sharp integral repair and a scoped abstract embedding

**PROVED scoped / FINITE-EXACT.** This is the final bounded addendum to
Session143, checked on 2026-11-03 UTC. All earlier scientific files are
unchanged. No conjecture closure, novelty priority, or unspecified cross-subject
transport is asserted. The square-of-squares problem remains **OPEN**.

For the regular five-tournament, write `A=P+P^2`, `P^5=I`. Then

    Q[A]=Q[P],             [Z[P]:Z[A]]=11.

Here `Z[A]` means the ordinary unital polynomial subring: it does not
silently include transpose, the all-ones matrix, or entrywise multiplication.
The quotient is an **additive module**, not a quotient ring.

## Why exactly eleven appears

The five columns `1,A,A^2,A^3,A^4` in the cyclic-shift basis have determinant
`-11` and Smith factors `(1,1,1,1,11)`. They span every polynomial in `A`
because `A^5=5A^2+5A+2`. Thus the additive index is exactly11, not an
approximation or a finite-height observation.

For `r` in `Z[P]`,

    r belongs to Z[A] iff r(1)=r(9) mod11.

The cyclotomic factor itself has **no** defect: for
`lambda=zeta5+zeta5^2`,

    zeta5=-lambda^3-2lambda^2-3lambda-2,
    Z[lambda]=Z[zeta5].

The additional congruence joins the trivial and cyclotomic factors: the
old norm5 condition is supplemented by norm11. This is a precise arithmetic
explanation for the denominator11 in the rational inverse, not a generic
meaning of that prime.

## All moduli, and the decisive hostile

Modulo every positive integer `m`, the generated subalgebra has additive
cokernel `Z/gcd(m,11)Z`. Hence ordinary polynomial generation is complete
exactly when11 is invertible modulo `m`. When11 divides `m`, the cokernel
has order11 even at `m=121,1331,...`; it does not grow with11-adic valuation.

Modulo11 the algebra splits and one eigenvalue collision identifies the
coordinates `P=1` and `P=9`. Modulo5 **all** eigenvalues coalesce, but `A`
still generates the full nonsemisimple algebra. Spectral-support counting
is therefore not a coefficient-ring-independent generation criterion.

## Full adjacency data repairs the loss integrally

The restriction above is essential. With `J` the all-ones matrix,

    J=I+A+A^T,
    P=5J-A^3-2A^2-3A-2I,
    P=A-A circ A^2,

where `circ` is entrywise/Hadamard multiplication. These integer matrix
identities hold over every commutative coefficient ring, including
characteristic11, and for each of the four orientations. Thus the complete
labeled graph data still recovers its cyclic generator. Adjoining `J` or
transpose, or permitting the stated entrywise operation, removes the
ordinary-polynomial obstruction.

No FRACTRAN state, transition or integrality predicate is transported by
this appearance of11. No Baby Monster embedding, Hecke map, octonion
multiplication or squarehood condition follows without its own explicit map.

[Proof and exact formulas](cyclic_integral_order/PROOF.md) ·
[Independent audit](independent_integral_order/REVIEW.md) ·
[Verification](VERIFICATION_CHECKPOINT3.md)


## Final existing-input corollary: an abstract embedding selects 5A

The explicitly verified puzzle subgroup `F20=C5:C4` embeds **abstractly**
into the cited Baby Monster subgroup's `(5:4) x HS` index-two base. This
refines the earlier no-embedding-constructed boundary: existence is now
proved from cited subgroup structure, but no Baby Monster generator words,
matrices, action identification or representation intertwiner are constructed.

The outer involution is **not** assumed to commute with the first factor.
Its normal index-two base has a characteristic normal order-five subgroup
`C`. A nonidentity element of `C` centralizes `C x HS`, of order221,760,000,
which exceeds the cited5B centralizer order6,000,000. It therefore lies in
**5A**. The bound `|N_B(C)| <= 4 |C_B(C)|` is attained at1,774,080,000,
identifying the cited normalizer without an ambiguous direct-product claim.

Over the complex numbers, `F20` has irreducible degrees `1,1,1,1,4`.
The cited4371-dimensional Baby Monster representation restricts as

    870 W4 + a0 lambda0 + a1 lambda1 + a2 lambda2 + a3 lambda3,
    a0+a1+a2+a3=891,    a1=a3.

Here `W4` is the deleted five-point affine permutation representation.
The individual linear multiplicities need further order-four data. In
particular, **870 four-dimensional constituents do not imply870 whole
five-dimensional permutation modules**: those additionally require870
trivial constituents. An explicit rational F20-module with the same C5
restriction and no trivial constituent proves the non-inference, without
claiming to be the actual Baby Monster restriction. A10-dimensional
ordinary complex irreducible of this particular F20 subgroup is impossible;
reducible ten-dimensional modules and other groups are not excluded.

[Abstract embedding and restriction proof](f20_embedding/report.md) ·
[Independent character/group audit](independent_f20_embedding/REVIEW.md)
