---
id: THM-2756
title: "Opposite-edge projectors, parity cancellation, and the integral clutch"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  Edge complementation is a central involution on the six-edge
  permutation module of K4.  Over Q its positive block is the three-matching
  module 1+[22], with kernel V4 and image S3, while its negative block is the
  faithful standard module [31].  Both block determinants equal the quartic
  sign, so the ambient six-edge determinant is their square and is trivial.
  The transposition/double-transposition trace collision is exactly
  (1+1)=(3-1).  Integrally the two eigensublattices have index 8, Smith form
  1^3,2^3, and quotient F2^3 carrying the matching-permutation action.  This
  is an exact finite representation/lattice theorem, not a Keller or LRC
  exclusion.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on: []
related:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2743-c2-c3-off-diagonal-projector-rank-and-s3-s4-compatibility-defect
  - THM-2746-c3-quotient-affine-v4-lifts-and-a4-projector-defect
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
script: 04-computation/s4_opposite_edge_projector_clutch_thm2756.py
output: 05-knowledge/results/s4_opposite_edge_projector_clutch_thm2756.out
script_sha256: 728cc67a83f00ea47e19176e8f87f57a11ce1a51d4e16e8e09dcac6a3410fd99
output_sha256: ecb97d57576212a7eab0ef55f19f3226dc14a50a71dda858623932ffcd665213
hash_basis: LF-normalized bytes
---

# THM-2756 -- opposite-edge projectors and the integral clutch

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2753 proves that the faithful action of `S4` on the six edges of `K4`
has trivial ambient sign, while the quotient action on the three perfect
matchings restores the quartic sign.  The mechanism is not merely that six is
larger than three.  The six-edge module contains **two copies of the same
determinant character**, separated by the canonical opposite-edge involution.
Their determinants multiply and erase one another.

There is also an integral residue of this rational splitting.  The two
eigensublattices do not span the edge lattice: their sum has index `2^3`, and
the missing quotient is exactly a binary coordinate on each of the three
matchings.  Thus the same object carries a literal

```text
binary clutch over a ternary matching base.                         (1)
```

This is the precise `2/3` structure proved below.  It is not a physical
Keller or LRC carrier without an additional realization map.

## 1. The central opposite-edge involution

Let `X={1,2,3,4}`, let `E=binom(X,2)`, and let

```text
M={12|34,13|24,14|23}.                                      (2)
```

Write `L=Z[E]` and `V=Q[E]`.  Define the opposition involution

```text
O(e)=X\e.                                                    (3)
```

For `g in S4`,

```text
g O(e)=g(X\e)=X\g(e)=O g(e),                                (4)
```

so `O` commutes with the sheet action.  The rational projectors

```text
P_+=(I+O)/2,                    P_-=(I-O)/2                  (5)
```

therefore give `S4`-submodules

```text
V=V_+ direct_sum V_-,             dim V_+=dim V_-=3.        (6)
```

For each matching `m={e,e'}`, choose an ordering and put

```text
p_m=e+e',                         n_m=e-e'.                  (7)
```

Then the `p_m` form a basis of `V_+` and the `n_m` form a basis of
`V_-`.  Changing the ordering only changes the sign of one `n_m`; none of
the representation statements depends on that gauge.

## 2. The two three-dimensional blocks

The positive basis in `(7)` is permuted exactly as the three matchings:

```text
V_+ = Q[M] = 1 direct_sum [22],
ker(S4 -> GL(V_+))=V4,
im(S4 -> GL(V_+))=S3.                                      (8)
```

Here the invariant line is spanned by `sum_m p_m`; its zero-sum complement
has character `(2,0,2,-1,0)` on the five `S4` classes and norm one, hence is
the irreducible `[22]` module.

The negative block is the standard sheet module:

```text
V_- isomorphic to Q[X]/Q(1,1,1,1) = [31].                  (9)
```

Indeed its trace at `g` is

```text
tr(g|V_-)=#{x in X:g(x)=x}-1,                              (10)
```

which is the standard character.  In particular `V_-` is faithful.  The
full class table is

| sheet class | size | `tr_+` | `det_+` | `tr_-` | `det_-` | `tr_E` | `det_E` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `1^4` | 1 | 3 | +1 | 3 | +1 | 6 | +1 |
| `2 1^2` | 6 | 1 | -1 | 1 | -1 | 2 | +1 |
| `2^2` | 3 | 3 | +1 | -1 | +1 | 2 | +1 |
| `3 1` | 8 | 0 | +1 | 0 | +1 | 0 | +1 |
| `4` | 6 | 1 | -1 | -1 | -1 | 0 | +1 |               (11)

The exact block diagonalization is checked for all `24` elements, not inferred
only from the class table.

## 3. Determinant cancellation and trace restoration

The determinant of the matching permutation block is its `S3` sign.  By the
same transposition-generator argument as in THM-2753, this is the quartic
sign.  The determinant of the standard block is also the quartic sign: the
four-dimensional sheet permutation module is the trivial line plus `[31]`,
and its determinant is the sheet sign.  Therefore

```text
det(g|V_+)=sgn(g)=det(g|V_-),
det(g|V)=sgn(g)^2=1.                                      (12)
```

Thus the ambient-even result of THM-2753 is a **two-copy cancellation**, not
the absence of the sign character from the edge representation.  Either
three-dimensional block carries the sign; forgetting the opposition grading
multiplies them.

The sharp transposition/double-transposition collision has the same anatomy.
For a transposition `tau` and a nontrivial `delta in V4`,

```text
(tr_+(tau),tr_-(tau))       =(1,1),
(tr_+(delta),tr_-(delta))   =(3,-1),
tr_E(tau)=tr_E(delta)       =2.                            (13)
```

So the matching block does not create new information from nothing: it stops
the `+2` and `-2` blockwise changes in `(13)` from cancelling.  The negative
block also separates the classes, with the opposite signed response.

## 4. The integral clutch

Define the integral eigensublattices

```text
L_+=L intersect V_+,                   L_-=L intersect V_-. (14)
```

On one opposite pair the columns `(1,1)` and `(1,-1)` have determinant `-2`.
Across the three pairs, the six columns `(p_m,n_m)` therefore have determinant
of absolute value `8`.  Equivalently their Smith form is

```text
diag(1,1,1,2,2,2),
[L:L_+ direct_sum L_-]=8.                                  (15)
```

The quotient has a canonical description.  Define

```text
kappa:L -> F2[M],
kappa(x)_m=x_e+x_e' mod 2                 for m={e,e'}.     (16)
```

It is surjective.  Moreover `kappa(x)=0` iff the two coordinates in every
opposite pair have the same parity, iff

```text
x=sum_m ((x_e+x_e')/2)p_m + ((x_e-x_e')/2)n_m
  belongs to L_+ direct_sum L_-.                           (17)
```

Consequently

```text
L/(L_+ direct_sum L_-) isomorphic to F2[M] isomorphic to F2^3. (18)
```

Equation `(4)` makes `(16)` `S4`-equivariant, so the quotient action is the
matching-permutation action, again with kernel `V4` and image `S3`.  The same
criterion can be stated directly in terms of the rational projectors:

```text
kappa(x)=0  iff  P_+x and P_-x are both integral.          (19)
```

This is the integral clutch hidden by the rational direct sum.  Modulo two,
the signs distinguishing `p_m` from `n_m` disappear; what remains is one
binary parity defect over each of the three matching directions.

## 5. What the `2/3` structure does and does not transfer

The result gives a canonical version of the user's modular-group heuristic:

- `3` is the number of nonzero `V4` directions and the size of the resolvent
  matching set `M`;
- `2` is the denominator of the opposition projectors and the torsion of the
  integral failure to split;
- the rational ternary quotient and the integral binary clutch occur on the
  same six-edge object, rather than as unrelated analogies.

For a quartic polynomial, `det(V_+)` is the cubic-resolvent discriminant
character, while `(12)` explains why the ungraded six-edge determinant cannot
see it.  The exact equality of the quartic and integral-resolvent
discriminants is still THM-2598; the present theorem proves only the
representation and lattice mechanism underneath its sign character.

THM-2606 and the octahedral LRC sidecars also use `L(K4)` and its opposite
matching.  A transfer would require a physical six-packet module on which the
opposition involution, its integral lattice, and the relevant current are all
equivariant.  No such physical identification is proved here.  In particular,
`kappa` is not automatically an endpoint current, relation address, tournament
orientation, or LRC owner label.

## 6. Exact verification

Run

```bash
python 04-computation/s4_opposite_edge_projector_clutch_thm2756.py
python -O 04-computation/s4_opposite_edge_projector_clutch_thm2756.py
```

Both executions byte-match the stored `16`-line transcript
`05-knowledge/results/s4_opposite_edge_projector_clutch_thm2756.out`.  The
companion uses explicit exceptions and no truth-bearing Python assertions.  It
enumerates all `24` sheet permutations, checks centrality of opposition and
the exact six-by-six block conjugacy, computes `(11)`, verifies both determinant
characters and the standard character, identifies the `[22]` norm and `V4`
kernel, proves determinant/index `8` and mod-two rank `3`, checks clutch
equivariance on every edge basis vector, and exhausts all `4^6` residue vectors
for the parity-kernel criterion.

## 7. Boundary ledger

```text
PROVED HERE (candidate):  central opposite-edge involution;
                          exact rational 3+3 block diagonalization;
                          E_+=1+[22] matching quotient with kernel V4;
                          E_-=[31] faithful standard module;
                          duplicate sign determinants and their cancellation;
                          exact blockwise trace anatomy of the class collision;
                          integral index 8 / Smith 1^3,2^3 clutch;
                          equivariant quotient F2[M].

NOT PROVED:               graph-quartic or Keller realization;
                          a new discriminant polynomial identity;
                          exclusion of A4 or S4 Keller monodromy;
                          equivalence with a six-vertex tournament;
                          physical LRC octahedral module/current;
                          endpoint, owner, or relation-address transport;
                          JC(2), DC(2), or LRC(14).                       (20)
```

QED (candidate).
