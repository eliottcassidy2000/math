---
id: THM-3042
title: "Subdirect graph-order common quotient and singleton-owner criterion"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  A subdirect graph order B inside a product S x D is exactly a fibre
  product over one common quotient Q=S/I=D/J.  Its singleton slice and the
  singleton slice of its conductor are both Ie, so the singleton idempotent
  belongs to B exactly when Q vanishes.  For a same-complement enlargement
  of the split monogenic quartic order in THM-3038 by generators
  (r_j,h_j), the exact ideal is I=(g(a),r_j-h_j(a)); thus extra graph
  generators may erase, reduce, or preserve the monogenic obstruction.
  This is a finite-order membership criterion, not an affine-owner or
  Keller exclusion.
source: codex-subdirect-graph-order-2026-08-01
depends_on:
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
related:
  - THM-2570-jelonek-cusp-cylinder-normalization-and-conductor
  - THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary
  - THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding
  - THM-3037-cusp-braid-s4-lift-dichotomy-and-common-sheet-owner-boundary
script: 04-computation/subdirect_graph_order_common_quotient_thm3042.py
output: 05-knowledge/results/subdirect_graph_order_common_quotient_thm3042.out
script_sha256: 72f69322a4bdd609cbd23b7b41395f8b6068d5818ceb84e8444723d6975abe14
output_sha256: cc13eedb294977d4ef33182c9ab4843db454e15e2fec4deff98aa0127213738b
hash_basis: LF-normalized bytes
---

# THM-3042 -- a full split graph order is governed by one common quotient

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## 1. Subdirect-product theorem

Let `S,D` be commutative unital rings and let

```text
B subset S x D                                             (1)
```

be a unital subring whose two coordinate projections are surjective.  Put

```text
I={s in S:(s,0) in B},       J={x in D:(0,x) in B}.       (2)
```

Then `I` and `J` are ideals and there is a unique ring isomorphism

```text
theta:S/I -> D/J                                           (3)
```

such that

```text
B=S x_theta D={(s,x):theta(s mod I)=x mod J}.             (4)
```

Thus a subdirect graph order is governed by one common quotient

```text
Q=S/I = D/J.                                               (5)
```

For the singleton idempotent `e=(1,0)` one has the exact slice identities

```text
B intersect Se = Ie,
f_B intersect Se = Ie,
f_B={m in S x D:m(S x D) subset B}.                       (6)
```

Consequently the following conditions are equivalent:

```text
e in B;
I=S;
J=D;
Q=0;
B=S x D.                                                   (7)
```

This is the ring-theoretic Goursat lemma with its load-bearing owner
coordinate exposed.  It applies to arbitrary full graph orders after taking
`D` to be their second projection; it does not assume that the complementary
order remains monogenic.

## 2. Proof

The two kernels in `(2)` are ideals: for example, if `s in I` and `u in S`,
choose `(u,x) in B` by surjectivity and multiply it by `(s,0)`.  The argument
for `J` is identical.  Surjectivity of the first projection then makes the rule

```text
s mod I |-> x mod J whenever (s,x) in B                    (8)
```

defined on every class of `S/I`.  If two elements of `B` have the same first
coordinate, their difference has the form `(0,x-x')`, so the right side of
`(8)` is well defined modulo `J`.  The same argument with the coordinates
reversed gives its inverse.  Addition, multiplication, and units are inherited
from the product, proving `(3)`.  Equation `(4)` follows in both directions
from the definitions of `I,J`.

The first equality in `(6)` is the definition of `I`.  If `s in I`, then
for every `(u,x) in S x D`,

```text
(s,0)(u,x)=(su,0) in B                                    (9)
```

because `I` is an ideal.  Hence `Ie` lies in the conductor.  Conversely the
conductor lies in `B` after multiplication by `1`, so its singleton slice is
contained in `B intersect Se=Ie`.  This proves `(6)`.  Finally `e in B` iff
`1 in I`; under the isomorphism `(3)` that is equivalent to both quotient
rings being zero, which proves `(7)`.

## 3. Exact enlargement formula for the split monogenic order

Adopt the hypotheses and notation of [THM-3038](THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary.md):

```text
g(Y) in R[Y] monic cubic,     a in R,     d=g(a),
C=R[Y]/(g),
A0=R[Y]/((Y-a)g)=R x_(R/(d)) C.                           (10)
```

Consider a same-complement graph-order enlargement

```text
B=A0[(r_1,hbar_1),...,(r_m,hbar_m)] subset R x C,         (11)
```

where `h_j in R[Y]` represents `hbar_j`.  Define the discrepancy ideal

```text
I_B=(d, delta_1,...,delta_m),
delta_j=r_j-h_j(a).                                       (12)
```

Changing `h_j` by a multiple of `g` changes `delta_j` by a multiple of `d`,
so `(12)` is independent of every representative choice.  Then

```text
B=R x_(R/I_B) C,
B intersect Re=I_B e,
e in B iff I_B=R.                                         (13)
```

Indeed `(d,0)` already belongs to `A0`, while

```text
(r_j,hbar_j)-(h_j(a),hbar_j)=(delta_j,0) in B.            (14)
```

Thus `I_B` is contained in the actual first kernel ideal.  Conversely every
generator in `(11)` satisfies the evaluation congruence modulo `I_B`, so the
whole generated algebra lies in the right side of `(13)`; its first kernel
is therefore contained in `I_B`.  Equality and `(13)` follow.

For completeness, equality of the two algebras in `(13)` is not being
inferred from the kernel alone.  If `(r,c)` lies in the displayed fibre
product, choose a representative `h` of `c`.  Then

```text
(r,c)=(h(a),c)+(r-h(a),0),                               (14a)
```

the first summand lies in `A0`, and compatibility says the second belongs to
the now-identified kernel `I_B e`.  Hence the reverse inclusion is literal.

For a local ring with `d` nonunit, `(13)` says that the singleton splits
exactly when at least one added discrepancy is a unit.  Over a DVR, with
`v(0)=infinity`, the remaining gluing length is

```text
length(R/I_B)=min(v(d),v(delta_1),...,v(delta_m)).         (15)
```

This is the promised full-order repair of the one-way boundary in THM-3038:
`d` a unit is sufficient, but a nonunit `d` is not an obstruction once the
full graph order supplies a unit discrepancy.

If additional generators enlarge the second projection beyond `C`, formula
`(12)` must not be used.  One instead sets `D=pr_2(B)` and computes the exact
kernel ideal `I` in `(2)` by module membership, Smith form over a DVR, or a
local Groebner/module calculation.  The equivalences `(6)--(7)` remain valid.

## 4. Sharp controls

The rank-two congruence order

```text
A_d={(r,c) in Z x Z:r=c mod d}                            (16)
```

is the split model with common quotient `Z/(d)`.  Adjoining one element
`b=(r,x)` gives

```text
A_d[b]=A_g,        g=gcd(d,r-x).                          (17)
```

Hence, with `d=4`:

- `b=(1,1)` has zero discrepancy and preserves the quotient `Z/4`;
- `b=(0,2)` shrinks it to `Z/2` but does not split the singleton;
- `b=(0,1)` has unit discrepancy and produces the full product.

The middle case is the sharp hostile to the false dichotomy "extra generators
either do nothing or produce an owner."  The last case is the sharp hostile
to using nonunit `g(a)` as a full-order exclusion.

There is an equally important geometric stopping boundary.  Even `Q=0`
constructs only a finite-order idempotent.  By THM-3038 the corresponding
normalization sheet belongs to the affine source only if every original
source coordinate is regular along its trace.  The source-coordinate pole
`2c/theta^2` in THM-2570 survives unchanged.  Therefore

```text
Q=0  does not imply  affine-source owner.                 (18)
```

## 5. Transfer contract and scope

| source | map | preserved | lost / undecided | required sidecar |
|---|---|---|---|---|
| arbitrary split full graph order | subdirect quotient `(5)` | exact finite gluing | affine-source incidence | source-coordinate regularity |
| monogenic order plus same-complement generators | discrepancy ideal `(12)` | complete singleton membership test | generators outside `C` | full projected order `D` |
| pointed THM-2993 star | THM-3038 cross-resultant | monogenic square `d^2` | enlargement discrepancies | graph-order presentation |

No quartic Keller map, forbidden `G1` degree, JC(2), DC(2), SFC(4), or
LRC(14) case is excluded.  The theorem identifies the exact finite-order
calculation still required by the quartic owner programme; it does not supply
that calculation for an unknown witness.

## 6. Exact companion

The companion exhausts the congruence subdirect products of
`(Z/12Z)^2`, checks their ring closure, projection kernels, common quotients,
singleton/conductor slices, and all `864` one-generator enlargement cases.
It also checks selected two-generator controls and the three sharp states in
`(17)`.  Every truth-bearing check uses an explicit runtime exception, so
ordinary and optimized Python execute the same verification.  Their
LF-normalized transcripts agree byte-for-byte with the stored output.
