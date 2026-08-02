---
id: THM-3064
title: "Pointed cubic-norm Keller decoder and inverse-different boundary"
status: >
  PROVED + VERIFIED-EXACT (independent audit pending).  After a fixed sheet
  of a split quartic f=(T-a)g has been marked and the branchwise
  primitive-element cofactor has been retained, the orbit/fixed physical
  Jacobian ratio is an element rho of the complementary cubic field.  The
  single descended scalar Norm(rho-1), equivalently an explicit cubic
  resultant quotient, vanishes exactly when the three conjugate orbit
  values equal the fixed value.  On a tame C3 divisor the complementary
  derivative generates the cubic different, so Keller constancy forces the
  cofactor to generate the inverse different, with exact ramified valuation
  -2.  But the fixed sheet is already separated by a unit cross-resultant:
  the THM-3038/3042 graph order is a direct product and supplies no
  orbit/fixed congruence.  An exact local family changes the normalized
  inverse-different unit, and hence rho, without changing the order or pole
  data.  This is an abstract local-order hostile, not a polynomial-map
  realization or a C3 Keller exclusion.
source: codex-jc-resolvent-bridge-2026-08-01
depends_on:
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
  - THM-3042-subdirect-graph-order-common-quotient-and-singleton-owner-criterion
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
related:
  - THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
script: 04-computation/pointed_cubic_norm_keller_decoder_thm3064.py
output: 05-knowledge/results/pointed_cubic_norm_keller_decoder_thm3064.out
script_sha256: 06d6332425bd3b4c4303fb301e6d788975fc95bf17b04a0ae409ec03dc642eec
output_sha256: 0b96d45307d6d9c5100e4e3bcbc115ac6d63d700a2ab881f785b9f7a5c101413
hash_basis: LF-normalized bytes
---

# THM-3064 -- the pointed cubic norm detects Keller equality

**PROVED + VERIFIED-EXACT (independent audit pending).**

## 1. Inheritance and the pointed ratio

[THM-3066](THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor.md)
proves that inserting four primitive-element cofactors through the natural
`K4` vertex gauge retains only their product.  The least-used surviving
sidecar is not another matching contraction: it is the actual fixed sheet
from [THM-3038](THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary.md)
together with the complementary cubic element carried by the cofactor.

Let `R` be an integral domain with fraction field `K`.  Take

```text
g(T) in R[T] monic cubic and separable irreducible over K,
a in R,                         d=g(a)!=0,
L=K[Y]/(g),                     y=Y mod g,
f(T)=(T-a)g(T).                                         (1)
```

Let `c_0 in K*` be the primitive-element cofactor on the fixed sheet and
`c in L*` its value on the transitive cubic factor.  The derivative stars are

```text
D_0=f'(a)=g(a)=d,
D_C=f'(y)=(y-a)g'(y).                                  (2)
```

Thus the physical Jacobian packet and its pointed orbit/fixed ratio are

```text
J_0=c_0d,                    J_C=c(y-a)g'(y),
rho=J_C/J_0 in L*.                                      (3)
```

Once the fixed sheet and the cofactor have been retained, `(3)` is the
cheapest unsymmetrized answer: `rho` is one element of the pointed
complementary cubic algebra.  The canonical one-scalar descent of the
**Keller equality test** is

```text
Delta_J=Norm_(L/K)(rho-1) in K.                         (4)
```

Because `L` is a field,

```text
Delta_J=0    iff    rho=1    iff    J_C=J_0 in L.       (5)
```

The last equality means that all three conjugate orbit values equal the
fixed value.  Irreducibility/transitivity is load bearing.  For a reduced
cubic product, a zero norm says that **at least one** component vanishes; it
does not test equality on all components.  In that scope one must retain the
element `rho-1` itself, or all coefficients of its regular representation.

## 2. Explicit resultant decoder

Write the cubic cofactor in the fraction algebra as

```text
c=h(y)/b(y),                 h,b in K[T],
b(y)!=0,                    J_0!=0.                     (6)
```

Since `g` is monic, field norms are resultants.  Equations `(3)--(4)` give
the exact rational formula

```text
Delta_J =
 Res_T(g(T), h(T)(T-a)g'(T)-J_0 b(T))
 ---------------------------------------------------- . (7)
 J_0^3 Res_T(g(T),b(T))
```

Every denominator in `(7)` is explicitly nonzero under `(1)` and `(6)`.
This is a fixed-sheet-pointed, `C3`-symmetric scalar.  It is not available
from the unpointed `S4/V4` matching packet because that quotient chooses
neither `a` nor the two cofactor projections.

The norm is applied to the **difference**, not separately to the ratio.
The distinction is sharp.  In `L=K(s)`, `s^3=t`, one has

```text
Tr_(L/K)(1+s)=3=Tr_(L/K)(1),             but 1+s!=1.    (8)
```

If `K` contains a primitive cube root `zeta`, then

```text
rho=(1+s)/(1+zeta s) !=1,
Norm_(L/K)(rho)=1.                                      (9)
```

Thus neither the normalized trace nor the ordinary norm of `rho` is a
Keller decoder.  The shifted norm `(4)` is.

## 3. What the graph-order common quotient says

There is a useful exact reformulation inside
[THM-3042](THM-3042-subdirect-graph-order-common-quotient-and-singleton-owner-criterion.md).
Let an `R`-subdirect graph order be

```text
B=R x_Q D,
I={r:(r,0) in B},       J={z:(0,z) in B},
Q=R/I=D/J.                                              (10)
```

For a unit `J_0 in R*` and a packet `(J_0,J_C) in B`, diagonal `R` lies in
`B`, so subtraction gives

```text
(0,J_C-J_0) in B,
J_C-J_0 in J,
rho-1 in J after division by J_0.                       (11)
```

Congruence modulo `J` is the complete relation supplied by the graph order.
It forces equality for every such packet only at the opposite extreme
`J=0`.  When the singleton idempotent splits, THM-3042 has `Q=0` and `J=D`,
so `(11)` is vacuous.  In other words, finite-order separation of the fixed
sheet removes rather than strengthens cross-sheet congruence.

On the tame `C3` skeleton of
[THM-3057](THM-3057-tame-quartic-inertia-clutch-index-resonance.md), the fixed
root has a distinct reduction from the cubic orbit.  Consequently

```text
d=g(a)=Norm(a-y) is a unit.                              (12)
```

THM-3038 therefore gives, already for the split monogenic order,

```text
R[T]/((T-a)g)=R x C,               C=R[T]/(g),
Q=0.                                                     (13)
```

Any same-complement graph-order enlargement is still the direct product.
Hence the THM-3038/3042 order coordinate supplies the marked factor needed
to *form* `(3)`, but it supplies no equation forcing `(3)` to be one.

## 4. The cofactor is the inverse-different coordinate

Now let `R` be a complete DVR with residue characteristic different from
three.  Suppose `L/K` is a totally tamely ramified cubic field, its integral
closure is the monogenic maximal order

```text
S=R[y],                                                   (14)
```

and the fixed root is separated as in `(12)`.  Normalize the ramified
valuation by `w(varpi_L)=1`, so `w(varpi_R)=3`.  The monogenic different is

```text
Different_(S/R)=(g'(y)).                                (15)
```

Since `y-a` is a unit, `(2)` and `(15)` imply

```text
(D_C)=Different_(S/R),                 w(D_C)=3-1=2.    (16)
```

Meanwhile `D_0=d` is a unit.  If the physical Jacobian is a unit on the
cubic branch, then `(3)` forces

```text
w(c)=-2,                       (c)=Different_(S/R)^(-1). (17)
```

For a Keller packet `J_0=J_C=kappa in R*`, this is exact:

```text
c_0=kappa/d,                   c=kappa/D_C.             (18)
```

Thus a tame `C3` Keller component does not ask for an integral cofactor.  It
asks for the canonical inverse-different pole, followed by one normalized
unit equality.  In particular:

```text
if w(c)>=-1, a unit Keller Jacobian on the C3 orbit is impossible.        (19)
```

Equation `(19)` is a genuine scoped exclusion lemma.  The missing theorem is
to derive such a lower bound from polynomial inverse-coordinate geometry.
Fixed-sheet affine regularity does not provide it: that predicate tests the
source-coordinate valuations at the fixed trace, not the cofactor at the
different cubic trace.

## 5. Exact split-order hostile

The smallest local model makes every boundary visible.  Over
`R=Q[[t]]`, take

```text
g(T)=T^3-t,                    f(T)=(T-1)(T^3-t),
S=R[s]/(s^3-t)=Q[[s]].                                  (20)
```

The cubic is transitive over `Q((t))` by Eisenstein at `t`, its local
inertia is tame `C3`, and

```text
d=g(1)=1-t in R*,              R[T]/(f)=R x S.          (21)
```

The pointed derivatives are

```text
D_0=1-t,                       D_C=3s^2(s-1),
w(D_0)=0,                      w(D_C)=2.                (22)
```

For any base unit `r in R*`, define the abstract cofactor packet

```text
c_0=1/(1-t),                   c_r=r/[3s^2(s-1)].       (23)
```

Then

```text
J_0=1,                         J_C=r,
rho=r,                         Delta_J=(r-1)^3.         (24)
```

Both `r=1` and `r=2` have the identical split graph order, fixed idempotent,
derivative valuations, and inverse-different cofactor valuation `-2`.
Their pointed ratios are respectively one and two.  Over `F_7`, the
normalized cubic derivative residue is `-3=4`; the normalized cofactor
residues are respectively `2` and `4`, giving Jacobian residues `1` and `2`.
The missing datum is exactly the normalized inverse-different unit.

The resultant form of `(24)` is literal.  With

```text
b(T)=3T^2(T-1),                  h(T)=r,                (25)
```

one has

```text
Res(g,b)=27t^2(t-1),
Res(g,(r-1)b)=27t^2(t-1)(r-1)^3.                       (26)
```

This hostile is an abstract local order plus branch-cofactor packet.  It
does **not** assert that either packet, especially `r=2`, is realized by a
polynomial map with the same target data.  That non-realization guard is
essential: physical realizability is precisely where the unresolved
Jacobian geometry lives.

## 6. Keller transport and the affine boundary

For an actual polynomial Keller map, once the cofactor in `(3)` is the true
chain-rule cofactor, the source identity

```text
det JF=kappa in K*                                      (27)
```

maps to `J_0=J_C=kappa` in every completed branch field.  Therefore `(27)`
forces `rho=1` and `Delta_J=0` without any boundary-regularity hypothesis.
This is faithful transport of the Keller identity, not a new exclusion.

THM-3038's affine-owner test answers a different question: whether every
original source coordinate is regular at the marked fixed trace.  Even when
that answer is positive, `(13)` has separated the complementary cubic
factor and imposes no condition on its normalized inverse-different unit.
The surviving proof obligation at a hypothetical `C3` Jelonek divisor is
therefore one of the following equivalent kinds:

1. compute the true branch cofactor and prove `Delta_J!=0`;
2. prove the cofactor pole is milder than `-2`, invoking `(19)`; or
3. prove that the required inverse-different principal part cannot arise
   from the Laurent expansions of polynomial inverse coordinates.

Another unpointed discriminant, matching product, or common-quotient
calculation cannot substitute for this obligation.

```text
PROVED HERE:       the pointed cubic norm/resultant Keller decoder;
                   its exact irreducible-field scope;
                   graph-order congruence (11) and C3 product collapse;
                   derivative=different and cofactor=inverse-different;
                   exact valuation -2 and scoped lower-bound exclusion;
                   exact split-order/unit-residue hostile.

NOT PROVED:        polynomial-map realization of the hostile packets;
                   cofactor integrality for an unknown Keller branch;
                   exclusion of a C3 Jelonek component, A4, S4, G1,
                   JC(2), or DC(2).                                     (28)
```

## 7. Exact companion

Run

```text
python3 04-computation/pointed_cubic_norm_keller_decoder_thm3064.py
python3 -O 04-computation/pointed_cubic_norm_keller_decoder_thm3064.py
```

Both modes must LF-byte-match the stored transcript.  The companion checks
the generic pointed derivative decomposition, exact singleton idempotent,
cubic norm/resultant quotient, tame different and inverse-different
valuations, residue-7 compensation, the `r=1,2` hostile, and the separate
trace/norm-one failures.  Every truth-bearing check uses explicit runtime
exceptions rather than Python assertions.
