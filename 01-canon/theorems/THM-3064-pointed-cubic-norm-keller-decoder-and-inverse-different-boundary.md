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
  orbit/fixed congruence.  THM-2621 identifies the cofactor as the reciprocal
  inverse-Jacobian numerator.  Its reduced degree makes the single shifted
  resultant equivalent to the full sheetwise Keller congruence; only the
  constant-field gate remains.  An exact split local pair
  changes the normalized inverse-different unit, and hence rho, without
  changing the order or pole data.  This is not a polynomial-map realization
  or a C3 Keller exclusion.
source: codex-jc-resolvent-bridge-2026-08-01
depends_on:
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
  - THM-3042-subdirect-graph-order-common-quotient-and-singleton-owner-criterion
  - THM-3057-tame-quartic-inertia-clutch-index-resonance
related:
  - THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding
  - THM-3059-quartic-twojet-even-jelonek-c3-escape-counterexample
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
script: 04-computation/pointed_cubic_norm_keller_decoder_thm3064.py
output: 05-knowledge/results/pointed_cubic_norm_keller_decoder_thm3064.out
script_sha256: fb82a95d8ccd17e4d509157353590e895e55a0fb90c14ba8d9f3b5fb82dcd696
output_sha256: 9b51211a8c6f24a29f79d74db4311242719e15b1ad84789411eb80fcf6686458
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

### 2.1 Exact specialization to the inverse-spectral pair

[THM-2621](THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger.md)
supplies the cofactor rather than leaving it abstract.  For its inverse pair

```text
f(T;u,v),                         y=b(x;u,v),
q(T)=f_v b_u-f_u b_v mod f,                              (7a)
```

coefficientwise differentiation at fixed `T` gives

```text
det partial(x,y)/partial(u,v)=q(x)/f_T(x),
J_phys(x)=f_T(x)/q(x).                                  (7b)
```

In a split completion `f=(T-a)g`, put

```text
q_0=q(a),                         q_C=q(y),
q_0!=0,                           Res(g,q)!=0.           (7c)
```

Thus `(3)` has `c_0=q_0^(-1)` and `c=q_C^(-1)`, and

```text
rho_q = q_0 (y-a)g'(y)/(d q(y)),                        (7d)

Delta_(f,b)=
 Res_T(g, q_0(T-a)g'(T)-d q(T))
 ------------------------------------------ .           (7e)
 d^3 Res_T(g,q(T))
```

The formula is independent of the representative of `q mod f`: both its
fixed evaluation and its cubic residue are unchanged.  Choose the unique
reduced representative with `deg q<=3`.  Then the shifted resultant is
stronger than a branchwise test:

```text
Delta_(f,b)=0
 iff d q(T)=q_0 f_T(T) in K[T]
 iff q == kappa^(-1)f_T mod f,             kappa=d/q_0. (7f)
```

Here is the clean proof, in the slightly more general `1+p` scope.  Let `g`
be irreducible of degree `p`, assume `deg q<=p` and
`q(a)d Res(g,q)!=0`, and put

```text
r(T)=q(a)(T-a)g'(T)-d q(T).                              (7g)
```

The polynomial `r` has degree at most `p`.  If its resultant with `g`
vanishes, irreducibility gives `r=lambda g`.  Evaluation at `T=a` gives
`-d q(a)=lambda d`, hence `lambda=-q(a)`.  Rearranging yields

```text
d q=q(a)[g+(T-a)g']=q(a)f_T,                            (7h)
```

which is `(7f)`; the converse is immediate.  Therefore a nonzero `(7e)` is
an exact obstruction to the THM-2621 congruence for a supplied inverse pair.
For a polynomial Keller map one retains one final scalar gate:

```text
kappa=d/q_0 must lie in the constant field C*.          (7i)
```

For a specified `kappa_0 in C*`, the exact two conditions are
`Delta_(f,b)=0` and `d=kappa_0 q_0`.  Even these algebraic conditions prove
the rational inverse-spectral congruence, not polynomial realization of that
pair on `A^2`.

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

## 5. Exact split-order supplied-pair hostile

The smallest rational inverse-spectral model makes every boundary visible.
Use target coordinates `(u,t)`, put `R=Q(u)[[t]]`, and take

```text
g(T)=T^3-t,                    f(T)=(T-u)(T^3-t),
S=Q(u)[[s]],                   s^3=t,
d=g(u)=u^3-t in R*.                                    (20)
```

The cubic is transitive by Eisenstein at `t`, its local inertia is tame
`C3`, and

```text
R[T]/(f)=R x S,                 e=g(T)/d=(1,0).         (21)
```

For any nonzero base scalar `r`, define the companion element

```text
b_r(T)=t e-(3u/r)T^2(1-e).                              (22)
```

It has the two exact branch descriptions

```text
fixed:       x=u,      y=b_r(u)=t;
cubic:       x=s,      y=b_r(s)=-(3u/r)s^2.             (23)
```

Their inverse Jacobians with respect to `(u,t)` are respectively `1` and
`1/r`.  Equivalently, direct coefficientwise calculation gives

```text
q_r=f_t(b_r)_u-f_u(b_r)_t
    ==g(T)+(T-u)g'(T)/r                    modulo f,     (24)

q_r(u)=d,                    q_r(y)=(y-u)g'(y)/r.
```

Consequently

```text
J_0=1,                        J_C=r,
rho_q=r,                      Delta_(f,b_r)=(r-1)^3.    (25)
```

At `r=1`, `(24)` is exactly `q_1=f_T mod f`, so the rational pair satisfies
the THM-2621 congruence with `kappa=1`.  At `r=2`, the fixed scalar gate still
holds but the pointed norm defect is one.  The two packets have the same
quartic, split graph order, fixed idempotent, derivative valuations, and
inverse-different cofactor valuation `-2`; only the normalized unit differs.

The resultant normalization in `(7e)` is literal:

```text
Res(g,q_r)=-27t^2 d/r^3,
Res(g,d(T-u)g'-d q_r)=-27t^2 d^4(r-1)^3/r^3.            (26)
```

After specializing the residue of `u` to one in `F_7`, the normalized cubic
derivative residue is `-3=4`; the normalized cofactor residues at `r=1,2`
are respectively `2,4`, giving Jacobian residues `1,2`.  The missing datum is
exactly the normalized inverse-different unit.

This is an exact rational inverse-spectral pair on a **split local cover**.
It is not the irreducible generic inverse quartic of a polynomial map, and it
does **not** assert a polynomial realization of `r=1` or `r=2` with the same
global target data.  That non-realization guard is essential: physical
global realization is precisely where the unresolved Jacobian geometry
lives.

## 6. Keller transport and the affine boundary

For an actual polynomial Keller map, once the cofactor in `(3)` is the true
chain-rule cofactor, the source identity

```text
det JF=kappa in C* subset K*                            (27)
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
                   the THM-2621 supplied-pair degree-p decoder;
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
cubic norm/resultant quotient, the supplied inverse-spectral companion and
its coefficientwise `q=f_t b_u-f_u b_t`, tame different and
inverse-different valuations, residue-7 compensation, the `r=1,2` hostile,
and the separate trace/norm-one failures.  Every truth-bearing check uses
explicit runtime exceptions rather than Python assertions.
