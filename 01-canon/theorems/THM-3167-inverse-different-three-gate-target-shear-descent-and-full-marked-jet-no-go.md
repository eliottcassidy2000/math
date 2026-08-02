---
id: THM-3167
title: "Inverse-different three-gate separation, target-shear descent, and full marked-pair finite-jet no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a supplied finite separable inverse-spectral pair
  (f,b), q=f_v b_u-f_u b_v and D=f_T, projective incidence [q]=[D] says
  exactly that the physical Jacobian is diagonal over the target field K;
  it does not say that the diagonal scalar lies in the geometric constant
  field k.  The full Keller predicate separates into the exact line gate,
  the constant-field gate, and the global connected polynomial-owner gate.
  The ratio q/D is invariant under polynomial target shears, so on an actual
  constant-response Faber source lift the line and constant-field gates are
  automatic and add no new passport or degree restriction; before a source
  lift, f,b,q are not defined by the abstract response passport.  An honest
  degree-four polynomial map has q=f_T/(4u) and Jacobian 4u, proving the
  constant-field gate is load bearing.  A second exact fixed-plus-cubic
  rational supplied pair realizes the determinant-one cofactor gauge by an
  actual companion b.  For every N it agrees with the diagonal pair in the
  complete marked inverse-spectral and target-differential jets through
  order N-1, while its pointed norm defect has order 3N.  Hence no invariant
  factoring through any fixed finite jet, including any order/tropical or
  tournament refinement built from such a jet, decides diagonal incidence.
  The hostile is a disconnected split/Laurent local cover, not a polynomial
  A2 map; global polynomial ownership remains a genuine possible
  obstruction.  No new Faber nonentry, C3 exclusion, or JC(2) theorem follows.
source: jc-inverse-different-three-gate-repair-2026-08-02
audit: >
  The primary companion verifies target-shear cancellation, the supplied
  companion's coefficientwise q modulo f, both branch scales, the exact
  shifted resultant, determinant-one product, sixteen full marked-pair jet
  hostiles, eighty target-differential jet checks, the polynomial
  constant-field hostile, and a punctured unit-Jacobian positive control.
  An independent branch/regular-representation audit evaluates the
  companion by CRT, computes the two forward Jacobians directly, reconstructs
  q from branch inverse determinants, obtains the shifted norm as a 3-by-3
  multiplication determinant, and checks thirty-two jet depths.  Fresh
  ordinary and optimized transcripts are byte-identical to the stored
  outputs.
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
related:
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
  - THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift
  - THM-3151-resonant-odd-bipole-equality-cell-nonentry-and-degree-floor
script: 04-computation/jc_inverse_different_three_gate_full_marked_jet_hostile_thm3167.py
output: 05-knowledge/results/jc_inverse_different_three_gate_full_marked_jet_hostile_thm3167.out
script_sha256: 20caa09e72aec73fa6d840186a07e093fce6015d6cb47f26686759ee5580d027
output_sha256: 2ad7549ca1293914fc4a2c82ebff460ff94d02a39fbc8cf07e3e8923ad1492a2
independent_script: 04-computation/jc_inverse_different_three_gate_full_marked_jet_hostile_independent_audit_thm3167.py
independent_output: 05-knowledge/results/jc_inverse_different_three_gate_full_marked_jet_hostile_independent_audit_thm3167.out
independent_script_sha256: 3abb4b14cc27d6ef5b2251248e2ed388ffe6c9b252ca8decf053bd6cd8336da0
independent_output_sha256: 0f1b44a8717c74aaa7ef8a24c65b9a8ecca210b09fc399b00e9b725707da3300
hash_basis: LF-normalized bytes
---

# THM-3167 -- inverse-different incidence has three separate gates

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

**Namespace correction.**  This theorem was first promoted as `THM-3165`;
a concurrent live-bank reservation created a duplicate YAML ID.  It is
canonically `THM-3167`; the mathematical statement is unchanged.

The closest proved mechanisms pull in opposite directions.  THM-3064 gives
the exact pointed cubic norm which detects equality of the fixed and cubic
physical Jacobian values, while THM-3066 and THM-3068 show that symmetric
products and coefficient-pole ledgers forget a normalized cofactor unit.
THM-3151 gives a strong degree floor inside a normalized constant-response
Faber chart.  The present theorem identifies the precise map between these
two languages and proves why the inverse-different line alone cannot sharpen
that Faber floor.

The corrected information hierarchy is

```text
exact inverse-different line:   sheetwise diagonality over K;
constant-field sidecar:         the diagonal scalar lies in k*;
polynomial-owner sidecar:       (f,b) comes from one connected polynomial A2 map.
```

These are three different predicates.  None may be silently renamed as
another.

## 1. Supplied inverse-spectral pair and the three gates

Let `k` be a characteristic-zero field, let

```text
K=k(u,v),
A=K[T]/(f),                  xi=T mod f,                 (1)
```

where `f` is monic and separable of degree `d`.  Let `b in K[T]` have degree
less than `d`, put `eta=b(xi)`, and take all coefficient derivatives below at
fixed `T`.  Define

```text
D=f_T(xi) in A*,
q=f_v(xi)b_u(xi)-f_u(xi)b_v(xi) in A*.                 (2)
```

THM-2621's chain-rule calculation is degree independent and gives

```text
det partial(xi,eta)/partial(u,v)=q/D,
J_phys=det partial(u,v)/partial(xi,eta)=D/q.            (3)
```

There are now three gates.

### Gate I: sheetwise diagonality

The derivative-normalized inverse-different line is

```text
K*D subset A*                     for q,
K*D^(-1) subset D^(-1)A          for c=q^(-1).          (4)
```

Equation `(3)` proves exactly

```text
q in K*D
 iff J_phys in K*1_A.                                      (5)
```

Thus `[q]=[D]` says that the physical Jacobian has the same value on every
generic sheet.  It is the smallest `K`-linear locus with this property: all
possible diagonal values `J in K*` give `q=J^(-1)D`, and conversely every
point of that line is diagonal.

### Gate II: the geometric constant field

A planar Keller pair requires more than `(5)`:

```text
J_phys=kappa*1_A,                  kappa in k*.          (6)
```

Consequently the exact full response condition for a supplied pair is

```text
q=kappa^(-1)D for some kappa in k*.                     (7)
```

The set in `(7)` is `k*D`, not the whole `K`-line in `(4)`.  When the
constant field of `K` under the target derivations is exactly `k`, Gate II
may equivalently be written

```text
d_(K/k)(D/q)=0.                                         (8)
```

The differential in `(8)` is evaluated only after Gate I has made `D/q` a
base scalar.  Projectivizing over `K` necessarily destroys this coordinate.

### Gate III: affine and global ownership

An arbitrary element of the line `(4)` is not yet the cofactor of the same
affine graph.  Locally one must retain the actual companion `b` and the exact
formula `(2)`.  Globally one must also prove that `(f,b)` arises from one
connected polynomial source `A2`, with polynomial target coordinates and no
omitted divisor or nonconstant source unit.

The implications are one-way:

```text
connected polynomial Keller map
  => supplied affine pair (f,b)
  => q from (2)
  => Gate I and Gate II.

abstract q on the line
  does not => supplied b;

supplied rational (f,b)
  does not => connected polynomial A2 realization.     (9)
```

This is the repaired typing of the inverse-different program.

## 2. Exact target-shear descent

Let `H in k[U]` and perform the polynomial target shear

```text
u'=u,                         v'=v+H(u).                (10)
```

Write

```text
f^H(T;u',v')=f(T;u',v'-H(u')),
b^H(T;u',v')=b(T;u',v'-H(u')).                         (11)
```

At fixed `T`, the derivative operators are

```text
partial_(v')=partial_v,
partial_(u')=partial_u-H'(u)partial_v.                 (12)
```

Therefore the two cross terms cancel:

```text
q^H
 =f_v(b_u-H'b_v)-(f_u-H'f_v)b_v
 =f_vb_u-f_ub_v
 =q,                                                    (13)
```

after the base-field identification in `(11)`.  Likewise
`(f^H)_T=f_T`.  Hence `q/D`, the line gate, and the scalar in Gate II all
descend through the exact target-shear quotient.

This is the direct connection to THM-2230.  Its full Faber gauge chooses one
representative of the target-shear response class; equation `(13)` proves
that it does not choose a new inverse-different incidence.

## 3. Consequence for the THM-3151 Faber/passport chart

An entrant in THM-3151 is already an actual normalized polynomial
constant-response source chart.  If its generic inverse degree is four and a
primitive inverse coordinate is chosen, `(3)` and the chart response
`Jac(P,Q)=kappa in k*` give

```text
q=kappa^(-1)f_T                                             (14)
```

identically.  Equation `(13)` shows that passing to the unique full Faber
representative does not change `(14)`.  Thus adding `[q]=[f_T]`, the pointed
norm equation, or the constant-field scalar to the hypotheses of THM-3151
cannot remove another pole passport or raise its proved floor `N>=4D`: those
conditions already hold for every source entrant.

Before chart entry, THM-2796's abstract balanced response retains the square
potential, Padé denominator, and pole passport, but it does not supply the
inverse minimal polynomial `f` or the affine companion `b`.  At that layer
`q` in `(2)` is not canonically defined.  The proposed bridge therefore has
an exact dichotomy:

```text
before an affine source lift:     [q]=[f_T] is undefined;
after a constant-response lift:   [q]=[f_T] is automatic.                (15)
```

This is a stopping reason for a passport-only degree argument, not a claim
that inverse-different methods are useless.  A viable new nonentry theorem
must derive an independent polynomial-owner or boundary restriction on the
**actual** `b` or `q` and then show that restriction is incompatible with
`(14)`.

The complete connection ledger is

```text
source:       polynomial source pair (P,Q), primitive xi, companion b;
map:          target-shear quotient followed by (f,b) -> q=f_vb_u-f_ub_v;
target:       the projective inverse-different class [q:D];
preserved:    constant Jacobian response and q/D;
destroyed:    none under target shear, but the abstract passport forgets f,b;
sidecar:      actual affine companion plus connected polynomial owner;
cheap test:   pointed shifted norm, then the constant-field differential. (16)
```

## 4. Polynomial degree-four hostile to the constant-field omission

The distinction between Gates I and II occurs in an honest polynomial map,
not only in an abstract algebra.  On `A2` with source coordinates
`(xi,eta)`, take

```text
u=xi^4,                         v=xi eta.               (17)
```

The map is generically finite and separable of degree four, and

```text
Jac_(xi,eta)(u,v)=4xi^4=4u.                             (18)
```

Its inverse quartic and reduced companion are

```text
f(T)=T^4-u,
b(T)=vT^3/u,                                            (19)
```

because `T b(T)=v` modulo `T^4-u`.  Coefficientwise differentiation gives

```text
q=f_vb_u-f_ub_v=T^3/u=f_T/(4u).                        (20)
```

Thus `[q]=[f_T]` holds exactly, but its diagonal scalar is `4u`, not an
element of `k*`.  This refutes

```text
[q]=[f_T]  =>  constant Jacobian.                       (21)
```

For a positive control with the same inverse quartic, on the punctured
source take `u=xi^4` and `v=eta/(4xi^3)`.  Then

```text
b(T)=4vT^3,                     q=f_T,
Jac=1.                                                   (22)
```

The positive control is Laurent, while the hostile `(17)` is polynomial.
Together they isolate exactly the constant-field coordinate rather than a
quartic-root artifact.

## 5. An actual supplied-companion determinant-one hostile

The finite-jet obstruction can be made inside the full marked pair, not only
inside an abstract cofactor module.  Let `k_0` contain a primitive cube root
of unity, put `k_hat=k_0((tau))`, let `u,t` be algebraically independent, and
work over

```text
K=k_hat(u,t),
L=K[s]/(s^3-t),
A=K x L.                                                (23)
```

The cubic is irreducible by the `t`-valuation.  In polynomial form put

```text
g(T)=T^3-t,
f(T)=(T-u)g(T),
d=g(u)=u^3-t,
e(T)=g(T)/d.                                            (24)
```

The idempotent `e` is one on the fixed factor and zero on the cubic factor.
For any nonzero `r` which is constant for `partial_u,partial_t`, define

```text
b_r(T)=r^3 t e(T)-(3u/r)T^2(1-e(T)).                   (25)
```

Its branch values are

```text
fixed:       X=u,             Y=r^3 t;
cubic:       X=s,             Y=-(3u/r)s^2.             (26)
```

Equivalently, the two forward source maps are

```text
fixed:       u=X,             t=Y/r^3;
cubic:       u=-rY/(3X^2),    t=X^3.                    (27)
```

Direct Jacobian calculation gives the physical packet

```text
J_0=r^(-3),                    J_C=r.                   (28)
```

In particular its four-sheet product is exactly one:

```text
J_0 Norm_(L/K)(J_C)=r^(-3)r^3=1.                       (29)
```

Let `D=f_T mod f`.  Either coefficientwise differentiation of `(25)`, or
the inverse determinants of `(27)` followed by Chinese remainder, gives

```text
q_r=f_t(b_r)_u-f_u(b_r)_t
   ==r^(-1)D+(r^3-r^(-1))g(T)                  modulo f. (30)
```

Indeed `(30)` evaluates to `r^3D_0` on the fixed factor and `r^(-1)D_C` on
the cubic factor, in agreement with `(28)` and `(3)`.

The pointed ratio and THM-3064 shifted norm are

```text
rho=J_C/J_0=r^4,
Delta=Norm_(L/K)(rho-1)=(r^4-1)^3.                     (31)
```

Thus

```text
q_r in K*D  iff  r^4=1.                                (32)
```

The exact resultant form of `(31)` is

```text
 Res_T(g, q_r(u)(T-u)g'(T)-d q_r(T))
 --------------------------------------------------
 d^3 Res_T(g,q_r(T))
   =(r^4-1)^3.                                          (33)
```

This is the determinant-one hostile of the inverse-different line, now with
the true affine formula `(2)` supplied by `(25)`.

## 6. Universal finite-jet no-go and the smallest exact survivor

Fix an arbitrary proposed jet depth `R>=0`.  Put

```text
N=R+1,                         r_N=1+tau^N.             (34)
```

Take `r=1` as the diagonal control.  Because all four scalar functions

```text
r_N^3, r_N^(-1), r_N^(-3), r_N                         (35)
```

are congruent to one modulo `tau^N`, equations `(25)` and `(30)` give

```text
f_N=f_1 exactly,
b_(r_N)==b_1 mod tau^N,
q_(r_N)==D mod tau^N A,
q_(r_N)^(-1)==D^(-1) mod tau^N D^(-1)A.                (36)
```

Since `tau` is constant for the target derivations, every fixed-order
`u,t` derivative of the two companions also agrees modulo `tau^N`.  The root
packet, derivative star, fixed idempotent, cubic field, every root-only order
clock, and the exact total product `(29)` are unchanged.  The branch
valuations of `q_r/D` and `D/q_r` are also unchanged.  Nevertheless

```text
Delta_N=((1+tau^N)^4-1)^3 !=0,
ord_tau(Delta_N)=3N.                                   (37)
```

Therefore the diagonal control and the non-diagonal hostile have identical
reductions of the **complete marked inverse-spectral differential pair**
through depth `R`.  It follows formally that every invariant which factors
through that truncation has the same value on them.  This includes any
valuation, associated-graded, finite leading-unit, order/tropical, or
tournament construction whose vertices and orientations are functions of
that bounded jet.  No choice of such a construction can decide `(32)` at a
uniform finite depth.

The quantifier is sharp.  The compatible full tower

```text
(f,b,q mod tau^n)_(n>=1)                               (38)
```

recovers the exact `tau`-adic object because the coefficient ring is
separated.  After affine ownership has supplied `q`, the least exact local
detector on a `1+3` field is the two-part packet

```text
Delta=Norm_(L/K)(J_C/J_0-1),
d_(K/k)J_0.                                             (39)
```

The first entry detects Gate I; after it vanishes, the second detects Gate
II.  Equivalently one may retain the exact class of `q/D` modulo diagonal
base scalars together with the surviving scalar's constant-field membership.
No bounded prefix of `(38)` replaces `(39)`.

## 7. Global failure boundary

The hostile `(23)--(37)` is an actual supplied rational inverse-spectral
pair, but it is not a polynomial Keller map on connected `A2`.  Its source is

```text
A2  disjoint-union  (G_m x A1),                         (40)
```

and the cubic forward coordinate `u=-rY/(3X^2)` is Laurent.  The cubic source
has the nonconstant unit `X`.  Hence a criterion using connectedness,
polynomial target regularity, or the no-nonconstant-unit owner sidecar may
still exclude it.  This is exactly Gate III, not a defect in the finite-jet
no-go.

The theorem proves

```text
PROVED:
  target-shear covariance of q and the inverse-different line;
  exact separation of sheet diagonality from constant-field response;
  the degree-four polynomial constant-field hostile;
  an actual supplied-companion determinant-one hostile;
  all-fixed-depth failure for complete marked differential jets;
  the tautology/undefined boundary at the Faber/passport interface.

NOT PROVED:
  that the hostile globalizes to a connected polynomial A2 map;
  that every conceivable global owner invariant factors through local jets;
  a new THM-3151 passport exclusion or degree floor;
  exclusion of C3, A4, S4, G1, JC(2), or DC(2).          (41)
```

The noncanonical branch candidate at commit `7862686ee3` motivated the
finite-jet construction, but it is not a proved-canon dependency.  The
present statement repairs its missing distinction between the `K`-diagonal
line and the geometric constant field, and upgrades its abstract cofactor
hostile to the supplied companion `(25)`.

## 8. Exact reproduction

Run

```text
python3 04-computation/jc_inverse_different_three_gate_full_marked_jet_hostile_thm3167.py
python3 -O 04-computation/jc_inverse_different_three_gate_full_marked_jet_hostile_thm3167.py
python3 04-computation/jc_inverse_different_three_gate_full_marked_jet_hostile_independent_audit_thm3167.py
python3 -O 04-computation/jc_inverse_different_three_gate_full_marked_jet_hostile_independent_audit_thm3167.py
```

Both modes must reproduce their stored outputs byte for byte.  The primary
path uses coefficientwise target derivatives, polynomial reduction, and the
shifted resultant.  The independent path uses direct forward Jacobians,
Chinese remainder branch evaluation, and a cubic regular-representation
determinant.  QED.
