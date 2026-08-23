---
id: THM-3846
title: "Every unimodular Russell arm jet has a formal Darboux lift, with an exact algebraization gate"
status: >
  PROOF CANDIDATE + VERIFIED-EXACT, PENDING INDEPENDENT HOSTILE AUDIT.  The
  completed triple arm of the THM-3785 Russell pseudo-plane has an explicit
  canonical coordinate s=e/[3(c^3+er)] with {z,s}=1.  Every polynomial arm
  immersion together with a unimodular normal Bezout row extends by a closed
  Catalan normal-coordinate series to an exact formal Darboux pair.  The
  completion is also the rational chart Frac(B)=k(s,z), so every nonzero
  normal Wronskian makes the displayed quadratic discriminant a simple
  linear prime and the lift is not even rational.  Vanishing Wronskian
  forces the arm map to be injective.  Hence every self-identifying arm
  required by THM-3843 has a formal lift but its canonical linear-normal
  resummation cannot algebraize.  Arbitrary higher-normal formal lifts and
  global polynomial Darboux pairs remain OPEN.
source: root / completed Russell pseudo-plane arm and algebraization lane, 2026-08-23
audit: >
  SELF-HOSTILE EXACT CANDIDATE.  The proof uses the relative formal implicit
  function theorem, checks both directions of the explicit (e,r)<->(s,z)
  coordinate change, derives the canonical bracket without dividing by z,
  and constructs the universal lift from a quadratic formal ODE.  The
  deterministic companion verifies the surface and Poisson identities,
  inverse coordinate formulas, generic Bezout/Jacobian calculation, Catalan
  truncations through order ten, the nodal Wronskian, and the Laurent-prime
  square obstruction.  Normal and optimized runs byte-match the frozen
  transcript.  Independent hostile audit remains.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3790-cubic-pseudoplane-arm-nodal-immersion-gate
related:
  - THM-3792-pure-first-normal-nodal-carriers-have-critical-points
  - THM-3812-nodal-arm-coefficient-second-normal-profile-nonentry
  - THM-3839-constant-tower-bichromatic-r2z-profile-nonentry
  - THM-3843-russell-arm-birational-immersion-and-forced-self-identification
script: 04-computation/jc2_russell_formal_arm_darboux_algebraization_thm3846.py
output: 05-knowledge/results/jc2_russell_formal_arm_darboux_algebraization_thm3846.out
script_sha256: b0d2a0c52461b413bb6763563306c53ea3569ced7dc6c20a2dc1fb6fa3aa3300
output_sha256: 5445e1a7a54f097cb03d06736910d0240bb48c3329b3db9ce21558fb1b1afca9
semantic_sha256: cb744130f854c584aa536d307aeba139b5cb2df9be96ea0e206781210b758e83
hash_basis: raw LF bytes
---

# THM-3846 -- the arm has no formal obstruction, only algebraization debt

**PROOF CANDIDATE + VERIFIED-EXACT, PENDING INDEPENDENT HOSTILE AUDIT.**
Let `k` be a field of characteristic zero, fix `c in k*`, and put

```text
B=k[r,z,e]/(r^2e-z^3+c^3r),                 I=(r,z).             (1)
```

Give `B` the THM-3785 Poisson bracket

```text
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3c^3+6re.              (2)
```

In the `I`-adic completion `Bhat`, define

```text
s=e/[3(c^3+er)].                                                (3)
```

The denominator is a unit along `I`.  Then

```text
Bhat = k[s][[z]],                         {z,s}=1,                (4)

e=3c^3s+9s^2z^3,
r=z^3/(c^3+3sz^3).                                             (5)
```

Thus the completed multiple-fibre arm is an ordinary formal symplectic
strip.  More strongly, every admissible first arm jet has an exact formal
Darboux lift, as follows.

Take `a,b,alpha,beta in k[s]` satisfying

```text
alpha b'-a' beta=1,                                             (6)
```

and put

```text
W=alpha beta'-alpha' beta.                                     (7)
```

There is a unique `Z in z k[s][[z]]` solving

```text
Z+(W/2)Z^2=z,                                                   (8)
```

namely

```text
Z=sum_(n>=1) Catalan_(n-1)(-W/2)^(n-1) z^n.                    (9)
```

Then

```text
Ahat=a(s)+alpha(s)Z,
Chat=b(s)+beta(s)Z                                             (10)
```

satisfy exactly

```text
{Ahat,Chat}=1,                                                  (11)
```

with arm restriction `(a,b)` and first normal row `(alpha,beta)`.

The construction has an exact rationality gate.  If `W!=0`, then

```text
Z=(-1+sqrt(1+2Wz))/W                                           (12)
```

in the formal completion, and therefore

```text
Z in Frac(B)  iff  1+2Wz is a square in Frac(B).                (13)
```

In fact `(5)` gives `Frac(B)=k(s,z)`.  A nonzero `W in k[s]` makes the
right side of `(13)` a degree-one polynomial in the transcendental `z`, so
it has a simple prime valuation and cannot be a square.  Therefore

```text
(Ahat,Chat) in Frac(B)^2  iff  Z in Frac(B)  iff  W=0.          (13a)
```

But `W=0` forces the arm map `(a,b)` to be injective.  Consequently every
self-identifying arm packet forced by THM-3843 has an exact formal lift of
the form `(10)`, but that lift is necessarily nonrational.

This is only a gate for the displayed lift: failure of `(13)` does not
exclude a different higher-normal formal lift or a global Darboux pair.

## 1. Exact completed-arm coordinates

The derivative with respect to `r` of the surface equation in `(1)` is
`2re+c^3`, which restricts to the unit `c^3` on `I`.  The relative formal
implicit-function theorem therefore gives

```text
Bhat ~= k[e][[z]],                                              (14)
```

with `r` the unique series beginning `z^3/c^3` that solves `(1)`.  In
particular the completed `I`-adic and `z`-adic topologies agree.

Put

```text
H=c^3+2er.
```

The surface relation gives, without a choice of square root,

```text
H^2=c^6+4ez^3,                    H=c^3 mod z.                   (15)
```

Since `H-c^3=2er` and `z^3=r(c^3+er)`, equation `(3)` is equivalently

```text
s=(H-c^3)/(6z^3).                                              (16)
```

The right side denotes its uniquely cancelled formal series; it makes no
division by `z` in `Bhat`.  Differentiating `(15)` at fixed `z` gives

```text
partial_e s=1/(3H).                                            (17)
```

Equation `(2)` now yields `{z,s}=3H partial_e(s)=1`.

Conversely, set

```text
e=3s(c^3+3sz^3),                 r=z^3/(c^3+3sz^3).             (18)
```

Direct substitution proves the surface equation and recovers `(3)`.  Since
`e=3c^3s mod z`, the change `e<->s` is a formal automorphism.  This proves
`(4)--(5)` with both directions controlled.

The coordinate `s` is deliberately local.  On the global surface the other
component of `z=0` is

```text
L'=V(z,c^3+er) ~= G_m,                                         (19)
```

and `e` is a unit at its generic point.  Hence `(3)` has a pole on `L'`.
The formal Darboux coordinate already carries a global pole-transfer debt.

## 2. Universal formal lift of a Bezout arm packet

First regard

```text
A_0=a(s)+alpha(s)Z,             C_0=b(s)+beta(s)Z               (20)
```

as functions of independent coordinates `(Z,s)`.  Their ordinary Jacobian
is

```text
J_(Z,s)(A_0,C_0)
 =alpha(b'+beta'Z)-(a'+alpha'Z)beta
 =1+WZ.                                                         (21)
```

Implicit differentiation of `(8)` gives

```text
partial_z Z=1/(1+WZ).                                          (22)
```

Because `(Z,s)` is obtained from `(z,s)` without changing `s`, equations
`(21)--(22)` multiply to one.  In the canonical bracket `(4)` this is
exactly `(11)`.  The initial condition `Z=z mod z^2` preserves the requested
arm value and first normal coefficient.

Existence and uniqueness in `(8)` follow either recursively or from the
formal implicit-function theorem at `(Z,z)=(0,0)`.  Lagrange inversion gives
the Catalan expansion `(9)`, which has coefficients in `k[s]` and uses only
characteristic zero.  Consequently every such arm packet has compatible
Darboux jets to arbitrary order: truncating sufficiently far gives elements
of `B` whose bracket is `1` to any prescribed `I`-adic order.  Precisely,
the canonical formula shows `{I^m,B} subset I^(m-1)`; approximate the formal
pair one order beyond the desired bracket precision.  The sparse
no-go theorems THM-3792, THM-3812, and THM-3839 are therefore
**algebraization/support obstructions**, not finite formal-neighborhood
obstructions.

Finally, if `W!=0`, the quadratic formula in the completed fraction field
gives `(12)`.  Since `W` lies in `k(s) subset Frac(B)`, one implication in
`(13)` is immediate.  Conversely, a square root in `Frac(B)` is integral
over the normal ring `k[s][[z]]`, hence belongs to it; changing its sign
chooses the root congruent to one modulo `z` and recovers `(12)`.  Moreover
`(6)` makes at least one of `alpha,beta` nonzero, so `(10)` shows that the
displayed pair lies in `Frac(B)` if and only if `Z` does.  When `W=0`,
equation `(8)` simply gives `Z=z`; the separate pole of `s` in `(19)` may
still prevent algebraization.  Equations `(3),(5)` also give both inclusions
in the global field identity

```text
Frac(B)=k(s,z).                                                 (22a)
```

For `W!=0`, the prime `z+1/(2W)` of `k(s)[z]` occurs to order one in
`1+2Wz`; this proves the nonsquare assertion and `(13a)`.

It remains to justify the injectivity claim when `W=0`.  If `alpha=0`,
`(6)` makes `a'` a nonzero constant, so `a` itself recovers `s`.  Otherwise

```text
0=W=alpha^2(beta/alpha)'                                      (22b)
```

in `k(s)`, hence `beta=lambda alpha` for some `lambda in k`.  Equation
`(6)` becomes

```text
alpha(b'-lambda a')=1.                                        (22c)
```

Both factors are polynomial units.  Thus

```text
b-lambda a=mu s+nu,                    mu in k*,                (22d)
```

which again recovers `s` from the target coordinates.  Therefore `W=0`
is incompatible with the noninjective arm demanded by THM-3843.

## 3. The minimal nodal packet crosses the square gate

Translate THM-3790's minimal nodal arm data to the canonical parameter
`e=3c^3s` on `z=0`:

```text
a=9c^6s^2,                    b=27c^9s^3-3c^3s,
alpha=-1/(3c^3),              beta=-3s/2.                       (23)
```

Then

```text
alpha b'-a'beta=1,             W=1/(2c^3),                      (24)
```

and `(8)` becomes

```text
Z+Z^2/(4c^3)=z,
Z=2c^3(sqrt(1+z/c^3)-1).                                      (25)
```

The square root in `(25)` is a genuine formal Hensel root, not a rational
function on `Y`.  Indeed `z+c^3` cuts out the irreducible Laurent divisor

```text
B/(z+c^3)
 ~= k[r,e]/(r^2e+c^3r+c^9)
 ~= k[r,r^-1].                                                  (26)
```

It occurs with multiplicity one, so

```text
ord_(z+c^3)(1+z/c^3)=1.                                       (27)
```

A square in `Frac(B)` has even valuation at every prime divisor.  Thus
`1+z/c^3` is not a square and `(13)` shows `Z notin Frac(B)`.  Since

```text
Ahat=9c^6s^2-Z/(3c^3),                                        (28)
```

while `s in Frac(B)`, even the first coordinate of this canonical formal
nodal lift lies outside `Frac(B)`.

This identifies the missing operation precisely.  Every allowable
self-identifying boundary and Bezout normal row is formally unobstructed, but
its canonical linear-normal symplectic resummation introduces a nonrational
quadratic Hensel sheet; the minimal node additionally displays the
opposite-arm pole `s` concretely.  A genuine counterexample on `Y` must use a
genuinely higher-normal seed whose resummation cancels these debts inside
`B`; repeating finite sparse rows without tracking the resummed normal
coordinate cannot decide that problem.  No such algebraization is
constructed here, and `JC(2)` remains open.  **QED, pending independent
hostile audit.**
