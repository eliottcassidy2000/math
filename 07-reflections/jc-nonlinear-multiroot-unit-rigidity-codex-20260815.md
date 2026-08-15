# Nonlinear multiroot unit rigidity: the generic observer has only one-root resonances

**Status:** derivation companion for
[THM-3424](../01-canon/theorems/THM-3424-nonlinear-monomial-fiber-unit-observer-rigidity.md),
which remains **RESERVED / PROVISIONAL + AUDIT REQUIRED**.  It classifies the
distinguished unit observer for every nonconstant `g`, but it does **not**
classify the full multiple-root Hamiltonian module.

## 1. Inheritance and connection contract

[THM-3418](../01-canon/theorems/THM-3418-one-monomial-nonlinear-fiber-keller-classification.md)
proves that a nonconstant top coefficient `g` admits no polynomial mate.
[THM-3419](../01-canon/theorems/THM-3419-generic-kummer-response-regular-sector-rank.md)
identifies the generic Hamiltonian response with `N` copies of the regular
Kummer packet, where `N=deg(rad(g))`.  The proved and independently audited
[THM-3422](../01-canon/theorems/THM-3422-one-root-nonlinear-integral-hamiltonian-response.md)
computes the complete integral module when `N=1`.  MISTAKE-374 is the closest
corrected near miss: generic vanishing of `[1]` loses integral vertical
torsion and is not a polynomial mate.

The remaining puzzle was whether a genuinely multiroot `g` could also make
the generic unit observer vanish.  It cannot.

| item | exact content |
|---|---|
| source | `theta=[1]` in the integral Hamiltonian cokernel for `P=ax+b+g(x)z^d` |
| target | its generic Kummer class and exact `K[P]` annihilator |
| map | localization `K[P] -> K(P)` followed by the weight-one generic-fiber equation |
| preserved | nonvanishing after localization implies zero integral annihilator |
| destroyed | localization kills the resonant one-root Prüfer position |
| required sidecars | fiber character, every root valuation, `x=infinity` degree, and `t=infinity` order |
| cheapest hostile | `d=2`, `g=x(x-1)^3`: the first infinity equation has a solution, but the full equation does not |

The concept board is deliberately small:

| object | invariant | operation | missing coordinate |
|---|---|---|---|
| generic Kummer packet | `N` equal sector ranks | project to the unit character | the actual class inside that rank-`N` sector |
| one-root Prüfer arm | annihilator depth `q` | localize in `P-beta` | whether another root can support the same vanishing |
| weight-one ODE | polynomial degree in `x` | expand at `t=infinity` | its leading polynomial equation |
| multiroot divisor | root count and multiplicities | local valuation at every puncture | global compatibility between punctures and infinity |

## 2. Candidate theorem

Let `K` have characteristic zero, let `d>=2`, and take

```text
P=ax+b+g(x)z^d,                  a in K*,                (1)
C_P=K[x,z]/D_P(K[x,z]),          D_P=Jac(P,-),
theta=[1] in C_P.
```

Assume `g` is nonconstant.  Then the generic class

```text
theta_gen=theta tensor_(K[P]) K(P)                       (2)
```

vanishes if and only if there are `c in K*`, `alpha in K`, and integers
`q>=1` such that

```text
g(x)=c(x-alpha)^(1+qd).                                  (3)
```

Combining this rigidity with the one-root calculation in
THM-3422 gives the exact integral annihilator

```text
Ann_(K[P])(theta)=
  ((P-(a alpha+b))^q)       in the case (3),
  (0)                       for every other nonconstant g. (4)
```

In the first line, `[1]` is nonzero vertical torsion; it is not a polynomial
mate.  In the second line, the distinguished class already survives on the
generic Kummer curve.  Thus a multiple-root `g` may have generic sector rank
`N>1`, but the unit occupies a nonzero vector in that sector and has no
polynomial annihilator.

This is an observer theorem, not a module theorem.  No direct sum of local
root arms, no higher-window multiroot classification, and no new `JC(2)`
case is asserted.

## 3. The exact generic-fiber equation

Replace `P` by `(P-b)/a` and `g` by `g/a`; this scales the Hamiltonian
differential by a unit and reduces to

```text
P=x+g(x)z^d.                                              (5)
```

Write `t` for the generic value of `P` and `F=K(t)`.  As in THM-3419,

```text
A=F[x,g^(-1),z]/(z^d-(t-x)/g).                           (6)
```

The Hamiltonian differential lowers fiber weight by one.  Therefore, if
`[1]` is zero in `A/D_P(A)`, the weight-one component of a primitive is
already a primitive.  Every weight-one element of `(6)` is uniquely

```text
Q=q(x,t)z,                    q in F[x,g^(-1)].           (7)
```

Since `D_P(t)=0`, direct differentiation at fixed `t` gives

```text
D_P(qz)
 =q+(t-x)((g'/g)q-d q').                                 (8)
```

Consequently `theta_gen=0` is equivalent to the rational ODE

```text
q+(t-x)((g'/g)q-d q')=1.                                 (9)
```

This typed equation is the bridge from the rank-only Kummer packet to the
actual unit vector.

## 4. Every root forces one simple zero

Extend scalars faithfully so that

```text
g=c product_(i=1)^N (x-alpha_i)^(e_i).                   (10)
```

At `y=x-alpha_i`, the factor `t-x` is a unit and

```text
g'/g=e_i/y+(regular).                                    (11)
```

If `q` had a pole `c_0 y^(-p)`, `p>0`, then the lowest term in the
parentheses in `(9)` would have the nonzero coefficient

```text
(e_i+d p)c_0 y^(-p-1),                                  (12)
```

which cannot be cancelled.  Thus `q` is regular at every root.  A nonzero
constant term would similarly create an uncancelled simple pole, so write
`ord_(alpha_i)(q)=n>=1`.  The first possible term on the left of `(9)` is

```text
(t-alpha_i)(e_i-d n)c_0 y^(n-1).                         (13)
```

Only `n=1` can supply the constant right side.  If `e_i=d`, that coefficient
also vanishes and no solution exists.  Hence every root forces exactly one
simple zero of `q`, and every allowed finite denominator has disappeared:

```text
q in F[x],                 rad(g) divides q.              (14)
```

This is the local half of the rigidity.  It is not sufficient by itself;
the hostile in Section 7 passes an analogous first asymptotic equation.

## 5. The two infinity locks

Put

```text
r=deg(g),                    m=deg_x(q).                  (15)
```

By `(14)`, `m>=N>=1`.  The coefficient of `x^m` at `x=infinity` in `(9)` is

```text
(1-r+d m) lc(q).                                           (16)
```

It must vanish, so

```text
r=d m+1.                                                  (17)
```

In particular `r=1 mod d`.

Now define the `t`-independent operator

```text
B(q)=(g'/g)q-dq'.                                        (18)
```

Equation `(9)` is

```text
q+(t-x)B(q)=1.                                           (19)
```

On polynomials divisible by `rad(g)`, `B` has no nonzero kernel under
`(17)`: a kernel element would satisfy `q^d=Cg`, forcing `d|r`, contrary to
`r=1 mod d`.

Expand the rational-in-`t` polynomial `q` at `t=infinity`.  If its leading
order is `t^L q_0(x)`, the injectivity of `B` and `(19)` force

```text
L=-1,                       B(q_0)=1.                    (20)
```

Every coefficient in the Laurent expansion still vanishes at all the
`alpha_i`, so `rad(g)|q_0`.  Let `n=deg(q_0)>=N`.  Multiplying `(20)` by `g`
gives

```text
g'q_0-dgq_0'=g.                                          (21)
```

If `n>1`, comparison of the degree `r+n-1` term forces

```text
r=d n,                                                    (22)
```

contradicting `(17)`.  Therefore `n=1`, and since `rad(g)|q_0`, necessarily

```text
N=1.                                                      (23)
```

The unique root descends to `K` in characteristic zero.  Equations
`(14),(17)` now say

```text
g=c(x-alpha)^e,            e=1+qd,            q=m>=1.    (24)
```

This proves necessity.

## 6. Closed converse and integral annihilator

Let `y=x-alpha`, `tau=t-alpha`, and suppose `e=1+qd`.  The THM-3422
polynomial primitive, divided by the invertible generic scalar `tau^q`,
becomes `Q=q_gen z`, where

```text
q_gen(y,tau)
 =tau^(-q) sum_(j=0)^(q-1)
   [binom(q-1,j)/(1+jd)] y^(q-j)(tau-y)^j.                (25)
```

Direct substitution gives

```text
q_gen+(t-x)((g'/g)q_gen-d q_gen')=1.                     (26)
```

Thus `(24)` is sufficient for generic vanishing.  Integrally, THM-3422
places `[1]` exactly `q` arrows before the unique zero arrow and supplies a
primitive for `(P-(a alpha+b))^q`.  Its minimality proves `(4)`.

Conversely, any nonzero integral annihilator would vanish after localization.
The generic necessity just proved therefore makes the annihilator zero in
every nonresonant or multiroot case.

## 7. Hostiles, boundaries, and exact replay

The sharp two-stage hostile is

```text
d=2,                    g=x(x-1)^3.                      (27)
```

The first `t=infinity` equation `(21)` does have the polynomial solution

```text
q_0=x(x-1).                                               (28)
```

Indeed `g'q_0-2gq_0'=g`.  But `deg(g)=4` violates the full infinity gate
`deg(g)=2m+1`, so no generic primitive exists.  This shows why either
infinity comparison alone is insufficient.

Further direct hostile systems include

```text
(d;e_1,e_2)=(2;3,4),             (3;2,5),
(d;e_1,e_2,e_3)=(4;3,5,5).                               (29)
```

All pass the first degree/root-count screen and fail the exact `K(t)` linear
system `(9)`.

The companion
[jc_nonlinear_multiroot_unit_rigidity_thm3424.py](../04-computation/jc_nonlinear_multiroot_unit_rigidity_thm3424.py)
uses exact SymPy arithmetic.  It checks closed primitives across translated
and rescaled one-root models, all local leading coefficients, the
incompatibility of the two infinity degrees, every ordered multiplicity
profile with `2<=d<=6`, `1<=N<=3`, and `1<=e_i<=5`, and the named hostiles.
Run it with

```bash
python3 04-computation/jc_nonlinear_multiroot_unit_rigidity_thm3424.py
```

Boundaries are exact:

- for constant nonzero `g`, THM-3419 gives `C_P=0`, so `[1]=0` rather than
  nonzero torsion;
- `g=0` is a separate non-Kummer zero-response boundary;
- at `d=1`, the conclusion matches THM-3348: generic vanishing occurs exactly
  for a single repeated root.  The weight-sector proof here assumes `d>=2`;
- characteristic zero is load bearing in `(12)--(13)` and in descent of the
  unique root;
- generic vanishing in `(3)` is still only localization-torsion.  By
  MISTAKE-374 it must not be restated as a polynomial mate.
