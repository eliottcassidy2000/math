---
id: THM-2808
title: "Three-pole e=2 Maxwell polynomial and finite accessory classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The remaining
  balanced-response chamber e=2,h=3 is finite in every degree.  For each
  ordered positive pole partition (a,b,c), the two non-pole critical points
  of x^a(x-1)^b(x-lambda)^c have equal critical value at exactly N-3 simple
  parameter values.  An explicit degree-(N-3) Maxwell polynomial gives all
  of them, every root is automatically admissible, and the same count is
  the marked noncrossing-chord Nielsen count.  The total
  orientation-preserving unmarked count has an exact parity-dependent
  cubic formula.  This is a response-layer theorem, not Keller-chart entry,
  JC(2), or DC(2).
source: root/jc-e2-three-pole-maxwell-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2799-one-pole-two-double-zero-chebyshev-response-classification
  - THM-2800-two-pole-two-double-zero-stieltjes-recurrence-and-first-nielsen-pair
  - THM-2805-general-two-pole-e2-maxwell-eliminant-and-nielsen-classification
script: 04-computation/jc2_three_pole_e2_maxwell_thm2808.py
output: 05-knowledge/results/jc2_three_pole_e2_maxwell_thm2808.out
independent_script: 04-computation/jc2_three_pole_e2_maxwell_independent_audit_thm2808.py
independent_output: 05-knowledge/results/jc2_three_pole_e2_maxwell_independent_audit_thm2808.out
script_sha256: 053958bfbc85baf06e40eb6d88d56ad23a2a03a8fbb04d9f8cf39cf9b5574f99
output_sha256: f8bb83a8acef3e3d541a94219506bc7dd6534e65672607715faa851a5967522d
independent_script_sha256: b0518f875330d04460decb5a149c575b804ff2ec4da95e5280e53658faaed499
independent_output_sha256: 6a657d27347d4bd72947d26e246cfb3109923eb79a3486bd55ce2d19450f2b2c
hash_basis: LF-normalized bytes
---

# THM-2808 -- three-pole e=2 Maxwell polynomial and finite accessory classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2796 leaves one apparent accessory parameter in the last `e=2`
balanced-response chamber.  It is not a positive-dimensional modulus.
It is the parameter of a one-variable Maxwell polynomial, and that
polynomial has exactly the degree predicted by the marked dessin atlas.

## 1. Exact critical-value reduction

Work over an algebraically closed field of characteristic zero.  In the
balanced passport notation of THM-2796, fix

```text
e=2,                    h=3,
pole partition=(a,b,c), a,b,c>=1,
N=a+b+c>=4,             s=N-4,                 r=N.    (1)
```

Order the three pole points and normalize them to `0,1,lambda`; normalize
the high point over the third branch value to infinity and
`F(infinity)=1`.  Put

```text
D=x^a(x-1)^b(x-lambda)^c.
```

The third fibre has partition `(N)`.  Therefore a monic numerator
`B=S E^2` must satisfy

```text
B-D=-v                                                   (2)
```

for a nonzero constant `v`.  Its two double zeros are consequently two
distinct non-pole critical points `gamma,delta` of `D`, and

```text
D(gamma)=D(delta)=v.                                    (3)
```

Conversely, any two distinct non-pole critical points satisfying `(3)`
give `(2)` with two double zeros.  Thus the whole chamber is the
equal-critical-value, or Maxwell, locus of `D`.

Remove the forced pole factors from `D'`.  The residual critical quadratic
is

```text
D'
 =x^(a-1)(x-1)^(b-1)(x-lambda)^(c-1) K,

K=Nx^2-Ux+a lambda,
U=(a+c)+(a+b)lambda.                                    (4)
```

Write

```text
E=K/N=x^2-(U/N)x+a lambda/N,
Delta=U^2-4Na lambda.                                   (5)
```

The roots of `E` are exactly `gamma,delta` whenever `Delta!=0`.

## 2. The exact degree-`N-3` Maxwell polynomial

Divide `D` by the monic quadratic `E`:

```text
D=H E+v_0(lambda)+R(lambda)x.                           (6)
```

Then

```text
D(gamma)-D(delta)=R(lambda)(gamma-delta).               (7)
```

The raw secant coefficient has the exact factorization

```text
R(lambda)=Delta(lambda) Q_(a,b,c)(lambda),              (8)
```

where

```text
Q_(a,b,c) in Q[lambda],          deg Q_(a,b,c)=N-3.      (9)
```

This defines the three-pole Maxwell polynomial without roots or
resultants.

### Proof of the factor and degree

Polynomial reduction by `E` shows `R in Q[lambda]`.  Explicitly, put

```text
sigma=U/N,                    pi=a lambda/N,
h_(-1)=0, h_0=1, h_1=sigma,
h_j=sigma h_(j-1)-pi h_(j-2).
```

Then

```text
x^m=h_(m-1)x-pi h_(m-2)                     modulo E,
```

so expanding `(x-1)^b(x-lambda)^c` gives `R` by a finite all-degree
recurrence.  In particular, the recurrence

```text
x^k=h_(k-1)x-p h_(k-2)       modulo E,
h_0=1, h_1=U/N,
h_j=(U/N)h_(j-1)-(a lambda/N)h_(j-2)                  (10)
```

bounds its parameter degree by `N-1`.

The quadratic `Delta` is squarefree as a polynomial in `lambda`, since

```text
Disc_lambda(Delta)=-16Nabc !=0.                         (11)
```

At a root of `Delta`, write `E=(x-u)^2`.  Differentiating `(6)` at `u`
gives

```text
R=D'(u)=0,
```

because `(4)` has the factor `E`.  Both roots of `Delta` are therefore
roots of `R`, proving divisibility in `(8)`.

The factor occurs exactly once.  Near a collision put

```text
E=(x-u)^2-w^2,               w^2=Delta/(4N^2),
P=x^(a-1)(x-1)^(b-1)(x-lambda)^(c-1).
```

Then `D'=NP E`, and the secant quotient gives

```text
R
 =[D(u+w)-D(u-w)]/(2w)
 =N/(2w) integral_(-w)^w P(u+t)(t^2-w^2) dt.           (12)
```

This is a polynomial identity, so the integral is only a compact notation
for termwise antidifferentiation.  It yields

```text
lim_(Delta->0) R/Delta=-P(u)/(6N) !=0.                  (13)
```

The collision roots do not put `u` at a pole, so the last quantity is
nonzero.  Hence `gcd(Q,Delta)=1`.

Finally let `d=a+b` and send `lambda` to infinity.  One critical point is

```text
gamma=(d/N)lambda+O(1),
```

and the other is `delta=a/d+O(lambda^(-1))`.  Equation `(7)` gives

```text
R(lambda)
 =[(d/N)^(d-1)(-c/N)^c] lambda^(N-1)
  +O(lambda^(N-2)).                                    (14)
```

Since the leading coefficient of `Delta` is `d^2`, this also gives

```text
[lambda^(N-3)]Q
 =(-1)^c d^(d-3)c^c/N^(N-1) !=0.                      (14a)
```

Thus `deg R=N-1`; removing the exact quadratic factor proves `(9)`.

## 3. Every root is simple and admissible

Neither `0` nor `1` is a root of `Q`.  At `lambda=0`, the critical
quadratic has roots `0,(a+c)/N`, and the latter has nonzero value under
`x^(a+c)(x-1)^b`.  At `lambda=1`, its roots are `1,a/N`, and the latter
has nonzero value under `x^a(x-1)^(b+c)`.  Since `Delta` is nonzero at
both parameters, `(7)--(8)` give

```text
Q(0)Q(1)!=0.                                           (15)
```

Let `lambda_0` be a root of `Q`.  By `(13)`, `Delta(lambda_0)!=0`, so the
two critical points are distinct.  They avoid all three poles, because

```text
K(0)=a lambda,
K(1)=b(1-lambda),
K(lambda)=c lambda(lambda-1).                          (16)
```

Their common critical value `v` is therefore nonzero.

There is also no hidden multiplicity in the parameter.  Write the parameter
dependence as `D(x;lambda)`, continue the two simple critical points locally,
and set

```text
J(lambda)
 =D(gamma(lambda);lambda)-D(delta(lambda);lambda).
```

At a zero of `J`, the critical-point derivatives drop out, and

```text
J'
 =-cv/(gamma-lambda)+cv/(delta-lambda)
 =cv(gamma-delta)/
   [(gamma-lambda)(delta-lambda)] !=0.                 (17)
```

By `(7)`, `J=R(gamma-delta)`.  Since `gamma-delta` and `Delta` are nonzero
at the root, `R=Delta Q` makes `Q'` nonzero there.  Therefore

```text
Q_(a,b,c) has exactly N-3 distinct roots.               (18)
```

For any such root, equation `(6)` has `R=0` and
`v_0=v`.  Define

```text
B=D-v,                    E=(x-gamma)(x-delta),
S=B/E^2,                  T=x(x-1)(x-lambda).           (19)
```

The criticality in `(3)` gives `E^2|B`.  Any repeated root of `B` is a
root of `D'`; the pole roots have `B=-v!=0`, while the only other critical
points are `gamma,delta`.  Moreover, from `D'=NPE`,

```text
D''(gamma)=NP(gamma)(gamma-delta) !=0,
D''(delta)=NP(delta)(delta-gamma) !=0.
```

Thus both roots have multiplicity exactly two in `B`; consequently `S` is
squarefree and disjoint from `ET`, and every converse gate is automatic.

Moreover,

```text
F=B/D,
F'/F
 =Nv/(SET).                                             (20)
```

Thus THM-2796 applies with `C=Nv`.  For a prescribed `kappa!=0`, one
explicit response is

```text
G=(Nv/(2kappa)) E/(DT),
V=(2kappa/(Nv))^2 S D T^2,
F=VG^2,
2VG'+V'G=2kappa.                                       (21)
```

Conversely, `(2)--(7)` show that every normalized response in `(1)` gives
a root of `Q`.  Hence `(18)--(21)` are the complete ordered
three-pole classification.

## 4. Marked Nielsen atlas

The same `N-3` has a direct dessin meaning.  Fix the full-cycle third
inertia and draw the two disjoint zero transpositions as chords of its
`N`-gon.  Crossing chords leave one cycle; noncrossing chords leave three.
If their four cyclic gaps are

```text
g_1,g_2,g_3,g_4>0,
```

then one noncrossing pairing has pole-cycle lengths

```text
g_1,                    g_3,                    g_2+g_4. (22)
```

For a marked pole cycle of length `p_j`, choose it as the last cycle in
`(22)` and split it into `g_2+g_4` in `p_j-1` positive ways.  The other
two gaps are forced by the other two marked pole lengths.  Quotienting by
the centralizer of the full cycle gives exactly

```text
sum_(j=1)^3 (p_j-1)=N-3                               (23)
```

marked Nielsen charts.  Their defects sum to `2N-2`, so they have genus
zero.  This independently explains the algebraic degree.

When pole multiplicities repeat, unmarked affine classes are the orbits of
the `N-3` roots under the exponent-preserving anharmonic subgroup.  For
example, exchanging the poles at `0` and `1` sends

```text
(a,b,c;lambda) -> (b,a,c;1-lambda).                    (24)
```

The theorem counts ordered normalized charts; it does not silently identify
the roots in `(24)`.

### 4.1 Total orientation-preserving unmarked three-pole count

There is nevertheless a closed count after all unordered three-part pole
partitions are combined.  Choose the four endpoints of the two zero chords.
There are two noncrossing pairings of four cyclically ordered points, hence

```text
2 binom(N,4)
```

marked-by-origin chord pairs.  Burnside's lemma for rotation by the full
cycle has only one possible nonidentity contribution.  When `N` is even,
the half-turn fixes the pair consisting of a chord and its opposite; there
are `N(N-2)/4` such pairs.  Two diameters cross and are not in this set.
No other rotation can preserve a two-element noncrossing chord set.
Therefore the total number of orientation-preserving unmarked `h=3`
Nielsen classes, modulo the centralizer of the fixed oriented full cycle,
is

```text
H_3(N)
 =(N-1)(N-2)(N-3)/12
  +1_(2|N)(N-2)/4.                                    (25)
```

The first values for `N=4,5,...` are

```text
1, 2, 6, 10, 19, 28, 44, 60, 85, 110, ...
```

This total quotient does not replace the ordered polynomial
`Q_(a,b,c)`: repeated pole parts can have smaller anharmonic orbits.

## 5. Consequence for the response frontier

The inequality `h<=e+1` from THM-2796 leaves only `h=1,2,3` when `e=2`.
THM-2799 handles one pole, THM-2800 and THM-2805 handle the two-pole
corridors, and the present result removes the final three-pole accessory
parameter.  All four results are proved and independently audited, so the
entire abstract `e=2` response layer is now a finite explicit atlas.

This is not yet a planar Jacobian theorem.  The response construction is
downstream of a particular nonsplit source-fibre chart; it neither proves
that an arbitrary Keller pair enters that chart nor supplies the missing
Faber-flux compatibility.  In particular, finiteness of the Maxwell atlas
is not degree descent.

## 6. Exact controls

Run

```bash
python 04-computation/jc2_three_pole_e2_maxwell_thm2808.py
python -O 04-computation/jc2_three_pole_e2_maxwell_thm2808.py
python 04-computation/jc2_three_pole_e2_maxwell_independent_audit_thm2808.py
python -O 04-computation/jc2_three_pole_e2_maxwell_independent_audit_thm2808.py
```

The companion uses exact rational polynomial arithmetic and checks every
unordered positive pole partition through `N=10`:

1. the raw degree `N-1` remainder coefficient and exact `Delta Q` quotient;
2. `deg Q=N-3`, squarefreeness, and separation from `0,1,Delta`;
3. the `E^2` division and every squarefree/disjoint converse gate;
4. the cleared response derivative;
5. the marked noncrossing-chord trace and count;
6. the orientation-preserving unmarked Burnside formula `(25)` through
   `N=20`; and
7. covariance under swapping the poles at `0` and `1`.

The independent companion uses no CAS.  Its custom exact `Q[lambda]`
engine checks all `363` ordered positive triples through `N=14`, comprising
`3003` simple Maxwell roots/marked charts, including the exact leading and
boundary formulas, gcd gates, `E^2` division, and response identity.  It
separately enumerates all products of two disjoint transpositions modulo the
full-cycle centralizer, labels the three pole cycles explicitly, and proves
set equality with the positive-gap atlas.  Its independent Burnside audit
finds only the identity and, in even degree, the half-turn fixed locus.
Normal, optimized, and stored outputs agree.  The finite controls support
but do not replace the all-degree proof.
