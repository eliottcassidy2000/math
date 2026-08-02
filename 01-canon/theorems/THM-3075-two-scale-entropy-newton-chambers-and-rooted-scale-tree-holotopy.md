---
id: THM-3075
title: "Two-scale entropy-Newton chambers and rooted scale-tree holotopy"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.  The two
  noncommuting ordered six-slot limits of THM-3069 open into two explicit
  simultaneous positive sectors.  A factorial-entropy invoice controls every
  term containing an outer lower-equation error, while THM-3073 exactly
  removes arbitrarily growing upper-triangular coefficients.  The resulting
  thin sectors have explicit slope thresholds, positive physical resultants,
  full carrier exponent ledgers, and sharp raw-convergence hostiles.  They are
  sufficient rooted scale-tree chambers, not a maximal fan or arbitrary
  six-slot theorem.
source: root-gmc-two-scale-tropical-fan-2026-08-01
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-3063-terminal-suspension-transverse-resultant-and-five-slot-tail-holotopy
  - THM-3069-one-normal-remote-terminal-suspension-and-physical-tropical-flag
  - THM-3073-upper-triangular-resultant-norm-and-torsion-character-death-barcode
related:
  - THM-2985-multiparameter-normal-map-and-arc-factor-separation
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3054-affine-moving-lower-tropical-beta-gamma-tail-holotopy
script: 04-computation/gmc_two_scale_entropy_newton_chambers_thm3075.py
output: 05-knowledge/results/gmc_two_scale_entropy_newton_chambers_thm3075.out
script_sha256: 1b9fb88fcdcbaeaa2a867f29b2e3e1b89f469f3e00ccc18f2e1c4fcdec08ba82
output_sha256: 4cfc2181e2e7e608261438bbcd99c9f810941094e1f31f6488851020aa013d8c
hash_basis: LF-normalized bytes
---

# THM-3075 -- two-scale entropy-Newton chambers and rooted scale-tree holotopy

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT REQUESTED.**

THM-3069 exhibits two different ordered six-slot cones.  In one, a remote
fourth low slot is frozen before a terminal pair moves; in the other, the
terminal pair is frozen before a sixth slot moves.  Ordered limits alone do
not imply that either bracketing occupies an open simultaneous-scale sector.

Both bracketings do occupy such a sector, but not by raw coefficientwise
convergence.  Some upper coefficients grow at every positive scale ratio.
The load-bearing repair is to combine a finite entropy invoice with
THM-3073's exact upper-triangular norm quotient.

## 1. The factorial entropy invoice

For one normalized factorial atom, put

```text
M_r(d_1,...,d_r)=
  L(f_(n+d_1)...f_(n+d_r))/L(f_n^r).                        (1)
```

If

```text
d_i=lambda_i T+O(1),       lambda_i>=0,                     (2)
```

then Stirling's formula gives

```text
log M_r=T Phi(lambda_1,...,lambda_r)+O(log T),
Phi(lambda)=(sum lambda_i)log(sum lambda_i)
            -sum lambda_i log lambda_i.                     (3)
```

Use the convention `0 log 0=0`.  With `u` slow factors of slope `x`, `v`
fast factors of slope one, and the remaining factors fixed, the exponential
weight is

```text
phi_(u,v)(x)=(ux+v)log(ux+v)-ux log x.                      (4)
```

For `0<x<=1` and `u+v<=r<=6`,

```text
0<=phi_(u,v)(x)-phi_(u,v)(0)
 <=r x(log(6/x)+1).                                        (5)
```

Indeed, for `u>0`,

```text
d phi/dx=u log(u+v/x)<=r log(6/x),                          (6)
```

and integration from zero proves `(5)`; `u=0` is constant.

The six-slot resultant uses forms of degrees `2,3,4,5,6`.  Its degree in the
coefficients of the degree-`r` form is

```text
mu_r=720/r=(360,240,180,144,120).                           (7)
```

After expanding every coefficient into factorial atoms, one resultant term
therefore contains

```text
sum_(r=2)^6 r mu_r=3600                                    (8)
```

underlying factorial factors, counted with multiplicity.  Define

```text
J(x)=3600 x(log(6/x)+1).                                   (9)
```

Every atomic resultant weight changes from the ordered boundary by at most
`J(x)` per unit fast scale.  The bank is finite, so all nonexponential costs
are absorbed by `poly(T)`.

## 2. Sector A: remote fourth low and faster terminal pair

Fix

```text
n>=1,       0<=a<b<c,       h>=1,
alpha,beta>0,               x=alpha/beta.                  (10)
```

Let integer offsets satisfy

```text
D_T=alpha T+O(1),           C_T=beta T+O(1),               (11)
```

and use the ordered support

```text
(a,b,c,D_T,C_T,C_T+h).                                     (12)
```

Write `S_3>0` for the arbitrary low-triple resultant of THM-2824,
`U_(r,L)` for the positive Gauss/Beta carrier of THM-3063/3069, and
`E_(56,C)` for THM-3063's exact binary normal resultant in degrees five and
six.

Put

```text
m_A(x)=min{
  x log(256/27),
  log(625/256)-J(x)
}.                                                          (13)
```

There is a unique `epsilon_A in (0,1)` satisfying

```text
J(epsilon_A)=log(625/256),                                 (14)
```

and

```text
epsilon_A=0.00001808125357979315927305509... .              (15)
```

Uniqueness follows from `J'(x)=3600 log(6/x)>0` on `(0,1]`.

For every fixed `0<x<epsilon_A`, the normalized physical six-slot resultant
satisfies

```text
R_A(T)=S_3^120 E_(56,C_T)^24
       U_(4,D_T)^180 U_(5,C_T)^144 U_(6,C_T)^120
       [1+O(poly(T) exp(-beta m_A(x)T))].                   (16)
```

In particular it is positive for every sufficiently large `T`.

### Proof of the triangular split

Use the determinant-one outer terminal-pair coordinates pivoted on the
**fixed** low offset `c`, not on the moving inner terminal `D_T`.  Permuting
the low coordinates is lawful, and any orientation change is killed by the
even six-slot resultant degree `6!=720`.

In this chart the lower equations of degrees `2,3,4` have the form

```text
H_r(D_T;y)+E_r(y,z),           r=2,3,4,                     (17)
```

where `E_r` is a genuine outer terminal-pair error.  Leave the upper equations
of degrees five and six completely arbitrary.  Expand the resultant according
to whether an `E_r` coefficient occurs.

The subpolynomial containing no lower error is exactly upper-triangular.
THM-3073 evaluates it, without requiring any upper coefficient to be small,
as

```text
S_4(D_T)^30 T_normal(C_T)^24.                               (18)
```

The outer one-variable covariance contributes the separate factors
`U_(5,C_T)^144 U_(6,C_T)^120`.  After restoring those standard normal
scalings,

```text
T_normal(C_T)=E_(56,C_T)
 [1+O(poly(C_T)(256/3125)^C_T)].                            (19)
```

Every other resultant term contains at least one lower error.  At the ordered
boundary, THM-3063's `k=6` contraction gives the outer gap `256/625`.  Moving
`D_T` to slope `x` can improve the exponential weight of any fully expanded
error term by at most `(9)`.  Thus the outer margin is at least

```text
log(625/256)-J(x).                                          (20)
```

Finally THM-3069 gives the inner exact physical sidecar

```text
S_4(D_T)=S_3^4 U_(4,D_T)^6
 [1+O(poly(D_T)(27/256)^D_T)].                              (21)
```

Raising `(21)` to the thirtieth power supplies the inner margin
`x log(256/27)`.  Combining `(18)--(21)` proves `(16)`.

## 3. Sector B: remote pair and faster sixth slot

Now use the support

```text
(a,b,c,C_T,C_T+h,R_T),
C_T=alpha T+O(1),        R_T=beta T+O(1).                   (22)
```

Let `E_(45,C)` be the exact degree-four/five binary normal resultant.  Put

```text
m_B(x)=min{
  x log(64/27),
  log(46656/3125)-J(x)
}.                                                          (23)
```

There is a unique `epsilon_B` with

```text
J(epsilon_B)=log(46656/3125),
epsilon_B=0.00006001388194673545523082595... .              (24)
```

For every fixed `0<x<epsilon_B`,

```text
R_B(T)=S_3^120 E_(45,C_T)^36
       U_(4,C_T)^180 U_(5,C_T)^144 U_(6,R_T)^120
       [1+O(poly(T) exp(-beta m_B(x)T))].                   (25)
```

Again the physical resultant is eventually positive.

Choose the outer one-normal chart pivoted on fixed `c`; the moving pair
`C_T,C_T+h` then belongs to the lower variables.  Split the degree-two through
degree-five equations into their pure lower parts and genuine sixth-slot
errors, leaving the degree-six upper equation arbitrary.  THM-3073 makes the
zero-lower-error part exactly

```text
R_5(C_T)^6 T_normal(R_T)^120,                               (26)
```

and outer covariance contributes `U_(6,R_T)^120`.  Here

```text
T_normal(R_T)=1+O(poly(R_T)(3125/46656)^R_T).               (27)
```

Every other term pays THM-3069's outer `k=6` gap `3125/46656`, reduced by at
most `J(x)`.  THM-3063 supplies

```text
R_5(C_T)=S_3^20 E_(45,C_T)^6
 U_(4,C_T)^30 U_(5,C_T)^24
 [1+O(poly(C_T)(27/64)^C_T)].                              (28)
```

Equations `(26)--(28)` prove `(25)` with the two margins in `(23)`.

## 4. Full positive exponent ledgers

THM-3063 gives

```text
E_(56,C)~(6^h-5^h)^30,
E_(45,C)~(5^h-4^h)^20.                                    (29)
```

Gauss multiplication gives

```text
U_(r,L)~kappa_(r,n) r^(rL)L^(-(r-1)/2),
kappa_(r,n)>0.                                              (30)
```

With

```text
K=S_3^120 kappa_(4,n)^180
  kappa_(5,n)^144 kappa_(6,n)^120>0,                        (31)
```

the two sectors have the explicit asymptotics

```text
R_A(T)~K (6^h-5^h)^720
  4^(720D_T) 30^(720C_T) D_T^-270 C_T^-588,                (32)

R_B(T)~K (5^h-4^h)^720
  20^(720C_T) 6^(720R_T) C_T^-558 R_T^-300.                (33)
```

Both total power exponents are `-858`.  The exponentially accurate formulas
are `(16)` and `(25)`, which retain the exact `E`.  Replacing `E` by `(29)`
leaves only a relative `1+O(T^-1)` expansion.

## 5. Why raw convergence fails

The triangular quotient is essential.  In sector A, a pure-slow degree-five
upper coefficient grows in the nested diagonal gauge at exponential rate

```text
5x log(5/4)>0                                               (34)
```

for every `x>0`.  In sector B, a pure-slow degree-six upper coefficient grows
at rate

```text
6x log(3/2)>0.                                              (35)
```

Both vanish under the normal restriction `y=0` and are erased exactly by
THM-3073.  Any proof claiming whole-system coefficientwise contraction is
false even inside the thin sectors above.

Among the `35` mixed upper atom types `(r,u,v)` with `u,v>=1` and `u+v<=r`,
the first nested-gauge walls are respectively

```text
(1+4x)log(1+4x)-4x log x-log 5-4x log 4=0,
x=0.34053358195087616137... ,                               (36)

(1+5x)log(1+5x)-5x log x-log 6-5x log 4=0,
x=0.24197406496062009851... .                               (37)
```

These are upper-atom walls, not proved resultant-sector boundaries.  Grouped
nilpotent cancellation can remove such a potential wall.  The full
entropy-Newton object is the finite signed atomic resultant bank; its chamber
test compares the finitely many analytic weights after exact grouping.

## 6. Rooted scale-tree holotopy

THM-3063 and THM-3069 are two grafting operations on a rooted scale tree.
Their ordered limits contract one internal scale edge.  Equations `(16)` and
`(25)` thicken the two boundary bracketings to honest open chambers: the
internal edge has a small but positive length `x`.

The chamber deformation preserves resultant multidegrees, scheme lengths,
the exposed sidecars, and positivity.  It forgets projective root addresses,
phase labels killed by THM-3073, the exact maximal wall, and all-order Hankel
data.  Crossing `(36)` or `(37)` is not covered; a larger chamber would need
an associated-graded quotient computation rather than the coarse invoice.

## 7. Exact evidence and scope

The exact companion verifies:

1. all five resultant coefficient degrees and the `3,600` factor invoice;
2. the `35` mixed upper atom types;
3. rigorous rational logarithm brackets of width `10^-17` for both slope
   thresholds and `10^-14` for both mixed upper walls;
4. the exact safe rays `x=1/100000` and `x=1/20000`;
5. every `S_3,E,U`, gap, exponential-base, and power exponent in `(16)`,
   `(25)`, `(32)`, and `(33)`; and
6. both pure-upper divergence hostiles.

Run

```text
python 04-computation/gmc_two_scale_entropy_newton_chambers_thm3075.py
python -O 04-computation/gmc_two_scale_entropy_newton_chambers_thm3075.py
```

Both modes must equal the stored transcript after LF normalization.

This theorem proves two thin simultaneous sectors with fixed lows, fixed gap,
and fixed positive slopes.  It does not find the maximal entropy-Newton
chambers, cross a wall, prove arbitrary six-slot SFC, an arbitrary multiscale
cone, a uniform threshold, one all-order Stieltjes tail, arbitrary-radial
GMC(2), NC2, LRC(14), JC(2), or DC(2).

**QED, pending independent audit and status promotion.**
