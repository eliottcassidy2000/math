---
id: THM-3211
title: "Uniform LRC channel limit, Bernoulli cubic, and sharp floor"
status: >
  PROVED + VERIFIED-EXACT.  In every admissible primitive cap-two channel
  P<Q<=2P, P+Q>=8 and every ordered lane of THM-3171's nine-edge reflected
  two-star, the common-dilation overlap has a limit.  The limit is an explicit
  periodic-Bernoulli-cubic expression and is at least 1/105.  Equality occurs
  exactly for primitive channel (3,5) and ordered lanes (1,2), (3,1), and
  (2,3).  The three equality numerators have exact all-g quadratic closed
  forms and approach 1/105 strictly from above.  A hostile lane has g=2 mass
  2030/280393 below 1/105 but limit 17/1680 above it, so the limit theorem
  does not replace finite heads, prove physical entry, or prove LRC(14).
audit: >
  A dependency-free exact floor-sum engine agrees with the Bernoulli formula
  on 27,342 lanes through Q=100 and exhausts the only four primitive channels
  with PQ<=30.  A hash-pinned independent THM-3171 affine-branch engine proves
  the three equality-ray formulas after checking every residue branch and
  finite head.  Ordinary, optimized, and stored certificate transcripts are
  byte-identical; both new scripts have no assert node.  The uniform limit,
  Fourier/Bernoulli identity, large-PQ bound, and small-channel reduction are
  proved analytically below rather than inferred from the finite audit.
source: root/frontier-synthesis-cont-2026-08-02
depends_on:
  - THM-3171-global-high-channel-cell90-floor-and-all-width-uniform-two-star-law
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
related:
  - THM-3135-directed-cycle-weak-order-lane-cover-and-reflected-h-boundary
engine_script: 04-computation/lrc_uniform_channel_limit_engine_thm3211.py
engine_script_sha256: 3fbcdef02c28644b0b585bdda5922edb2892323c64b7a1246a4430eb687bc57f
script: 04-computation/lrc_uniform_channel_limit_bernoulli_certificate_thm3211.py
output: 05-knowledge/results/lrc_uniform_channel_limit_bernoulli_certificate_thm3211.out
script_sha256: e279ac590748c37083d151b018d38b713646f158d69070cba246f2075e5a7b13
output_sha256: ae53bf9e1367f1dfbcb04305a2b89a8ddaf00da1c0cb3f988b448f102f4c7f9f
independent_engine_commit: 75d0c078d2c204b5fd37051e4fb2d2e1b64f286e
independent_engine_sha256: d73273a4cf4b88bea2890e001166d96cb07dd9b61f3a248ff1538ec44579796a
hash_basis: LF-normalized bytes
---

# THM-3211 -- uniform LRC channel limits and the Bernoulli cubic

**PROVED + VERIFIED-EXACT.**

## 1. Universe and nondegenerate cross-coordinate

Fix THM-3171's reflected cell-90 body

```text
H=(1,2,3,4,6,12)
```

and the nine unordered label edges

```text
{1,2}, {1,3}, {1,4}, {1,6}, {1,12},
       {2,3}, {2,4}, {2,6}, {2,12}.                       (1)
```

Both orientations `(e,f)` of every edge are allowed.  Let

```text
P<Q<=2P,                  gcd(P,Q)=1,                  P+Q>=8, (2)
p=gP,                     q=gQ,                         g>=1.
```

Put `R=90e mod 168`, `S=90f mod 168`, with residues in `{0,...,167}`, and

```text
z_g=168Pg-e,              w_g=168Qg-f,
C=Qe-Pf,                  D=QR-PS.                         (3)
```

Let `A_e(p)` be the reflected arc union from THM-3171/3200 and set

```text
I_g=mu(A_e(gP) intersect A_f(gQ)),
N_g=z_g w_g I_g.                                           (4)
```

The integer `C` never vanishes in this universe.  Indeed, `C=0` would make
`P/Q=e/f`.  Among the increasing ratios on `(1)`, the only reduced ratios in
`[1/2,1)` are `1/2` and `2/3`; their primitive channels have sums `3` and `5`,
contrary to `(2)`.  Repetitions such as `2/4` reduce to the same excluded
`1/2` channel.  Thus the period parameter

```text
M=|C|                                                     (5)
```

is always a positive integer.

## 2. The two-scale torus limit

Let

```text
chi(y)=1 if ||y||<=1/14, and 0 otherwise,                 (6)
```

viewed as a one-periodic function.  The exact reflected-arc indicator is,
away from its measure-zero boundary,

```text
1_(A_e(gP))(t)=chi(gPt-(R+et)/168).                       (7)
```

The maps `t -> ({gt},t)` equidistribute on the two-torus as `g->infinity`.
For continuous periodic test functions this follows by checking Fourier
characters; nonconstant characters have oscillatory integral tending to
zero.  Approximating `(6)` from above and below, whose boundary has measure
zero, extends the statement to the product in `(7)`.  Therefore

```text
lim_(g->infinity) I_g = L(P,Q;e,f),                       (8)
```

where

```text
L(P,Q;e,f)
 = integral_[0,1]^2
   chi(Px-(R+es)/168) chi(Qx-(S+fs)/168) dx ds.            (9)
```

This proves existence uniformly in the sense of one theorem for every
channel and lane.  It does not assert a channel-independent convergence rate.

## 3. Absolutely convergent Fourier formula

The Fourier coefficients of `chi` are

```text
h(0)=1/7,
h(n)=sin(pi n/7)/(pi n),                    n!=0,          (10)
sinc(y)=sin(y)/y,                           sinc(0)=1.
```

Integrating first in `x` forces the frequency relation `nP+mQ=0`.  Since
`P,Q` are coprime, write `(n,m)=(kQ,-kP)`.  Integrating the remaining affine
phase in `s` gives the absolutely convergent series

```text
L=1/49
 +2 sum_(k>=1) h(kP)h(kQ) sinc(pi k C/168)
      cos(2pi k(D+C/2)/168).                              (11)
```

Absolute convergence follows from `|h(kP)h(kQ)|<=1/(pi^2 k^2 P Q)`.  Thus no
pointwise Fourier convergence at an interval endpoint is being assumed.

## 4. Exact periodic-Bernoulli cubic

For `t={x}` define the continuous periodic Bernoulli cubic

```text
Bbar_3(x)=t^3-(3/2)t^2+(1/2)t.                            (12)
```

Set

```text
a=P/14,               b=Q/14,
u=(D+C)/168,           v=-D/168,                          (13)
```

and

```text
Psi=
 Bbar_3(u+a-b)+Bbar_3(u-a+b)
+Bbar_3(v+a-b)+Bbar_3(v-a+b)
-Bbar_3(u+a+b)-Bbar_3(u-a-b)
-Bbar_3(v+a+b)-Bbar_3(v-a-b).                            (14)
```

Then the limit is the exact rational number

```text
L(P,Q;e,f)=1/49+28 Psi/(P Q C).                           (15)
```

To prove `(15)`, integrate the final phase in `(11)` before pairing signs.
The nonconstant part becomes

```text
168/(pi^3 P Q C) sum_(k>=1)
 sin(2pi k a) sin(2pi k b)
 [sin(2pi k u)+sin(2pi k v)]/k^3.                        (16)
```

Apply the four-sine product identity separately at `u` and `v`, followed by

```text
sum_(k>=1) sin(2pi k x)/k^3=(2pi^3/3) Bbar_3(x).          (17)
```

The eight terms are exactly `(14)`, with coefficient `28/(PQC)`.

## 5. Uniform sharp floor and equality classification

Taking absolute values in `(11)` and using `|sinc|<=1` gives

```text
|L-1/49|
 <=2/(pi^2 P Q) sum_(k>=1) 1/k^2
 =1/(3P Q).                                               (18)
```

If `PQ>=31`, then

```text
L>=1/49-1/(3PQ)>=1/105+1/7595>1/105.                     (19)
```

Under `(2)`, the complete list with `PQ<=30` is

```text
(P,Q)=(3,5), (4,5), (4,7), (5,6).                        (20)
```

Direct substitution into the rational formula `(15)` across all eighteen
ordered lanes gives

| primitive channel | minimum limit | equality lanes |
|---|---:|---|
| `(3,5)` | `1/105` | `(1,2), (3,1), (2,3)` |
| `(4,5)` | `1/70` | none |
| `(4,7)` | `1/49` | none |
| `(5,6)` | `2/105` | none |

Equations `(19)--(20)` prove the global statement

```text
L(P,Q;e,f)>=1/105,                                        (21)
```

with equality **exactly** in the three displayed primitive-`3:5` lanes.
This is a sharp floor for the ray limits, not for every finite `I_g`.

## 6. Exact closed forms on the equality rays

**VERIFIED-EXACT for every integer `g>=1`.**  On the three equality lanes,
the cleared numerator `(4)` is

```text
(P,Q;e,f)=(3,5;1,2):       N_g=4032g^2+ 96g,
(P,Q;e,f)=(3,5;3,1):       N_g=4032g^2+744g,
(P,Q;e,f)=(3,5;2,3):       N_g=4032g^2+ 48g.              (22)
```

The certificate proves every affine branch is stable on each residue class
modulo `M`, checks the resulting quadratic at four tail points, and checks
every preceding positive head.  It imports the independently published
THM-3171 engine by an immutable commit and verifies its byte hash before
execution.

Although the limits equal `1/105`, the finite masses are strictly larger:

```text
105N_g-z_gw_g =
 11928g-2,                 (1,2),
 81144g-3,                 (3,1),
  8232g-6,                 (2,3),                         (23)
```

and every expression in `(23)` is positive for `g>=1`.  Thus `(22)` gives
closed-form, constant-time evaluation of the sharp sequences and determines
the direction of approach.

## 7. Refined recurrence structure

THM-3200 proves that on each residue `r mod M`,

```text
N_(r+Mh)=A_r(h),                    deg A_r<=2             (24)
```

eventually.  The existence of the common limit `(8)` now identifies the top
coefficient on **every** residue:

```text
[h^2]A_r = 168^2 P Q M^2 L(P,Q;e,f).                      (25)
```

Indeed, the denominator `z_(r+Mh)w_(r+Mh)` has `h^2` coefficient
`168^2PQM^2`, and its ratio with `(24)` tends to `L`.  Equivalently, in the
global variable `g`, every residue polynomial has common leading coefficient

```text
168^2 P Q L(P,Q;e,f).                                    (26)
```

Consequently

```text
N_g-168^2PQL(P,Q;e,f)g^2                                 (27)
```

is an eventual degree-at-most-one quasipolynomial and is annihilated by

```text
(E^M-1)^2.                                                (28)
```

Thus all nontrivial roots-of-unity modes have polynomial degree at most one;
the quadratic growth is a single global mode.  This sharpens the universal
`(E^M-1)^3` numerator recurrence of THM-3200 without making the normalized
mass C-finite.

## 8. Failure boundary and next test

The lane

```text
(P,Q;e,f;g)=(3,5;6,1;2)                                  (29)
```

is the required hostile control:

```text
I_2=2030/280393 < 1/105,
lim I_g=17/1680 > 1/105.                                  (30)
```

Hence a proof may use `(21)` to classify tails or equality modes, but may not
replace a finite physical head by its limit.  The finite period scout in the
certificate for `(P,Q;e,f)=(Q-1,Q;12,1)`, `5<=Q<=30`, is **FINITE-EXACT only**;
it does not prove that minimal periods are unbounded.

The next lawful experiment is to classify the signed `1/g` correction (and
its residue dependence) for all eighteen lanes.  That can separate rays that
approach `(21)` from above from finite hostile heads approaching from below,
while retaining THM-3171's exact head certificate.  Nothing here establishes
physical survivor entry, another reflected cell, the rung, or `LRC(14)`.

## Exact reproduction

Run

```text
python3 04-computation/lrc_uniform_channel_limit_bernoulli_certificate_thm3211.py
python3 -O 04-computation/lrc_uniform_channel_limit_bernoulli_certificate_thm3211.py
```

Both modes must reproduce the declared output byte for byte.  The new engine
uses exact integer floor sums and rational arithmetic.  The certificate uses
explicit exceptions rather than `assert`, hash-pins its independent historical
engine, and labels its finite period scout separately from the proved theorem.
QED for sections 1--5 and 7; section 6 is verified-exact in its stated infinite
integer universe.
