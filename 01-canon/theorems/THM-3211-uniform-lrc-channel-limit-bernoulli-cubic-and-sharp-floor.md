---
id: THM-3211
title: "Uniform LRC channel limit, Bernoulli cubic, sharp floor, and boundary cocycle"
status: >
  PROVED + VERIFIED-EXACT.  In every admissible primitive cap-two channel
  P<Q<=2P, P+Q>=8 and every ordered lane of THM-3171's nine-edge reflected
  two-star, the common-dilation overlap has a limit.  The limit is an explicit
  periodic-Bernoulli-cubic expression and is at least 1/105.  Equality occurs
  exactly for primitive channel (3,5) and ordered lanes (1,2), (3,1), and
  (2,3).  The three equality numerators have exact all-g quadratic closed
  forms and approach 1/105 strictly from above.  Every ray has an exact signed
  1/g correction given by an endpoint overlap barycenter; all residue
  dependence in the eventual numerator begins in the constant term.  Every
  ordered lane realizes both correction signs.  The correction is an exact
  periodic-Bernoulli-quadratic phase coboundary.  Around a complete reflected
  cell orbit it sums to zero, the bulk sums to T/49, and the exact finite
  aggregate differs from T/49 by at most T h_g^2/(3n_gm_g)=O(g^-2).
  A hostile lane approaches
  17/1680 from below and has g=2 mass 2030/280393 below 1/105, so neither the
  limit nor its correction replaces finite heads, proves physical entry, or
  proves LRC(14).
audit: >
  A dependency-free exact floor-sum engine agrees with the Bernoulli formula
  on 27,342 lanes through Q=100 and exhausts the only four primitive channels
  with PQ<=30.  A hash-pinned independent THM-3171 affine-branch engine proves
  the three equality-ray formulas after checking every residue branch and
  finite head.  Ordinary, optimized, and stored certificate transcripts are
  byte-identical; all three maintained scripts have no assert node.  The uniform limit,
  Fourier/Bernoulli identity, large-PQ bound, and small-channel reduction are
  proved analytically below rather than inferred from the finite audit.  A
  second exact certificate checks the signed correction on all 27,342 rays
  through Q=100 and the common residue-linear term on 190,786 states through
  Q=30; its uniform statement follows from the corrector proof and THM-3200.
  A dependency-free boundary-cocycle companion checks the endpoint Bbar_2
  identity, closed-path telescoping, reflection palindrome, universal T/49
  bulk, exact closed-geodesic reduction, and the all-g aggregate bound on
  246,186 phase states through Q=25.  Its zero-correction and negative-
  aggregate hostiles prevent cellwise or one-sided overclaims.
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
correction_script: 04-computation/lrc_signed_one_over_g_correction_certificate_thm3211.py
correction_output: 05-knowledge/results/lrc_signed_one_over_g_correction_certificate_thm3211.out
correction_script_sha256: 275610333bda349a997c9a73c85709f50ee265f9a1746d037838fa2f773e1be6
correction_output_sha256: 2b89b1ec33a457d15fe8441d16fe74705e0440d693feabe148387eba646f5218
cocycle_script: 04-computation/lrc_boundary_cocycle_complete_orbit_thm3211.py
cocycle_output: 05-knowledge/results/lrc_boundary_cocycle_complete_orbit_thm3211.out
cocycle_script_sha256: e02886ce5f0bc3e9171809bc1cf70939817397b724bcd0602bfdd9e14a62c93e
cocycle_output_sha256: e32f996441b7e487d137afad1a5125b41a70b4bb1aa8bd99f403d847d9274a83
independent_engine_commit: 75d0c078d2c204b5fd37051e4fb2d2e1b64f286e
independent_engine_sha256: d73273a4cf4b88bea2890e001166d96cb07dd9b61f3a248ff1538ec44579796a
hash_basis: LF-normalized bytes
---

# THM-3211 -- uniform LRC channel limits, Bernoulli bulk, and boundary cocycle

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

## 8. Signed first correction and constant-only residue modes

For rational phases `alpha,beta`, define the centered static overlap
barycenter

```text
B_(P,Q)(alpha,beta)
 =integral_0^1 (x-1/2) chi(Px-alpha)chi(Qx-beta) dx.       (29)
```

Then every fixed admissible channel and ordered lane has the exact first
correction

```text
c(P,Q;e,f):=lim_(g->infinity) g(I_g-L)
 =B_(P,Q)((R+e)/168,(S+f)/168)
  -B_(P,Q)(R/168,S/168),                                  (30)

I_g=L+c/g+O_(P,Q,e,f)(g^-2).                              (31)
```

To prove `(30)`, put

```text
F(x,s)=chi(Px-(R+es)/168)chi(Qx-(S+fs)/168),
bar F(s)=integral_0^1 F(x,s)dx.                            (32)
```

Let `G` be the one-periodic, x-mean-zero primitive with
`partial_x G=F-bar F`.  The finite affine interval arrangement makes `G`
Lipschitz and piecewise polynomial; its almost-everywhere `s` derivative is
bounded and Riemann integrable.  The a.e. chain rule, periodicity, and integer
`g` give

```text
g(I_g-L)=G(0,1)-G(0,0)
          -integral_0^1 partial_s G(gt,t)dt.               (33)
```

Two-scale equidistribution sends the last integral to its torus mean, which
is zero because `G` has x-mean zero.  With the explicit normalization

```text
G(x,s)=integral_0^x(F(u,s)-bar F(s))du
       -integral_0^1 integral_0^v(F(u,s)-bar F(s))du dv,
```

one has `G(0,s)=B_(P,Q)((R+es)/168,(S+fs)/168)`.  This proves
`(30)`.

Now set

```text
d_2=168^2 P Q,                    d_1=-168(Pf+Qe).         (34)
```

THM-3200 makes `N_g` eventually quadratic on every residue modulo `M`.
Comparing its ratio with `(30)` forces the linear coefficient on every
residue to equal `d_2c+d_1L`.  Hence rational constants `kappa_r` exist such
that, for all sufficiently large `g`,

```text
N_g=d_2L g^2+(d_2c+d_1L)g+kappa_(g mod M).                (35)
```

Dividing `(35)` by `z_gw_g=d_2g^2+d_1g+ef` proves the rate `(31)`.  It also
sharpens `(28)` to

```text
(E-1)^2(E^M-1)N_g=0                                      (36)
```

eventually: every nontrivial root-of-unity mode has degree zero, not one.
Equivalently, subtracting the displayed global quadratic and linear terms
leaves an eventually `M`-periodic sequence.

Formula `(29)` is terminating and rational.  Intersect the two finite
interval combs in `[0,1]`; each overlap component `[ell,r]` contributes

```text
((r-1/2)^2-(ell-1/2)^2)/2.                               (37)
```

The three equality rays and the canonical hostile have

```text
c(3,5;1,2)=71/264600,      c(3,5;3,1)=23/12600,
c(3,5;2,3)=1/5400,         c(3,5;6,1)=-8213/1411200.      (38)
```

Thus the equality rays approach `1/105` from above, while the hostile
approaches `17/1680` from below.  Exact witnesses of both signs occur in each
of the eighteen ordered lanes.  Therefore lane and orientation alone do not
classify the correction; the phase-dependent barycenter `(30)` is the missing
sidecar.  The certificate's absence of a zero through `Q<=100` is
**FINITE-EXACT only**, not a global nonvanishing theorem.

## 9. Boundary cocycle and complete reflected-cell orbit

The first correction has an exact operation law.  For a phase
`phi=(alpha,beta)` in the two-torus, put

```text
S_phi={x in R/Z: chi(Px-alpha)chi(Qx-beta)=1},
B(phi)=integral_0^1 ({x}-1/2)1_(S_phi)(x)dx.              (BC1)
```

Write the signed BV boundary of the reduced interval representation as

```text
partial 1_(S_phi)=sum_(xi in partial S_phi)sigma_xi delta_xi,
sigma=+1 at an entry and -1 at an exit,                    (BC2)
```

and let

```text
Bbar_2(x)={x}^2-{x}+1/6.                                  (BC3)
```

Then

```text
B(phi)=-1/2 sum_(xi in partial S_phi)sigma_xi Bbar_2(xi). (BC4)
```

Indeed, every reduced component `[ell,r]` contributes
`(Bbar_2(r)-Bbar_2(ell))/2`; splitting at zero changes nothing because
`Bbar_2(0)=Bbar_2(1)`.  Thus `(30)` is the exact phase coboundary

```text
C(phi_0,phi_1)=B(phi_1)-B(phi_0).                         (BC5)
```

It is additive on polygonal phase paths and sums to zero on closed loops.
For fixed `P,Q`, the finite endpoint-order/activation atlas makes `B`
piecewise quadratic.  Hence the zero/sign locus of a fixed phase-step
correction is a finite exact semialgebraic problem, rather than a grid
nonvanishing question.

### 9.1 Closed reflected-cell path

Now allow the reflected-cell address to run.  With `L_0=168`, fixed labels
`e,f`, define

```text
d=gcd(L_0,e,f),              T=L_0/d,
phi_j=(je/L_0,jf/L_0),       delta=(e/L_0,f/L_0).          (BC6)
```

The cell-`j` overlap is

```text
I_(g,j)=integral_0^1
 chi(gPu-e(j+u)/L_0)chi(gQu-f(j+u)/L_0)du.               (BC7)
```

The corrector proof gives

```text
I_(g,j)=L_j+c_j/g+O_j(g^-2),
L_j=integral_0^1 A(phi_j+s delta)ds,
c_j=B(phi_(j+1))-B(phi_j).                               (BC8)
```

Therefore every consecutive block telescopes and the complete orbit obeys

```text
sum_(j=a)^(b-1)c_j=B(phi_b)-B(phi_a),
sum_(j=0)^(T-1)c_j=0.                                    (BC9)
```

The substitution `x->1-x` gives `B(-phi)=-B(phi)`.  Since
`phi_(T-j)=-phi_j`, the correction word is a reflection palindrome:

```text
c_(T-1-j)=c_j.                                           (BC10)
```

Thus any nonzero complete correction orbit contains both signs; there can be
no uniform positive `1/g` improvement around the whole reflected cell.

### 9.2 Universal bulk and all-`g` aggregate bound

Assume the standing nondegeneracy `Qe-Pf!=0`, and set

```text
a=e/d,                     b=f/d,
n_g=gPT-a,                 m_g=gQT-b,
h_g=gcd(n_g,m_g),          K=|Qe-Pf|/d.                  (BC11)
```

The full orbit has universal bulk

```text
sum_(j=0)^(T-1)L_j=T/49.                                 (BC12)
```

To prove this, the integer torus map
`(x,u)->(Px-au,Qx-bu)` has nonzero determinant `(Qe-Pf)/d`, so it pushes Haar
measure to Haar measure; each target interval has mass `1/7`.

Concatenating the finite cells also gives an exact closed geodesic:

```text
J_g:=sum_(j=0)^(T-1)I_(g,j)
 =T integral_0^1 chi(n_gu)chi(m_gu)du.                   (BC13)
```

Moreover `Qn_g-Pm_g=(Pf-Qe)/d`, so `h_g|K`.  Reducing the two frequencies by
`h_g` and applying the absolutely convergent Fourier bound from Section 5
gives, for every `g>=1`,

```text
|J_g-T/49|
 <=T h_g^2/(3n_gm_g)
 <=T K^2/(3n_gm_g).                                      (BC14)
```

The complete orbit therefore cancels its entire `1/g` term and has an
explicit two-sided `O(g^-2)` error, including finite heads.

This aggregate law is not cellwise safety.  The exact control

```text
(P,Q;e,f;j)=(5,6;1,12;56)
```

has `c_j=0` and `L_j=4/189`, but its finite errors at `g=1,2,6,9` have signs
`-,+,-,-`.  The complete-orbit hostile `(3,5;1,12;g=2)` has

```text
J_2-T/49=-10/979811.                                     (BC15)
```

Thus zero first correction does not kill the periodic constant, and the
aggregate bound is not a lower bound.  It also forgets the cell owner, so it
does not replace exact finite heads, prove physical survivor entry, settle
the rung, or prove `LRC(14)`.

The exact companion checks `(BC4),(BC9)--(BC15)` on `246,186` rational phase
states through `Q<=25`, all eighteen lanes and their full `84`/`168`-cell
orbits.  Its finite census has `133,412` positive, `111,024` negative, and
four zero cell corrections; this census is not an all-channel zero
classification.

## 10. Failure boundary and next test

The lane

```text
(P,Q;e,f;g)=(3,5;6,1;2)                                  (39)
```

is the required hostile control:

```text
I_2=2030/280393 < 1/105,
lim I_g=17/1680 > 1/105.                                  (40)
```

Hence a proof may use `(21)` to classify tails or equality modes, but may not
replace a finite physical head by its limit.  The finite period scout in the
certificate for `(P,Q;e,f)=(Q-1,Q;12,1)`, `5<=Q<=30`, is **FINITE-EXACT only**;
it does not prove that minimal periods are unbounded.

The signed correction, its chamber type, and its complete-orbit cancellation
are now classified.  The next lawful experiment is the next
Euler--Maclaurin rung: identify the residue-periodic `g^-2` coefficient as a
signed endpoint-owner functional, or exhibit its minimal obstruction.  Exact
cell owners and finite heads must remain attached while testing partial-cell
telescoping or a different reflected cell.  Nothing here establishes physical
survivor entry, the rung, or `LRC(14)`.

## Exact reproduction

Run

```text
python3 04-computation/lrc_uniform_channel_limit_bernoulli_certificate_thm3211.py
python3 -O 04-computation/lrc_uniform_channel_limit_bernoulli_certificate_thm3211.py
python3 04-computation/lrc_signed_one_over_g_correction_certificate_thm3211.py
python3 -O 04-computation/lrc_signed_one_over_g_correction_certificate_thm3211.py
python3 04-computation/lrc_boundary_cocycle_complete_orbit_thm3211.py
python3 -O 04-computation/lrc_boundary_cocycle_complete_orbit_thm3211.py
```

All three mode pairs must reproduce their declared outputs byte for byte.  The engine
uses exact integer floor sums and rational arithmetic.  The certificate uses
explicit exceptions rather than `assert`, hash-pins its independent historical
engine, and labels its finite period scout separately from the proved theorem.
The correction companion imports the same exact engine and uses no floating
point or random choices.  The cocycle companion is dependency-free and checks
the exact endpoint, closed-path, reflection, bulk, aggregate, and hostile
identities over its declared finite universe.  QED for sections 1--5 and
7--10; section 6 is
verified-exact in its stated infinite integer universe.
