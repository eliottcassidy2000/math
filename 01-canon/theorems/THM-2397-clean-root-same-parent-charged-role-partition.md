---
id: THM-2397
title: "Clean-root same-parent actual-role charged correlation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On a fixed THM-2392 clean-root owner/status/translate cell, each of
  the five actual guard/lower-q roles has nonzero same-parent
  correlation with q_* on every one of the twelve nonzero root-target
  colours. Uniformly, every named role has a negative Hermitian-pair
  floor rho/1014, while some actual role/colour has the sharper signed
  floor rho/845. Refining by at most 13 two-root translates makes one fixed
  literal lower-q role quantitatively live on all twelve colours with
  coefficient magnitude greater than 4rho/371293. THM-2396 makes
  these statements uniform throughout the last septimal lane.
  Every nonnegative root-constant parent filter preserves the signed
  pair role by role, and every rational such filter preserves all
  twelve colours. BV mixing makes the common-parent intersection with
  any fixed positive rational terminal word nonzero at every
  sufficiently large 13^k clock. The exclusive q_* word is exactly the
  literal single-factor deletion support at the common unshifted
  indicator boundary; the owner-resolved cell carries a same-parent
  F_7 x F_13 coefficient, distinct from the external product of
  separately integrated currents. No lawful target-shifted current,
  canonical-full-endpoint, row-exclusion, ledger-decrement, or LRC(14)
  conclusion is asserted.
source: codex-2026-07-26-clean-root-charged-role
depends_on:
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2396-common-core-forty-nine-orbit-word-incompatibility
related:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2383-polarized-complete-subcube-gram-tomography
  - THM-2398-prime-cyclic-rational-restoration-dichotomy
script: 04-computation/lrc14_clean_root_charged_role_thm2397.py
output: 05-knowledge/results/lrc14_clean_root_charged_role_thm2397.out
script_sha256: 0593f15559dbedd37cc9465b07ccdef8e01261abeefadc018147c10d985235c2
output_sha256: 33fa800491a8a8bb0d280e208705425a375e098f1ba9ba5328d4c002e98fa05d
hash_basis: working-tree bytes (LF)
---

# THM-2397 -- a clean root has an actual charged role neighbour

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2392 gives a positive parent set carrying a fixed singleton or
adjacent two-root exclusive `q_*` word. THM-2396 makes that parent mass
uniform. This theorem keeps the five other **actual** scalar roles,
rather than first turning them into disjoint derived selectors.

Two Parseval identities and one prime-cyclotomic rigidity argument do
different jobs:

```text
each actual role is disjoint from the exclusive q_* word
  -> every named role has a negative Hermitian target pair;

the rational circular cross-correlation has a zero base shift
  -> every named role is nonzero on all twelve target colours;

sum of the five actual role words
  -> one role/colour has a sharper global charged floor.         (1)
```

The root-target alignment is therefore already present at one common
parent. The remaining problem is terminal transport, not disjoint
target support at the clean root.

## 1. Fixed clean-root cell

Retain THM-2392's clean-parent set `S`. For `y in S`, write

```text
x_r=(y+r)/13,                         r in F_13.       (2)
```

Let `q_*` be the distinguished top septimal ordinary label. Partition
`S` by

```text
one deterministic low-blocker owner:      2 choices;

q_* singleton/adjacent status:             2 choices;

exclusive-root translate:                 13 choices.  (3)
```

On one fixed cell `Y`,

```text
rho=mu(Y)>=delta/52.                                   (4)
```

The exclusive `q_*` mask `A` is constant on `Y`. It is:

```text
singleton, if q_* belongs to the unique double pair;

two adjacent roots, if q_* is outside the double pair. (5)
```

Put `W=1-A`. With the normalized root-target transform

```text
a_k=(1/13)sum_r A(r) zeta^(-kr),

w_k=(1/13)sum_r W(r) zeta^(-kr),          zeta^13=1,  (6)
```

THM-2392 proves, for every `k!=0`,

```text
w_k=-a_k!=0,                                          (7)

sum_(k!=0)|a_k|^2
 =12/169 in the singleton status,
 =22/169 in the adjacent status.                      (8)
```

## 2. Every actual role has a charged Hermitian pair

Let

```text
L={H} union {q_i:q_i!=q_*}                            (9)
```

be the five actual remaining roles. For `i in L`, let
`C_i(y,r)` be its root-activity indicator and write

```text
c_(i,k)(y)
 =(1/13)sum_r C_i(y,r)zeta^(-kr).                    (10)
```

The root count is fixed:

```text
n_H=4,                    n_(q_i)=2.                 (11)
```

Moreover `A C_i=0` pointwise, because `A` consists precisely of the
roots at which `q_*` is exclusive. Normalized Parseval gives, for every
`y in Y`,

```text
sum_(k in F_13) a_k conjugate(c_(i,k)(y))=0.
```

Since

```text
a_0=|A|/13,                    c_(i,0)=n_i/13,
```

one has the exact rolewise identity

```text
sum_(k!=0) a_k conjugate(c_(i,k)(y))
 =-|A|n_i/169.                                      (12)
```

The word `A` has an exact physical deletion typing. Let
`D_*(y,r)` be the actual `q_*` danger indicator on the root fibre.
Then

```text
A(y,r)
 =D_*(y,r) product_(i in L)(1-C_i(y,r)).            (12a)
```

Indeed, `A` consists exactly of the `q_*`-dangerous roots at which
every other guard/lower-`q` role is safe. If

```text
B(y,r)=product_(i in L)(1-C_i(y,r)),

M_*(y,r)=1-D_*(y,r),
```

then

```text
B-BM_*=BD_*=A.                                      (12b)
```

Thus `A` is the literal single-factor deletion **support** obtained by
inserting the actual `q_*`-safe factor after the other five safe
factors, not an artificial least-role selector.

The fixed low-blocker owner and the deepest-safe blocker are
root-constant on this fibre. Indeed, for every blocker speed
`c_j=13C_j`,

```text
1_(D_(c_j))(x_r)
 =1_(D_(C_j))(y),                                   (12c)
```

because `C_j r` is integral. Multiplying (12b) by those fixed Boolean
factors, and later by the root-constant delayed word from Section 4,
preserves the deletion identity at the **common unshifted
indicator/bare-endpoint boundary**. This does not supply THM-2370's
lawfully co-shifted packet, target quotient, diagonal zero, or
Poisson--Abel typing. Its boundary remains load-bearing:
target-shift-covariant realization and replacement of the common bare
right endpoint by the canonical fully masked owner endpoint are not
performed here.

Define the same-parent mixed coefficient

```text
G_i(k)
 =integral_Y a_k conjugate(c_(i,k)(y))dy.            (13)
```

Integrating (12) yields

```text
sum_(k!=0)G_i(k)=-rho |A|n_i/169<0.                  (14)
```

Thus **every one of the five named roles** has some nonzero colour
`k` with

```text
Re G_i(k)<=-rho |A|n_i/(169*12)
          <=-rho/1014.                              (15)
```

All masks are real, so the conjugate colour `-k` has the same negative
real part. Every actual role therefore shares at least one Hermitian
pair of root-target colours with `q_*`.

## 3. Rational cross-correlation forces all twelve colours

The fixed cell `Y` is a finite Boolean combination of rational comb
intervals. For each actual role `i`, define its normalized circular
cross-correlation with the exclusive word by

```text
K_i(s)
 =(1/13)integral_Y sum_r A(r+s)C_i(y,r)dy,
                                             s in F_13.          (15a)
```

Every `K_i(s)` is a nonnegative rational number. Exclusivity and the
fixed root count give

```text
K_i(0)=0,

sum_s K_i(s)=rho |A|n_i/13>0.                       (15b)
```

Let

```text
P_i(X)=sum_(s=0)^12 K_i(s)X^s in Q[X].              (15c)
```

If `P_i(zeta^k)=0` for one `k!=0`, then `zeta^k` is primitive and its
minimal polynomial over `Q` is

```text
Phi_13(X)=1+X+...+X^12.
```

Since `deg P_i<=12`, either `P_i=0` or `P_i` is a rational multiple of
`Phi_13`. The first alternative contradicts the positive sum in
(15b); the second makes every coefficient equal, contradicting
`K_i(0)=0`.

With the normalized transform

```text
Khat_i(k)=(1/13)sum_s K_i(s)zeta^(-ks),
```

the substitution `u=r+s` gives exactly

```text
Khat_i(k)
 =integral_Y a_k conjugate(c_(i,k)(y))dy
 =G_i(k),

P_i(zeta^(-k))=13G_i(k).                             (15d)
```

Therefore

```text
G_i(k)!=0
 for every actual role i and every k!=0.             (15e)
```

This is stronger than merely finding one Hermitian pair. Its load-bearing
input is rationality: an arbitrary real nonnegative correlation kernel
with zero base shift can have a vanishing nontrivial Fourier coefficient.
The argument is the prime-cyclic all-or-flat mechanism; no sign is
claimed colour by colour.

## 4. Root-constant parent filtering is hereditary

The signed conclusion does not require the whole parent cell. Let
`f:Y->[0,infinity)` be any nonzero integrable weight which depends on
the parent coordinate `y` but not on the root sheet `r`, and put

```text
rho_f=integral_Y f(y)dy,

G_i^f(k)
 =integral_Y f(y)a_k conjugate(c_(i,k)(y))dy.         (15f)
```

Multiplying the pointwise identity (12) by `f(y)` gives

```text
sum_(k!=0)G_i^f(k)
 =-rho_f |A|n_i/169.                                (15g)
```

Thus every actual role still has a negative Hermitian pair with

```text
-Re G_i^f(k)>=rho_f/1014.                           (15h)
```

The whole-bank identities below are pointwise as well, so one filtered
role/colour has the sharper signed floor `rho_f/845`. If

```text
Q_f(k)=rho_f a_k,

R_i^f(k)=integral_Y f(y)c_(i,k)(y)dy,
```

then

```text
Q_f(k)conjugate(R_i^f(k))=rho_f G_i^f(k),            (15i)
```

and the corresponding aggregate floors are `rho_f^2/1014` and
`rho_f^2/845`.

If `f` is the indicator of a finite rational parent cell, or more
generally a nonnegative rational step weight on such cells, then the
weighted correlations in (15a) remain rational. Section 3 therefore
also gives

```text
G_i^f(k)!=0
 for every actual role i and every k!=0.             (15j)
```

Consequently nonnegative root-constant filtering cannot destroy the
clean-root alignment. A surviving obstruction must use root-sheet
dependence, a change of endpoint gauge/pushforward, or signed mixing.

There is an exact delayed-clock corollary. Let `Q` be a nonnegative
rational circle word and let `R=13^k`, `k>=1`, as in THM-2365. On the
root fibre (2),

```text
Q(Rx_r)
 =Q(13^(k-1)(y+r))
 =Q(13^(k-1)y),                                    (15k)
```

because `Q` is one-periodic. The delayed word is therefore
root-constant. Put

```text
f_R(y)=Q(13^(k-1)y).
```

If

```text
rho_R=integral_Y f_R(y)dy>0,                        (15l)
```

then every literal actual role retains its signed pair and all twelve
nonzero target colours after multiplication by the delayed word. Thus
the delayed word cannot spectrally cancel a positive common-parent
survivor.

In fact positive survival is automatic for every sufficiently large
clock. Both `1_Y` and `Q` have bounded variation. Put

```text
V_Y=Var(1_Y),                    V_Q=Var(Q).
```

The standard two-BV Fourier estimate, now at frequency
`M=13^(k-1)`, gives

```text
|rho_R-rho mu(Q)|
 <=V_Y V_Q/(12*13^(k-1)).                           (15m)
```

Indeed, the nonzero Fourier sum is
`sum_(n!=0) Qhat(n) 1_Yhat(-nM)` and
`sum_(n!=0)n^(-2)=pi^2/3`. Therefore

```text
13^(k-1)>V_Y V_Q/(12rho mu(Q))
  -> rho_R>0.                                       (15n)
```

So a fixed positive delayed terminal word cannot even erase the
carrier once the clock is sufficiently large.

## 5. The whole actual-role bank gives a sharper pair

The sum of the five actual role masks has an exact two-case form.

If `q_*` belongs to the double pair, `A` is its exclusive singleton
and every `W` root has exactly one incidence among the other roles:

```text
sum_i C_i=W.                                         (16)
```

If `q_*` is outside the double pair, let `D_y` be the singleton mask
of the unique double root. Then

```text
sum_i C_i=W+D_y,                    A D_y=0.          (17)
```

Write `d_(y,k)` for the normalized transform of `D_y`. In the adjacent
case, normalized Parseval and disjointness give

```text
sum_(k!=0)a_k conjugate(d_(y,k))
 =-a_0 d_(y,0)
 =-2/169.                                            (18)
```

Consequently

```text
sum_i sum_(k!=0) Re G_i(k)
 =-12rho/169,              q_* in the double pair;

 =-24rho/169,              q_* outside the double.   (19)
```

There are `5*12=60` role/colour pairs. Hence some **actual named**
role and nonzero colour satisfy

```text
-Re G_i(k)>=rho/845                                  (20)
```

uniformly. In the adjacent status the stronger floor is

```text
-Re G_i(k)>=2rho/845.                                (21)
```

The chosen negative colour automatically brings its conjugate, so the
role in (20) is live on at least two colours.

The singleton status is stronger colour by colour. Equation (16) gives

```text
sum_i G_i(k)=-rho|a_k|^2<0
```

for every `k!=0`. Choosing a witness on each of the six conjugate
colour pairs and pigeonholing among five roles shows that one actual
role is live on at least four distinct colours.

The corresponding colourwise negativity is false in the adjacent
status. For example, with

```text
A={0,1},                    D={3},                    (22)
```

the numerator of `169` times the real all-role sum at `k=5` is

```text
-2-2cos(10pi/13)+cos(4pi/13)+cos(6pi/13)

 =-2+2cos(3pi/13)+cos(4pi/13)+cos(6pi/13).
```

Using `cos x>=1-x^2/2`, positivity of `cos(6pi/13)`, and
`pi<22/7`, this is strictly greater than

```text
1-17pi^2/169
 >1-17(22/7)^2/169
 =53/8281>0.                                        (22a)
```

Thus (19), not a colourwise sign assertion, is the correct raw-role
statement there.

## 6. Aggregate currents and the same charged pair

Put

```text
Q(k)=integral_Y a_k dy=rho a_k,

R_i(k)=integral_Y c_(i,k)(y)dy.                      (23)
```

Because `A` is fixed on `Y`,

```text
Q(k)conjugate(R_i(k))=rho G_i(k).                    (24)
```

Therefore the same `(i,k)` in (15), (20), or (21) gives an aggregate
same-parent charged product. In particular:

```text
every actual role:
  some -Re[Q(k)conjugate(R_i(k))]>=rho^2/1014;

some actual role/colour:
  -Re[Q(k)conjugate(R_i(k))]>=rho^2/845.             (25)
```

This is the exact same-target support overlap tested abstractly by
THM-2380, on the clean-root `F_13` line and in one common parent gauge.
It proves the coefficient product is nonzero; it does not manufacture
THM-2380's separate physical translated-correlation measurement.

## 7. A 13-translate quantitative all-colour refinement

Section 3 is qualitative and applies to all five actual roles without
refining `Y`. A small finite refinement gives one literal lower-`q`
role a simultaneous numerical floor on all twelve colours.

Fix any one actual ordinary role `q_i!=q_*`. Its mask has exactly two
roots on every clean fibre. Since `q_i` is a thirteen-unit, multiplication
by `q_i` permutes the root labels. The two active points are adjacent in
the `q_i r` ordering, so their difference in the original root coordinate
is the fixed nonzero class `+-q_i^(-1)`. As `y` varies, the mask therefore
ranges over only the thirteen cyclic translates of one fixed two-root
pattern. Partition `Y` by that translate. One cell `Z` has

```text
rho_Z=mu(Z)>=rho/13.                                 (26)
```

Write `F` for its fixed two-root mask and `f_k` for the normalized
transform. Both `A` and `F` are nonempty, proper, and disjoint. For
every `k!=0`, prime-cyclic irreducibility already gives
`a_k f_k!=0`. More quantitatively,

```text
min_(j!=0)|1+zeta^j|=2sin(pi/26)>2/13,               (27)
```

where the strict inequality follows from the chord bound
`sin x>2x/pi` on `(0,pi/2)`. Thus a singleton `A` is even cheaper,
while in the two-root case

```text
|a_k f_k|>4/28561.
```

Consequently the **same fixed literal lower-q role and same fixed
two-root pattern** satisfy, simultaneously for all `k!=0`,

```text
|integral_Z a_k conjugate(f_k)dy|
 =rho_Z|a_k f_k|
 >4rho/371293.                                      (28)
```

No least-role selector, double-root address, or full-word grouping is
used.

## 8. Uniform last-lane and common-core constants

THM-2396 and the complementary THM-2393 branches give, throughout the
last septimal lane,

```text
delta>=1/26754,

rho>=1/1391208.                                      (29)
```

The universal actual-role floors are therefore

```text
every actual role:
  -Re G_i(k)>=1/1410684912,
  -Re[Q(k)conjugate(R_i(k))]
     >=1/1962556135053696;

some actual role/colour:
  -Re G_i(k)>=1/1175570760,
  -Re[Q(k)conjugate(R_i(k))]
     >=1/1635463445878080.                           (30)
```

On the literal common-core chain, THM-2396 gives

```text
delta>=66/4459,

rho>=33/115934.                                     (31)
```

Hence

```text
every actual role:
  -Re G_i(k)>=11/39185692,
  -Re[Q(k)conjugate(R_i(k))]
     >=363/4542954016328;

some actual role/colour:
  -Re G_i(k)>=33/97964230,
  -Re[Q(k)conjugate(R_i(k))]
     >=1089/11357385040820.                         (32)
```

The 13-translate refinement in Section 7 gives one fixed literal
lower-`q` role, simultaneously on all twelve colours,

```text
|G_i^Z(k)|>1/129136447986             universally,

|G_i^Z(k)|>66/21522741331             in the common core.
                                                        (33)
```

All fractions are exact and reduced.

## 9. Owner-resolved septimal coefficients

Instead use THM-2392's owner-resolved cell `Y_(7x13)`. It fixes the
septimal excess address `d`, the root status, and the exclusive-root
translate, and in the common core has mass

```text
rho_o>=delta/338>=33/753571.                         (34)
```

For every `ell in F_7`, its fixed septimal coefficient is

```text
j_ell=(1/7)zeta_7^(-ell d).                          (35)
```

The actual-role argument gives some `(i,k)` with

```text
-Re G_i(k)>=33/636767495,

-Re[Q(k)conjugate(R_i(k))]
 >=1089/479849517974645.                            (36)
```

There are two different products:

```text
j_ell G_i(k)
 =integral_Y j_ell a_k conjugate(c_(i,k)(y))dy       (37)
```

is the natural derived joint same-parent `F_7 x F_13` tensor
coefficient. It is nonzero for every `ell` and has magnitude at least

```text
33/4457372465.                                      (38)
```

By contrast,

```text
j_ell Q(k)conjugate(R_i(k))                          (39)
```

is an external product of separately integrated currents. Its magnitude
is at least

```text
1089/3358946625822515.                              (40)
```

Neither (37) nor (39) is automatically one physical scalar
`91`-unit terminal mode.

Role by role, the corresponding owner-resolved floors are

```text
-Re G_i(k)>=11/254706998,

-Re[Q(k)conjugate(R_i(k))]
 >=363/191939807189858,                             (41)
```

and the external septimal product has magnitude at least

```text
363/1343578650329006.                               (42)
```

Refining this owner-resolved cell by the at most `13` translates from
Section 7 gives one fixed literal lower-`q` role and one fixed pattern
such that every `k!=0` satisfies

```text
|G_i^Z(k)|>132/279795637303.                        (42a)
```

The fixed septimal coefficient then gives, for every `ell`,

```text
|j_ell G_i^Z(k)|>132/1958569461121.                 (42b)
```

This is a natural joint same-parent tensor coefficient, not the
external product in (39).

## 10. Scope and next target

The theorem proves

```text
uniform positive clean-root mass
  -> every actual guard/lower-q role shares all twelve
     root-target colours with q_*
  -> one fixed literal lower-q role has a simultaneous
     quantitative all-colour floor after at most 13 cells.       (43)
```

The `C_i` are literal guard/lower-`q` factor indicators; no derived
least-role selector is used. The correlations `G_i` are
coefficient-derived same-parent groupings. Rationality and primality
are essential for the unrefined all-colour conclusion.

The theorem identifies `A` with literal single-factor deletion support
at the common unshifted indicator boundary, but it does not construct
the lawfully co-shifted THM-2370 packet, replace the bare right endpoint
by the canonical fully masked owner endpoint, or make (37)--(42)
separately supplied physical pair-twist observables. Thus THM-2380's
disjoint-support hostile is impossible at the clean root, and the
`13^k` delayed word preserves a positive all-colour survivor at every
sufficiently large clock. Only target-shift-covariant current
realization/owner routing, canonical-right-endpoint replacement,
additional root-sheet-dependent filtration, or signed mixing remains.

The next sharp target is:

```text
realize the literal deletion support as a lawfully co-shifted current,
then transport it from the common bare right endpoint to the canonical
fully masked owner endpoint without changing its root-target phase.      (44)
```

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 11. Exact companion

The dependency-free exact companion:

- verifies the two actual-role sum laws (16)--(19);
- checks the rolewise disjointness/Parseval identity for counts
  `4,2,2,2,2`;
- checks the exact single-factor deletion identity (12b);
- exhausts every disjoint size-`1/2` exclusive mask and size-`2/4`
  actual-role mask, checking the rational cross-correlation has zero
  base coefficient, positive sum, and is non-flat;
- checks all `13` translates of each fixed two-root gap and the
  simultaneous all-colour
  quantitative floors;
- checks the adjacent raw colour-sign hostile;
- checks the hereditary signed floors under a nonnegative
  root-constant parent weight;
- checks the arithmetic `1/12` constant in the two-BV estimate;
- verifies the `5*12=60` and Hermitian-pair pigeonholes;
- checks every rational floor in (29)--(42b); and
- keeps the joint integral (37) distinct from the external product
  (39).

Run

```bash
python3 04-computation/lrc14_clean_root_charged_role_thm2397.py
python3 -O 04-computation/lrc14_clean_root_charged_role_thm2397.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_clean_root_charged_role_thm2397.out
```

after LF normalization. Every assertion remains active under optimized
Python.
