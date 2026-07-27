---
id: THM-2513
title: "Anchored moment spectrum, ternary packing, and delayed self-joining"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every rational
  C_7-by-C_13 table with the positive delta owner anchor and a nonflat
  owner row has all 72 primitive mixed colours either in the table itself
  or in its entrywise square. In the only first-moment failure, the table
  is delta plus six replicas; the square has decisive rectangles 2aw_s
  and exact interaction energy (24a^2/7)sum_s(w_s-wbar)^2. Both sides of
  the disjunction are sharp; one nonflat target column already gives six
  nonzero first/square rectangle pairs. For the lawful THM-2449 table this removes
  the replica branch as a spectral-cancellation obstruction: either the
  live first-moment cut bundle of THM-2512 survives, or the complete bundle
  survives on the independent-pair table. The two channels form a fixed
  nonzero vector local system and admit a cancellation-free ternary scalar
  tag with all 144 primitive table modes and 10,368 tagged cut modes nonzero.
  The square has both its canonical independent-pair realization and an exact
  one-copy nonnegative scalar-reweighted realization that preserves the deep
  diagonal-zero lift. For the lawful BV density, delayed self-joinings on one
  circle converge entrywise to the square at O(13^-L); in the replica branch,
  every sufficiently delayed joining has all 72 mixed colours, an exact owner
  anchor, a Boolean-before-summing realization, and the old-sheet deep lift.
  The second packet is not a same-time typed owner/arrival factor; no owner-loop
  current, row exclusion, or LRC(14) closure follows.
source: codex-2026-07-27-anchored-moment-dichotomy
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2466-delayed-word-simultaneous-drift-service-retention
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
related:
  - THM-2456-two-root-replica-uniform-offset-boundary
  - THM-2511-affine-cut-quadratic-root-service-and-pair-space-boundary
script: 04-computation/lrc14_anchored_moment_spectrum_thm2513.py
output: 05-knowledge/results/lrc14_anchored_moment_spectrum_thm2513.out
script_sha256: 6328c92f6ce5be3826316eb83f93b01be09a29106078bbdb6492b774f2ff3928
output_sha256: 1dcdbe295bd832d2051130d1849bd2ffae886b97724b9a1c5756b3eb63ed7313
hash_basis: working-tree bytes (LF)
---

# THM-2513 -- anchored moments admit a branch-free spectral package

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2449 and THM-2512 leave one exact zero branch: a positive owner delta
plus six identical target replicas.  That branch is invisible only to the
first moment.  Squaring the response masses creates a cross term between the
owner delta and the replica profile, and that cross term has a complete mixed
spectrum.

The result is an exact two-level alternative:

```text
lawful anchored response A
  -> mixed spectrum in A
  or
  -> mixed spectrum in the independent-pair table A o A.        (1)
```

The canonical meaning of the second arrow is an independent-pair moment.
There is nevertheless an exact one-copy **weighted** realization, proved in
Section 7, and a delayed Boolean self-joining on one circle.  The remaining
distinction is same-time Boolean owner/arrival provenance and sheet typing,
not mere scalar integration.  A fixed two-channel vector, and then a ternary
scalar tag, remove even the need to choose a branch.

## 1. Anchored tables and their two moments

Let

```text
A=(A_(ell,s))_(ell in F_7,s in F_13) in Q^(7 x 13)             (2)
```

obey the positive owner anchor

```text
A_(ell,0)=a 1_(ell=0),                    a>0.                  (3)
```

Assume that the owner row

```text
s -> A_(0,s)                                                    (4)
```

is nonconstant.  Nonnegativity is not needed for the algebraic theorem,
although every lawful THM-2449 table is nonnegative.

Write `I(A)` for the doubly centred ANOVA interaction of THM-2449, and put

```text
B=A o A,                 B_(ell,s)=A_(ell,s)^2.                 (5)
```

Then

```text
I(A)!=0                         or I(B)!=0.                     (6)
```

Because both tables are rational, the coprime Galois all-or-all law turns
(6) into the stronger spectral statement:

```text
every one of the 72 mixed coefficients of A is nonzero,

or

every one of the 72 mixed coefficients of B is nonzero.        (7)
```

The alternatives need not be exclusive; a generic table fires both.

## 2. The square breaks the replica boundary

If `I(A)!=0`, the first line of (7) is THM-2449.  Suppose instead that

```text
I(A)=0.                                                         (8)
```

The anchored ANOVA classification gives a rational vector `w`, with
`w_0=0`, such that

```text
A_(0,s)=a+w_s,
A_(ell,s)=w_s                         for every ell!=0.          (9)
```

For nonnegative `A`, every `w_s>=0`, but the proof below only needs
rationality.  For each `ell!=0`, the anchored rectangle of `B` is

```text
B_(0,s)-B_(ell,s)-B_(0,0)+B_(ell,0)

 =(a+w_s)^2-w_s^2-a^2+0

 =2a w_s.                                                       (10)
```

The owner row in (9) is nonconstant, so `w` is nonconstant.  Since `w_0=0`,
some `w_s!=0`; since `a!=0`, equation (10) is nonzero.  The rectangle/ANOVA
equivalence proves `I(B)!=0`, and rational Galois transitivity proves the
second line of (7). QED.

There is a more local proof that exposes the underlying geometry.  Fix one
target column `s` and one `ell!=0`, and write

```text
x=A_(0,s),                         y=A_(ell,s),
Delta_1=x-y-a,                     Delta_2=x^2-y^2-a^2.        (S1)
```

These are the anchored rectangles of `A` and `B`.  Since `a!=0`,

```text
(Delta_1,Delta_2)=(0,0)       iff       (x,y)=(a,0).           (S2)
```

Indeed, `Delta_1=0` gives `x=y+a`, after which
`Delta_2=2ay`.  Choose one `s_*` with `A_(0,s_*)!=a`, which exists by
owner-row nonflatness.  Then the same target column supplies

```text
(Delta_1(ell,s_*),Delta_2(ell,s_*))!=(0,0)
                         for every ell!=0.                     (S3)
```

Thus degree two is already a local six-source certificate before ANOVA,
Galois propagation, or cut transforms.  Geometrically, the parabola
`u -> (u,u^2)` separates every anchored secant from the anchor point.  The
ternary tag of Section 6 turns each rational pair in (S3) into two nonzero
scalar rectangle witnesses without choosing a source leg or moment branch.

The nonflatness hypothesis is sharp.  If `w=0`, then

```text
A_(0,s)=a for all s,              A_(ell,s)=0 for ell!=0,       (11)
```

and both `A` and `B` have zero interaction.  This is exactly the flat owner
row excluded by (4).

## 3. Exact pair-space energy in the replica branch

The square interaction can be written completely.  In (9), set

```text
wbar=(1/13)sum_s w_s.
```

The common `w_s^2` part of all seven rows is a column main effect.  The extra
owner-row part is

```text
f_s=a^2+2a w_s.                                                 (12)
```

After ANOVA centring,

```text
I(B)_(0,s)=(6/7)(f_s-fbar),
I(B)_(ell,s)=-(1/7)(f_s-fbar)             for ell!=0.           (13)
```

Therefore

```text
sum_(ell,s)|I(B)_(ell,s)|^2
 =(24a^2/7)sum_s(w_s-wbar)^2.                                  (14)
```

THM-2449's normalized Parseval identity gives the corresponding complete
mixed spectral energy:

```text
sum_(kappa!=0,b!=0)|Bhat(kappa,b)|^2
 =(24a^2/637)sum_s(w_s-wbar)^2>0.                              (15)
```

Thus the second-moment branch has an exact positive invoice, not only an
existence argument.  There is no denominator-free numerical floor for
arbitrary rational `A`.

## 4. Both halves of the disjunction are necessary

Neither moment can replace the pair in (7).

### First moment can fail while the square succeeds

Take any nonconstant rational `w` with `w_0=0` and use (9).  Then

```text
I(A)=0,                         I(A o A)!=0                      (16)
```

by (10).

### The square can fail while the first moment succeeds

Let `a=3`.  At `s=0`, put

```text
A_(0,0)=3,                     A_(ell,0)=0 for ell!=0.
```

At every `s!=0`, put

```text
A_(0,s)=5,                     A_(ell,s)=4 for ell!=0.          (17)
```

The first-moment anchored rectangle is

```text
5-4-3+0=-2,
```

so `I(A)!=0`.  But the square has anchor `9` and, away from `s=0`, owner
entry `25` and replica entries `16`:

```text
(A o A)_(0,s)=9+16,
(A o A)_(ell,s)=16.                                         (18)
```

Thus `A o A` is itself an exact replica table and `I(A o A)=0`.  Both tables
in (16)--(18) are nonnegative and have nonflat owner row.  The first-or-second
form of (7) is sharp in both directions.

## 5. Consequence for the lawful live table

For the delayed lawful table `A^R` of THM-2449, the anchor (3) is positive
and the owner row (4) is nonflat.  Hence (7) applies at every clock after the
eventual-overlap threshold.

There are two exact branches.

```text
non-replica:
  I(A^R)!=0;
  THM-2512 puts all 5,184 primitive cut coefficients on the live
  first-moment response table;

replica:
  I(A^R)=0;
  I(A^R o A^R)!=0;
  the same ANOVA-transpose/cut algebra puts all 5,184 primitive
  coefficients on the square table, canonically on pair-space and also
  on the one-copy scalar-reweighted lift of Section 7.           (19)
```

Thus the delta-plus-six-replicas branch is no longer a spectral-cancellation
obstruction at moment degree at most two.  No higher colour census is needed
to decide it.  What remains is to turn the weighted second line of (19) into
a lawful Boolean ancestry current with the required owner/arrival semantics.

The entrywise square also retains an anchored table:

```text
B_(ell,0)=a^2 1_(ell=0).                                       (20)
```

Affine **index relabellings** act compatibly on `B`, but `B` need not inherit
THM-2449's linear `M+C/R` clock law: entrywise squaring produces terms through
`R^(-2)`.  No linear-clock conclusion is claimed.

## 6. Branch-free vector bundle and the minimal ternary tag

The alternative in (7) can be retained without choosing a branch.  For each
mixed character put

```text
Vhat(kappa,b)=(Ahat(kappa,b),Bhat(kappa,b)).                    (T1)
```

Equation (7) says, more invariantly,

```text
Vhat(kappa,b)!=0                    for all kappa,b!=0.         (T2)
```

This is a fixed `Q(zeta_91)^2`-valued local system.  Applying the same
THM-2508 cut factor to both coordinates gives

```text
Psi^V_(tau,a_0)(alpha,beta)
 =(Psi^A_(tau,a_0)(alpha,beta),Psi^B_(tau,a_0)(alpha,beta))
 !=(0,0)                                                        (T3)
```

for every one of the `5,184` primitive cut quadruples.  Thus no adaptive
first-versus-second-moment choice is needed, even before scalarization.

There is also a canonical cancellation-free scalar package.  Let
`omega=zeta_3` and add a bookkeeping coordinate `t in F_3`:

```text
C_(0,ell,s)=A_(ell,s),
C_(1,ell,s)=B_(ell,s),
C_(2,ell,s)=0.                                                  (T4)
```

For `gamma!=0`, the normalized three-dimensional transform is

```text
Chat(gamma,kappa,b)
 =1/3 (Ahat(kappa,b)+omega^gamma Bhat(kappa,b)).                (T5)
```

Put `F=Q(zeta_91)`.  Since `3` and `91` are coprime,

```text
[F(omega):F]=phi(273)/phi(91)=2,                               (T6)
```

so `{1,omega}` is an `F`-basis.  If `gamma=1`, then
`x+omega y=0` forces `x=y=0`; if `gamma=2`, use

```text
x+omega^2 y=(x-y)-omega y.                                    (T7)
```

Together with (T2), this proves

```text
Chat(gamma,kappa,b)!=0
 for all gamma,kappa,b!=0:             2*6*12=144 modes.       (T8)
```

Fourier-transforming the cut bank also in `t` packages (T3) into

```text
2*5,184=10,368
```

nonzero tagged primitive cut coefficients.

The ternary tag is the smallest prime cyclic tag covered by the
cancellation-free degree criterion below.  The raw nonnegative binary
placement has only the nontrivial coefficient `x-y`, and this placement fails
inside the anchored class: take `a=1` and let `A` be supported only at
`(ell,s)=(0,0)`.  Then `B=A`, every mixed coefficient is `1/91`, but the
binary-tag difference vanishes identically.  The following general
channel-packing lemma explains (T5)--(T8).

> If `F=Q(zeta_N)`, `p` is prime with `p` not dividing `N`, `r<p`, and
> `x_0,...,x_(r-1) in F` are not all zero, then
> `sum_(j<r)x_j zeta_p^(gamma j)` is nonzero for every
> `gamma in F_p^*`.

Indeed, `zeta_p^gamma` has degree `p-1` over `F`, whereas a nonzero
polynomial with those `r` coefficients has degree at most `r-1<p-1`.

The `F_3` coordinate is an algebraic selector, not a lawful new LRC time,
residue, or root.  In particular, (T8) does not supply modulus-`273` root
service.  It records the branch-free two-channel object faithfully while
preserving the same-time typing boundary below.

## 7. Weighted one-copy lift and delayed Boolean self-joining

Suppose one response entry is represented by a nonnegative density

```text
A_(ell,s)=integral_T F_(ell,s)(x) dx.                           (21)
```

Then its square is

```text
B_(ell,s)
 =integral_(T x T) F_(ell,s)(x)F_(ell,s)(y) dxdy.               (22)
```

Equation (22) is the canonical independent-pair coefficient.  It is not the
same-time diagonal quantity

```text
integral_T F_(ell,s)(x)^2 dx,                                  (23)
```

but it does have the exact one-copy scalar-reweighted realization

```text
B_(ell,s)=integral_T A_(ell,s)F_(ell,s)(x) dx.                 (24)
```

Thus pair-space is not an obstruction to scalar integration.  When `F` is
nonnegative rational BV, the density in (24) is too.  It is generally not
Boolean: the global mass `A_(ell,s)` is a weight, not a local ancestry
factor.

For the lawful THM-2449 array, write its integrand as

```text
H_ell(r,s,t)=integral_T h_(ell,r,s,t)(x) dx,
A_(ell,s)=sum_r H_ell(r,s,0).                                  (25)
```

The reweighted array

```text
H^[2]_ell(r,s,t)=A_(ell,s)H_ell(r,s,t)                         (26)
```

satisfies

```text
sum_r H^[2]_ell(r,s,0)=B_(ell,s),
H^[2]_ell(t,s,t)=0.                                            (27)
```

Therefore THM-2449's deep transform applies verbatim:

```text
J^[2](kappa,0,b)=Bhat(kappa,b)/13,
sum_alpha J^[2](kappa,alpha,b)=0.                              (W1)
```

Whenever `Bhat(kappa,b)!=0`, some `alpha!=0` and `tau` give a nonzero
deep coefficient.  In the replica branch this holds for every one of the
`72` mixed `(kappa,b)`.  Hence the second moment already has an exact
one-copy nonnegative **weighted** deep lift.

There is also a genuine Boolean-before-summing realization on one circle,
at two times.  Put

```text
F_(ell,s)(x)=sum_r h_(ell,r,s,0)(x),

D^L_(ell,s)
 =integral_T F_(ell,s)(x)F_(ell,s)(13^L x) dx.                 (W2)
```

Every `F_(ell,s)` is nonnegative rational BV.  The exact two-BV covariance
estimate used in THM-2466 and THM-2478 gives

```text
|D^L_(ell,s)-A_(ell,s)^2|
 <=Var(F_(ell,s))^2/(12*13^L).                                (W3)
```

Thus `D^L -> B` simultaneously over all `91` entries, on either prescribed
parity class of `L`.  Every `D^L` is rational.  There is an explicit
branch threshold.  Put `V_(ell,s)=Var(F_(ell,s))`.  ANOVA centring is an
orthogonal projection, so

```text
||I(D^L)-I(B)||_2
 <=(1/(12*13^L))(sum_(ell,s)V_(ell,s)^4)^(1/2).                (W3a)
```

In the replica branch, (14) says

```text
||I(B)||_2
 =a ((24/7)sum_s(w_s-wbar)^2)^(1/2).                          (W3b)
```

Consequently `I(D^L)!=0` whenever

```text
13^L
 >(sum_(ell,s)V_(ell,s)^4)^(1/2)
   /(12a ((24/7)sum_s(w_s-wbar)^2)^(1/2)).                    (W3c)
```

The denominator is positive by nonflatness.  Take the next `L` in either
desired parity class.  Rational Galois transitivity then gives

```text
Dhat^L(kappa,b)!=0                   for all kappa,b!=0.        (W4)
```

The anchor is retained exactly.  For `ell!=0`, nonnegativity and
`A_(ell,0)=0` force `F_(ell,0)=0` almost everywhere, hence
`D^L_(ell,0)=0`.  Meanwhile `D^L_(0,0)->a^2>0`, so it is positive for every
sufficiently large `L`.

To expose the Boolean provenance and retain the old deep sheet, refine (W2):

```text
mathcal H^L_ell(r,r',s,t)
 =integral_T h_(ell,r,s,t)(x)h_(ell,r',s,0)(13^L x) dx,

H^L_ell(r,s,t)=sum_(r') mathcal H^L_ell(r,r',s,t).             (W5)
```

Then `sum_r H^L_ell(r,s,0)=D^L_(ell,s)` and the old pointwise diagonal zero
gives `mathcal H^L_ell(t,r',s,t)=0`, hence `H^L_ell(t,s,t)=0`.  Every refined
entry is the mass of a literal Boolean two-time intersection before either
DFT.  More precisely, the same covariance estimate gives

```text
|H^L_ell(r,s,t)-H^[2]_ell(r,s,t)|
 <=Var(h_(ell,r,s,t))Var(F_(ell,s))/(12*13^L).                 (W5a)
```

The convergence is simultaneous on the finite deep array.  Choose one
nonzero deep coefficient of `H^[2]` for each of the `72` mixed colours in
the replica branch; all chosen coefficients remain nonzero at one common
sufficiently large `L`.  The deep lift (W1) also repeats directly with
`D^L`.  This closes the independent-pair-to-one-circle descent: at a
sufficiently delayed same-circle Boolean stalk built from lawful packet
factors, every second-moment mixed colour and a corresponding old-sheet deep
lift survive exactly.

The scope boundary has moved but not disappeared.  The second factor in
(W2)--(W5) is a full packet at `13^Lx`; it is not a same-time owner/arrival
factor on the first packet.  Its own local-time source/word skew is lawful,
but identifying those labels with the first packet's global-time labels can
rebase ancestry sheets; under that identification source orientation depends
on parity because `13^L=(-1)^L mod 7`.  The even subsequence is
source-covariant as written.  On the odd subsequence one may replace the
future `ell` by `-ell`; in the replica branch its mean is unchanged because
all six nonowner rows coincide.  No equality of source and arrival
projections, owner-loop positivity, common relation address, or terminal
phase follows.

Moreover, the two moments do not determine overlap with an independently
typed same-time probe.  On a four-point probability space, fix

```text
g=(1,1,0,0),
f_same=(1,1,0,0),
f_disjoint=(0,0,1,1).                                         (W6)
```

Both responses have mean `1/2` and square-mass `1/4`, but their same-time
overlaps with the fixed probe are

```text
E(g f_same)=1/2,                  E(g f_disjoint)=0.            (W7)
```

Thus `A` and `B` alone cannot select a typed same-time intersection.  The
missing theorem must use the actual ancestry factors in (W5), not only their
first two scalar moments.

THM-2513 changes the frontier from

```text
perhaps every mixed colour cancels
```

to

```text
a complete mixed spectrum and old-sheet deep lift exist by degree two;
identify the delayed second sheet with a typed owner/arrival fibre.        (W8)
```

No scalar row is removed and LRC(14) remains open.

## 8. Exact companion and independent audit

Run

```bash
python3 04-computation/lrc14_anchored_moment_spectrum_thm2513.py
python3 -O 04-computation/lrc14_anchored_moment_spectrum_thm2513.py
```

Both executions reproduce the stored transcript byte-for-byte.  The exact
companion:

- exhausts all `4,096` binary replica profiles, finding `4,095` nonflat
  second-moment successes and the unique flat boundary;
- checks (10) and the quantitative identity (14) on every profile;
- verifies one sharp hostile in each direction of Section 4;
- checks the two ternary basis formulas, the exact anchored binary-tag
  hostile, and the `144`/`10,368` mode counts;
- checks `60,000` single-target parabola-secant witnesses, including all six
  source legs at each selected target;
- checks `273` exact initial-arc delayed self-joining entries at both
  parities, their BV invoices, anchors, and nonzero interactions;
- tests `10,000` deterministic anchored tables; and
- reproduces the scalar-reweighted identity and typed-probe separation
  (W6)--(W7).

The independent audit rederived the anchored replica normal form, rectangle
identity, both sharp examples, Galois all-or-all step, and moment scope.
It separately checked the Fourier normalization, degree-two cyclotomic
extension, binary-tag hostile, general `r<p` channel lemma, and the warning
that the tag is bookkeeping rather than physical root service.  A second
audit found the exact scalar reweighting (24), rederived the delayed
self-joining from the BV covariance bound, checked anchor/deep-diagonal
inheritance, and isolated the remaining two-time typing/rebase boundary.  It
confirmed that nonflatness and the same-time qualification are load-bearing
and found no unresolved proof defect. QED.
