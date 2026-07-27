---
id: THM-2539
title: "Diagonal cubic owner-clock boundary current"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  On the full-support branch of THM-2517, for
  every live row and every sufficiently large even delay, the diagonal
  cubic scheduler has all 72 mixed source-target means simultaneously.
  Its two late response factors and all seven scheduler factors are exactly
  invariant on the physical mod-13 predecessor fibre, and the literal
  phase-zero owner occurs at slot d=-ell.  The visible source character
  kappa is therefore the opposite clock character -kappa.  Common guard
  support and a consecutive-root Vandermonde force at least three nonzero
  root colours beyond the mean on one seed channel; Galois covariance turns
  them into three distinct target-aligned slopes, uniform over all 72
  channels for that fixed row and delay.  After refining the three response
  sums, the carrier is a labelled direct sum of literal ten-factor Boolean
  atoms.  Applying the THM-2537 positive boundary factorisation atomwise
  gives two positive Boolean packet systems whose signed difference has 216
  nonzero source-target-root incidences and whose occupied endpoints (tail
  for +, head for -) retain the old terminal/deep carrier.  The aggregate
  cubic density and its character
  sums are not one Boolean event.  If one future phase has zero mean, every
  diagonal scheduler cell is null.  A two-atom anchored hostile shows that
  either fixed positive orientation, and in particular the canonical
  selected-source packet, can lose every mixed mean.  No common branch,
  row-independent slopes, inherited ancestry-clock identification,
  semantic arrival, scalar row exclusion, or LRC(14) follows.  The hostile is
  an exact labelled-index obstruction, not a realization of a live O9c row.
source: codex-2026-07-27-diagonal-cubic-owner-clock-boundary
depends_on:
  - THM-2449-coprime-owner-anova-and-delta-replica-boundary
  - THM-2517-cubic-anchored-spectrum-flt3-and-three-time-boolean-lift
  - THM-2533-owner-weighted-phase-and-mixed-gain-radon-ladders
  - THM-2537-cayley-wall-scalarization-and-positive-selector-intertwiner
related:
  - THM-2521-k13-drift-k14-potential-module-bridge
  - THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient
  - THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse
  - THM-2536-deep-comb-selector-flow-target-drift-and-centered-toothpick-boundary
  - THM-2540-weighted-live-event-kakeya-flux-and-transverse-gain-boundary-refinement
script: 04-computation/lrc14_diagonal_cubic_owner_clock_boundary_thm2539.py
output: 05-knowledge/results/lrc14_diagonal_cubic_owner_clock_boundary_thm2539.out
script_sha256: 5ec56c1966ee043c46b83fcf3d1daae1e12f5ff549fb9d2c66986fc80f4c5e50
output_sha256: 2cd425b6e3c7608446066694a275169f0d9a1f4021eff7b942c6cb5322c5a62e
hash_basis: working-tree bytes (LF)
---

# THM-2539 -- the diagonal cubic gives a clock-marked boundary local system

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem composes three already proved mechanisms without erasing their
labels:

```text
THM-2517 diagonal cubic scheduler
  -> THM-2533 physical root fibre and guard Vandermonde
  -> THM-2537 positive boundary pair.                         (1)
```

The output is a source/target-character-active **labelled local system of
Boolean boundary atoms**.  It is not one Boolean event after character
summation.  The theorem applies only to THM-2517's full-support phase branch.
Those two qualifications are load-bearing.

## 1. Scope and exact quantifiers

Fix one of the live THM-2449 response rows.  Write

```text
F_(ell,s)(x)=sum_r h_(ell,r,s,0)(x),

A_(ell,s)=integral_T F_(ell,s)(x) dx,                         (2)
```

where `ell in F_7` is the visible source row, `s in F_13` is the lawful
first-target label, and every `h` is a rational Boolean packet.  The table
`A` has its positive delta anchor and a nonflat owner row.

Retain the seven lawfully shifted future owner--word blocks

```text
G_gamma,                    q_gamma=integral_T G_gamma,
q_0>0,                      gamma in F_7.                     (3)
```

Assume the **full-support branch**

```text
q_gamma>0                         for every gamma in F_7,      (4)

P=product_gamma q_gamma>0.                                  (5)
```

For an even positive integer `L`, put

```text
N=13^L,
D_d=(3+rep(d))L,                         d in F_7,

P_c^L(x)=product_(d in F_7)G_(c+d)(13^(D_d)x),               (6)

R^L_(ell,s)(x)
 =F_(ell,s)(x)F_(ell,s)(Nx)F_(ell,s)(N^2x),

Q^L_(ell,s)(x)=R^L_(ell,s)(x)P_ell^L(x),
K^L_(ell,s)=integral_T Q^L_(ell,s)(x)dx.                     (7)
```

THM-2517 equation (O9c) gives, simultaneously on the finite table,

```text
K^L_(ell,s) -> P A_(ell,s)^3                                (8)
```

along the even subsequence.  The Fermat-cubic rectangle of THM-2517 is
nonzero, so there exists an even `L_0` such that for every even `L>=L_0`,

```text
Khat^L(kappa,b)!=0
       for every kappa in F_7^* and b in F_13^*.             (9)
```

All `72` assertions in (9) hold at the **same** `L`.  This follows either
from uniform convergence of the finite table followed by rational Galois
transitivity, or directly from THM-2517's stated eventual conclusion.
Because there are only `165` live rows, one may take the maximum of their
thresholds and obtain one common `L_0` over all rows satisfying (4).
Nothing here says that every live row satisfies (4).

## 2. Exact predecessor-root invariance and the owner clock

Disintegrate the circle as

```text
x=(z+u)/13,                  z in T, u in F_13.              (10)
```

For every integer `e>=1`,

```text
13^e(x+v/13)=13^e x+v13^(e-1) mod 1=13^e x mod 1.           (11)
```

Consequently the two late response factors and all seven scheduler factors
in (7) are exactly root-invariant.  More explicitly,

```text
Q^L_(ell,s)((z+u)/13)
 =F_(ell,s)((z+u)/13) g^L_(ell,s)(z),                        (12)

g^L_(ell,s)(z)
 =F_(ell,s)(13^(L-1)z)F_(ell,s)(13^(2L-1)z)
   product_d G_(ell+d)(13^(D_d-1)z).                         (13)
```

The only varying predecessor-root coordinate is the first response factor.
This is stronger than asymptotic mixing: (12)--(13) are exact for every
allowed `L`.

The Latin scheduler contains phase `gamma=ell+d` at slot `d`.  Its literal
phase-zero future owner--word block therefore occurs exactly once, at

```text
d=-ell.                                                       (14)
```

If the clock mark is retained, the graph (14) converts characters exactly:

```text
zeta_7^(kappa ell)=zeta_7^(-kappa d).                         (15)
```

Thus visible source colour `kappa` is clock colour `-kappa`.  Equation (15)
does not identify `ell` with the old inverse-ancestry or deep label.  It is
the exact character law on THM-2517's visible row/clock anti-diagonal.

## 3. Common guard support forces three nonzero root colours

Let `E_a` be THM-2533's physical last-root projector and define the mixed
function channel

```text
Z^L_(kappa,b)(x)
 =1/91 sum_(ell,s)
   Q^L_(ell,s)(x)zeta_7^(kappa ell)zeta_13^(bs).              (16)
```

Equation (9) says

```text
integral_T Z^L_(kappa,b)=Khat^L(kappa,b)!=0.                 (17)
```

In particular `E_0Z^L_(kappa,b)!=0`.  Every first response factor obeys the
common guard law of THM-2533,

```text
F_(ell,s)=F_(ell,s)1_(C_H),
```

so (7) and (16) have the same support law.  On a generic predecessor fibre,
the guard forbids three or four cyclically consecutive roots after the unit
relabeling by `H`.

Put

```text
mathcal A_(kappa,b)
 ={a in F_13:E_aZ^L_(kappa,b)!=0}.                           (18)
```

If `#mathcal A_(kappa,b)<=3`, choose that many consecutive forbidden roots.
Reconstruction on those roots gives a square Vandermonde matrix

```text
(zeta_13^(au))_(u forbidden,a in mathcal A_(kappa,b)),       (19)
```

which is invertible.  Every active projector would therefore vanish almost
everywhere, contradicting (17).  Hence

```text
0 in mathcal A_(kappa,b),
#mathcal A_(kappa,b)>=4.                                    (20)
```

Fix one seed `(kappa_0,b_0)` and choose three distinct

```text
a_1,a_2,a_3 in mathcal A_(kappa_0,b_0) minus {0},
lambda_i=a_i/b_0 in F_13^*.                                 (21)
```

The entries of every `Q^L_(ell,s)` are rational step functions.  Therefore
the independent cyclotomic Galois action gives

```text
sigma_(v_7,v_13)(E_aZ^L_(kappa,b))
 =E_(v_13a)Z^L_(v_7kappa,v_13b).                            (22)
```

The action is transitive on the `72` mixed channels.  Equations (21)--(22)
give the three **channel-uniform gain diagonals**

```text
E_(lambda_i b)Z^L_(kappa,b)!=0

for every kappa,b!=0 and i in {1,2,3}.                       (23)
```

For a fixed live row and fixed `L`, the same three `lambda_i` work on all
`72` channels.  They may depend on the row and on `L`; neither a common
triple over all live rows nor stability as `L` varies is asserted.

## 4. The honest Boolean refinement

The aggregate density `Q^L_(ell,s)` is generally **not Boolean**.  Indeed,
`F_(ell,s)` is a finite sum of Boolean deep packets.  Expand all three sums.
For the old deep extension define

```text
q^L_(ell,r_0,r_1,r_2,s,t)(x)
 =h_(ell,r_0,s,t)(x)
  h_(ell,r_1,s,0)(Nx)
  h_(ell,r_2,s,0)(N^2x)
  P_ell^L(x).                                                (24)
```

Every `q^L` is a literal ten-factor Boolean intersection: three response
packets and seven composite handoff blocks.  At `t=0`, summing the three
deep labels gives `Q^L_(ell,s)`.  The inherited deep diagonal is exact:

```text
q^L_(ell,t,r_1,r_2,s,t)=0.                                  (25)
```

After (10), write a refined atom as

```text
q^L((z+u)/13)=e_u(z)g(z),                                   (26)
```

where `e_u` is the first response bit and `g` is the root-independent
product of the other nine factors.  In particular `g` contains the actual
phase-zero owner--word block at `d=-ell`.

For `tau!=0`, apply THM-2537 atomwise:

```text
partial^+_tau q(z,u)=g(z)e_u(z)(1-e_(u+tau)(z)),

partial^-_tau q(z,u)=g(z)(1-e_u(z))e_(u+tau)(z).             (27)
```

Both packets are nonnegative Boolean events on the labelled refined space,
and

```text
partial^-_tau q-partial^+_tau q
 =(P_tau-I)q
 =(I+P_tau)C_tau q.                                         (28)
```

Sum (28) over the three deep labels at `t=0`, then take the source/target
character transform.  Call the resulting signed current

```text
mathcal J^L_(kappa,b;tau)
 =(P_tau-I)Z^L_(kappa,b).                                   (29)
```

For `tau=lambda_i`, equations (23) and (29) give

```text
E_(lambda_i b)mathcal J^L_(kappa,b;lambda_i)
 =(zeta_13^(lambda_i^2 b)-1)
   E_(lambda_i b)Z^L_(kappa,b)
 !=0.                                                       (30)
```

Thus every full-support live row and every sufficiently large even `L` has

```text
3*6*12=216                                                   (31)
```

nonzero source-target-root boundary incidences.  The two positive systems
in (27), retained as an ordered pair, are not both zero on any incidence in
(30).  At least one orientation contributes, but the contributing
orientation may depend on `(kappa,b,i)`.

The character contraction in (16), (29), or (30) is complex and is not a
Boolean event.  The exact statement is that it is the contraction of a
labelled direct sum of positive Boolean atoms.  Forgetting the refinement
may add overlapping indicators with multiplicity.

## 5. Exact sidecar ledger

The coordinates have different types:

```text
ell       visible source row in F_7;
s         lawful first-target label in F_13;
u         physical predecessor root in F_13;
a         Fourier character dual to u;
d=-ell    future owner slot in F_7;
r_0,t     old first-response deep chart;
r_1,r_2   the two later response-root sidecars.              (32)
```

They are not interchangeable.  In particular, `a`, `s`, and `r_0` are
three distinct mod-13 coordinates.

Every packet in (27) retains the external labels `(ell,s,d=-ell)`, the two
late response factors, and all seven scheduler blocks.  Hence both
orientations retain the literal phase-zero owner at the row-dependent clock
and the seven shifted terminal-word times.  The occupied-to-empty tail
`partial^+` is a subset of the original refined atom, so it additionally
retains:

```text
the common guard-safe first response;
the old terminal word and internal clock;
the old (r_0,t) deep chart and diagonal zero (25);
the future sidecars r_1,r_2.                                 (33)
```

For `partial^-`, the stored tail is empty in the first response; its
occupied head retains (33).  Keeping the ordered edge retains that
incidence, but it does not turn the empty tail into a carrier event.  The
occupied endpoint is guard-safe.  The empty endpoint need not be a guard
danger, a target-active terminal component, or a chronological arrival.

## 6. The zero-phase branch annihilates the diagonal carrier

Suppose instead that

```text
q_lambda=0
```

for some phase `lambda`.  Since `G_lambda` is nonnegative Boolean,

```text
G_lambda=0 almost everywhere.                               (34)
```

Every Latin cell `P_c^L` uses every phase exactly once, so every cell
contains the null factor (34).  Therefore

```text
P_c^L=0 almost everywhere for every c,
Q^L_(ell,s)=0,
mathcal J^L_(kappa,b;tau)=0.                                (35)
```

This is exact, not an asymptotic loss.  THM-2517's separate zero-phase
construction still makes a covariant rectangle, but only its owner row has
the actual phase-zero owner at the fixed common epoch.  It does not supply
the full anti-diagonal clock graph (14).  Hence (31) is conditional on (4)
and cannot be stated uniformly over all `165` rows without another branch
argument or a proof that (4) always holds.

## 7. Canonical selection and either fixed orientation can erase the mixed mean

The signed pair in (28) is sharp.  A two-atom construction already prevents
promotion to one fixed positive orientation or to the canonical selected
source packet.

Fix a slope `tau=1` and two adjacent roots `0,1` in a common guard-safe run.
Work on one rational base cylinder and use labelled Boolean atoms.  At the
anchor target `s=0`, give the owner row one singleton mask `{1}` and every
nonowner row no atom.  At every `s!=0`, give each nonowner row one singleton
`{1}`, while the owner row has two labelled atoms, the adjacent pair `{0,1}`
and the singleton `{1}`.  The aggregate occupancy-mass table, after a common
positive scaling, is

```text
K_(0,0)=1,               K_(ell,0)=0 for ell!=0,

K_(0,s)=3,               K_(ell,s)=1 for ell!=0,s!=0.        (36)
```

Every anchored rectangle is

```text
3-1-1=1,                                                     (37)
```

so rational Galois propagation makes all `72` mixed means nonzero.

For both `{1}` and `{0,1}`, the canonical occupied-to-empty `tau=1` wall
has tail `1`.  The selected-source, equivalently the positive `partial^+`
packet in this construction, has table

```text
S_(0,0)=1,               S_(ell,0)=0 for ell!=0,

S_(0,s)=2,               S_(ell,s)=1 for ell!=0,s!=0.        (38)
```

This is exactly delta plus six replicas: its anchored rectangles are
`2-1-1=0`, so all `72` mixed means vanish.  Reflection gives the analogous
hostile in which the fixed `partial^-` orientation vanishes.  Thus (30)
guarantees the ordered positive pair and its signed difference, not one
orientation chosen uniformly in advance.  The nonlinear canonical selector
does not preserve the mixed spectrum.

The hostile uses rational Boolean atoms, the delta anchor, nonnegativity,
and common guard support.  The external anti-diagonal clock label `d=-ell`
can be retained without changing (36)--(38).  Multiplication by an arbitrary
actual row-dependent scheduler cell `P_ell` need not preserve those mass
tables, so no such claim is made.  The hostile is not a realization of every
lawful terminal/deep continuation coordinate or of an actual live O9c table.
It proves that the labelled-index, guard, Galois, and boundary mechanisms of
Sections 1--4 alone cannot yield the stronger conclusion.

## 8. Exact referee and stopping boundary

Run

```bash
python3 04-computation/lrc14_diagonal_cubic_owner_clock_boundary_thm2539.py
python3 -O 04-computation/lrc14_diagonal_cubic_owner_clock_boundary_thm2539.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_diagonal_cubic_owner_clock_boundary_thm2539.out
```

The dependency-free referee uses exact integers, `Fraction`, and arithmetic
in `F_547`, which contains primitive seventh and thirteenth roots.  Universal
late-factor root invariance is the symbolic identity (11); the referee samples
it at `L=2,4,6` on two rational grids.  It also checks the Latin anti-diagonal
and source/clock character law, every guard Vandermonde of size at most three,
all `72` mixed means of (36), Galois-uniform root diagonals and `216` boundary
incidences on a positive control, atomwise positive boundary factorisation,
zero-phase annihilation, the non-Boolean aggregate warning, and both fixed-
orientation/canonical-selector hostiles.  The hostile vanishing in
characteristic zero follows from the exactly checked single-root support and
delta-plus-six-replicas integer table; it is not inferred from a zero modulo
`547`.

The theorem closes the following conditional algebraic seam:

```text
full phase support
  -> visible source/clock anti-diagonal
  -> all mixed source-target means
  -> three Galois-uniform physical-root slopes
  -> an owner-clock-attached positive Boolean boundary pair. (39)
```

THM-2540 proves a separate product-torus/weighted-live-event Kakeya
refinement.  It is not a dependency or evidence for (39).  Conversely, the
present O9c theorem does not claim THM-2540's broader uniform live-event scope.

It does not prove full phase support, one Boolean character-summed event, a
row-independent slope triple, equality of the visible source and inherited
deep/ancestry coordinates, a target-active empty head, chronological owner
arrival, scalar row exclusion, or LRC(14).

**QED.**

The independent audits reconstructed the O9c limit, predecessor-root module,
guard Vandermonde, independent `C_7 x C_13` Galois action, atomwise positive
boundary identity, and both orientation hostiles.  Fresh normal and optimized
runs byte-matched the stored transcript.  The audits also enforced the sharp
scope above: the hostile carries an external clock label rather than an
arbitrary live scheduler factor, and the occupied carrier endpoint is the
tail for `partial^+` but the head for `partial^-`.
