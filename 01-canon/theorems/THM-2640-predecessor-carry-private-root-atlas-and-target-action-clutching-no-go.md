---
id: THM-2640
title: "Predecessor-carry private-root atlas and target-action clutching no-go"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Refine THM-2623 before delayed integration by the physical predecessor
  carry c=floor(13^6 x) mod 13 and retain the future half-digit
  kappa.  If d=2h+kappa, b=floor(d/13), and epsilon is one on the left
  half-tooth and zero on the right, every nonzero refined row has the
  singleton physical root r=2c+b+epsilon mod 13.  The thirteen carry cells
  partition the old numerator exactly, the global primitive content remains
  26, and the guard-safe/danger cospan gives every one of the 84 base cells
  all twelve nonzero private unit roots.  At fixed c it realizes exactly the
  full geometrically possible set {2c,2c+1,2c+2}\{0}.  This is a maximally
  saturated static local-section atlas, not a principal transition: the
  carry selector is target-neutral, while THM-2365 moves r, so preserving
  the graph would require the unsupported translation c->c+7 delta.  No
  common transition branch, holonomy trivialization, row exclusion, or
  LRC(14) conclusion follows.
source: carry-transition-cell-2026-07-28-private-carry-atlas
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
related:
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2630-old-wall-affine-clutching-and-successor-sector-no-go
  - THM-2637-derangement-character-fixed-branch-holotopy-principle
  - THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate
script: 04-computation/lrc14_predecessor_carry_private_root_atlas_thm2640.py
output: 05-knowledge/results/lrc14_predecessor_carry_private_root_atlas_thm2640.out
script_sha256: a28b03a5903256c1c1c294ea5af389c7991fc0a5ad6908f0f25a5b0cc6e71abf
output_sha256: b3c1f5106a4f7c908435862e4a5824758128b8d42f833f0a7701b7d36b69a940
hash_basis: LF-normalized bytes
---

# THM-2640 -- the missing carry creates private roots, not their transition

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2623 closes the delayed guard support puncture and repairs four nonlinear
unit holes, but every row refined only by the future half-digit still carries
eleven or twelve physical roots.  Its exact diagnosis is that the future
fractional coordinate forgets the integer predecessor carry.  Retaining that
carry produces the strongest possible local repair: every allowed carry/root
branch becomes private and unit.  The repair still does not commute with the
lawful target action.  This separates two notions that support-only arguments
repeatedly conflate:

```text
complete static local sections       -- proved here;
one equivariant transition section   -- still absent.             (1)
```

## 1. The physical carry selector

Keep all notation and both delayed guard sectors of THM-2623.  In particular,

```text
R=13^6,                    S=R/13=13^5,
c3=2*13^5=2S.                                             (2)
```

For a physical point `x`, write

```text
Sx=N+u,              0<=u<1.
```

Then

```text
c=floor(13u)=floor(Rx) mod 13,
y={13u}={Rx}.                                                (3)
```

Thus the thirteen half-open Boolean cells

```text
P_c(x)=1_[c/13,(c+1)/13)({Sx})                           (4)
```

are a physical partition before delayed integration.  They do not arise by
splitting a coefficient after the nonlinear unit test.  The endpoint
convention is immaterial to all integrals, and the companion works throughout
with the canonical half-open interval convention.

Retain also THM-2623's future half-digit

```text
h=floor(13y),
kappa=floor(26y)-2h in {0,1},
d=2h+kappa in {0,...,25}.                                (5)
```

For each rail, delayed guard sector, left/right deep half-tooth, owner clock,
and nonzero deep probe, intersect the common physical packet with `(4)` and
`(5)` before integrating.  Denote the resulting raw integer numerator by

```text
A_(e,sector,epsilon,c,h,kappa;ell5,r),                    (6)
```

where `epsilon=1` is the left half and `epsilon=0` the right half.

## 2. Exact one-step prefix descent

The carry refinement can be evaluated without materializing `13^6` pullback
intervals.  Let `Q_d` be either delayed guard-sector word restricted to the
half-digit `d`, let `L_d` be its integer-grid length, and let `Phi_d` be its
prefix function on the canonical `T` grid.  In the variable `u={Sx}`, its
carry-`c` copy is

```text
Q_(c,d)(u)=1_[c/13,(c+1)/13)(u) Q_d({13u}).               (7)
```

Fix an evaluation endpoint `u`, put `c=floor(13u)`, and index the vector of
carry copies by `c'`.  The prefix of the `c'`-copy at `u` has the exact
triangular form

```text
L_d/13,             if c'<c;
Phi_d(13u-c)/13,    if c'=c;
0,                  if c'>c.                              (8)
```

All displayed divisions are integral on the resolved grid.  The usual
floor/prefix identity at scale `S` has denominator `S*T`, while THM-2623's
canonical raw numerator has denominator `R*T=13*S*T`.  Therefore the raw
carry numerator is thirteen times the descended prefix numerator.  With this
normalization, the exact partition is

```text
sum_(c in F13) A_(...,c,h,kappa;ell5,r)
 =A^old_(...,h,kappa;ell5,r).                              (9)
```

The companion checks `(9)` at every one of `1,415,232` selected half-digit
slots.  Across the complete carry-refined carrier it obtains

```text
positive raw entries:       230,042,
global gcd/content:                26.                    (10)
```

Thus fixing `c` before the unit test neither reprimitivizes a row nor changes
the integral lattice.

## 3. The exact singleton law

Write `Rx=n+y`.  On carry cell `c`, `n=c mod 13`; on half-digit `d`,

```text
floor(2y)=floor(d/13)=:b in {0,1}.                        (11)
```

The absolute high-speed deep digit is therefore

```text
a=floor(13{c3 x})
 =2c+b                         mod 13.                    (12)
```

The left half of the translated danger tooth has `r=a+1`, while the right
half has `r=a`.  Hence every raw row in `(6)` obeys

```text
r=2c+floor((2h+kappa)/13)+epsilon       mod 13.           (13)
```

If the right side is zero, the retained nonzero-root bank has the literal
zero row.  Otherwise it is a private physical root.  This is a pointwise
support identity, before cyclotomic reduction or summation over clocks.  The
companion checks all

```text
18,398,016                                                    (14)
```

candidate `(configuration,r)` entries.  Of these, `1,415,232` are the
possible on-graph slots selected by `(13)`, while the other `16,982,784`
are off-graph slots; every off-graph slot is exactly zero.  Thus `(14)` is
the total singleton-law test count, not the number of excluded slots.

## 4. The private unit atlas is maximally saturated

Normalize every nonzero entry in the full carry carrier by the one global
content `26`.  At fixed

```text
(e,sector,epsilon,c,h,kappa),
```

let `r` be the private root in `(13)` and define the seven-clock polynomial
exactly as in THM-2623:

```text
Y_ell5=(A_ell5/26) r^(-1) mod 13,
Y(z)=sum_(ell5=0)^6 Y_ell5 z^ell5 in F13[z]/(Phi_7).      (15)
```

Call the row unit when `Y` is a unit.  Unite the two rails and both guard
sectors over one base cell `(s,ell4)`.  Exact computation gives

```text
private nonzero unit-root set = F13^*       on all 84 cells. (16)
```

There is no hidden loss at fixed carry.  Formula `(13)` allows, as `h`,
`kappa`, and the edge vary, exactly

```text
R_c^geom={2c,2c+1,2c+2}\{0}.                             (17)
```

The unit atlas realizes equality with `(17)` for every carry and every base
cell.  Therefore

```text
|R_c^unit|=2 on 252 cell/carry pairs,
|R_c^unit|=3 on 840 cell/carry pairs.                     (18)
```

The `252=84*3` two-root cases are exactly those for which zero occurs in the
three-term geometric set.  Across all carries their union is `(16)`.  There
are `3,024=84*36` fixed `(cell,c,r)` private-unit witnesses, one for every
geometrically possible nonzero branch.

The delayed label `q=h` is also nearly complete:

```text
|Q^unit|=13 on 80 base cells,
|Q^unit|=12 on  4 base cells.                             (19)
```

The four missing pairs are exactly

```text
((2,2),11),       ((6,5),2),       ((7,2),11),       ((11,5),2).
```

The fixed
`(c,h)` edge census sharply distinguishes complete local digit transitions
from the aggregate coverage `(16)`.  Among all

```text
84*13^2=14,196
```

labelled digit edges, the private-unit root-set size histogram is

```text
size 0:  2,050,       size 1:    170,
size 2: 11,136,       size 3:    840.                    (19a)
```

Exactly `2,052` fixed edges fail to realize their full geometric root set:
`2,050` are empty, and the only two partial failures are

```text
(cell,c,h)=((6,3),8,12):   {5} instead of {4,5},
(cell,c,h)=((6,3),9,12):   {7} instead of {6,7}.          (19b)
```

The number of surviving directed edges `c->h` per base cell has histogram

```text
119^4,             142^2,             145^2,             146^76. (19c)
```

The `80` graphs with `142`, `145`, or `146` edges are strongly connected,
with diameters respectively `2`, `3`, and `2`.  The remaining four are
exactly the four cells in `(19)` with a missing future label; they have
`119` edges, minimum indegree zero, and are not strongly connected.  Thus
the symbolic carry/future correspondence is extremely dense, but it is not
a complete same-labelled transition atlas.

The guard cospan is load-bearing at the middle carry.  On the guard-safe
sector, every lane at `c=6` has zero unit rows.  The guard-danger sector
supplies the complementary seam lanes:

```text
left/kappa=1:   134 unit rows,
right/kappa=0:  134 unit rows,                            (20)
```

with the other two middle-carry danger lanes zero.  Thus the complement is
not merely filling the old endpoint digits `0,12`; it glues the unique
doubling-carry seam.

## 5. Uniformly forgetting the carry is charged, not invariant

There is one natural attempt to avoid choosing a carry: sum all thirteen
private rows uniformly.  Equation `(9)` shows that this is lawful, but it is
exactly the old THM-2623 half-digit root profile.  For every one of its
`20,778` nonempty fixed-clock instances, write

```text
A(r)=sum_c A(c,r),                  r in F13.             (20a)
```

The exact support census is

```text
A(0)=0,                  |supp A| in {11,12}.             (20b)
```

Thus the uniform carry sum is never translation-invariant.  It is in fact
maximally charged.  If a nonzero character Fourier coefficient of the
rational integer vector `A` vanished, its degree-at-most-twelve coefficient
polynomial would vanish at a primitive thirteenth root.  Irreducibility of
`Phi_13` would force all thirteen coefficients to be equal, contradicting
`(20b)`.  Hence

```text
Ahat(k)!=0                    for every k in F13,          (20c)
```

including `k=0` because `A` is nonnegative and nonzero.  Abstractly, the
thirteen circulant translates of `A` therefore form an invertible signed
regular-representation decoder.

This is the sharp positive survivor, but not yet a physical multisection.
The abstract translates in the circulant inverse only relabel `r`; the lawful
THM-2365 orbit simultaneously moves the present `t`-factor.  Canon supplies
no identification of those translated vectors with the same frozen physical
packet, and the inverse need not be positive.  This is the same signed-versus-
physical boundary isolated independently by THM-2624.

## 6. Why private rows still do not give the target action

The selector `(4)` is target-neutral.  It is a delayed speed-`13^5` observer,
so every THM-2365 target phase is multiplied by a multiple of thirteen and
disappears.  In contrast, the lawful target action contains the generator

```text
(r,s,t)->(r+delta,s,t+delta).                             (21)
```

The root in `(13)` is the absolute deep-probe label in

```text
Delta_r(x)=d(c3 x-r/13);                                  (22)
```

it is not the invariant difference `r-t`.  Consequently `(21)` keeps `c`
fixed and moves a point off the graph `(13)` whenever `delta!=0`.

Algebraically one could preserve `(13)` by adjoining the rule

```text
c->c+7 delta,                                             (23)
```

because `7=2^(-1) mod 13`.  But `(23)` is not the action of the physical
carry observer in `(4)`.  Installing it would add an auxiliary shifted probe
coordinate and change the packet/gauge.  Static private rows therefore do
not supply a principal `C13` bibundle or a branch of the already selected
target transition.

There is a useful near miss which makes the failure completely physical.
Choose the integer representative `delta in {1,...,12}` and translate the
base circle by the canonical integer-linear lift

```text
tau_delta=7 delta/R.                                      (23a)
```

Then `R tau_delta=7 delta` is integral, so the delayed prefix is unchanged,
`c` moves to `c+7 delta`, and

```text
c3 tau_delta=14 delta/13=delta/13       mod 1.            (23b)
```

Thus `(23a)` implements exactly the desired `c->c+7 delta`,
`r->r+delta` covariance on the carry/root/delayed-prefix subpacket.  It is
not a symmetry of the retained present carrier.  The first factor already
fails: the speed-one guard-safe set `C_H={||x||>1/7}` satisfies

```text
mu(C_H triangle (C_H-tau_delta))
   =2 tau_delta=14 delta/R>0.                             (23c)
```

Here `0<tau_delta<1/7`, so `(23c)` is the exact circular-interval
symmetric-difference formula, not an endpoint estimate.  The companion
records `14/4826809` at `delta=1`.  Other low-speed present factors move as
well.  More precisely, a speed `v` receives the phase

```text
v tau_delta=7 delta v/R.
```

Because `delta` and `7` are thirteenth-units, this lies on the lawful
`1/13` phase grid exactly when `13^5` divides `v`.  Among the nine present
speeds of the canonical row this selects only `c3=2*13^5`; the guard, five
ordinary units, `c1`, and `c2` all acquire unresolved finer phases.  Hence
even this canonical small physical lift of the formal clutch changes the
packet on positive measure; the obstruction is not merely a choice of
notation for `c`.  It is not the only lift: THM-2657 classifies all global
translation lifts as `tau=k/R` with `k=7 delta mod 13`.  Every such nonzero
translation also moves the speed-one guard, whose rotational stabilizer is
trivial.  Formula `(23c)` is the exact measure only for the displayed
integer-linear representatives, for which `0<tau_delta<1/7`.

This is exactly the hypothesis boundary in THM-2637.  Its derangement
principle kills nonzero carry holonomy only after **one common transition**
has a fixed branch.  Equations `(16)`--`(18)` provide separate static branch
choices at each base cell; `(21)` fails to transport them.  No holonomy
conclusion may be drawn.

## 7. Consecutive digits and the next two-clock test

There is nevertheless a canonical symbolic adjacency.  At the next scale,

```text
c_next=floor(13R x) mod 13
      =floor(13{Rx})
      =h.                                                 (24)
```

Thus a fixed `(c,h)` row is a genuine consecutive base-thirteen digit edge.
The current theorem proves unit/private survival only in the clock-`R` chart.
To turn `(24)` into a transition, one must prove that the clock-`13R` packet
has a unit row on the **same physical component and same retained labels**.
THM-2624's two signed clock reconstructions may live on disjoint carriers and
does not supply that identification.

The cheapest exact follow-up is therefore:

```text
clock R private chart + clock 13R private chart
 -> intersect at c_next=h on one physical component
 -> test the resulting carry holonomy.                    (25)
```

If direct positive matching still fails, THM-2639 suggests a genuinely
different test: eliminate the carry between the first two clock-compatibility
polynomials and compute their resultant/Bezout ideal.  Its Gaussian-moment
two-rung mechanism proves that permanent failure of a private *linear*
decoder need not prevent a nonlinear certificate.  That is a next experiment,
not a consequence or dependency here.

## 8. Exact boundary

The gain/loss ledger is:

| item | exact status |
|---|---|
| physical predecessor carry | retained before integration |
| primitive lattice | unchanged, global content `26` |
| private root | exact formula `(13)` |
| unit-root coverage | all twelve roots on all `84` cells |
| fixed-c coverage | all and only `(17)` |
| target action | does not preserve the graph |
| consecutive-clock identity | `c_next=h`, but no same-component second chart |
| holonomy/LRC exit | absent |

In particular, this theorem proves no canonical owner/root transition,
same-component adjacent-clock branch, equivariant carry translation,
holonomy trivialization, row exclusion, or LRC(14).

The exact companion runs as

```bash
python 04-computation/lrc14_predecessor_carry_private_root_atlas_thm2640.py
python -O 04-computation/lrc14_predecessor_carry_private_root_atlas_thm2640.py
```

and both executions must byte-match

```text
05-knowledge/results/lrc14_predecessor_carry_private_root_atlas_thm2640.out.
```

An independent hostile audit replayed both normal and optimized executions,
reproducing all `1,415,232` carry partitions, all `18,398,016` singleton-law
tests (split into `1,415,232` possible on-graph slots and `16,982,784`
off-graph zero slots), the complete root/digit-edge censuses, and the stored
transcript byte for byte.  It independently rederived `(3)`, `(9)`, `(13)`,
the consecutive-digit identity `(24)`, the `Phi_13` argument in `(20c)`, and
the slope-seven physical near miss `(23a)`--`(23c)`.  In particular, it
checked on the canonical speed tuple

```text
(1,14,27,40,53,66,13,13^3,2*13^5)
```

that `v(7 delta/R)` lies on the `1/13` phase grid exactly when `13^5|v`,
so only `c3` survives.  The LF-normalized script/output hashes were reproduced
as `a28b03a5...e71abf` and `b3c1f510...9a940`.  A subsequent scope audit
identified and repaired the total-versus-off-graph count, the overloaded
prefix variable, and the nonuniqueness of the physical lift; none changes
the private-root theorem or its exact computation.

QED.
