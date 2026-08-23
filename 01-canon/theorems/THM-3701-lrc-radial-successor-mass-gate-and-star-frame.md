---
id: THM-3701
title: "LRC radial successor-mass gate, six-chart compression, and star frame"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On every THM-3670 owner/word stratum, any one of the twelve radial
  inequalities A_k != S_0 or B_l != S_0 already makes one lawful successor
  table nonconstant.  THM-2365 then forces positive H-drift, a nonzero deep
  colour and nonzero target at some exact 91-unit frequency triangle.  Thus
  the thirteen-mass failure locus is the one-dimensional all-equal scalar
  line, not THM-3670's two-dimensional forward-mask locus.  Six charts along
  any derangement are coordinate-minimal, the radial map is the K_(1,12)
  star incidence frame with spectrum 0,1^11,13, and a symmetrized sixty-row
  three-site frame has spectrum 0,21^5,29^5,45,65.  Exact two- and three-value
  variance tariffs improve THM-3672's strongest control-row drift floor to
  2.435186817e-12.  A strict typed bare-owner divisibility chain realizes
  constant successor mass on every packet, proving that the scalar-cover or
  delayed-word input remains load-bearing.  The all-equal mass line is not
  THM-3679's scalar address diagonal; no LRC(14) conclusion is claimed.
source: codex-lrc14-20260822 / THM-2365 radial composition
audit: >
  PASS AFTER ONE TYPING CAUTION -- two independent session audits reproduced
  the rank-12 gate, legal-packet counts, sharp variance constants, six-chart
  minimum, control invoices and bare-owner hostile.  One audit noted that the
  B radial difference appears in THM-3674's displayed three-site notation as
  a=0, kappa=8 with a factor two; direct normalization gives the same genuine
  two-site kappa=2 tariff.  The sixty-row forward/reverse spectrum is therefore
  labeled an algebraic successor-statistic frame rather than a source-role
  interchange.  Both normal and optimized companions match the pinned output.
depends_on:
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
  - THM-3666-lrc-owner-pivot-dual-pair-swap-twist-basis
  - THM-3670-lrc-successor-mass-transfer-and-thirteen-number-gate
  - THM-3672-lrc-successor-mass-all-packet-positive-control
  - THM-3674-sharp-successor-variance-drift-and-target-energy-tariff
related:
  - THM-3665-lrc-support-minimal-three-twist-target-detector
  - THM-3679-lrc-p-adic-scalar-seam-and-total-speed-checksum
script: 04-computation/lrc_radial_successor_mass_star_frame_thm3701.py
output: 05-knowledge/results/lrc_radial_successor_mass_star_frame_thm3701.out
script_sha256: fadcdfa6ff2c4c1c4f48434ad6268c440fdacc8c3ea5fdb18f9de6a41ec69f0e
output_sha256: 7b0043197502a6b38795ff4cab502f85acc4009bd5e04f4157bf8dc95b383d87
hash_basis: raw LF bytes
---

# THM-3701 -- one radial step shrinks the successor hostile to a scalar line

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  THM-3670 used one support-minimal three-site mask on
each legal graft chart.  Its rank-11 system is correct, but its
two-dimensional nullspace is only the nullspace of that chosen mask bank.
THM-2365 needs merely a nonconstant successor marginal.  Comparing either
shifted value directly with the centre supplies a strictly stronger gate.

## 1. The twelve radial tests

Fix exactly the THM-3670 scope: one owner/word stratum whose successor
blocker is one of the two target blockers.  For distinct grafts `k,l` its
lawful successor table is

```text
S_(k,l)(s,t)=sum_r H_(k,l)(r,s,t)
            =integral_T F_(s,t)(x)(2-d(13cx)) dx.   (1)
```

The thirteen shared rational masses are

```text
S_0=S_(k,l)(0,0),
A_k=S_(k,l)(1,0),
B_l=S_(k,l)(0,1).                                  (2)
```

They are independent of the unused graft and omitted label for the reasons
proved in THM-3670.  For every fixed `k`, choose any `l!=k` and an omitted
label outside `{k,l}`.  If

```text
A_k!=S_0,                                          (3)
```

then the two entries `S_(k,l)(1,0)` and `S_(k,l)(0,0)` of that one lawful
table differ.  Hence its successor marginal is nonconstant.  THM-2365 gives

```text
D_H>0
 -> some exact gcd(m,91)=1 and X have
    a nonzero deep-colour, nonzero-target fibre.    (4)
```

The identical argument with `S_(k,l)(0,1)` proves (4) whenever

```text
B_l!=S_0.                                          (5)
```

There are five choices of the other graft and four choices of the omitted
label.  Thus one nonzero radial difference fires at least twenty of the 120
owner packets.

The simultaneous failure locus of all twelve radial screens is exactly

```text
A_0=...=A_5=B_0=...=B_5=S_0.                       (6)
```

It is a one-dimensional scalar line in the thirteen-dimensional mass space.
On a positive selected stratum it is the corresponding positive ray.  This
strictly improves THM-3670's two-dimensional forward-mask locus

```text
A_0=...=A_5=A,
B_0=...=B_5=B,
S_0+A=2B.                                          (7)
```

Equation (6) says only that these thirteen sampled masses coincide.  It does
not say that any of the 169 successor tables is constant; an unsampled entry
may still force (4).

## 2. Six charts are enough and coordinate-minimal

One legal chart `(k,l)` exposes the triple

```text
(S_0,A_k,B_l).                                     (8)
```

Choose a derangement `pi` of the six graft labels and inspect only

```text
(k,pi(k)),                   k=0,...,5.             (9)
```

Every `A_k` and every `B_l` occurs exactly once.  If all six triples in (9)
are constant, (6) follows; otherwise one chart has nonconstant successor
mass and fires (4).  There are

```text
!6=265                                             (10)
```

such chart banks.  Six is minimal for this coordinate-edge strategy: one
chart exposes only one `A` label and one `B` label, while all six labels of
each type must be compared with `S_0`.  Any one of the four legal omitted
labels can realize each chart.

## 3. The radial star frame

Write

```text
x=(S_0,A_0,...,A_5,B_0,...,B_5) in C^13             (11)
```

and let `R:C^13->C^12` be the twelve radial differences

```text
R x=(A_0-S_0,...,A_5-S_0,B_0-S_0,...,B_5-S_0).     (12)
```

This is the oriented incidence matrix of the star `K_(1,12)`.  Its Gram
spectrum is

```text
spec(R^*R)={0^1,1^11,13^1}.                        (13)
```

Indeed, the scalar vector spans the kernel; the eleven leaf vectors of sum
zero have eigenvalue one; and `(-12,1,...,1)` has eigenvalue thirteen.  If
`L_sc=C(1,...,1)`, (13) gives the sharp frame inequality

```text
dist(x,L_sc)^2
 <=sum_k|A_k-S_0|^2+sum_l|B_l-S_0|^2
 <=13 dist(x,L_sc)^2.                              (14)
```

The middle expression combines entries exposed by different packet charts;
it is a cross-packet diagnostic, not the variance of one successor table.
Section 4 passes it to one legal chart before applying THM-3674.

For comparison, retain THM-3670's thirty forward statistics and add their
algebraic reverse orientations:

```text
f_(k,l)=S_0+A_k-2B_l,
g_(k,l)=S_0-2A_k+B_l,             k!=l.             (15)
```

The sixty-row frame has rank twelve and exact spectrum

```text
{0^1,21^5,29^5,45^1,65^1}.                        (16)
```

The standard `S_6` blocks are `[[25,4],[4,25]]`, with eigenvalues 21 and
29.  On the three-dimensional constant-coordinate block, the transverse
vectors

```text
(0,-1^6,1^6),             (-12,1^6,1^6)            (17)
```

have eigenvalues 45 and 65.  Consequently

```text
21 dist(x,L_sc)^2
 <=sum_(k!=l)(|f_(k,l)|^2+|g_(k,l)|^2)
 <=65 dist(x,L_sc)^2.                              (18)
```

This is an algebraic successor-statistic frame.  It does not exchange the
selected source or deepest target role.

## 4. Sharp two- and three-value energy tariffs

Let `S:F_13^2->C` be one lawful successor table and fix two of its entries
whose difference is `delta`.  Among all 169-entry tables with those values,
the normalized variance satisfies sharply

```text
Var(S)>=|delta|^2/(2*13^2).                         (19)
```

This is the two-point variance identity: their least squared deviation is
attained when the other entries equal their midpoint.  THM-3674 therefore
gives

```text
D_H>=|delta|^2/52728,
E_dt>=|delta|^2/685464.                             (20)
```

For the three entries in (8), put

```text
Q_(k,l)=|S_0-A_k|^2+|A_k-B_l|^2+|B_l-S_0|^2.      (21)
```

The least total squared deviation of these three fixed values is
`Q_(k,l)/3`, attained at their mean.  Hence

```text
Var(S_(k,l))>=Q_(k,l)/(3*13^2),

D_H>=Q_(k,l)/(3*13^3*12),
E_dt>=Q_(k,l)/(3*13^4*12).                         (22)
```

These constants are sharp at the finite-table level.  THM-3674's
nonnegative lawful extremizers show that positivity and the diagonal zero do
not improve the variance-to-drift tariff uniformly.

Put `R_star=||Rx||^2`, the middle term of (14).  On any derangement bank
from (9),

```text
sum_k Q_(k,pi(k))
 =R_star+sum_k|A_k-B_(pi(k))|^2
 >=R_star.                                         (22a)
```

Therefore one of those six legal charts satisfies

```text
Q_(k,pi(k))>=R_star/6>=dist(x,L_sc)^2/6.
```

Composing with (22), some one packet in every minimal bank obeys

```text
D_H>=R_star/474552
    >=dist(x,L_sc)^2/474552,

E_dt>=R_star/6169176
     >=dist(x,L_sc)^2/6169176.                     (22b)
```

The factor `1/6` is sharp: take every `A_k` and `B_l` equal to the same
offset from `S_0`.  Thus (22b) is a typed cross-packet-to-one-packet
transfer, not a sum of unrelated packet energies.

Once (4) extracts a nonconstant exact complex target-twist profile, a bank
of the two coordinate edge masks

```text
delta_0-delta_alpha,              delta_0-delta_beta              (23)
```

also detects it.  Their common kernel consists of constants because
`alpha,beta` generate the target plane.  THM-3666 realizes each surviving
edge as one physical blocker/graft pair shift.  Thus some two-site physical
pair-swap current difference is nonzero, with the direction chosen after
the target is known.  This does not contradict THM-3665: no *single* fixed
two-site mask detects every nonzero target character.

## 5. The typed control gets a stronger invoice

On THM-3672's typed non-cover row, all twelve radial differences are nonzero.
The largest is

```text
A_5-S_0=-249094179/931187992280.                    (24)
```

Equation (20) gives

```text
D_H>=20682636670561347
     /15240344288762454200781478400
    =1.357097732090724...*10^-12,

E_dt>=20682636670561347
      /198124475753911904610159219200
     =1.043921332377480...*10^-13.                  (25)
```

This is `1.258799250...` times THM-3672's previous strongest three-site
invoice.  The joint chart `(k,l)=(5,0)` is stronger still.  Equation (22)
gives

```text
D_H>=360926324521774869394441
     /148212992112761067316289860457462400
    =2.435186817139354...*10^-12,

E_dt>=360926324521774869394441
      /1926768897465893875111768185947011200
     =1.873220628568734...*10^-13.                  (26)
```

The drift floor in (26) is `2.258799250...` times the previous invoice.
These remain control-row quantities, not uniform covering-row floors.

## 6. A strict all-packet bare-owner hostile

The scalar line (6) cannot be excluded using only local present factors,
correct packet rank, distinct unit residues and a strict blocker-depth
profile.  Put

```text
(H,q_1,...,q_5)=(1,7,49,343,2401,16807),

c_1=13*7*q_5=1529437,
c_2=7*13*c_1=139178767,
c_3=7*13^3*c_2=2140430257693.                      (27)
```

The six unit residues modulo thirteen are

```text
(1,7,10,5,9,11),                                  (28)
```

and the blocker valuations are `(1,2,5)`.  Every adjacent quotient in the
displayed divisibility chain is a multiple of seven.  THM-2367 (19) therefore
peels every arbitrarily phase-shifted lower `d`, `g`, or `u_2` factor at its
Haar mean, for all 120 omitted/graft packets with the fixed owner `c_1` and
targets `c_2,c_3`.  Every resulting bare-owner tensor is

```text
H_E(r,s,t)=A J(r-t),

A=(5/7)(6/7)^5(1/7)(6/7)=233280/5764801,

J(0)=0,
J(1)=J(12)=1/13,
J(u)=1/7 otherwise.                                (29)
```

Consequently every packet and every `(s,t)` has the same successor mass

```text
sum_r H_E(r,s,t)=33592320/524596891.                (30)
```

This is genuine simultaneous all-packet recirculation, not six unrelated
pair hostiles.  It is deliberately not a scalar cover: its exact uncovered
mass is

```text
(5/7)(6/7)^8=8398080/40353607>0.                   (31)
```

There is also no delayed word in (29).  Thus a valid exclusion of (6) must
use the global scalar cover, a delayed-word correction, or another
nonlacunarity sidecar; purely local packet geometry cannot suffice.

## 7. Exact remaining interface and scope

The covering-row obligation inherited from THM-3670 is now strictly
sharper:

```text
on a genuine scalar cover with an eligible low target,
exclude S_0=A_0=...=A_5=B_0=...=B_5,               (32)

or find any unsampled successor-table variation.                  (33)
```

The scalar **mass** line in (32) consists of thirteen equal rational
integrals.  It is not THM-3679's scalar **address** diagonal `c(1,...,1)`.
The latter is killed at a row-dependent 13-adic depth by a total-speed
checksum; no address-to-mass intertwiner currently identifies these two
objects.

The conclusion (4) selects its own exact `(X,m)` and nonzero target.  It does
not retain THM-2327's preselected triangle, impose all-nine-coordinate
`91`-unit support, preserve visible height or terminal phase, exclude (32)
on a genuine cover, or prove LRC(14).

The dependency-free companion pins THM-3672's exact ledger, verifies both
Gram spectra and eigenbases, exhausts 2,627 sparse frame controls, counts all
265 derangement banks, checks the exact tariffs, and reproduces the hostile
arithmetic.  Normal and optimized runs are byte-identical to the stored
transcript.

```bash
python -B 04-computation/lrc_radial_successor_mass_star_frame_thm3701.py
python -B -O 04-computation/lrc_radial_successor_mass_star_frame_thm3701.py
```

**QED.**
