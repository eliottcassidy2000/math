---
id: THM-2682
title: "Central-arrival clock trap and three-event dilation nilpotence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  THM-2616 arrival-six raw nonnegative carrier, use the only currently typed
  physical handoff D(x)={13x}.  A positive two-event atom fibre product
  exists, but every three-event product
  E_0(x) intersect E_1(Dx) intersect E_2(D^2x) is empty, uniformly over all
  clocks, sources, rails, sharp states, predecessor/carry labels, and later
  restrictions.  Indeed E_1 and E_2 force z=Dx and Dz into the central digit
  interval I_6.  On that return cylinder,
  z-1/2=(Dz-1/2)/13, so the consecutive shallow seven-clock labels agree;
  the first stored incidence edge requires them to differ.  In event-product
  convention the support nilpotence length is exactly three (a two-factor
  product survives and every three-factor product vanishes).  This closes
  iterated D-chronology on the arrival-six carrier, not alternate carriers,
  boundary handoffs, or LRC(14).
source: root-2026-07-28-central-arrival-clock-trap
depends_on:
  - THM-2616-cross-time-target-future-diagonal-and-principal-action-no-go
  - THM-2680-dilation-reversed-two-edge-clock-fibre-products-and-source-drift-boundary
related:
  - THM-2624-two-clock-root-tomography-and-disjoint-carrier-holotopy-boundary
  - THM-2635-half-tooth-opposite-graph-unit-section-and-reversed-digit-closure
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2670-sharp-graph-clock-incidence-atlas-and-physical-gluing-boundary
script: 04-computation/lrc14_dilation_three_edge_nilpotence_probe.py
output: 05-knowledge/results/lrc14_dilation_three_edge_nilpotence_probe.out
script_sha256: 59475b5b4fe236960252e533f728d2d02a774089fb2e72f482f7d8b7adae2578
output_sha256: f7dd4dc9c9bedcc9a47f7abf53548e73272106296d42fa9ba9d5c38db0c90b66
hash_basis: LF-normalized bytes
---

# THM-2682 -- the central arrival kills the third event

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2680 supplies the first correctly typed physical composition in the
post-THM-2616 clock program.  It proves that a two-event dilation fibre
product can be positive, despite the failure of the earlier formal matrix
products.  Iterating that handoff once more reveals a much simpler and
uniform obstruction: the fixed arrival digit and the shallow seven-clock
cannot make two consecutive nonzero moves.

The proof uses only two half-open intervals and the scale identity `C1=13`.
It occurs before the sharp graph, delayed digits, source transport, or
Bockstein arithmetic.

## 1. Typed three-event object

Let

```text
D(x)={13x}.                                                (1)
```

Following THM-2680, write a chronological clock path as

```text
ell_0 -> ell_1 -> ell_2 -> ell_3,                         (2)
```

with adjacent labels distinct, and write the event at time `i` as

```text
E_i=E^(s_i)_(owner=ell_(i+1), shallow=ell_i;
             alpha_i,j_i,h_i,epsilon_i,kappa_i).          (3)
```

The dilation covariance laws are

```text
owner(E_i)=shallow(E_(i+1)),
h_i=j_(i+1).                                              (4)
```

Thus the first genuine three-event object is

```text
C_3=E_0(x) intersect E_1(Dx) intersect E_2(D^2x).         (5)
```

The sources `s_0,s_1,s_2` and every rail/private label may drift
independently.  This only enlarges the candidate universe, so an obstruction
before those labels is uniform for every later transport rule.

## 2. Arrival-six return cylinder

Every one of THM-2616's `162` retained middle rails is contained in

```text
I_6=[6/13,7/13).                                          (6)
```

Put

```text
z=Dx,                  w=Dz=D^2x.                         (7)
```

The rail factor in `E_1(Dx)` forces `z in I_6`, and the rail factor in
`E_2(D^2x)` forces `w in I_6`.  Since no wrap occurs on `I_6`,

```text
w=13z-6.                                                  (8)
```

Equivalently,

```text
z-1/2=(w-1/2)/13.                                        (9)
```

The exact return cylinder is

```text
I_6 intersect D^(-1)(I_6)=[84/169,85/169),               (10)
```

of mass `1/169`.

## 3. The shallow-clock diagonal trap

The shallow clock in the inherited carrier is the seven-clock of
`{C1*x}` with

```text
C1=13.                                                    (11)
```

Therefore the shallow label `ell_0` of `E_0(x)` is the depth-one phase-clock
label of `z=Dx`, while the shallow label `ell_1` of `E_1(z)` is the
depth-one phase-clock label of `w=Dz`.  This one-index shift is essential;
the companion checks all seven exact shallow-to-next-phase interval
identities.  Owner/shallow covariance in (4) also says that `ell_1` is the
owner label of `E_0`.

Inside `I_6`, the half-open seven-clock partition has only two cells:

```text
clock(z)=3  below 1/2,
clock(z)=4  at or above 1/2.                              (12)
```

Equation (9) preserves the sign relative to `1/2`, including the half-open
boundary convention.  Hence

```text
clock(z)=clock(w),                                        (13)
```

and the exact labelled support on (10) is

```text
(ell_0,ell_1)=(3,3) with mass 1/338,
(ell_0,ell_1)=(4,4) with mass 1/338.                      (14)
```

But the first stored incidence edge in (3) is

```text
shallow(E_0)=ell_0 -> owner(E_0)=ell_1,                  (15)
```

and the clock graph retains only nonzero steps, so `ell_0!=ell_1`.
Equations (14) and (15) are incompatible.  Therefore

```text
C_3=empty                                                 (16)
```

for every legal choice of all remaining labels.

## 4. Exact scope census and sharpness

The exact companion checks the interval proof and then exhausts the complete
label universe inherited at this stage:

```text
7*6^3 = 1,512                    clock quadruples;
1,512*12^3 = 2,612,736           clock/source candidates;
18,740,124                       rail-labelled triples;
1,596,358,722,816                nominal predecessor/carry-atom triples. (17)
```

All are eliminated by the same source-independent arrival/clock gate.  The
last number is explicitly a nominal labelled-atom census, not a count of
distinct geometric components or points.

THM-2680's exact positive two-event numerator certificate is replayed:

```text
126816337986097204341478787325120 > 0.                   (18)
```

Thus one event is not the maximum: a genuine two-factor support product
exists.  Equations (16) and (18) give the exact convention-independent
statement

```text
some two-event D-chain exists;
no three-event D-chain exists.                            (19)
```

We call this **support-product nilpotence length three**, counting event
factors.  If instead one encodes a handoff as a binary transition relation,
then (19) says that relation is nonempty but its relational square is empty;
this prevents an off-by-one ambiguity in the word "index".  Every longer
event product contains a three-event prefix and is therefore empty.

## 5. Mechanism and exact frontier

The obstruction is not a numerical coincidence in a large atom table.  It is
the central-digit fixed-point identity (9): the affine return branch fixes
`1/2` and preserves its two sides, while the seven-clock boundary at `1/2`
records exactly those sides.  The proposed nonzero clock move demands the
information that this return branch cannot change.

Consequently the canonical THM-2616 arrival-six carrier cannot support an
iterated chronology under the dilation handoff (1).  Restricting it further
by a Bockstein unit, rainbow chart, endpoint owner, sharp graph, Perron sheet,
source law, or carry cocycle cannot restore a support already empty in (16).
This closes that entire order of downstream refinements.

The theorem does **not** exclude:

```text
* a carrier using another physical arrival stratum;
* a lawful handoff other than D(x)={13x};
* a boundary correspondence not represented by ordinary pullback;
* a carrier whose edge permits a zero clock step;
* returning to the independent AP, Euler/phase, peel, or rank-box routes. (20)
```

No global source transport, Bockstein unit, endpoint section, scalar cover,
row exclusion, or `LRC(14)` conclusion follows.  Rather, (20) is now the
minimal bypass invoice: a successful chronology proof must change at least
one of arrival stratum, handoff, boundary semantics, or edge grammar.

## 6. Reproduction and hostile audit

Run

```bash
python3 04-computation/lrc14_dilation_three_edge_nilpotence_probe.py
python3 -O 04-computation/lrc14_dilation_three_edge_nilpotence_probe.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_dilation_three_edge_nilpotence_probe.out.
```

The companion uses exact integer interval arithmetic.  It verifies all
`1,050` weighted rail pieces stay in (6), all seven owner/shallow dilation
covariances, all seven shallow-to-next-phase identities, the exact cylinder
(10), the two phase masses (14), the full census (17), and the positive
control (18).  The independent hostile audit rederived the scale indexing
directly from `C1=13`, checked which two event rails force (6), reproduced
the interval proof, and replayed normal and optimized executions.  Both modes
match the stored output and hashes above.

QED.
