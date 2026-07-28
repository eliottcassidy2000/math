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
  product survives and every three-factor product vanishes).  A general
  p=2q-1 calculation explains the trap: a same-arrival return has
  off-diagonal nearest-clock support exactly when its arrival digit is even
  and noncentral.  In the inherited LRC tensor the only arrival labels with
  rail support are 0,6,12, so only endpoints 0,12 remain as phase-level
  evasion labels.  This label intersection is not a same-x endpoint chain;
  alternate physical carriers, boundary handoffs, and LRC(14) remain open.
source: root-2026-07-28-central-arrival-clock-trap
depends_on:
  - THM-2584-b-word-depth-five-absolute-deep-root-tensor
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
secondary_script: 04-computation/lrc_central_arrival_clock_return_classification.py
secondary_output: 05-knowledge/results/lrc_central_arrival_clock_return_classification.out
secondary_script_sha256: 1d794246547b539dee725bc947deec41ed2d6138c3448b948eb0e480ca2c42e5
secondary_output_sha256: d762f0cc55066f3ed225daf1e7bdd0535bfce1f9a1536cf68614cad7a411583c
tertiary_script: 04-computation/lrc14_alternate_arrival_rail_phase_filter.py
tertiary_output: 05-knowledge/results/lrc14_alternate_arrival_rail_phase_filter.out
tertiary_script_sha256: c910fa1bc086eba65320c374e1effb22ec46adfa3cf79b93a74bebec992b5831
tertiary_output_sha256: b16575f69accb88a1536d9d93539725f86ef6754514350faa81841391beed542
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

## 6. Why the central digit is exceptional for every `p=2q-1`

The clock trap is not special to the numbers `13` and `7`.  Let `q>=2`, put

```text
p=2q-1,
I_m=[m/p,(m+1)/p),                    0<=m<p,
c_q(x)=floor(q{x}+1/2) mod q.                            (21)
```

The internal boundaries of the nearest `q`-clock are

```text
b_r=(2r+1)/(2q)=(2r+1)/(p+1).                           (22)
```

Since

```text
floor(p b_r)=2r,                                        (23)
```

the clock boundaries occur in exactly the even `p`-digits.  On `I_m`, the
same-arrival return condition `D(z)={pz} in I_m` cuts out the single interval

```text
J_m=[m(p+1)/p^2,(m(p+1)+1)/p^2),                        (24)
```

of mass `1/p^2`.  If `m` is odd, `I_m` contains no clock boundary, hence the
pair `(c_q(z),c_q(Dz))` is diagonal throughout `J_m`.

Let `m=2r`.  Apart from the central case `m=q-1`, the short interval `J_m`
lies wholly in one clock cell while its image traverses `I_m` and crosses
`b_r`.  Pulling the two lengths on either side of `b_r` back through this
slope-`p` branch gives

```text
m<q-1:
  (r,r)       (p-m)/(2q p^2),
  (r,r+1)     (m+1)/(2q p^2);

m>q-1:
  (r+1,r)     (p-m)/(2q p^2),
  (r+1,r+1)   (m+1)/(2q p^2),                           (25)
```

with labels modulo `q`.  All displayed masses are positive.  Therefore

```text
same-arrival off-diagonal support
  iff m is even and m!=q-1.                              (26)
```

If `q` is odd, the central digit `m=q-1` is even and its boundary is exactly
`1/2`.  Both `z` and `Dz` cross that boundary, but

```text
z-1/2=(Dz-1/2)/p,                                       (27)
```

so the two halves remain diagonal, each with mass `1/(2p^2)`.  It is the
unique boundary-straddling digit with diagonal return.  If `q` is even, the
central digit is odd and lies in one clock cell; its return is again diagonal,
but there is no boundary-straddling claim.  Thus central diagonality holds for
every `q`, whereas the stronger fixed-boundary description uses odd `q`.

For `p=13,q=7`, (26) gives the fixed-arrival evasion labels

```text
{0,2,4,8,10,12};                                        (28)
```

the nearest to six are `4,8`.  Allowing the two arrival digits to vary makes
`5->6` and `7->6` the nearest phase-level evasions, each of off-diagonal mass
`1/338`.  These are statements about the clock interval alone.  They do not
assert that an inherited positive rail exists on those labels.

The central inverse branch reverses the same diagonal relation.  The raw map
`rho o D`, with `rho(x)={-x}`, looks off-diagonal for `q=7`, but that
comparison forgets the reflected current clock.  With its typed label
`c_q(rho z)`, its positive-length support is diagonal again.  This statement
concerns exactly the forward branch, central inverse branch, and typed
reflected branch just named; it is not a quantification over every conceivable
transporter.  Null endpoints may change under reflection because of the
half-open convention.

## 7. The inherited rail tensor leaves only endpoint labels

The exact THM-2584 tensor answers the next, separate question: which arrival
labels in (28) even have upstream rail support?  Rebuilding all

```text
13 arrivals * 2 inherited deep labels * 12 sources * 7 clocks
  =2,184 cells                                            (29)
```

from the nonnegative step profile gives

```text
(arrival,deep) positive-cell counts:
(0,0):81,  (6,0):81,  (6,12):81,  (12,12):81;            (30)
all other pairs:0.
```

No deep label outside `{0,12}` has positive mass anywhere in the full tensor,
and the rebuilt arrival-six bank agrees piecewise with all `162` rails used by
THM-2616.  Hence the exact rail **arrival-label inventory** is

```text
{0,6,12}.                                                (31)
```

Intersecting (31) with (28), as label sets only, leaves

```text
{0,12}.                                                  (32)
```

The endpoint banks project to two `81`-element subsets of the `84` static
`(source,clock)` labels, with overlap `78` and union `84`.  Their private
labels are

```text
m=0 private:  (6,1),(6,2),(6,3),
m=12 private: (7,4),(7,5),(7,6).                         (33)
```

Canonical half-open reflection, equal to literal reflection away from null
endpoints, acts by

```text
(m,t,s,ell,[a,b),w)
  ->(12-m,12-t,13-s,-ell,[T-b,T-a),w)                    (34)
```

and matches all `81` positive endpoint cell profiles, comprising `537`
constituent intervals on each leg.  The large totals in the companion are
tensor numerators over the common THM-2584 denominator
`6939029398456584394954868880`, not unnormalized measures.

As a useful hostile marginal, restrict the phase atlas to arrival labels
`(0,6,12)`.  In units of `1/2366`, its off-diagonal mass matrix is

```text
      next 0  6 12
now 0      1 14 14
    6     14  0 14
   12     14 14  1.                                     (35)
```

Thus every projected phase arrow except `6->6` looks available.  Equations
(29)--(35) deliberately do **not** multiply a positive rail marginal by a
positive phase marginal: the two witnesses may be different points.  In
particular, the `78`-cell overlap in (33) is only a projected static-label
overlap; the endpoint rails themselves occupy disjoint base-time digit
intervals.  The next lawful object is their actual pulled-back same-`x` rail
intersection with typed owner/shallow clocks.

## 8. Reproduction and hostile audit

Run

```bash
python3 04-computation/lrc14_dilation_three_edge_nilpotence_probe.py
python3 -O 04-computation/lrc14_dilation_three_edge_nilpotence_probe.py
python3 04-computation/lrc_central_arrival_clock_return_classification.py
python3 -O 04-computation/lrc_central_arrival_clock_return_classification.py
python3 04-computation/lrc14_alternate_arrival_rail_phase_filter.py
python3 -O 04-computation/lrc14_alternate_arrival_rail_phase_filter.py
```

Each normal/optimized pair must byte-match its corresponding transcript:

```text
05-knowledge/results/lrc14_dilation_three_edge_nilpotence_probe.out,
05-knowledge/results/lrc_central_arrival_clock_return_classification.out,
05-knowledge/results/lrc14_alternate_arrival_rail_phase_filter.out.
```

The companion uses exact integer interval arithmetic.  It verifies all
`1,050` weighted rail pieces stay in (6), all seven owner/shallow dilation
covariances, all seven shallow-to-next-phase identities, the exact cylinder
(10), the two phase masses (14), the full census (17), and the positive
control (18).  The independent hostile audit rederived the scale indexing
directly from `C1=13`, checked which two event rails force (6), reproduced
the interval proof, and replayed normal and optimized executions.  The general
classification is proved by (22)--(27); its exact `q=2..64` sweep is a hostile
implementation control, not the source of the all-`q` quantifier.  A second
independent audit reconstructed the THM-2584 rail census, all `537` reflected
endpoint intervals, and the phase atlas, and caught the numerator/measure and
literal/a.e. reflection distinctions now made explicit above.  All modes
match the stored outputs and hashes above.

QED.
