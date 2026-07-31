---
id: THM-2744
title: "Relative-present unit repair and root-zero overlap clutch"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the canonical
  LRC(14) source-one deepest-semantic fibre, the relative complement of the
  frozen one-target present packet repairs the missing primitive-unit side.
  Edge-preserving scalar transport sends right root 12 to the forbidden
  right-root-zero label, but the right-root-zero and left-root-one physical
  half-tooth charts overlap openly.  Translation by 7/13^6 and recharting on
  that overlap give exact equality of the two recomputed raw coefficient
  vectors.  Their root-normalized classes differ by the sign -1.  Ten of
  fourteen rails work without an auxiliary cut and all fourteen on their equal-weight
  loci; rail 8 carries a strict E3-to-D^6-Q_(3,{1,2}) witness and the identical
  vector (0,a,a,a,a,a,a), which is a private unit at roots 12 and 1.  This is
  not a global action, endpoint current, row exclusion, or LRC(14) proof.
source: root/relative-present-root-zero-overlap-clutch-2026-07-28
audit: coordinate-first-audit/root-zero-overlap-2026-07-28; semantic-present-wall-hostile-audit-2026-07-28 (independent counterexample, direct complement and polynomial-gcd unit reconstruction, count/overlap/scope audit, and normal/-O/hash replay: ACCEPT)
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2640-predecessor-carry-private-root-atlas-and-target-action-clutching-no-go
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2672-slope-seven-carry-nerve-exact-eleven-simplex-and-root-zero-cap
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
  - THM-2707-full-physical-lift-fibre-common-simplex-and-packet-scc
  - THM-2712-semantic-following-congruence-lock-and-address-coboundary-descent
  - THM-2720-unshifted-deepest-source-present-packet-global-disjointness
  - THM-2742-full-two-target-present-sheet-deepest-source-semantic-current
related:
  - THM-2635-half-tooth-opposite-graph-unit-section-and-reversed-digit-closure
  - THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary
  - THM-2727-fixed-rail-target-deck-frobenius-realization-no-go
  - THM-2749-fully-marked-root-zero-clutch-and-target-character-profile
script: 04-computation/lrc14_root_zero_overlap_clutch_20260728.py
output: 05-knowledge/results/lrc14_root_zero_overlap_clutch_20260728.out
script_sha256: e10fa7c9a5a238461ef422ea314dc334f7e65ec1787cf65d4e4bea12b96aefb8
output_sha256: ba9d0a67dfede0b64cf97ff55af7e86c9bb46c962c669d33015df7e574e8e91e
secondary_script: 04-computation/lrc14_relative_present_semantic_lift_probe_20260728.py
secondary_output: 05-knowledge/results/lrc14_relative_present_semantic_lift_probe_20260728.out
secondary_script_sha256: f16754bd38ae0dfa0d7d91cc404b4447dbf359635101aa7b4223363f8064352f
secondary_output_sha256: 861e920b945aafde3bc24c6089ba630763035f0919d738b77a0525b91ba6fa74
hash_basis: LF-normalized bytes
---

# THM-2744 -- relative present has an open root-zero chart clutch

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The one-target present wall of THM-2720 is a projection loss, and THM-2742
repairs its missing target coordinate.  Rebuilding the physical sidecars
creates a second apparent wall: an edge-preserving lift reaches private-root
label zero.  That label is forbidden, but its physical half-tooth is not
disjoint from the adjacent lawful chart.  The overlap supplies a genuine
partial clutch and corrects the false implication recorded in MISTAKE-310.

## 1. Relative-present coefficient and endpoint bank

Keep the canonical typed row and source-one fibre of THM-2707.  Broaden its
fixed rail to all fourteen source-one rails, retain delayed sector zero,
`h=6`, `kappa=1`, and clocks `1,...,6`, and require the semantic record

```text
E3 -> D^6 -> Q_(3,{1,2}).                              (1)
```

There are exactly `12,848` such fork points.  Their address residues are

```text
residue 7: 6404,                 residue 8: 6444.       (2)
```

Every point has one source-one rail and one strict private half-tooth and is
outside all thirteen frozen one-target present labels.

For one rail and each of the seven clocks, let `Present_(ell,7)` denote the
single frozen one-target packet with label `7`, and set

```text
V_old  = L_(rho,ell)(rail intersect H_rho intersect Present_(ell,7)),
V_free = L_(rho,ell)(rail intersect H_rho),
V_rel  = L_(rho,ell)(rail intersect H_rho
                     intersect Present_(ell,7)^c).      (3)
```

Exact delayed-carry linearity, before determinants, gives

```text
V_free=V_old+V_rel.                                    (4)
```

This complement is the complement of the **fixed label `7`**, not the
complement of the union of all thirteen one-target labels.  Separately, the
`12,848` semantic midpoints above happen to lie outside that entire union.
The independent audit constructed `[0,T)\Present_(ell,7)` directly as its gap
union and applied the nonnegative functional, without using subtraction; all
`14*2*7=196` coefficients agreed with `(3)--(4)`.

The raw gcd of each complete coefficient bank is
`2,122,120=26*81,620`; unit tests remain in THM-2640's inherited
content-`26` lattice.  The exact unit-rail table is

```text
                 residue 7                 residue 8
V_old            all 14                    none
V_free           all except rail 13        all except rail 3
V_rel            all 14                    all except rail 3. (5)
```

Thus `V_rel` leaves `6404+6178=12,582` unit endpoints.  It has
`79,127,824` directed cross-residue scalar pairs and `2,166,112` directed
reverse-clock pairs.  Every endpoint retains its inherited open cylinder and
has between nine and eleven lawful nonzero second-target labels `t`; none has
`t=0` at its midpoint.  The full-target assertion for all `12,582` endpoints
is midpoint-pointwise only; stability of a common `(s,t)` on the whole open
cylinder is proved below only for the displayed rail-`8` pair.

Here “unit endpoint” means that the endpoint inherits the unit of its
rail/residue **aggregate coefficient vector**.  It is not a pointwise
transition amplitude.  Likewise `79,127,824` counts all directed raw
cross-address translations; it does not assert a common rail, reverse clock,
target label, coefficient equality, or current for each pair.  The separate
`2,166,112` count imposes only the reverse-clock equality condition.

## 2. The sharp edge-preserving no-go

The two endpoint classes have the uniform private data

```text
residue 7: carry 12, right edge, root 12,
residue 8: carry  6, left  edge, root  1.               (6)
```

A forward scalar numerator has `2k=1 mod13`; its edge-preserving root action
sends right root `12` to right root `0`.  The reverse has `2k=-1` and sends
left root `1` to left root `0`.  Hence the unrestricted labelled
edge-preserving **nonzero-root label relation** is empty.  The precursor
transcript's line `private_root_transport=EMPTY` is retained only as this
fixed-edge label statement; MISTAKE-310 forbids reading it as physical-support
emptiness.

It is false to infer that the physical support is empty.  In the THM-2640
convention the open half-tooth charts are

```text
H^L_r=(14r-13,14r)/182,        H^R_r=(14r,14r+13)/182. (7)
```

They satisfy, uniformly in `r mod13`,

```text
H^R_r intersect H^L_(r+1)=(14r+1,14r+13)/182.          (8)
```

This width-`12/182` overlap is the coordinate omitted by the label-level
argument.

## 3. The canonical physical rechart

Put

```text
R=13^6,                         tau=7/R.                (9)
```

Since `c3=2*13^5`, translation by `tau` moves the deep `182`-coordinate by
exactly `14`.  On the present fibre its source and target coordinates are

```text
z =169+56447/742586,
z'=  1+56447/742586.                                  (10)
```

Consequently

```text
H^R_12 intersect T_tau^(-1)H^L_1=(169,181)/182,
T_tau((169,181)/182)=(1,13)/182=H^R_0 intersect H^L_1. (11)
```

The circle-coordinate margin is

```text
56447/100360982066072,                                 (12)
```

exactly `56,447` times the inherited open-cylinder radius.  Thus the whole
common cylinder, not just one midpoint, survives the rechart.  The physical
point is unchanged at the target; only its chart changes from
`(right,root0)` to `(left,root1)`.

## 4. A semantic and full-target witness

The adjacent address pair

```text
n=6715:
q =47850889647341/100360982066072,

n=6716:
q'=47851035194197/100360982066072=q+7/R              (13)
```

has record `(1)` at both endpoints.  Both points lie on rail `8`, with
metadata `(1,4,12)`, in the same weighted segment; their reverse clock is
`(shallow,owner)=(1,4)` and their delayed coordinate agrees.

Their lawful target banks are

```text
source s=(0,1,2,3,8,9,10,11,12),
source t=(3,4,5,6,7,8,9,10,11),

target s=(0,1,2,3,8,9,10,11,12),
target t=(3,4,5,6,7,8,9,10,11,12).                    (14)
```

In particular `(s,t)=(0,3)` is lawful and stable on both endpoint
cylinders.  The independent hostile margin is `1,541,619` times the inherited
open-cylinder radius.  This whole-cylinder full-target statement is asserted
only for the pair `(13)`.  The extra target label `t=12` is retained as a real
sidecar difference; the two banks are not silently identified.

## 5. Coefficient intertwining before the determinant

For a weighted rail profile `rho_j` and present complement `P_ell^c`, pull
the target back to the source coordinate and restrict to

```text
support(rho_j) intersect T_tau^(-1)support(rho_j)
 intersect P_ell^c intersect T_tau^(-1)P_ell^c
 intersect H^R_12 intersect T_tau^(-1)H^L_1.           (15)
```

Apply the delayed-carry functional to `(15)` with source carry `12`.
Translate the restricted target profile and independently apply it with
target carry `6`.  The delayed phase is unchanged because

```text
{R(x+7/R)}={Rx}.                                      (16)
```

On rail `8`, source and target weights agree on the whole overlap.  The two
independently computed seven-clock vectors are

```text
A^-=A^+=(0,a,a,a,a,a,a),
a=5359949020444386606638400.                          (17)
```

After division by `26`, root normalization, and reduction modulo `13`, their
classes in `F_13[z]/(Phi_7)` are the nonzero constants

```text
source root 12: 11,                 target root 1: 2. (18)
```

Since `2=-11 mod13`, the raw equality `(17)` is **not** equality of the two
root-normalized rows: those classes differ by `-1`.  Both determinants are
nevertheless nonzero.  Calling `(17)` a normalized intertwiner would require
an additional typed root gauge `-1`, which is not constructed here.  In
characteristic zero, if
`P(X)=a(X+...+X^6)`, then

```text
P(zeta_7^k)=-a,                    k=1,...,6,           (19)
```

so all six nontrivial clock characters survive.  These are clock characters,
not THM-2742's deepest-target characters.

## 6. Exact uniformity and hostile controls

The overlap law `(8)` holds for all thirteen roots.  Without any auxiliary
weight cut, source and target vectors agree on exactly the ten rails

```text
0,4,5,6,7,8,9,10,11,12.                              (20)
```

Restricting each rail to its exact equal-weight locus gives a positive vector
of shape `(0,a_j,...,a_j)` and private units at roots `12` and `1` on all
fourteen rails.  The four added cuts are lawful derived subcarriers, not a
canonical global action.  Rail `8` is stronger: no such cut is needed and it
contains `(13)`.  Semantic/common-full-target/whole-cylinder verification is
proved only for that rail-`8` pair, not uniformly over either rail census.

Four controls are sharp:

- on the unrestricted half-tooths, the rail-`8` source and target relative
  rows have support sizes seven and six and do not intertwine;
- the one-unit outer flanks `(0,1)/182` and `(13,14)/182` admit no opposite-
  edge rechart;
- the physical identity changes chart labels, so it does not contradict the
  edge-preserving no-go of Section 2.
- the root-normalized endpoint classes are negatives, so the raw equality
  does not silently supply a normalized-row action.

The proved chain is therefore

```text
semantic scalar lift + relative-present unit endpoint
 + open adjacent-tooth overlap
 -> partial physical chart isomorphism
 -> equal recomputed raw coefficient vector
 -> private units at roots 12 and 1.                  (21)
```

It does not give a global `C_13` action, endpoint transition amplitude,
target/physical diagonal, row exclusion, or LRC(14).  The residual ledger is
still `165`.  THM-2749 now supplies the formerly next decisive refinement:
it places a common `(s,t)` and the full semantic factor in `(1)` **inside**
`(15)`, proves a uniform all-fourteen unit bank, and computes a target-label
character window.  The next debt is the fixed-triangle THM-2334 transplant
that retains an exact relation address and endpoint allocation.

## 7. Exact reproduction and audit

Run

```bash
python3 04-computation/lrc14_relative_present_semantic_lift_probe_20260728.py
python3 -O 04-computation/lrc14_relative_present_semantic_lift_probe_20260728.py
python3 04-computation/lrc14_root_zero_overlap_clutch_20260728.py
python3 -O 04-computation/lrc14_root_zero_overlap_clutch_20260728.py
```

against the two declared outputs.  Normal and optimized modes byte-match
their stored transcripts, all four declared hashes agree, and neither script
contains truth-bearing `assert` nodes.  The overlap companion pins all six
direct dependencies before reconstructing the thirteen chart overlaps,
fourteen rail vectors, strict witness, target labels, cylinder margins,
determinants, and hostile controls.

MISTAKE-312 records the cross-platform evidence repair: those six dependency
pins declare LF images and are now checked after CRLF-to-LF normalization.
The repair changes no transcript value or mathematical conclusion; normal and
optimized Windows replays again match the stored output byte-for-byte.

The independent audit began by trying to prove the former zero-root no-go and
found the open overlap `(8)` instead.  It independently reconstructed the
forward and reverse translations, caught and repaired an overstatement that
the two `t`-label banks were identical, recomputed every rail before unit
testing, replayed both modes and stored bytes, and verified the exact hashes.
The second hostile audit directly reconstructed the fixed-label-`7`
complement on all `196` rail/residue/clock cells and obtained exactly the same
coefficients.  It independently used Euclidean polynomial gcd with `Phi_7`
over `F_13`, rather than the determinant helper, to reproduce the unit table;
it also rederived every endpoint, address, reverse-clock, overlap, margin, and
normalized-sign count.  The two audits found no remaining hypothesis gap.

## 8. Boundary ledger

```text
PROVED HERE:       fixed-label-7 complement identity and content;
                   aggregate unit table and endpoint/address/clock counts;
                   fixed-edge zero-label law, not physical emptiness;
                   all-root adjacent-chart overlap;
                   one strict rail-8 semantic/full-target cylinder pair;
                   equality of its restricted raw vectors and two unit tests;
                   natural 10/14 and derived equal-weight 14/14 rail censuses.

DESTROYED:         edge label and root label under rechart;
                   equality of unrestricted half-tooth rows;
                   normalized-row equality without a typed -1 gauge.

SUPPLIED LATER:    THM-2749 places semantic/full-target factors inside the
                   integrated vector on a two-sided common section.

NOT CONSTRUCTED:   a global clutch functor or C13 action;
                   a fixed THM-2334 triangle/address transplant;
                   pointwise endpoint amplitude or spectral endpoint pair;
                   endpoint current, row exclusion, or LRC(14).           (22)
```

QED.
