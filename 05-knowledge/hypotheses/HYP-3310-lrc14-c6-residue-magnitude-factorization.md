---
id: HYP-3310
title: LRC14 C6 residue-magnitude factorization
status: SYNTHESIS / exact CRT-C6 scout plus proof-route split; not an LRC14 proof
source: user prompt plus codex-2026-06-28 integration of S256b, HYP-3265, HYP-3266, and HYP-3300
tangent: T1360
technique: LTI-360
tournament_technique: LTT-260
script: 04-computation/lrc14_c6_residue_magnitude_factorization_codex_20260628.py
result: 05-knowledge/results/lrc14_c6_residue_magnitude_factorization_codex_20260628.out
reflection: 07-reflections/lrc14-c6-residue-magnitude-factorization-codex-20260628.md
related:
  - HYP-3310
  - HYP-3300
  - HYP-3266
  - HYP-3265
  - HYP-3259
  - HYP-3258
  - HYP-3257
  - HYP-3256
  - HYP-3255
  - HYP-3254
  - HYP-3253
  - HYP-3250
  - HYP-3248
  - HYP-3246
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-2909
  - THM-523
  - OPEN-Q-108
---

# HYP-3310: LRC14 C6 Residue-Magnitude Factorization

## Claim

The AP/Goddyn-Wong six-contact picture should be read through a two-layer
factorization:

```text
14 = 2 * 7
Z/14 residue data = 7-adic residue skeleton x 2-adic parity/magnitude gate
```

The binding runners are exactly the units

```text
(Z/14)* = {1,3,5,9,11,13}.
```

They form the full cyclotomic `C6` Galois skeleton.  Quotienting by antipodal
conjugation `+/-1` gives the `C3` action on the three binding slots:

```text
(1,13) -> (3,11) -> (5,9) -> (1,13)        under multiplication by 3 mod 14.
```

The covering runners are the nonunits:

```text
even cover residues = {2,4,6,8,10,12}
ramified apex cover = {7}
```

The apex `7` is not a binding runner.  It is the ramified prime in the
covering layer: CRT label `(1 mod 2, 0 mod 7)`, chi_7 value `0`, and no unit
inverse.

This joins the user's field cue with the current contact graph:

```text
Gal(Q(zeta_7)/Q) ~= (Z/7)* ~= (Z/14)* ~= C6 = C2 x C3.
```

Here `C2` is conjugation/antipodal pairing, the `C3` quotient acts on binding
slots, and the order-3 subgroup fixes the quadratic field `Q(sqrt(-7))`.
Thus the `C3` slot action and `Q(sqrt(-7))` are complementary projections of
the same cyclotomic package, not two independent explanations.

## Exact Scout Readout

The script records:

```text
units_mod14 = (1, 3, 5, 9, 11, 13)
nonunits_mod14 = (2, 4, 6, 7, 8, 10, 12)
even_cover = (2, 4, 6, 8, 10, 12)
apex7 = 7
C6 generator 3 powers = (1, 3, 9, 13, 11, 5)
C3 binder orbit = ((1, 13), (3, 11), (5, 9))
```

At the six AP unit touch points `a/14`, the binding pair is
`{+a^{-1}, -a^{-1}} mod 14`:

```text
1/14, 13/14 -> (1,13)
3/14, 11/14 -> (5,9)
5/14,  9/14 -> (3,11)
```

Modulo `7`, each binding pair has one quadratic-residue side and one
nonresidue side:

```text
(1,13) -> (1,6) chi=(+,-)
(3,11) -> (3,4) chi=(-,+)
(5,9)  -> (5,2) chi=(-,+)
```

This is the precise residue skeleton behind the AP safety function touching
`1/14` at exactly six points in three mirror-symmetric antipodal pairs about
`t=1/2`.

## The 12->24 Hinge

The Goddyn-Wong row replaces `12` by `24`.  The exact hinge ledger is:

```text
12 -> 24
12 mod 14 = 12, 24 mod 14 = 10
CRT: (0,5) -> (0,3)
v2: 2 -> 3
odd part: 3 -> 3
unit-contact killer: false
```

So the hinge is genuinely a 2-adic doubling in the integer/magnitude layer,
but it is not a residue-preserving mod-14 move.  It stays in the even covering
branch and raises `v2` by one while changing the 7-residue.  This is exactly
why the proof packet must retain a magnitude sidecar: the residue skeleton
organizes the binders, but the census distinguishing AP, GW, loose `12->26`,
and positive near-misses is magnitude-level.

## Proof Route

The creative proof route suggested by this factorization is:

```text
1. Prove one binding-pair lemma, such as the HYP-2909 pair.
2. Transport it by C3 to all three binding slots.
3. Split the nonunit covering layer into even cover residues and apex-7.
4. Prove a covering-floor/flex theorem on that layer.
5. Classify the 2-adic magnitude hinge, with 12->24 as the live equality move.
6. Glue the layers using the HYP-3300 observability/Morse packet.
```

The resulting theorem target is not "C3 proves LRC14."  It is:

```text
unit C6/C3 skeleton proves the binding/contact fragment
+ covering branch proves floor/off-grid positivity
+ magnitude hinge theorem isolates the only equality flex
+ observability/Morse glue prevents illegal quotient forgetting.
```

This also explains why `Q(sqrt(-7))` is useful but nonterminal.  It organizes
the `7`-adic residue and chi_7 splitting; it does not by itself see the
2-adic magnitude layer where `12->24`, `12->26`, and `12->36` separate.

## Assumption Challenge

The explored vertex sets were:

```text
runners
gaps
fixed circle sections
section boundaries
wall-crossing events
residue/CRT fibers
cover arcs
Fourier/cyclotomic modes
matroid contact circuits
proof obligations / sidecar columns
```

The chosen Tournament Analysis vertices are proof obligations and sidecar
columns.  The quotient preserves the LRC predicate only as:

```text
unit-contact certificate
+ covering-floor/off-grid witness route
+ 2-adic magnitude/hinge sidecar
+ apex ramification status
+ finite chamber or named residual exit.
```

It destroys exact heights, off-grid safety profile, endpoint owners, and
covering-component data unless those are retained.

## Tournament Analysis

Pairwise observable:

```text
retained LRC predicate coordinates minus destroyed sidecars
```

Switch/gauge:

```text
higher weighted retained payload; ties by fewer destroyed sidecars and priority
```

Exact fingerprint:

```text
vertices = 9
score_hist = {-26:1, 13:1, 21:1, 28:1, 33:1, 41:1, 47:1, 57:1, 127:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1,1,1,1,1,1]
hamiltonian_path_count = 1
priority_hamiltonian_path =
  observability_morse_glue
  -> c6_unit_group
  -> seven_adic_residue_skeleton
  -> two_adic_magnitude_layer
  -> jacobsthal_doubling_hinge_12_24
  -> ramified_apex7_cover
  -> c3_binding_slot_orbit
  -> c2_qr_nqr_conjugation
  -> raw_runner_partition
```

The path is a controlled-forgetting warning.  The `C6` skeleton is the best
local binding carrier after the global observability/Morse glue, but it must
be paired with the 2-adic magnitude layer and the ramified-apex covering
branch before it can support a full LRC14 proof route.

## Next Pull

Turn this into a finite theorem interface:

```text
Binding theorem:
  one complement pair + C3 transport gives all six unit contacts.

Covering theorem:
  even-cover and apex-7 rows have strict off-grid floor unless they are in the
  named AP/GW equality chamber.

Magnitude hinge theorem:
  the only integer equality flex in the covering manifold is 12->24, while
  other same-layer moves open strict mass or route to known floor packets.

Observability theorem:
  no quotient may forget whether a row changed residue, v2 magnitude, apex
  ramification, unit-contact graph, endpoint owner, or off-grid floor route.
```

This is the natural successor to HYP-3265's contact graph and HYP-3300's
observability/Morse proof-angle scaffold.

## Rebase Integration: HYP-3266

The incoming HYP-3266 formal/analytic obligation ledger makes this packet more
targeted.  HYP-3310 should feed three named obligations:

```text
O15 full tight-locus rigidity:
  C6/C3 unit-contact transport plus blind magnitude/covering sidecars.

O12 Part A / off-grid bulk survivor positivity:
  covering-floor route after the unit skeleton is killed or stops being global.

O16 Qsqrt(-7) signed-floor reorganization:
  7-adic chi_7 organization, explicitly guarded by the 2-adic magnitude layer.
```

Thus HYP-3310 is not a replacement for the obligation ledger.  It is the
coordinate dictionary that says which residue/magnitude/ramification fields
must be present before O15 or O16 can be used in a proof and which floor route
must handle the O12 analytic side.
