---
id: HYP-2615
title: LRC(14) signed-mass sequence spine
status: OPEN
source: codex-2026-06-19-S13
depends_on:
  - THM-538
  - HYP-2598
  - HYP-2608
  - HYP-2614
related:
  - HYP-2612
  - HYP-2613
  - MISTAKE-078
  - OPEN-Q-108
---

# HYP-2615 - LRC(14) Signed-Mass Sequence Spine

## Claim

The HYP-2608(a) support-6 tail should be organized around a small set of
integer and fractional sequences, not around a single absolute-volume estimate.
The recurring pattern is:

```text
large absolute mass on boundary/cusp faces
+ exact support floor
+ residue-addressed reciprocal coimage
= tiny signed mass
```

The useful proof object is the coimage of the support-6 relation lattice under
the residue character map modulo `7`, together with finite deletion of
low-height additive walls.  Boundary faces can have huge integer relation
counts, but their signed reciprocal sums are tiny after the seven-sector
kernel is applied.

## Sequence Evidence

Script:

- `04-computation/lrc14_signed_mass_sequence_spine_codex.py`
- output: `05-knowledge/results/lrc14_signed_mass_sequence_spine_codex.out`

The script records five linked spines.

### 1. Support Floor

THM-538 gives the exact-support sequence

```text
support size: 0 1 2 3 4 5 6
live?       : 0 0 0 0 0 0 1
```

for every ambient dimension checked in the proof route.  This is the hard
cancellation: all genuine support `<=5` relation terms vanish.  The remaining
tail is intrinsically six-body.

Guardrail: this is not tested by simply plugging zero coordinates into the raw
ambient `K` function, since that includes boundary/main pieces.  The proof
object is the exact-support coefficient after signed inclusion-exclusion.

### 2. Residue-Constant Decay

For six-support coefficients,

```text
K(n_1,...,n_6,0,...) = C_d(n mod 7)/(n_1...n_6).
```

The normalized residue constants decay with ambient dimension:

```text
d=6..13 max |C_d|:
0.643084862, 0.321542431, 0.168426988, 0.091869266,
0.0519342586, 0.0303106082, 0.0182058608, 0.0112237893
```

Relative to the blunt `64*c1^6 = 7.35714264`, the ratios run from about
`0.0874` at `d=6` down to `0.00153` at `d=13`.

The one-coordinate marginals remain nonzero:

```text
0.384780076, 0.193862704, 0.107446515, 0.0637696096,
0.0397214435, 0.0256286869, 0.0169882801, 0.011508678
```

So cancellation is not free coordinate-marginal cancellation.  It must use the
relation hyperplane `sum e_i n_i = 0`, residue addresses, and finite wall
deletion.

### 3. Cusp Ratio Sequence

Parsing the exact S12 ledger gives the abs/signed ratios:

```text
final support-6 ratio sequence:
29.0667, 217.015, 58.2626, 1118.76, 447.629

one-face boundary ratio sequence:
202.386, 3484.64, 13.3239, 34.8034, 112.732
```

The corresponding one-face boundary relation counts are huge:

```text
5,291,542; 1,825,394; 353,652; 783,346; 996,976.
```

This is the user's "large absolute mass but tiny signed mass" signal in a
form future agents can reuse.  The integer side is the boundary relation
count; the fractional side is the small signed/absolute coimage error.

### 4. Universal-Center Survivor Sequence

HYP-2598 contributes the small-speed denominator sequence

```text
survivor(s) = C(7,s) + C(9,s) - C(5,s)
            = 1, 11, 47, 109, 156, 146, 91, 37, 9, 1, 0, 0, 0, 0.
```

The mixed residual complement is

```text
0, 2, 31, 177, 559, 1141, 1625, 1679, 1278, 714, 286, 78, 13, 1.
```

This separates the fixed denominator-cap skeleton from the colored-resonance
tail.

### 5. Empty-Window Certificate Sequence

HYP-2608 gives the AP-clearing degree sequence

```text
k=8..12 degree: 4, 3, 3, 2, 1.
```

with AP margins

```text
431/24696, 1769/140140, 12937/321048,
29749/917280, 11917/543312.
```

This is the region-side companion to the support floor: as `k` grows, fewer
moments are needed, while the true analytic residual is pushed into the
support-6 tail.

## Category / Number-Theory Reading

The absolute count lives before quotienting: it counts all lattice points on
boundary faces and low-height walls.  The signed mass lives after a coimage:

```text
relation lattice -> residue addresses mod 7 -> signed reciprocal coefficient.
```

In categorical language, the useful map is not the inclusion of all lattice
points but the coimage after killing lower support faces.  In number-theory
language, the coefficient `C_d` is a finite character/Dedekind coefficient
over `(Z/7Z)^6`, and the relation hyperplane supplies the cancellation.  The
prime `7` is the half-denominator gate for `n=14`; the resonance is visible
only after the `2n` problem has collapsed to the seven-sector signed quotient.

## Tournament Analysis

Vertices are sequence families, not runners:

- support floor;
- residue-constant decay;
- cusp signed ratio;
- empty-window degree drop;
- universal-center survivors;
- low-height wall ledger;
- coimage reciprocal tail.

The script's quotient tournament is transitive:

```text
support_floor
> universal_center_survivors
> empty_window_degree_drop
> cusp_signed_ratio
> residue_constant_decay
> low_height_wall_ledger
> coimage_reciprocal_tail.
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
SCC sizes = [1,1,1,1,1,1,1]
Hamiltonian path count = 1
```

This preserves the analytic predicate "support-6 correction below the cap
margin" and the recurrence across runner count / ambient dimension.  It
destroys witness-time geometry, so it is a proof-route quotient, not a witness
constructor.

## Proof Target

HYP-2615 does not prove LRC(14).  It sharpens the live sub-problem:

```text
finite low-height wall ledger
+ residue-addressed signed reciprocal hyperplane tail
+ no-scale cluster quotient
< cap margin.
```

The next useful move is to prove a summation-by-parts / cotangent-Dedekind
bound for the coimage reciprocal tail after deleting the explicit low-height
walls.  The evidence says the absolute cusp mass is a ghost; the signed coimage
is the object to bound.
