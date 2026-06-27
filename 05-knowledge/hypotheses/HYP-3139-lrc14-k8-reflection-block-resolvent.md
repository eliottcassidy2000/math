---
id: HYP-3139
title: LRC14 k=8 reflection-block resolvent splits the hard bounded-core node into an inner shell page, center coupling, antisymmetric oscillation, and boundary leakage
status: EVIDENCE / exact finite block scout and proof-target refinement; not a proof
source: codex-2026-06-27-S273
tangent: T1204
technique: LTI-265
tournament_technique: LTT-163
script: 04-computation/lrc14_k8_reflection_block_resolvent_codex_s273.py
result: 05-knowledge/results/lrc14_k8_reflection_block_resolvent_codex_s273.out
extends:
  - HYP-3132
  - HYP-3122
  - HYP-3085
related:
  - HYP-3133
  - HYP-3129
  - HYP-3128
  - HYP-3131
  - THM-577
  - OPEN-Q-108
---

# HYP-3139: LRC14 k=8 Reflection-Block Resolvent

## Claim

The HYP-3132 De Moivre/biquadratic reduction of the hard `k=8` bounded-core
row is visible directly in the exact pairwise sector co-emptiness matrix.  After
removing the pinned sector `0`, sectors `1..5` split under reflection
`s->6-s` into:

```text
inner reflected 2x2 shell page on (1+5, 2+4)
+ fixed-center coupling through sector 3
+ antisymmetric 2x2 oscillation on (1-5, 2-4)
+ boundary-sector-6 leakage vector.
```

The two roots of the HYP-3132 biquadratic are exactly the two inner reflected
shells.  Thus the remaining `k=8` proof target can be made finite: prove the
inner 2x2 shell bound, add center and boundary ceilings, then attach the known
`phi4` sign for the `S3/S4` correction.

## Exact Evidence

The scout
[`04-computation/lrc14_k8_reflection_block_resolvent_codex_s273.py`](/Users/e/Documents/GitHub/math/04-computation/lrc14_k8_reflection_block_resolvent_codex_s273.py)
computes exact rational pairwise co-emptiness blocks for consecutive cores
`k=8,9,10`.  For `consec_8` it finds:

```text
core_reflection_exact=True
S2_core_1..5=874/735
S2_boundary_with_sector6=442/735
S2_total=188/105
core_fraction_of_S2=0.664134

reflection_symmetric_3x3_block_on_(1+5,2+4,3):
  [11/28, 5/21, 32/245]
  [5/21, 272/735, 5/49]
  [64/245, 10/49, 117/490]

inner_shell_2x2_page_on_(1+5,2+4):
  [11/28, 5/21]
  [5/21, 272/735]

center_shell_coupling=[32/245, 5/49, 64/245, 10/49]
boundary_sector6_vector=[149/1470, 311/2940, 117/980, 37/294, 73/490]
```

The shell polynomial enumeration around center `3` gives exactly one reflected
four-root fold avoiding the boundary shell:

```text
shells=(1, 2) roots_t=(1, 2, 4, 5) y_poly=y**2 - 5*y + 4 discr=9
```

The two other reflected folds use boundary shell `3`; they remain solvable but
are not the bounded-core `k=8` node.

## Why This Matters

HYP-3085 localized the covering-moment bound to a low-order moment functional
led by pairwise `S2`, with a reflection-symmetric `3x3` Perron block and a
`-9S3+6S4` correction still to control.  HYP-3122 identified `k=8` as the
unique negative `kappa4` / `phi4` row.  HYP-3132 reframed that row as a solvable
De Moivre biquadratic resolvent.  HYP-3139 makes the bridge exact at matrix
level:

- HYP-3132's biquadratic is the inner shell page, not the whole matrix.
- The Perron mass lives in the full symmetric `3x3` block.
- The antisymmetric block is a nonmaximizing oscillation.
- Sector `6` is boundary leakage and must remain a sidecar.

This is a sharper proof obligation than "bound a quartic cumulant" and safer
than "the quartic is solvable."  It says precisely what a quotient may forget
and what it must retain.

## Tournament Analysis

Vertices are proof carriers, not runners or arcs:

```text
exact_reflection_3x3_core_block
hyp3132_inner_biquadratic_fold
hyp3129_spec_resonance_certificate
hyp3133_a000568_middle_shadow
worpitzky_path_moment_vocabulary
antisymmetric_2x2_nonmax_block
boundary_sector6_leak_vector
raw_numeric_spectrum
```

The pairwise observable is majority retention over status preservation, route
preservation, exact arithmetic, folding of the `k=8` obligation, SPEC-floor
connection, A000568 connection, formalization readiness, and scalar guardrail.
The tie Hamiltonian path is the carrier order above.

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
selected_path=exact_reflection_3x3_core_block
  -> hyp3132_inner_biquadratic_fold
  -> hyp3129_spec_resonance_certificate
  -> boundary_sector6_leak_vector
  -> hyp3133_a000568_middle_shadow
  -> antisymmetric_2x2_nonmax_block
  -> worpitzky_path_moment_vocabulary
  -> raw_numeric_spectrum
```

## Guardrails

- This is not a proof of LRC14.  It is an exact reduction of the current hard
  row into smaller finite obligations.
- The inner 2x2 page alone does not preserve the LRC predicate; it destroys
  center coupling and sector-6 boundary leakage unless those are retained.
- A000568/HYP-3133 remains a middle-shadow diagnostic for SPEC rows, not a
  zero-free engine.
- Worpitzky/Eulerian moments are useful expansion vocabulary, not the
  load-bearing certificate here.
- Raw floating eigenvalues are telemetry; the proof target is exact rational
  matrix inequality plus the HYP-3122 `phi4` sign.

## Next Tasks

1. Prove an exact upper bound for the inner reflected 2x2 shell page in the
   `k=8` bounded-core row.
2. Bound the fixed-center coupling and sector-6 leakage without losing the
   HYP-3132 inner fold.
3. Combine the resulting `S2` ceiling with the known `-9S3+6S4` / `phi4`
   sign from HYP-3122.
4. Translate the exact rational block inequalities into the HYP-3107/THM-577
   Lean-facing coverage-extremality packet.
5. Use HYP-3133's A000568 shadow only as a stratifier for the HYP-3129 finite
   SPEC constant chase, especially when a row tries to mimic the inner-shell
   quotient without preserving center/boundary data.
