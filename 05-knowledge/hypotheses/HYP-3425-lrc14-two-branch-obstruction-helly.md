---
id: HYP-3425
title: LRC14 two-branch obstruction / Helly certificate for the two-adic floor
status: EVIDENCE / exact finite-ruler interval audit; not an LRC14 proof
source: codex-2026-06-28
tangent: T1386
technique: LTI-386
tournament_technique: LTT-286
script: 04-computation/lrc14_two_branch_obstruction_helly_codex_20260628.py
result: 05-knowledge/results/lrc14_two_branch_obstruction_helly_codex_20260628.out
reflection: 07-reflections/lrc14-two-branch-obstruction-helly-codex-20260628.md
related:
  - HYP-3424
  - HYP-3423
  - HYP-3422
  - HYP-3421
  - HYP-3420
  - HYP-3419
  - HYP-3418
  - HYP-3417
  - HYP-3415
  - HYP-3129
  - HYP-2963
  - THM-523
  - OPEN-Q-108
---

# HYP-3425: LRC14 Two-Branch Obstruction / Helly Certificate

## Claim

Downstream of HYP-3424's covering-floor duality transfer, HYP-3422's two-adic
relocation lemma can be sharpened from an overlap slogan to a finite-ruler
obstruction theorem.  Split a covering row as

```text
S = O union 2E,     u = 2t.
```

Then relocation succeeds exactly when

```text
E_safe cap (branch0_good union branch1_good) != empty.
```

HYP-3425 rewrites the failure set:

```text
relocation_good = E_safe minus (B0_odd cap B1_odd)
```

where `B0_odd` is the union of odd near-integer bad intervals for
`||o*u/2|| < 1/14`, and `B1_odd` is the union of odd near-half bad intervals
for `||o*u/2|| > 3/7`.  Thus a true failure must cover every even-safe
component by a two-color odd bad core.  This is a one-dimensional Helly /
interval-piercing target, not a raw resonance or topology claim.

## Exact Readout

Script:

```text
04-computation/lrc14_two_branch_obstruction_helly_codex_20260628.py
```

Stored result:

```text
05-knowledge/results/lrc14_two_branch_obstruction_helly_codex_20260628.out
```

The exact audit covers `62` rows: the HYP-3422 curated rows, the canonical
`{1..11,13,84m}` tower through `m=12`, six two-tail probes, and forty
deterministic primitive covering rows.

```text
positive two-branch good union:       62/62
selected relocation score >= 1/14:    62/62
smallest good-union measure:          1/105 (covering_AP_with_84)
smallest surviving component gap:     1/118692 (random_covering_16)
max nonempty odd-pair obstructions:   56 (random_covering_12)
```

The tight canonical row

```text
S = {1,2,3,4,5,6,7,8,9,10,11,13,84}
```

has

```text
E_safe measure        = 107/245
bad two-color core    = 314/735
good union            = 1/105
surviving components  = 4
component gaps        = 1/588 .. 3/980
selected branch       = 1
selected t            = 2293/3920
selected score        = 59/784
```

So the most visible resonant one-tail row is not merely "safe off-grid"; it has
four explicit surviving finite-ruler windows after the two-color obstruction
core is removed.

## Canonical Tower Signal

For `S_m = {1..11,13,84m}`, the audit through `m=12` keeps positive
two-branch good measure on every row.  The number of even-safe components grows
with `m`, but only a small number of components survive the two-color odd
obstruction.  This is useful: the proof target should control the surviving
component windows, not estimate the whole safe mass too coarsely.

## Proof Target

The next finite lemma should be stated as:

```text
For every primitive covering 13-row S=O union 2E,
E_safe is not contained in B0_odd cap B1_odd.
```

Equivalently, at least one component of `E_safe` has positive length outside
the two-color odd bad core.  A proof can try:

1. A finite-ruler bound on every component of `E_safe`.
2. A Helly/interval-piercing argument showing the near-integer and near-half
   odd bad intervals cannot jointly cover every component.
3. A two-adic descent step on `E` when even speeds remain.
4. HYP-3417/HYP-3419 owner-current labels only for named exceptional packets,
   especially the `2:g2` even-cover shadow.

## Guardrail

Do not replace this target by apex-7 topology, raw resonance transparency, or
a scalar `Rprime` statement too early.  HYP-3424 moves the covering floor into
the dual/positive-negative transfer lane; HYP-3423 says topology is q-uniform
and cannot close the q-specific magnitude/floor branch by itself.  HYP-3421
says resonant speeds are transparent at the full off-grid optimum.  HYP-3425
says what must be proved inside that transparency: the two-color odd
obstruction core must fail to cover the even-safe interval family.

## Tournament Analysis

Tournament vertices are proof obligations, not runners, raw resonant speeds,
or constants.  The observable compares predicate retention, finite exactness,
interval/Helly shape, two-adic induction, owner-current glue, `Rprime` glue,
and failure guardrails.

```text
score_hist={21:1, 55:2, 56:1, 60:2, 61:1}
directed_3cycles=0
hamiltonian_path=
  two_branch_bad_core_identity
  -> component_gap_helly_certificate
  -> two_adic_descent_induction
  -> owner_current_exception_router
  -> canonical_84m_surviving_windows
  -> signed_SPEC_Rprime_floor
  -> raw_resonance_transparency_slogan
```

The leading carrier is the exact identity; the second is the component-gap
certificate.  That is the proof route to push next.
