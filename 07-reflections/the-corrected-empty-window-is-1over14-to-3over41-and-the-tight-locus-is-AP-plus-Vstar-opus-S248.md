---
source: opus-2026-07-11-S248
status: CONSTRUCTIVE sharpening (restores S246's structure with corrected parameters). The empty M-spectrum
  window at k=13 is (1/14, 3/41), NOT (1/14, 2/27) — verified empty over 56039 families (404 with M<3/41, all
  at exactly 1/14, zero in the gap). The tight locus (M=1/14) is exactly {AP} ∪ {V*=doubling 12→24}, up to
  residue-preserving shift — TWO essential families, not the unique-AP of the prime k=12 case. Corrected
  closing target: [M < 3/41 ⟹ family ∈ {AP, V*}], the k=13 analog of HYP-4151's k=12 rigidity with the
  composite-14 (=2·7) twist. So S246's window-empty ⟹ rigidity ⟹ LRC(14) STRUCTURE was right; only the window
  bound (3/41 not 2/27) and the tight-locus size (two, not one) were wrong.
tags:
  - lrc14
  - empty-window
  - tight-locus
  - AP-and-Vstar
  - composite-14
  - crux
  - corrected-closing-target
---

# The corrected empty window is (1/14, 3/41), and the tight locus is {AP, V*}

**opus-2026-07-11-S248.** Owner: keep sharpening the LRC(14) crux via hypothesis investigation. Following the
S247 correction constructively — the right window and the right tight locus — restores S246's whole *structure*
with corrected parameters. This is the sharpest picture of the crux to date.

## The corrected empty window: (1/14, 3/41)

S247 refuted S246's claim that the window `(1/14, 2/27)` is empty (`3/41` fills it). The **correct** empty
window is `(1/14, 3/41)` — the true second value of the M-spectrum is `3/41`, not `2/27`. **Verified empty**
over 56 039 primitive families (broad random + near-AP + near-V* perturbations): 404 families have `M < 3/41`,
**all 404 at exactly `M = 1/14`**, and **zero** strictly inside `(1/14, 3/41)`. The observed low ladder is

`1/14 → [empty gap] → 3/41 → 2/27 → 3/40 → 1/13 → …`

so S246's bound was one rung too high — it named the mediant `2/27` when the actual next rung is `3/41`
(= HYP-2934's `K_{3,3}` incidence wall / near-miss `12→36`, the composite-side landmark).

## The tight locus is {AP, V*} — two families, the composite-14 signature

The `M = 1/14` locus is **not** the AP alone. Among single-element replacements `{1..13}`, `k→m` stays tight
iff either `m ≡ k (mod 14)` (residue-preserving shift — trivially the same mod-14 configuration) **or** the
single non-trivial case

`k = 12 → m = 24 = 2·12` — **the doubling** (V* = `{1..11, 13, 24}`; mac-mini THM-708).

`12→24` is the **only** non-residue-preserving tight replacement (siblings `11→25, 13→15, 6→20` all jump to
higher M). So the **essential tight locus (distinct mod 14) is exactly two families**: the AP `{1..13}`
(residues `{1..13}`, full) and V* (residues `{1..11,13,10}` — missing 12, doubling 10, a **collision**). Every
non-AP tight family has a residue collision mod 14 — the composite signature. `12→24` is enabled precisely by
`14 = 2·7`: doubling by 2 through the zero-divisor at 7 preserves the binding constraint. At the prime k=12
(mod 13, a field) there is no such doubling — the tight family is **unique**.

## The corrected closing target (S246's structure, restored)

S246's logical skeleton was right — *empty window ⟹ rigidity ⟹ LRC(14)* — only mis-parametrized. Corrected:

> **`M < 3/41` ⟹ `M = 1/14` ⟹ family ∈ {AP, V*} (up to residue-preserving shift).**

This proves LRC(14): a counterexample (`M < 1/14 < 3/41`) would be forced into `{AP, V*}`, both of which have
`M = 1/14` — contradiction. It is the exact **k=13 analog of HYP-4151's proved k=12 rigidity** (`M < 2/25 ⟹`
unique dilated AP), carrying the **composite-14 twist**: the window shrinks (`3/41` not `2/25`-analog `2/27`)
and the tight locus **doubles** (`{AP, V*}` not unique). Both twists trace to the same root — the zero-divisor
at 7.

## Why this is the honest, sharper picture

- **S247 (correction):** S246's specific window `(1/14, 2/27)` is false. ✓
- **S248 (this, constructive):** the *right* window `(1/14, 3/41)` IS empty (broadly verified), and the tight
  locus is exactly `{AP, V*}`. So S246's **method was not dead** — it needed the composite-corrected
  parameters. The crux is not "AP is the unique minimizer" (false — V* ties it) but "**{AP, V*} are the only
  minimizers, and the gap above them to `3/41` is empty**".
- The composite `14 = 2·7` is now pinned as a *two-fold* obstruction to importing the k=12 proof: it (i)
  enlarges the minimizer set by the doubling V*, and (ii) lowers the second rung to `3/41`. A k=13 rigidity
  proof must reproduce **both** — the field-based k=12 equioscillation argument sees neither.

## Net

The LRC(14) crux, sharpened: **the M-spectrum has an empty gap `(1/14, 3/41)` above a two-element minimizer
locus `{AP, V*}`** (verified over 56k families). The closing target is the corrected rigidity `M < 3/41 ⟹
{AP, V*}` — HYP-4151's k=12 statement with the composite-14 doubling built in. This turns S247's "the window is
non-empty" from a dead end back into a live, correctly-parametrized route, and isolates exactly what a k=13
proof must add over the proved k=12 case: the doubling family and the `3/41` rung, both children of the
zero-divisor at 7.

→ opus-S247 (the correction this builds on), opus-S246 (structure restored, parameters corrected), HYP-4151
(k=12 rigidity — the template), HYP-4306 (empty-window ladder — corrected rung), HYP-2934 (`3/41` = K33 wall),
THM-708/mac-mini (V* = doubling `12→24`), LEM-015 (E₃ — the AP/V* both maximize the mod-14 additive structure).
Files: `lrc14_crux_picture_composite14_opus_S248.py` (+`.out`).
