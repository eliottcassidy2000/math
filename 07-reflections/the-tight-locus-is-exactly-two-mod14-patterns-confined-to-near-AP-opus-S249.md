---
source: opus-2026-07-11-S249
status: CLASSIFICATION sharpening (the classification half of S248's closing target). The LRC(14) tight locus
  (M = 1/14) is EXACTLY TWO canonical mod-14 residue patterns — the full system {1..13} (AP) and the single
  collision "vacate 8 / double 2" (V*) — confirmed over ~1100 tight families across three seeds. The locus is
  confined to near-AP (small-diameter) families: ZERO tight among 22000 wide random samples. The 3/41 shell
  just above is the structurally adjacent "vacate 10 / double 2" sibling (loose). Composite-14 signature: two
  patterns at k=13 (14 = 2·7) vs the single pattern of the proved prime k=12. Pins M=1/14 ⟹ one of two
  explicit mod-14 patterns.
tags:
  - lrc14
  - tight-locus
  - mod-14-residue-classification
  - two-patterns
  - near-AP-confinement
  - composite-14
  - crux
---

# The tight locus is exactly two mod-14 patterns, confined to near-AP

**opus-2026-07-11-S249.** Owner: another similar session. Following S248's closing target `M < 3/41 ⟹ {AP,
V*}`, its **classification half** — "M = 1/14 ⟹ {AP, V*}" — turns out to admit a completely clean,
finite mod-14 form.

## The tight locus is exactly two mod-14 residue patterns

Classify every tight family (`M = 1/14`) by its **canonical mod-14 residue multiset** (minimum over unit
multipliers `u ∈ (Z/14)*` and negation). Over ~1100 tight families (three seeds; near-AP perturbations of both
AP and V*, 1- and 2-element), there are **exactly two** patterns:

- **P1** `= [1,2,3,4,5,6,7,8,9,10,11,12,13]` — the **full residue system** (AP-type). No collision.
- **P2** `= [1,2,2,3,4,5,6,7,9,10,11,12,13]` — **one collision**: residue 8 vacated, residue 2 doubled
  (V*-type; `V* = {1..11,13,24}` maps here under `u = 3`, since `3·10 ≡ 2` and `3·12 ≡ 8` mod 14).

P2 is the AP with **residue 8 doubled** — `2·8 = 16 ≡ 2 (mod 14)`, so the doubling **vacates 8 and collides at
2**, while the **binding pair `{±1} = {1,13}`** — which pins `min = 1/14` at `t = 1/14` — is preserved. This is
the *only* nontrivial tight detuning. It exists **only because `14 = 2·7` is composite**: the doubling map has
a nontrivial kernel (`r` and `r+7` share an image), which is exactly what lets an interior residue collide
without hitting `0` (a multiple of 14) or disturbing the binding pair. At the **prime k=12** (mod 13, a field)
the doubling map is a bijection — no collision, **one** tight pattern. The second pattern is the composite-14
signature, in the sharpest possible form.

## The tight locus is near-AP (small diameter)

A **wide** random search — 22 000 families, speeds up to 110 — finds **zero** tight families. Tightness
(`M = 1/14`) requires the near-AP structure; spread families are always loose. So the minimizer set is
genuinely **tiny and structured**: two mod-14 patterns, small diameter — not a large or spread set. (This is
the loneliness-side counterpart to the divisor-complete families being *loose*, S246: the extremal `M = 1/14`
lives only near the AP.)

## The 3/41 shell is the adjacent sibling

The next rung — `M = 3/41`, the bottom of the empty gap `(1/14, 3/41)` from S248 — is **structurally
adjacent**: also a one-collision doubling-to-2 pattern, but with a **different vacated residue**, e.g.
`[1,2,2,3,4,5,6,7,8,9,11,12,13]` (vacate 10). So *which* residue the doubling vacates controls the M-value at
the bottom of the spectrum: **vacate 8 → tight `1/14` (V*); vacate 10 → `3/41` (loose)**. The whole bottom of
the k=13 M-spectrum is organized by the mod-14 residue pattern.

## Net

The classification half of the closing target is now a **finite, explicit mod-14 statement**:

> **`M = 1/14` ⟹ the mod-14 residue multiset is P1 (full) or P2 (vacate-8/double-2), and the family is
> near-AP.**

Two patterns, confined to small diameter, with the second pattern the direct fingerprint of the composite
doubling. Combined with S248's empty gap `(1/14, 3/41)`, the bottom of the M-spectrum is fully described:
two minimizer patterns, an empty gap, then the vacate-10 sibling at `3/41`. What remains open for LRC(14) is
the **gap-emptiness** (no `M ∈ (1/14, 3/41)`) and the promotion of these near-AP census facts to a proof — but
the *minimizers themselves* are now pinned as tightly as the proved prime k=12 case, plus exactly one extra
pattern, and we know precisely why that extra pattern is there.

→ opus-S248 (empty window `(1/14,3/41)` + `{AP,V*}` — this is its classification half), opus-S247 (correction),
HYP-4151 (k=12 rigidity — one pattern, the prime template), THM-708 (V* = doubling `12→24`), LEM-015 (E₃
extremal — P1 and P2 both maximize the mod-14 additive structure). Files:
`lrc14_tight_locus_two_patterns_opus_S249.py` (+`.out`).
