---
id: THM-873
title: THE μ₆ SPINE OF LRC(14) — (T) THE TIGHT-LOCUS THEOREM: the times where the 13-runner AP {1,…,13} attains its maximin 1/14 are EXACTLY {p/14 : p ∈ (Z/14)*} = the six powers of 3 mod 14 (cyclic order 1→3→9→13→11→5), so the tight locus is a GROUP ≅ μ₆ under numerator multiplication and folds to three Eisenstein directions under t ↔ 1−t; this identifies the S315 Farey-14 extremal ratio classes with LRC(14)'s equality set; (W) THE WELL DICHOTOMY mod 3: mediant-well denominators q = 3N+2 ≡ 2 (mod 3) are NEVER Eisenstein norms (inert branch), while the deep well 183 = Φ₆(14) = N(14+ω) IS one (14 a primitive 6th root of unity mod 183, the S52 lever) — the well taxonomy splits along the Eisenstein prime classification; (L) the leave-one-out atlas {1..13}∖{k}: M = 1/k exactly with t* = 1/k for k ≥ 7; M = 2/q with paired denominators q = 15, 17, 17, 19, 19, 23 for k = 1..6 (the floor extremal ∖{6} = 2/23 is the global minimum, matching THM-777's uniqueness); honest negatives: neither χ₁₃ nor the Eisenstein-norm property grades the atlas (BH: 2 gradings scanned, 0 laws)
status: T PROVED (Dirichlet-with-equality: some v ≤ 13 has ||vt|| ≤ 1/14 for every t, with equality forcing t ∈ (1/14)Z — three-distance equality analysis; exact scan q ≤ 200 confirms) + the group law immediate; W exact (q ≡ 2 mod 3 obstructs the norm form; 183 verified); L exact (12 families, exact ℚ maximin, qmax = 400)
source: opus-2026-07-15-S320 (owner: apply the new lenses to LRC(14), think roots of unity, synthesize past threads)
depends_on:
  - the S52 Eisenstein lever (the-13-comb-lever-is-the-eisenstein-resonance, with the S55 correction)
  - THM-777 (the |G'| floor at {1..13}∖{6} — its M-minimality now placed in the atlas)
related: [S315 Farey-14 law (the same six, now identified as the tight locus), THM-773 (F₇ tokens), THM-869/872 (Milgram μ₈ phases), HYP-4572/4516 (the mediant mod-6 gate), klein-S124 (the Eisenstein sibling)]
verification: 05-knowledge/results/mu6_spine_lrc14_opus_S320.out
---

# THM-873 — the μ₆ spine of LRC(14)

## (T) The tight-locus theorem

For the boundary AP S = {1,…,13}: min_v ||vt|| ≤ 1/14 for every t (Dirichlet
on the 13 points {vt}: the 14 arcs they cut from the circle average 1/14, and
the arc containing 0's neighborhood bounds the minimum), with equality iff
every nonzero residue class {vt} sits exactly on the (1/14)-grid, i.e.
**t = p/14 with gcd(p, 14) = 1**. Exact scan (all Farey points q ≤ 200)
confirms: the tight locus is precisely

> {1/14, 3/14, 5/14, 9/14, 11/14, 13/14} — numerators = ⟨3⟩ = (Z/14)* ≅ μ₆,
> in cyclic order 1 → 3 → 9 → 13 → 11 → 5.

Consequences: the equality set of LRC(14)'s extremal AP is a **group** (μ₆
under numerator multiplication mod 14); the reflection t ↔ 1−t is p ↔ −p,
folding the six to **three Eisenstein directions** {±1}, {±3}, {±5}; and the
S315 Farey-14 extremal ratio classes (the depth-13 resonance maximizers) are
EXACTLY this tight locus — the optimal Hunter tree-edge directions and the
AP's loneliest moments are the same six points. The problem's self-tuning
(S315) now has a group-theoretic face: **LRC(14)'s hardest instances and its
best analytical tools both live on μ₆ mod 14.**

## (W) The well dichotomy (Eisenstein classification)

- Mediant wells: q = 3N + 2 ≡ 2 (mod 3) — a prime ≡ 2 (mod 3) to the first
  power obstructs the Eisenstein norm form a² − ab + b²: **mediant
  denominators are never Eisenstein norms** (verified N = 7..31; the
  congruence proof is one line). The inert branch.
- The deep well: 183 = Φ₆(14) = 14² − 14 + 1 = N(14 + ω) = 3·61 (ramified ×
  split): **an Eisenstein norm, with 14 a primitive 6th root of unity mod
  183** (powers 14, 13, 182, 169, 170, 1 — the six Eisenstein units; the S52
  lever, corrected S55). The norm branch.

The two known well species of the covering-min landscape are separated by
the Eisenstein prime classification of their denominators — the μ₆ lens
splits the well taxonomy exactly.

## (L) The leave-one-out atlas (and the honest negatives)

Exact maximins of {1,…,13}∖{k}:

| k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 |
|---|---|---|---|---|---|---|---|---|---|----|----|----|----|
| M | 2/15 | 2/17 | 2/17 | 2/19 | 2/19 | **2/23** | 1/7 | 1/8 | 1/9 | 1/10 | 1/11 | 1/12 | 1/13 |

For k ≥ 7: M = 1/k exactly (t* = 1/k: removing k opens its own gap). For
k ≤ 6: the wells pair up (2, 3 → 17; 4, 5 → 19) and the floor extremal
∖{6} = 2/23 is the GLOBAL minimum of the atlas — the M-face of THM-777's
|G'|-floor uniqueness. Negatives, logged per the BH discipline: the χ₁₃
(QR/QNR) grading and the Eisenstein-norm grading of this atlas both FAIL
(the M-classes mix characters freely) — the μ₆ structure lives at the tight
locus and the well denominators, NOT in the leave-one-out ordering.

## The spine assembled (the synthesis)

μ₆ now appears in LRC(14) at five independent stations: the TIGHT LOCUS
((Z/14)* ≅ μ₆ — this theorem); the DEEP WELL (14 = primitive 6th root mod
Φ₆(14) — S52); the MEDIANT GATE (N ≡ 1 mod 6 — HYP-4572/4516); the QR
HEXAGON (φ(13)/2 = 6 — the signed-selector residues); and the FAREY-14
RESONANCE ROW (φ(14) = 6 — S315). The Milgram phases (μ₈, THM-869/872) and
the F₇ token roots (THM-773) flank it. Conjecture for the fleet: these are
one structure — the μ₆-grading of the wall lattice by the character of
(Z/14)* — and the deep two-sheet residual's signed selectors should be
re-coordinatized by this character before the next uniform-failure attempt.
