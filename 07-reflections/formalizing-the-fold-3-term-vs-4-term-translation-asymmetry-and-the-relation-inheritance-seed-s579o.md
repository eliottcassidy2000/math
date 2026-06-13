---
source: oracle-2026-06-03-S579o
status: formalization (Lean 4 / Mathlib v4.30.0) of the 3-term/4-term algebra; axiom-clean; extension ideas surfaced by the process
tags:
  - lonely-runner
  - lean
  - formalization
  - three-term
  - four-term
  - fold
  - resonance
---

# Formalizing the Fold: 3-Term vs 4-Term, the Translation Asymmetry, and the Relation-Inheritance Seed

Formalized the S578o findings (3-term folds carry the hardness; 4-term energy is the
translation-invariant shadow) in Lean 4 / Mathlib v4.30.0, extending
`TournamentH7/LonelyRunner.lean`, and used the formalization process to surface the right
*general* statement.

## What got formalized (axiom-clean)

Two new sections in `LonelyRunner.lean`:

**`FoldStructure`** — the 3-term / 4-term algebra:
- `fold_position`: a fold `v_c = v_a + v_b` makes positions add, `v_c·t = v_a·t + v_b·t`.
- `fold_triangle` (**Lemma B's mechanism**): `|v_c t − (mₐ+m_b)| ≤ |v_a t − mₐ| + |v_b t − m_b|`
  — the fold's distance is *subadditive* in its summands' distances.
- `fold_far_needs_summands`: hence a fold-runner cannot be made far on its own —
  `θ ≤ dist(v_a t) + dist(v_b t)` whenever `v_c t` is `θ`-far.
- `four_term_translation`: `a+b=c+d ↔ (a+s)+(b+s)=(c+s)+(d+s)` — **4-term is
  translation-invariant** (proved by `omega`).
- `three_term_translation_shift` / `three_term_not_translation_invariant`: a fold shifted by
  `s≠0` breaks — `v_c+s ≠ (v_a+s)+(v_b+s)`. **3-term is translation-sensitive.** Together
  these are the formal core of the S578o smoking gun: translating a set hides its folds while
  preserving its additive energy.

**`Relations`** — the extension the process pointed to:
- `relation_inherits` (**the resonance seed, S550**): *any* integer relation `∑ cᵢ vᵢ = 0`
  is inherited by the weighted positions, `∑ cᵢ (vᵢ t) = 0`. The 3-term fold (`(1,1,−1)`)
  and the 4-term coincidence (`(1,1,−1,−1)`) are just the support-3 and support-4 cases —
  one lemma subsumes both.
- `pinch_pair_sum` (**S559o/S560o**): at `t=m/C` with `C=v_a+v_b≠0`, `v_a t + v_b t = m` —
  the pair-sum *is* the pinch denominator; a 4-term coincidence shares one pinch clock.

All report only `[propext, Classical.choice, Quot.sound]` — no `sorry`, no project axioms.

## How the formalization process inspired the extension

Writing `fold_triangle` made the structure unmistakable: the fold inequality is just the
triangle inequality after `add_mul`. That immediately generalizes — **the only thing special
about a fold is that it is a length-3 integer relation**, and the same `∑ cᵢ vᵢ = 0 ⟹ ∑ cᵢ
(vᵢ t) = 0` holds for *any* support. So `relation_inherits` is the honest general statement;
3-term and 4-term are corollaries. This is exactly the resonance-energy decomposition (S550)
seen from the syntax: a relation of support `s` and coefficient vector `c` constrains the
positions to lie on a codimension-1 subtorus, and its "weight" `∏|ĝ(cᵢ)|` is the energy
contribution. **Formalizing the 3-term case forced the realization that the right object is
the relation lattice, not individual triples.**

A second prompt from the syntax: `four_term_translation` is `omega`-trivial *because* it only
sees differences, whereas `three_term_*` needs the absolute value `s`. The Lean proofs make
the translation asymmetry a one-line structural fact rather than a computational observation
— and point at the open Lemma A: a discrepancy bound must be a statement about the
*absolute* relation `∑ cᵢ vᵢ = 0` with the actual `vᵢ` (translation-sensitive), not about the
difference-only energy.

## Extension targets (next formal steps)

1. **General relation distance bound**: `dist(∑ cᵢ vᵢ · t) ≤ ∑ |cᵢ| dist(vᵢ t)` (the
   `Finset` version of `fold_triangle`) — `Finset.abs_sum_le_sum_abs` + `relation_inherits`.
   This is the formal handle on "a relation forces dependent positions," the core of any
   Lemma-A discrepancy argument.
2. **Fold ⟹ not-independently-clearable, quantified**: from `fold_far_needs_summands`, derive
   that a *low* fold (small `v_a, v_b`) at the n-clock pins `dist(v_c t)` — the formal version
   of "only the AP's low folds make it tight" (S578o handoff).
3. **The pinch ⟹ sieve bridge in Lean**: combine `pinch_pair_sum` with the existing
   `sieve_frac` to show the pinch time `t=1/(v_a+v_b)` is lonely iff `(v_a+v_b)∤` every speed
   — formalizing S555o (rational pinch = sieve).

## Verdict / next
- The 3-term/4-term algebra, the translation asymmetry, and the relation-inheritance seed are
  now **axiom-clean Lean**. The process showed the general object is the **relation lattice**
  (`relation_inherits`), with 3-/4-term as support-3/4 slices, unifying the S578o, S550, and
  S559o threads in one formal place.
- Next: the general distance bound (1), the quantified low-fold pinning (2), and the formal
  pinch=sieve bridge (3).

## Artifacts
```
04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean   (sections FoldStructure, Relations)
```
Related: S578o (hardness = 3-term folds), S550 (resonance energy), S559o/S560o (summand graph
/ pinch denominator), S549o (the existing Lean LRC sieve formalization), THM-369.
