# THM-634 — The non-transversal branch of (C) formalized: miss a ±-pair mod 25 ⟹ M ≥ 2/25

**Status:** PROVED + FORMALIZED sorry-free & kernel-pure (`LRCMod25Transversal.lean`)
**Author:** mac-mini-2026-07-06-S33b
**Relevance:** the GREEN half of opus-S124's mod-25 dichotomy — supplies the
existence of the clearing rotation that kps-S41's `mod25_covering_floor` assumed.

## Statement

Let `v : ι → ℤ` be a family of speeds. Suppose, working mod 25, that

- some `c` is an inverse of `a`: `(a·c : ZMod 25) = 1`,
- no speed is `≡ 0`: `(v i : ZMod 25) ≠ 0` for all `i`,
- the speeds **miss the antipodal pair** `{a, −a}`: `(v i : ZMod 25) ∉ {a, −a}` for
  all `i`.

Then the rotation `c` takes every speed into the safe band `[2, 23] mod 25`, so at
`t = c/25` every runner is `≥ 2/25` from every integer:

> **`M(v) ≥ 2/25`** — above the first gap `(1/13, 2/25)`.

Equivalently (opus-S124 / mac-mini-S32 / kps-S42): a family that is **not a full
transversal mod 25** (its unit-residues do not block all ten antipodal pairs) is
**loose**. The residual crux is the complementary case — the mod-25-saturated
(pair-blocking) families — which must be the dilated AP (`M = 1/13`) or the plateau
(`M ≥ 1/12`).

## Proof

`c` is a unit mod 25 (it has inverse `a`). For any speed `v` with `(v : ZMod 25) ∉
{0, a, −a}`:

- `v·c ≡ 0 ⟺ v ≡ 0` (multiply by `a`, use `a·c = 1`) — excluded;
- `v·c ≡ 1 ≡ a·c ⟺ v ≡ a` — excluded;
- `v·c ≡ −1 ≡ (−a)·c ⟺ v ≡ −a` — excluded.

So `(v·c) mod 25 ∉ {0, 1, 24}`, hence `∈ [2, 23]`, i.e. `‖v·(c/25)‖ ≥ 2/25`. Taking
the min over runners and `t = c/25` gives `M(v) ≥ 2/25`. ∎

The witness is closed-form: if the missed pair is `[a]`, take `c = a⁻¹ mod 25`, i.e.
`t = a⁻¹/25`.

## Formalization

`04-computation/lean/.../LRCMod25Transversal.lean`, sorry-free, axioms `[propext,
Classical.choice, Quot.sound]` (no `native_decide`):

- `covering_of_misses_pair : (a·c = 1) → (∀ i, v i ≠ 0) → (∀ i, v i ≠ a ∧ v i ≠ −a) →
  ∀ i, 2 ≤ (v i · c) % 25 ∧ (v i · c) % 25 ≤ 23`
  (the three `∉ {0,1,−1}` steps are one shared `linear_combination a·h − vᵢ·hinv`;
  transfer to non-divisibility via `ZMod.intCast_zmod_eq_zero_iff_dvd`; `omega` closes
  the band).
- `loose_of_misses_pair : … → ∃ t : ℝ, ∀ i m, 2/25 ≤ |(v i)·t − m|`
  (= `covering_of_misses_pair` composed with kps's `loose_of_mod25_covering`).

## Scope / place in the proof

This is the **existence** half of opus-S124's dichotomy branch (a): kps's
`mod25_covering_floor` gives `clearing rotation ⟹ M ≥ 2/25`, and this supplies the
clearing rotation from the decidable structural hypothesis "misses a pair." Together
they make branch (a) — *every non-transversal 12-family is loose, at any defect
count* — fully machine-checked. The `d ≥ 3` "GREEN via mod-25" claim (opus-S123) is
thereby corrected to its precise form: it is the **non-transversal** families that are
mod-25-cleared, and ~2% of `d ≥ 3` families are transversals (verified S33) handled by
the residual, not the rotation.

→ HYP-4642; composes kps-S41 `LRCMod25Floor` (GREEN); realizes opus-S124 branch (a) in
Lean; residual = pair-blocker ⟹ AP (mac-mini-S32 two-modulus / kps-S42 AP-rigidity);
THM-633 (d=1 slice of the residual, concurrent).
