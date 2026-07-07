# The non-transversal branch of the crux is now formalized — and it corrects the "d≥3 GREEN" claim

*mac-mini-2026-07-06-S33b (HYP-4642 / THM-634). Owner: finish off the crux and then
formalize it; pull often and integrate. Verified: `LRCMod25Transversal.lean` (GREEN,
kernel-pure), `lrc_d1_ladder_reconcile_macmini_S33.out`.*

## What got finished

The crux (C) = "(G): no 12-family has `M ∈ (1/13, 2/25)`" splits, along opus-S124's
mod-25 dividing line, into:

- **(a) non-transversal (cleared): `M ≥ 2/25`** — the family misses an antipodal pair
  `{a,−a}` mod 25, so the rotation `c = a⁻¹` clears every speed off `{0,±1}`.
- **(b) full transversal (saturated): `M = 1/13` (AP) or `M ≥ 1/12` (plateau)** — the
  residual near-AP moat.

Branch (a) was "GREEN via kps" in the fleet's telling, but there was a **hole in the
Lean chain**: kps-S41's `mod25_covering_floor` takes a clearing rotation `c` as a
*hypothesis* and concludes `M ≥ 2/25`; nothing produced `c`. This session closes that
hole. `LRCMod25Transversal.lean` (kernel-pure, `[propext, Classical.choice,
Quot.sound]`) proves:

> `loose_of_misses_pair`: if the speeds miss the pair `{a,−a}` mod 25 (and none is
> `≡ 0`), then `∃ t, ∀ i m, 2/25 ≤ ‖vᵢ t − m‖`.

with the explicit `t = a⁻¹/25`. The proof is three lines of one shared
`linear_combination a·h − vᵢ·hinv` (ruling out `vᵢc ∈ {0,1,−1}` from `a·c = 1` and
`vᵢ ∉ {0,a,−a}`), a divisibility transfer, and `omega`. So **branch (a) — every
non-transversal 12-family is loose, at any defect count — is now fully machine-checked.**

## The correction it forces

opus-S123 filed the mod-25 clearance under "**d ≥ 3** GREEN." That mis-locates it. The
mod-25 rotation clears a family **iff it is not a full transversal** — a condition
orthogonal to defect count. Verified this session: among 4,755 structured `d ≥ 3`
families, **95 (~2%) are full transversals** and are *not* cleared by any rotation;
they are loose for other reasons (0 in-gap), but not by kps's lemma. So the honest
partition is **transversal / non-transversal**, not `d < 3 / d ≥ 3` — which is exactly
kps-S43's "the crux is defect-agnostic." The Lean lemma is stated on the correct
hypothesis (misses a pair), so it is immune to the mis-slicing.

## Where the crux stands

Two independent same-session mac-mini runs plus opus and kps have converged on one
picture, now partly formalized:

| branch | statement | status |
|---|---|---|
| (a) not a transversal | `M ≥ 2/25` | **FORMALIZED** (`loose_of_misses_pair` + kps `mod25_covering_floor`) |
| (b) transversal, `d=0` | dilated AP, `M = 1/13` | boundary (not in open gap) |
| (b) transversal, `d=1` | `{1..11,x}`, `M ≥ 2/25` | **FORMALIZED** (THM-633, concurrent) |
| (b) transversal, `d=2` + plateau | `M ≥ 1/12` | residual (opus-S124/S115) |

So the only mathematically-open piece is the **saturated `d ≥ 2` plateau**: a mod-25
pair-blocking family with longest-AP `≤ 10` has `M ≥ 1/12`. Everything else in (C) is
either formalized or a boundary value. Because `1/12 > 2/25`, closing that plateau
closes (C), hence LRC(14) via the proof map.

## Integration note (shared working directory)

The concurrent mac-mini-S33 shares this checkout; its `git add -A` (THM-633 commit)
swept up my then-uncommitted `LRCMod25Transversal.lean`, so both formalizations landed
together — a happy accident of the shared tree, but a reminder to commit narrowly and
often. (Also: THM-633 reused `HYP-4632`, which collides with my S32b two-modulus entry;
flagged for reconcile.)

## Net

- **Non-transversal branch of (C) FORMALIZED** (kernel-pure): miss a pair mod 25 ⟹
  `M ≥ 2/25`, closing the existence hole in kps's mod-25 certificate.
- **"d ≥ 3 GREEN" corrected** to "non-transversal GREEN" (defect-agnostic; ~2% of `d≥3`
  are transversals) — matches kps-S43.
- **Residual reduced** to the saturated `d ≥ 2` plateau (`M ≥ 1/12`); d=0 boundary,
  d=1 formalized (THM-633), non-transversal formalized (THM-634).

## Pointers

- Lean: `LRCMod25Transversal.lean` (`covering_of_misses_pair`, `loose_of_misses_pair`).
  Canon: `THM-634-non-transversal-branch-formalized.md`. Compute:
  `lrc_d1_ladder_reconcile_macmini_S33.py`.
- Composes: kps `LRCMod25Floor` (HYP-4567). Realizes: opus-S124 branch (a) (HYP-4566).
  Extends: my S32 pair-blocking/two-modulus (HYP-4622/4632). Reconciles: opus-S123
  d-strata (HYP-4556), kps-S43 defect-agnostic (HYP-4587). Sibling: THM-633 (d=1).
