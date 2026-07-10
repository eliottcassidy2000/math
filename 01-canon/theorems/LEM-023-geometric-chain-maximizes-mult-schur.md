---
id: LEM-023
title: The geometric power-chain maximizes multiplicative Schur triples — among k integers ≥ 2, M₃(S) = #{(a,b)∈S×S : a·b∈S} ≤ C(k,2), with equality iff S = {a, a², …, aᵏ} (a ≥ 2); the multiplicative dual of LEM-015, transported through x ↦ log x
status: PROVED (Lean, sorry-free, kernel-pure [propext, Classical.choice, Quot.sound]); `LRCMultRigidity.multCount_eq_choose_iff_geometric`
source: kind-pasteur-2026-07-09-S127 — the multiplicative rank-1 geometric-chain rigidity klein's character frame sees
depends_on:
  - LEM-015                        # the additive twin (E₃ ≤ C(k,2), equality iff dilated interval)
  - LRCSchurRigidity.dilated_of_closedUnderDiff   # reused verbatim on the exponent set
related:
  - THM-682(d)                     # doublings are the sole W₀-carriers; the a=2 case {2,4,…,2ᵏ} is the doubling chain
  - HYP-5835 (klein)               # the multiplicative character frame at prime rulers that resonates on geometric progressions
  - the-final-rung-lives-on-the-diagonal-dyadic-carrier-free-triangles-kps-S127   # the diagonal split placing this in the endgame
---

# LEM-023 — the geometric power-chain maximizes multiplicative Schur triples

## Statement

For a finite set `S` of integers `≥ 2` with `|S| = k`, the **multiplicative Schur count**

  `M₃(S) = #{(a,b) ∈ S × S : a·b ∈ S}`

satisfies `M₃(S) ≤ C(k,2)`, with **equality iff `S` is a geometric power-chain** `{a, a², …, aᵏ}` for
some `a ≥ 2`.  Equivalently `M₃(S) = C(k,2) ⟺ MultClosedUnderQuot S` (every `x < y` in `S` has `x ∣ y`
and `y/x ∈ S`).

This is the exact multiplicative dual of **LEM-015** (`E₃(S) ≤ C(k,2)`, equality iff dilated interval
`c·{1,…,k}`).  Under `x ↦ log_a x` the multiplicative structure of `S` becomes the additive structure of
its exponent set, and the extremal is the exp-image of the additive dilated interval `{1,2,…,k}`.

## Proof (Lean, `LRCMultRigidity.lean`)

- **The count `⟺` quotient-closure** (`multCount_eq_choose_iff_multClosed`): the injection
  `(a,b) ↦ {a, a·b}` sends multiplicative Schur pairs into the 2-subsets of `S`, injectively (since
  `a < a·b` when `b ≥ 2`), and is a bijection exactly when every 2-subset `{x,y}` (`x<y`) is hit, i.e.
  `(x, y/x)` is a pair, i.e. `x ∣ y` and `y/x ∈ S`.  Mirror of `schurCount_eq_choose_iff_closedUnderDiff`.
- **Quotient-closure `⟺` geometric chain** (`geometric_of_multClosed` / `multClosed_of_geometric`):
  every element is a power of `min S` (strong induction: `y > min ⟹ y/min ∈ S`, smaller, a power by IH),
  so the exponent set `E = S.image (log_a)` is well-defined, closed under differences, and contains `1`;
  the **additive rigidity `dilated_of_closedUnderDiff` applied to `E`** gives `E = {1,…,k}`, whence
  `S = {a¹, …, aᵏ}`.  The additive theorem is reused *verbatim* — the two rigidities are one theorem in
  two coordinates. ∎

## Why it matters (the endgame)

The `a = 2` instance `{2, 4, 8, …, 2ᵏ}` is the pure **doubling chain** — a geometric progression of
ratio 2, a rank-1 GAP in the multiplicative group.  By THM-682(d) the **doublings are the sole carriers
of the exact-load `W₀`** in the LRC(14) final rung, and (kps-S127 diagonal split) they are exactly the
diagonal `(a,a)` of the additive `E₃`.  So LEM-023 is the extremal characterization of the 2-adic carrier
that the final rung lives on: max multiplicative-Schur content ⟺ a geometric (2-adic, for `a=2`) chain —
the object klein's multiplicative character frame (HYP-5835) resonates on at prime rulers.  Together with
LEM-015 it completes the **structural/extremal characterization of the residual class on both Freiman
axes** (additive dilated interval + multiplicative geometric chain).

*Verification: `04-computation/lean/TournamentH7/TournamentH7/LRCMultRigidity.lean` (build green, 8476
jobs; axioms = the standard three). Numerics: `lrc14_e3_diagonal_split_kps_S127` confirms the doubling
diagonal; the extremizers `{2,4,8,16}`, `{3,9,27,81}` at k=4.*
