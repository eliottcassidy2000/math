---
source: claudebox-2026-06-03-S587
status: REFLECTION + RESULT — the recursive fractal of translated APs; Riesz-product
  factorization; a scale-invariant loneliness gap (the digit gap is a fixed point).
tags: [LRC, fractal, self-similar, translated-AP, Riesz-product, lacunary, additive-energy,
  scale-invariance, sum-free, HYP-2120, HYP-2125]
---

# The loneliness gap is a fixed point of the fractal of translated APs

**Prompt (user):** consider the recursive fractal concept of translated APs.

S585 found the flip: translate an AP above its own diameter and every 3-term relation dies while
the additive energy survives. The natural question the word *recursive* asks: do it again, inside
itself. Use a sum-free translated-AP as the **digit set** of a base-`b` numeral system with
`b > 2·max`, so no two digits ever carry. The set `S_d` of all `d`-digit numbers is then a
self-similar Cantor-like set — `|digits|^d` speeds — and "no carry" turns every additive question
into a digit-wise one.

## What factorizes (rep theory) and what tensors (Haskell)

Because there are no carries, the character polynomial *factors over scales*:
`p_{S_d}(x) = Π_{i<d} D(x^{b^i})` with `D` the single-digit polynomial — a **Riesz product**, a
lacunary spectrum. The functional reading is the same fact as a fold: `S_d` is the base-`b`
evaluation of `replicateM d digits`, the `d`-fold Cartesian power. And the additive invariants
inherit the product: additive energy is *exactly multiplicative*, `E(S_d) = E(digits)^d` (verified
to the integer: 19, 361, 6859, 130321). Energy is a tensor power; 3-term relations are digit-wise,
so sum-free digits stay sum-free at every scale, and the 3-term count of a non-sum-free digit set
grows self-similarly (`{1,2,3}` gives `(3^d+1)/2`). Even the difference set is a fractal —
`|S−S| = |D−D|^d`, dimension `log 5 / log 3 ≈ 1.465`.

## The surprise: the gap is a fixed point

I expected the gap to drift. It doesn't. `G(S_d)` sits on `G(digits)` to within a percent for every
`d` — the digit gap is an **attractor of the recursion**. The reason is the lacunarity that the
no-carry base forces: at the witness time `t*` that makes the digits lonely, the finer scales `b^i
t*` spin so fast (geometric gaps) that the fine-scale runners are essentially independent of the
coarse witness — they land away from 0 on their own, by Lemma A at each scale, so the coarse witness
survives one more level with only a small correction. Dilation-invariance of `G` says the same thing
from the other side: each scale is a dilated copy, and `G` cannot tell them apart, so the coarsest
one — the digits — sets the value.

That gives a clean engineering fact for the conjecture: since `δ = 1/(|S_d|+1)` collapses while
`G(S_d)` holds at `G(digits)`, the margin `G − δ` *grows* with depth. **The fractal of translated
APs manufactures arbitrarily large lonely-runner-easy sets with a fixed, uniform gap.** The safe
(Lemma A) regime is not just nonempty — it has exponentially large, explicitly self-similar members.

## Where the hardness actually lives

The mirror image is the lesson. The conjecture's tight cases — `{1,…,k}` — are exactly the sets
that are *maximally non-fractal*: dense, 3-term-closed, all sums landing back inside. The fractal
construction takes that same `{1,…,k}` as a digit set and *dilutes* its tightness away: only the
`d=1` base is tight; every deeper level is safe, because the set outgrows its danger. Hardness is a
boundary phenomenon — between the lacunary/fractal sets where the scales decouple (safe) and the
dense self-summing sets where they don't (hard). The translated-AP flip was the first crack of light
through that boundary; its fractal is a whole family living on the safe side.

## The transcending pattern

Twice now the structure has been *multiplicative across scales*: the relation lattice's theta sum
(S585) and now the character polynomial's Riesz product. A no-carry numeral system turns addition
into a product of independent digit-problems, and the loneliness gap — being dilation-blind —
reads off only the coarsest factor. The conjecture, in this corner, is a statement about lacunary
independence: when the scales of a speed set decouple, loneliness is free; the difficulty is exactly
the failure to decouple, the few low-degree fusions that tie the scales together.

**Artifacts:** `04-computation/lrc_fractal_translated_ap_s587.py` (+`.out`); formal
`Math/LonelyRunner/FractalSumFree.lean` (`8443da1`); new **HYP-2125**. Builds on HYP-2120/S585.
