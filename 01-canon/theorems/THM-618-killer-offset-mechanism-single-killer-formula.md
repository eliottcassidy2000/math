---
id: THM-618
title: The killer-offset mechanism — covering forces a runner divisible by n−1 that sits at 0 at the small-base optimum t=1/(n−1), so the hiding spot offsets to where the killer binds, giving the covering-min = 1/(n−1) − 1/((n−1)Φ₆) = n/Φ₆. Exact formula for the single-killer ladder: M({1,…,n−2,X}) = X/((n−1)(X+1)) for (n−1)n | X (a runner-1 vs killer equioscillation), monotone increasing in X, minimized at the smallest killer X=(n−1)n=182 → 14/183.
status: PROVED (the single-killer formula, by explicit equioscillation + exact-M check; the mechanism, elementary). VERIFIED exactly for X=182,364,…,2366.
source: mac-mini-2026-07-04-S43
depends_on:
  - THM-523   # q-witness: {1,…,n−2} hides at 1/(n−1) (the small-base optimum); covering forces the (n−1)|a killer
related:
  - HYP-4060  # kps: primitive covering-min = 14/183, deep well unique extremizer
  - HYP-4076  # (klein) the Ostrowski ladder; this is the {1..n−2}-base ladder with a geometric derivation
  - HYP-4082  # klein: one-swap covering stratum = 13 formula-closable ladders floored by 14/183
external: LRC(14); covering-min n/Φ₆(n).
---

# THM-618 — the killer-offset mechanism and the single-killer formula

## The mechanism (why the covering-min is `14/183`)
A covering family covers `q=n−1`, so it has a **killer** `a` with `(n−1)∣a`. The small runners `{1,…,n−2}`
have their max-min at `t=1/(n−1)` (`M({1,…,n−2})=1/(n−1)`, the q-witness). But there the killer sits at
`‖a/(n−1)‖ = ‖(a/(n−1))·1‖ = 0` (an integer) — **the killer blocks the small-base optimum.** So the family
cannot hide at `1/(n−1)`; the hiding spot **offsets** to `t* = 1/(n−1) − δ`, where the killer just becomes
safe (`‖a t*‖` rises from `0` to `M`). The covering-min is that offset value:
```
   covering-min = 1/(n−1) − 1/((n−1)Φ₆) = (Φ₆−1)/((n−1)Φ₆) = n(n−1)/((n−1)Φ₆) = n/Φ₆ = 14/183.
```
The offset `δ = 1/((n−1)Φ₆) = 1/2379` is one `1/((n−1)Φ₆)` unit — the arithmetic shadow of the killer.
(Verified: every covering family has an `(n−1)∣a` killer and min-dist `0` at `t=1/(n−1)`; deep well hides at
`14/183`, killer `182` binding, offset `1/2379`.)

## The single-killer formula (PROVED)
> For `X` a multiple of `(n−1)n` (a single killer covering both `n−1` and `n`),
> **`M({1,…,n−2, X}) = X/((n−1)(X+1))`.**

**Proof.** Take `t = 1/(n−1) − ε` with `0<ε` small. Runner `1`: `‖t‖ = 1/(n−1)−ε`. Killer `X` (`(n−1)∣X`):
`‖Xt‖ = ‖X/(n−1) − Xε‖ = ‖−Xε‖ = Xε`. Runner `k` (`2≤k≤n−2`): `‖kt‖ = ‖k/(n−1) − kε‖`, which for small `ε`
is `≥ 1/(n−1) − (n−2)ε > ` the binding value below (the extreme small runner `k=n−2` gives `‖−1/(n−1)−(n−2)ε‖
= 1/(n−1)+(n−2)ε`, larger still). So the min is `min(1/(n−1)−ε,\ Xε)`, maximized where they meet:
`1/(n−1)−ε = Xε ⟹ ε = 1/((n−1)(X+1))`, giving `M = 1/(n−1) − 1/((n−1)(X+1)) = X/((n−1)(X+1))`. This is a
2-point equioscillation (runner `1` ↔ killer `X`), hence the global max. ∎ *(Exact-M verified `X=182,…,2366`.)*

## Consequence: the covering-min over the single-killer `{1,…,n−2}`-ladder is the smallest killer
`X/((n−1)(X+1))` is **strictly increasing in `X`** (derivative `1/((n−1)(X+1)²)>0`), so it is minimized at
the **smallest** admissible killer `X = (n−1)n = 182` (the least common multiple of `n−1,n`), giving
`182/(13·183) = 14/183`. Larger killers (or split killers `{…,n−1,n}`: `M({1..11,13,14})=1/12`, looser) only
increase `M`. So on this ladder the covering-min is exactly the deep well, and the value `14/183` is the
**minimal-killer / base-optimum-minus-one-offset** point.

## Scope (honest)
This proves the covering-min on the single-killer `{1,…,n−2}`-base ladder (the extremizer's own ladder). The
*full* covering-min needs the other strata too — split killers, non-`{1..n−2}` bases, multi-swaps — which are
klein's 13 ladders (HYP-4082) + kps's per-rung residue formulas, and are all `≥ 14/183` (looser). The
killer-offset mechanism gives the **geometric why** (`base-optimum − killer-offset`) uniformly, and a clean
equioscillation proof for the extremal ladder; it does not by itself close the other strata.
