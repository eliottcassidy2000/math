# The compact case is equivalent to LRC(14), and descent is the wrong tool for it

*boxeph-2026-07-18-S113. Owner: prove the compact case `ρ<13 covering ⟹ M ≥ 1/13`. Honest outcome: this
is the **sole residual** of LRC(14) (S86) — equivalent to the full conjecture — and it is **sharp**
(boundary families at `M=1/13` exactly), so no crude bound reaches it. The natural route (the descent
recursion THM-1010) is **provably too weak** for compact families (loses a factor ~2 at `ρ≈1`). The
dilated-AP-core compact families are handled by THM-1013; the rest is the crux. Not proved. Verified S113
computation.*

## The compact case is the whole conjecture

The LRC(14) reduction map (S86): non-covering (sieve) + `≥2` outliers (THM-726) + single-killer
(THM-724/THM-1007) + **compact (`ρ = v_max/v_2nd < 13`) ⟹ `M ≥ 1/13`**, the *sole residual*. So

> **`ρ<13 covering ⟹ M ≥ 1/13`  is equivalent to LRC(14).** Proving it proves the conjecture.

Empirically robust: over ~100+ compact covering families (range `[1,26]`), `min M = 3/31 ≈ 0.097 > 1/13`,
all with `q* ≤ 47` (shallow binding), **zero** with `M<1/13` — exactly the definitions' "compact ⟹
shallow, `M~0.10–0.14`".

## It is sharp — no slack for a crude bound

The bound is *attained*: `{2·\{1,…,12\}, 13\} = \{2,4,…,24,13\}` is compact (`ρ=24/22≈1.09`), covering,
with **`M = 1/13` exactly** (`M(core) = 1/12` drops to `1/13` when the resonant `v_max=24` is added). So
`ρ<13 ⟹ M≥1/13` is tight — any proof must be sharp at these boundary families, ruling out an
approximate/gap argument (unlike, one might hope, the `1/12`-gap of S110–S111, itself already refuted).

## Descent is the wrong tool (the concrete negative)

The natural attack is the descent recursion **THM-1010**: `M(V) ≥ ρ·M(core)/(ρ+1)`. For **compact**
families `ρ ≈ 1` (`v_max ≈ v_2nd`), so `ρ/(ρ+1) ≈ 1/2` and the bound is only `≈ M(core)/2`. But the
*actual* `M(V) ≈ M(core)` (removing `v_max ≈ v_2nd` barely changes `M`). So descent **loses a factor ~2**
exactly where compactness lives. Verified: descent proves `M≥1/13` for only **5 of 15** compact families;
the failures have `M(core)` up to `0.148` yet descent LB only `0.075 < 1/13`, while the true `M ≈ 0.14`.
Descent is sharp for large `ρ` (single-killer, where it gave THM-1008's `ρ≥13 ⟹ M≥1/14`) and useless for
small `ρ`. The compact case needs a different mechanism.

## What is provable (the dilated-AP sub-case)

The boundary and near-boundary compact families with a **dilated-AP core** are handled by **THM-1013**
(dilated sieve, kernel-pure Lean): if the speeds are `\{d,2d,…,12d\}` (plus a killer), they avoid `13d·ℤ`
by `≥ d`, so `t = 1/(13d)` is `1/13`-lonely — `M ≥ 1/13`. So:

> **dilated-AP-core compact ⟹ M ≥ 1/13** (PROVED, THM-1013).

The residual is the **non-(dilated-AP)-core compact** families — which is the crux, and (like THM-724's
single-killer residual, S112) bottoms on the near-dilated-core rigidity: a non-AP compact core near the
tight locus would need a resonant `v_max`, and controlling that is the inverse theorem.

## Net (honest)

- **Confirmed:** the compact case `ρ<13 ⟹ M≥1/13` **is** LRC(14) (the sole residual, S86), sharp at
  boundary families (`M=1/13` attained), hence not reachable by any crude/gap bound.
- **New negative:** the descent recursion (THM-1010) is provably too weak for compact families — it loses
  a factor ~2 at `ρ≈1` (5/15 proven), so it is the wrong tool. Descent belongs to the single-killer
  (large-`ρ`) regime.
- **Provable fragment:** dilated-AP-core compact ⟹ `M≥1/13` (THM-1013). The rest is the crux.
- **Not proved:** the compact case entails the full conjecture; I did not close it.

So both cases that could carry INV — single-killer (S112) and compact (S113) — bottom on the same
near-dilated-core rigidity: single-killer via THM-724's residual, compact via the non-AP-core residual.
That is one wall, the inverse theorem, seen from two sides. LRC(14) rests there.

Cross-links:
[[INV-val14-is-the-single-killer-case-essentially-done-not-the-open-compact-crux-boxeph-S112]],
THM-1010/1008 (descent), THM-1013 (dilated sieve), THM-724 (single-killer),
[[the-route-B-crux-is-the-open-inverse-theorem-what-covering-gives-and-why-maximality-cannot-finish-boxeph-S101]].
