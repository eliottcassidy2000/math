# The razor's-edge frame: LRC(14) is an EXACT identity (L=0, M=1/14, perfect cancellation at the 6 units mod 14) plus a ONE-SIGNED variance — every quantity is an exact critical value + a non-negative perturbation, not a range

*opus-2026-06-29. Owner: think of everything in terms of perfect cancellation, razor's edge — exact
values plus possible variances, not ranges. This reframes the whole proof: the AP is an exact critical
point (a cyclotomic identity), and the conjecture is that the variance around it is one-signed.*

## The razor's edge is an EXACT identity (verified)
The AP `{1,…,13}` is not "near" `M=1/14` — it is `M=1/14` EXACTLY, with:
> **`L(AP) = 0` EXACTLY** (lonely measure zero — the danger arcs COVER the circle minus measure-zero
> points), and the razor's edge is **exactly the 6 units mod 14**: `min_i‖i·τ‖ = 1/14` precisely at
> `τ ∈ {1,3,5,9,11,13}/14` — the primitive 14th roots, `φ(14)=6`.
The same `6` as the safe `Z₇` cyclotomic modes (`ρ_j=(0,1,1,1,1,1,1)`): `φ(14)=φ(2)φ(7)=1·6`, the units
`k` with `gcd(k,14)=1` are exactly where `i·k ≡ ±1 (14)` makes `min‖·‖=1/14`. **The razor's edge is a
cyclotomic locus** — the perfect cancellation lives on `Q(ζ_14)`, gentlest because `7` is Heegner.

## Perfect cancellation = an exact algebraic identity
`L(S) = Σ_k (−1)^k E_k = E_τ[(1−1)^{X(τ)}] = E_τ[0^{X}] = P(X=0)`. At the AP this is **EXACTLY 0** because
`X(τ)≥1` almost everywhere (the danger covers the circle): the inclusion–exclusion telescopes to zero,
the moments `E_k(AP)` cancel perfectly. **The AP is the configuration where the alternating energy series
sums EXACTLY to 0** — a razor's-edge identity, not a small number.

## Everything is an exact value + a one-signed variance (the reframe)
Stop reading the structures as ranges (`L ≥ floor`, `M ∈ [1/14, …]`). Read them as **exact critical value
+ variance**:
| quantity | EXACT value (razor's edge) | VARIANCE (covering, EXACT rational) |
|---|---|---|
| peak `M` | `1/14` | `M−1/14`: `{1..11,13,84}→9/1246`, `{2..14}→3/56` |
| lonely measure `L` | `0` (perfect cancellation) | `L`: `563/105105`, `3/49` |
| witness `τ` | the 6 units `k/14` (`φ(14)=6`) | the Farey-jump off the units |
| `Σ(−1)^k E_k` | `0` (exact telescoping) | the energy-deficit `> 0` |
| first moment `E_1` | `13/7` (set-independent) | — (none) |
| cyclotomic `ρ_j` | `(0,1,1,1,1,1,1)` | resonance off the 7-scale |
| metagraph `E[H]` | `n!/2^{n−1}` (exact mean) | `Var(H) = (2/n)E[H]² (1+o(1))` |
| OCF `H` | `1` (exact transitive `α₀`) | `2Σα_k` (the cycles) |
> **The mean / `α₀` / the unit / `1/14` is the EXACT value; the higher moments / cycles / energy /
> looseness is the VARIANCE.** The AP sets the variance to its critical value (`L=0`, perfect
> cancellation); every other set is the AP **plus a variance**.

## The conjecture, in this frame: the variance is ONE-SIGNED
LRC(14) `⟺` `M(S) = 1/14 + δ(S)` with **`δ(S) ≥ 0` for every set** (`δ=0` only at the razor's edge,
the AP and its dilates). Equivalently `L(S) = 0 + (a non-negative variance)`. **The AP is the exact
MINIMUM; the proof is that the variance never goes negative — no perturbation overshoots THROUGH the
razor's edge to `M < 1/14`.** This is a second-order (Hessian/variance) statement around an exact
critical point, NOT a range bound:
- The covering perturbation (forced mult-of-14) is one specific direction; its variance is the exact
  rational deficit (`9/1246`, etc.), always `> 0` (verified).
- The "variance" is literally the deviation of the alternating energy sum from its perfect-cancellation
  value `0` — and the covering ENERGY-DEFICIT (`E_k(AP) − E_k(covering) > 0`, verified `k=3,4`) is what
  makes `δ > 0`.

## Why this frame is the right one (it dissolves the old obstructions)
- **The union bound is the EXACT value `E_1=13/7`** — it is the mean, carrying ZERO variance information;
  of course it is "vacuous," it is not the variance.
- **The truncations alternate because they are partial sums of the EXACT cancellation** — the razor's
  edge `Σ(−1)^k E_k=0` is reached only in the limit (conditional convergence, THM-504); a partial sum is
  the exact value minus a tail, never the value.
- **The floor is the set-independent exact part (the cyclotomic `ρ_j`, the 6 units); the cap is the
  variance** (the set-dependent energy-deficit). Floor = exact value; cap = variance.
The whole proof is: **certify that the variance around the exact razor's-edge identity is non-negative**,
i.e. the perfect-cancellation `Σ(−1)^k E_k=0` of the AP is a strict minimum of `L`, and the covering
energy-deficit is one-signed positive. This is the Riesz-product resummation (THM-515) read as a
variance/Hessian certificate at the cyclotomic critical point.

## Status
- **Verified (opus, EXACT):** `L(AP)=0`, `M(AP)=1/14`; razor's edge = the 6 units mod 14
  (`{1,3,5,9,11,13}/14`, `φ(14)=6`); covering variances exact rationals (`9/1246`, `3/56`, `563/105105`,
  `3/49`), all `> 0`.
- **Reframe (owner+opus):** every structure = exact critical value + one-signed variance; LRC(14) =
  the variance `δ(S) ≥ 0` around the exact razor's-edge identity; the perfect cancellation
  `Σ(−1)^k E_k=0` is the AP's exact minimum, lives on `Q(ζ_14)`.
- **Target (sharpened):** the variance/Hessian certificate — the AP's perfect cancellation is a strict
  minimum and the covering energy-deficit is one-signed (the Riesz resummation as a positivity at the
  cyclotomic critical point).

Related: the additive-energy/THM-504-wall reflection (`L=Σ(−1)^k E_k`, the energy-deficit), the
Z₇-cyclotomic-SOS-floor (the 6 safe modes = the 6 units), the master-reframe (bulk/tail = exact/variance),
THM-515 (Riesz program), the peak-cut-deficit (`M=1/14+`cut-deficit), the concentration `CV²~2/n` (exact
mean + variance), OPEN-Q-108.
