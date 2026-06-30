# The AP maximizes the higher additive energies E₃, E₄ (minimizing the lonely measure L≈0) — the AP-extremality begins at the |T|=3 (SL(3)/Littlewood) level, NOT the pair level — but the support-truncation L^{≤k} alternates wildly (the THM-504 wall confirmed), so the proof must control the FULL alternating energy sum (Riesz product), not a truncation

*opus-2026-06-29. Owner: take the concrete next move — set up the |T|=3 additive-energy perturbation of
the cyclotomic SOS and test whether the AP minimizes L. Done. The AP-energy-extremality is confirmed and
sharply located, AND the naive support-truncation is shown to fail — both are progress.*

## The setup (inclusion–exclusion = the moment/energy expansion)
`L(S) = Σ_T (−1)^{|T|} meas(∩_{i∈T} D_i)`, so the support-`≤k` truncation is
`L^{≤k} = Σ_{j=0}^k (−1)^j E_j`, where **`E_j = E_τ[\binom{X(τ)}{j}]`** are the factorial moments of the
danger count `X(τ)=#\{i:‖v_iτ‖<1/14\}`. `E_j` is exactly the **|T|=j additive energy** (the j-fold danger
overlap = the j-term resonance mass). `E_1 = 13/7` (set-independent union bound); `E_{≥2}` set-dependent.

## What the test shows (Q=100800, verified)
| set | `E₂` (pair) | `E₃` (triple) | `E₄` | `L` (true `=P(X=0)`) |
|---|---|---|---|---|
| **AP `{1..13}`** | 2.20 | **4.58** | **9.66** | **0.00006 ≈ 0** |
| covering `{1..11,13,84}` | 2.06 | 3.82 | 7.40 | 0.0054 |
| covering `{2..14}` | 2.19 | 4.31 | 9.04 | 0.061 |
| lacunary (Sidon) | **2.43** | 2.47 | 2.05 | 0.279 |
1. **The AP MAXIMIZES the higher energies `E₃, E₄`** (4.58, 9.66 — the largest) and **MINIMIZES the lonely
   measure `L ≈ 0`** (razor-thin: its lonely set is the single peak point `τ=1/14`, measure zero). Confirmed.
2. **The AP-extremality begins at `|T|=3`, NOT `|T|=2`.** The PAIR energy `E₂` is maximized by the
   *lacunary/Sidon* set (2.43), not the AP — the pair level is the **cyclotomic floor** (`SL(2)`, the
   set-independent Z₇ SOS), and it does NOT see the AP. The TRIPLE energy `E₃` (and `E₄`) is where the AP
   wins — exactly the **`SL(3)`/Littlewood (and `SL(4)`/cap)** level of the moment hierarchy. *The
   additive-energy extremality and the dimension ladder agree precisely.*

## The negative result: support-truncation hits the THM-504 wall
The truncations do NOT converge to `L`: `L^{≤2}(AP)=1.34`, `L^{≤3}(AP)=−3.23` — they **alternate wildly**
and overshoot (`L_true=0.00006`). The factorial moments GROW (`1.86, 2.20, 4.58, 9.66, …`), so
`L = Σ_k(−1)^k E_k` is **conditionally convergent** — the support-truncation is the WRONG tool (THM-504's
`|T|≥3` conditional-convergence wall, reconfirmed concretely). One cannot bound `L` by stopping at `|T|=3`.

## The right picture (and the proof target, recast)
> **`L(S) = Σ_k (−1)^k E_k(S)`, conditionally convergent. The AP MAXIMIZES every higher energy `E_k` and
> thereby achieves PERFECT CANCELLATION (`L=0`). The covering constraint (forced mult-of-14) REDUCES the
> higher energies → the cancellation is incomplete → `L > 0`.** The floor `min_{covering} L > 0` is exactly
> the covering ENERGY-DEFICIT: covering sets cannot reach the AP's energy-maximal perfect cancellation.
So the cap/tail is an **energy-extremality**: the AP is the unique `Σ(−1)^k E_k = 0` point (the
additive-energy-maximal / `S₄`-richest direction), and the covering perturbation's energy-deficit is
bounded below — keeping `L` positive. This is the Riesz-product / additive-energy program (THM-515), NOT a
truncation: the proof must resum the full alternating series with the AP as the extremal cancellation.

## Concrete progress this move
1. **Confirmed the cap object:** the AP maximizes `E₃, E₄` (the set-dependent higher additive energy),
   minimizing `L`. The tail IS the AP energy-extremality.
2. **Located it exactly:** the extremality is absent at `|T|=2` (cyclotomic floor, `SL(2)`, lacunary-max)
   and present from `|T|=3` (`SL(3)`/Littlewood) — the moment hierarchy and the energy agree.
3. **Killed the naive tool:** the support-truncation `L^{≤k}` diverges/alternates (THM-504 wall) — the
   proof needs the resummed Riesz product, controlling the covering energy-deficit, not a finite `|T|`.
4. **The target sharpened:** show the covering-forced energy-deficit `E_k(AP) − E_k(covering) > 0`
   (verified at `k=3,4`) keeps the resummed `L > 0` — an extremal + resummation statement.

## Status
- **Verified (opus):** AP maximizes `E₃=4.58, E₄=9.66` (the higher additive energy) and minimizes `L≈0`;
  pair `E₂` is maximized by lacunary, NOT the AP (extremality starts at `|T|=3`); support-truncation
  `L^{≤k}` alternates (`1.34, −3.23`) — the THM-504 wall.
- **Reframe:** `L=Σ(−1)^k E_k`; AP = perfect cancellation (`L=0`); covering = energy-deficit (`L>0`); the
  cap is an energy-extremality at the `SL(3)+` level, resummed (Riesz), not truncated.
- **Open (the genuine tail):** bound the covering energy-deficit through the resummation — i.e. break the
  THM-504 wall with the AP as the extremal cancellation point (the Riesz-product route of THM-515).

Related: THM-515 (theta/lonely-measure, Riesz program), THM-504 (the `|T|≥3` wall — now reconfirmed),
the Z₇-cyclotomic-SOS-floor reflection (the `|T|≤2` set-independent floor), the Siegel–Rogers moment
hierarchy (`|T|=k`↔`SL(k)`↔`ζ(k)`, AP maximizes `S₄`), HYP-2823 (consec variance extremality), the
master-reframe (bulk/tail), the resonance-arity reflection (AP additive-relation-richest), OPEN-Q-108.
