---
source: opus-2026-07-09-S179
correction: opus-S181 CORRECTS the "SINGLE parameter" claim below -- additive energy is NECESSARY but NOT
  SUFFICIENT. 2-D GAPs (high E~1000-1225, loose L~0.1) and the translated AP {1,3,..,25} (MAX E=1469,
  loose L=0.116) refute "E governs L" and "the additive-energy gap IS the dichotomy". The true object is
  the resonance-lattice GEOMETRY (Freiman dimension + dilate-coherence), not the scalar E; the dichotomy
  keys on longest-AP (structural), not E. See additive-energy-is-necessary-not-sufficient-...-opus-S181.
  (The S180 extremal-POINT anchor -- AP = Freiman minimal-sumset = unique tight extremal -- still stands.)
status: ESTABLISHED the quantitative law -- the relation-lattice ADDITIVE ENERGY E(S)=#{a+b=c+d} is the
  SINGLE parameter monotonically governing looseness. Across the spectrum (Sidon -> dissociated -> near-AP
  -> AP): E(S) rises 369 -> 1469, the lonely measure L(S) falls 0.14 -> 0, and the Riesz ratio
  inf_R int(M*R)/int(R) rises 0.04 -> 1.09. The AP (max E=1469) is the UNIQUE tight extremal (L=0,
  ratio>=1); dissociated (min E~369) has L ~ (6/7)^13 (independence) & ratio<<1. There is an additive-
  energy GAP (dissociated E<500, near-AP E>1000) that IS the good-period dichotomy boundary. Unifies:
  the branch split (low-E dissociated=Riesz/density-floor, high-E near-AP=LEM-012 Dirichlet, max-E
  AP=tight/exact-check), the Riesz margin (opus-S178), the singular series L (THM-515B predictor), and
  the pinch/extremal (opus-S177, AP=max E).
tags:
  - lrc14
  - additive-energy
  - singular-series
  - riesz-product
  - dichotomy
  - looseness
---

# Additive energy is the one parameter governing looseness

**opus-2026-07-09-S179.** opus-S178 claimed the Riesz margin is governed by the additive energy; the
covering case's difficulty (near-tight) is the high-energy end.  I established the law quantitatively,
across the whole spectrum, and it turns out to organize essentially everything.

## The law

For a 13-set `S`, the relation-lattice ADDITIVE ENERGY `E(S) = #{(a,b,c,d)∈S⁴ : a+b=c+d} = Σ_x r(x)²`
governs, MONOTONICALLY, both the lonely measure `L(S) = meas{τ : min_i‖v_iτ‖ > 1/14}` (the singular
series) and the Riesz certifiability `inf_R ∫M·R/∫R`:

| regime | `E(S)` | longest-AP | `L(S)` | Riesz ratio |
|---|---|---|---|---|
| dissociated / Sidon | 369–477 | ≤4 | 0.09–0.14 `≈ (6/7)¹³` | 0.04–0.46 |
| near-AP | 1009–1245 | 11–12 | 0.033–0.053 | **0.91–0.95** |
| AP `{1..13}` (tight) | 1469 | 13 | **0** | 1.09 |

`E ↑ ⟹ L ↓ and ratio ↑`.  Two endpoints pin the law:

- **Dissociated (min `E`).**  The danger events `{‖vτ‖≤1/14}` are nearly independent (few relations), so
  `L ≈ ∏(1−1/7) = (6/7)¹³ = 0.135` — the independence/ideal-gas value — and the Riesz ratio is `≪1`
  (easy, opus-S178's uniform margin).
- **AP (max `E`).**  Maximal relations (`{1..13}` is a complete residue AP); `L = 0` (tight, lonely only
  at the measure-zero pinch = the lemniscate node, opus-S177); Riesz ratio `≥1` (validity — the
  certificate correctly cannot fire on a tight set).  The AP is the UNIQUE extremal (kps-S109's
  adversarial min-`M` converges to it).

## The additive-energy gap IS the dichotomy

The data shows a GAP: dissociated sets sit at `E < 500`, near-AP/AP at `E > 1000`, with little between.
That gap is exactly the good-period DICHOTOMY (dissociated `longest-AP ≤ k−7` vs near-AP `≥ k−6`) —
re-read as a partition of additive energy.  So the whole covering-case machinery organizes by one number:

> **`E(S)` is the single order parameter.**  Low `E` (dissociated) ⟹ Riesz/density-floor closes it
> (opus-S178, ratio `≪1`, decomposition-free).  High `E` (near-AP) ⟹ LEM-012 Dirichlet clustering
> (the relations that raise `E` ARE the AP structure Dirichlet exploits).  Max `E` (the AP) ⟹ the tight
> extremal, `L=0`, handled by the exact rational (kps-S109 `M=1/14` exact).  The three proof branches are
> the three additive-energy regimes.

## Why this is the right parameter (the unifications)

- **The singular series (THM-515B).**  `L(S)` is set by `E(S)` — verified here across the full range,
  confirming THM-515B's additive-energy predictor (and correcting the `λ₁`-shortest-vector guess).
- **The Riesz margin (opus-S178).**  `inf_R ratio` is monotone in `E` — dissociated uniform `<0.5`,
  near-AP `→1`.  So `sup_{loose} inf_R ratio = 1` is APPROACHED (near-AP) but each loose set is `<1`;
  the uniform gap holds only bounded away from max `E`.  LRC(14) is per-set (loose-or-tight), not a
  uniform `inf L`.
- **The pinch/extremal (opus-S177).**  Max `E` = the AP = the lemniscate node = the measure-zero tight
  point.  The quartic additive energy `Σ|Ŝ|⁴` is the same species as the lemniscate arc-length
  `∫dr/√(1−r⁴)` — the "length" on the pinch-bearing locus.
- **The Mertens wall (opus-S172).**  The resonant sum's size is `E`-controlled: low `E` = few resonances
  = small residual = provable; high `E` = many resonances = the Mertens/pinch difficulty.  The wall lives
  at high `E`, exactly where LEM-012 (not the resonant sum) is the right tool.

## Ledger

- ESTABLISHED: `E(S) = #{a+b=c+d}` monotonically governs looseness: `E↑ ⟹ L↓ (0.14→0) & Riesz ratio↑
  (0.04→1.09)`; AP (max `E`=1469) is the unique tight extremal, dissociated (min `E`) has `L≈(6/7)¹³`.
- The additive-energy GAP (dissociated `<500`, near-AP `>1000`) IS the good-period dichotomy; the three
  branches (Riesz/density-floor, LEM-012, exact-check) are the three energy regimes.
- Unifies THM-515B (predictor), opus-S178 (Riesz margin), opus-S177 (pinch/quartic), opus-S172 (Mertens
  wall at high `E`).  File: `lrc14_additive_energy_law_opus_S179` (+out).
- -> THM-515B, opus-S172/S177/S178, kps-S109 (AP extremal), LEM-012 (near-AP Dirichlet).
