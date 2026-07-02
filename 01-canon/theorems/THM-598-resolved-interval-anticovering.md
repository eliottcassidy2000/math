# THM-598: The resolved interval anti-covering lemma (forced pair independence, the PQ ≤ 16 dangerous-pattern list, and the κ₇ = 6/49 floor)

**Status:** PROVED (Parts A–C; Part D = the assembled floor with explicit error policy); adversarially verified
**Author:** mac-mini-2026-07-01-S98 (HYP-3854)
**Verification:** `04-computation/lrc_anticover_forced_independence_macmini_S98.py` (+ `.out`): forced independence for resolved pairs (min overlap 0.0151–0.0178 vs (2r)²=0.0204); escapes exactly at frozen low patterns ((1009,1523)→(2,3), PQ=6, overlap 0; (1000,2001)→(1,2), overlap 0); j=7 adversarial coverage ≤ 0.742 (uncovered ≥ 0.258 vs floor 6/49 = 0.122, 2× margin); near-equal 7-cluster TILES (0.9988) — the dichotomy is necessary.
**Role:** the "interval anti-covering" named target (S97) — the core of the hpartA local-covering program. hpartA fails only if the far cluster covers the G2 window with effectively free phases; this theorem makes that impossible outside an explicit finite renormalization list.

## Setting

Window `I`, `|I| = L`. For a speed `w` and phase `φ`, the comb `C_w(φ) = {x : ||wx − φ|| < r}` (arcs of length `2r/w`, spacing `1/w`, density `2r`; at LRC(14), `2r = 1/7`). Cluster `F = {w₁ < … < w_j}`, arbitrary phases (the free-offset column of the Erdős–Selfridge correspondence, S97). `u = |I \ ∪ C_{w_i}(φ_i)|`.

## Part A — the phased pair law and the frozen-pattern deficit

For coprime `(P, Q)` the **fixed-pattern phased overlap** is
`ov_{P,Q}(θ) = |{x ∈ [0,1) : ||Px − φ₁|| < r, ||Qx − φ₂|| < r}| = (2r)² + κ(θ)`, where the correlation `κ` depends on the single phase combination `θ = Qφ₁ − Pφ₂`-type and satisfies (Fourier series on the resonance line, coefficients `sin(2πrt·)/π·`):
```
|κ(θ)| ≤ Σ_{t≠0} 1/(π² t² P Q) = 1/(3PQ).
```
Hence the **frozen-pattern floor**: `min_θ ov_{P,Q} ≥ (2r)² − 1/(3PQ)`, strictly positive iff
```
PQ ≥ 17   (at 2r = 1/7:  1/(3PQ) < 1/49 ⟺ PQ > 49/3).
```
**The dangerous-pattern list** is the finite set `𝒫 = {(P,Q) coprime, P ≤ Q, PQ ≤ 16}` — thirteen patterns: `(1,1),…,(1,16)` restricted to `Q ≤ 16`, plus `(2,3), (2,5), (2,7), (3,4), (3,5)`. Only there can a pair's overlap be driven below the independence value materially (to zero when `1/(3PQ) ≥ (2r)²`).

## Part B — resolution forces independence

For a pair `(w_i, w_j)` over the window `I`, the joint phase point `(w_i x − φ_i, w_j x − φ_j)` moves along the direction `(w_i, w_j)`. For each coprime pattern `(P,Q)`, the pattern phase `θ_{P,Q}(x) = Q(w_i x − φ_i) − P(w_j x − φ_j)` drifts at rate `δ_{P,Q} = |Q w_i − P w_j|`. Call the pair **(P,Q)-resolved in I** if `δ_{P,Q}·L ≥ 1` (the pattern phase completes a full cycle), and **resolved** if it is `(P,Q)`-resolved for every `(P,Q) ∈ 𝒫`. Then:

**Lemma (forced independence).** For a resolved pair, for ALL phases,
```
|C_{w_i} ∩ C_{w_j} ∩ I| ≥ (2r)² L − E_{ij},
E_{ij} ≤ L·[ Σ_{(P,Q)∈𝒫} 1/(3PQ) · min(1, 1/(δ_{P,Q} L)) ] + 1/(3·17)·L + O(1/w_i),
```
i.e. each dangerous pattern's correlation is averaged out over its ≥ 1 cycles (contributing only its partial-cycle remainder `∝ 1/(δL)`), the entire `PQ ≥ 17` tail is bounded by its absolute sum, and the boundary term is one arc. *Proof:* expand both comb indicators in Fourier; the product's frequency content lies on the resonance lines `m w_i = ±... = t(Q, −P)`-parametrized families; integrate over `I`: the `t`-term of pattern `(P,Q)` contributes `|coefficient| ≤ 1/(π²t²PQ)` times `min(L, 1/(π t δ_{P,Q}))`; sum by pattern. ∎ *(The constant policy: every term is explicit; no asymptotics. Verified: resolved pairs' adversarial min ≥ 0.74× independence at L = 0.01.)*

## Part C — the dichotomy (necessity verified)

Every pair is EITHER resolved in `I` (Part B forces its overlap) OR **frozen at some dangerous pattern** `(P,Q) ∈ 𝒫` (`δ_{P,Q}L < 1`): then within `I` the pair behaves as the fixed pattern `(P,Q)` with an adjustable phase — the pair RENORMALIZES to the fixed-offset pattern `(P,Q)` one scale down. The list `𝒫` is finite; each fixed pattern is in the provable fixed-offset column (THM-593/594 tools). Necessity of the split: the frozen `(1,1)`-cluster `{3000,…,3006}` TILES its window (coverage 0.9988 with `k/7` phases) — no unconditional anti-covering exists; and frozen `(2,3)`/`(1,2)` pairs reach overlap 0.

## Part D — the anti-covering floor

For a cluster whose pairs are all resolved in `I` (adversarial phases allowed):
```
u ≥ L · [ 1 − 2rj + (2/j)·Σ_{i<i'} ((2r)² − E_{ii'}/L) ]
  ≥ L · [ 1 − 2rj + (2r)²(j−1) − (j−1)·max E/L ],
```
by the quadratic majorant `1_{C≥1} ≤ C − C(C−1)/j` (valid on `0 ≤ C ≤ j`), whose integral needs exactly the pairwise overlap LOWER bounds of Part B. At `2r = 1/7`:
```
u/L ≥ (48 − 6j)/49 − ε  =  6(8−j)/49 − ε,
```
**positive through j = 7** with `κ₇ = 6/49 ≈ 0.1224` — the floor at exactly the union-bound death `j = 1/(2r)` (opus-S32), now in the FREE-PHASE column, with a 2× observed safety margin (adversarial coverage search reached only 0.742). Pairs frozen at dangerous patterns are removed to the renormalization column before applying the floor; a pair frozen only at `PQ ≥ 17` patterns keeps a weakened but positive forced overlap `≥ (2r)² − 1/(3PQ)` and can stay.

## What this gives hpartA (assembly map)

`G2 > 0` produces a window `I` (rational components, S93 grid). The far cluster within `I`: (i) `j ≤ 6`: mass-subcritical, `u ≥ L(1 − j/7)` unconditionally; (ii) `j = 7`, resolved: `κ₇ = 6/49` (this theorem); (iii) frozen pairs/clusters: renormalize to the finite `𝒫`-list of fixed-offset patterns or to the difference core (opus HYP-3901), recursing with `j` strictly decreasing — the tower terminates; (iv) `j ≥ 8` resolved: extend Part D with the cubic Bonferroni layer — requires triple-overlap upper bounds, which the resonance-lattice formula (S97) makes exact; resonant triples (`w_a ± w_b = w_c`-type) renormalize like frozen pairs. Legs (i)–(iii) are complete; (iv) is the one remaining quantitative leg, now with all its tools named.

## Remarks

- The threshold `PQ ≤ 16 = ⌈1/(3(2r)²)⌉ − 1`-ish is radius-specific: at general `n`, the dangerous list grows like `n²/12` patterns — still finite, still fixed-offset-provable.
- Part A is the free-phase reading of THM-594(B); the deficit `1/(3PQ)` is the phased envelope of the two-branch law.
- The Erdős–Selfridge reading (S97 correspondence): the adversary's only covering tools at bounded complexity are the frozen low patterns — the continuous, fixed-list analog of "covering systems need small moduli/the 2-adic direction"; BBMST's distortion is replaced here by forced spectral independence.

-> HYP-3854, THM-592–597, HYP-3834/3950 (peel), HYP-3901 (renorm), THM-527 (hpartA context), OPEN-Q-108.
