# HYP-2200 — Four lenses on the covering-depth object: Helly order, Vitali wall, Collatz two-block, partition-function sibling

**Session:** claudebox-2026-06-03-S618. **Consolidates:** HYP-2195 (covering-depth master object), HYP-2185 (apex
sheaf), HYP-2175 (Collatz twin), HYP-2180 (altitude). **Prompt (user):** Helly = "how many orders of overlap you
must keep," Vitali = "finite moments cannot decide p₀," Collatz/two-block = "correlated residue where density is
blind," OCF/partition functions = "the sibling world where independence polynomials already play the
depth-distribution role."

## The single object
The covering-depth distribution `{p_k}` is governed by the **overlap hierarchy** `S_k = Σ_{|I|=k} meas(∩_{i∈I} A_i)`
(the factorial moments / inclusion-exclusion terms of the forbidden arcs `A_i`). The lonely measure is the full
alternating sum
`p₀ = Σ_{k=0}^n (−1)^k S_k`  (verified exactly: incl-excl = direct, all configs).
The four lenses are four readings of this one hierarchy.

## Lens 1 — Helly: how many overlap-orders you must keep (verified)
The Bonferroni truncations `T_m = Σ_{k≤m}(−1)^k S_k` bracket `p₀` (odd m = lower bound, even m = upper). The **order-3
lower bound** `T₃ = 1 − S₁ + S₂ − S₃` is a **rigorous loneliness certificate**: `T₃ > 0 ⟹ p₀ > 0`. Verified: gap
configs are decided at order 3 (`{1,4,6,9}` T₃=+0.067, `{2,3,7,8}` +0.057, `{2,5,9,11}` +0.070) — the **circular-arc
Helly number 3**. The collapse family always has `T₃ ≤ 0` (the resonance lives at the order-3 triples — sum
relations `a+b=c`). **Refinement (honest):** order 3 is sufficient but not necessary off-boundary — a thin boundary
layer exists (`{1,5,8,11,13}`: T₃=−0.030 yet p₀=0.092), where order 3 fails and higher orders/construction are
needed. So: keep 3 overlap-orders to decide loneliness EXCEPT in a boundary layer hugging collapse.

## Lens 2 — Vitali: finite moments cannot decide p₀ at the boundary (verified)
At collapse (`p₀ = 0`) NO finite truncation `T_m`, `m < n`, equals 0; the alternating series closes to exactly 0
only at the top order `m = n` (`{1,3,4,7}`: T = +1,−.6,+.281,−.057,+.000). Finitely many moments give Bonferroni
brackets that pinch `p₀` only in the limit — at the collapse boundary the measure functional abdicates, and `p₀=0`
cannot be separated from `p₀=0`-with-a-witness (tight, LRC holds) vs `p₀=0`-empty (LRC fails). The Vitali wall =
measure blind, construction required. (HYP-2150's measure-vs-construction boundary, made quantitative.)

## Lens 3 — Collatz two-block: the multiplicative twin of the Vitali wall
On the multiplicative face (HYP-2175), "almost all" / density sees density-1 sets; the structured residue (cycles,
the 2-adic apex/2-block) is density-zero and **invisible to density** — exactly the Vitali wall transported across
the additive↔multiplicative dictionary. Density (first moment / measure) is blind to the correlated residue on both
faces; the residue is the collapse family (additive) = the cycle/apex (multiplicative).

## Lens 4 — partition-function sibling (verified + formalized)
The depth generating function `P(z) = ∫ z^{depth(t)} dt = ∫ ∏ᵢ(1+(z−1)1_{Aᵢ}) dt` is a **hard-core partition
function**; `p₀ = P(0)`. Its integrand factors over non-interacting blocks (formalized `depthGF_union`) — the
independence-polynomial / hard-core-gas multiplicativity. **Sharpening:** at the LRC gap `δ = 1/(n+1)` the pairwise
arc-dependency graph is the COMPLETE graph (every arc pair correlated), so the crude independence polynomial is
degenerate `1+nλ` (Shearer λ*=1/(n−1) for K_n) — exactly as the first moment is vacuous (HYP-2195). The
discriminating object is therefore the FULL `P(z)` / all overlap orders, not the pairwise sibling. All four lenses
converge: **you need the whole hierarchy.**

## Formalized (math-lean, sorry-free) — `Math/LonelyRunner/DepthGenerating.lean`
`depthGF` (the partition function), `depthGF_one` (normalization P(1)=1), `depthGF_union` (factorization over
disjoint blocks — lens 4), `prod_one_sub_eq_inclusion_exclusion` (the overlap hierarchy — lenses 1,2),
`depthGF_zero` (`p₀ = P(0)` = the inclusion-exclusion alternating sum).

## The consolidated statement
`p₀ = Σ(−1)^k S_k` is the master identity; Helly says order 3 decides off-boundary, Vitali says no finite order
decides on the boundary, Collatz is the boundary's multiplicative twin, and the partition function packages the
whole hierarchy whose zeros decide `p₀ > 0`. The collapse family (additive chains = resonances = apex sheaf H⁰=∅)
is exactly the boundary where order 3 fails, the truncations don't close, density goes blind, and the partition
function's value at 0 vanishes — one boundary, four names.

## Open
- Prove `T₃ > 0` (order-3 Bonferroni) for an explicit large class (the off-boundary loneliness certificate).
- Characterize the boundary layer where order 3 fails but `p₀ > 0` (between Helly-decidable and collapse).
- Lee–Yang / Shearer zeros of `P(z)` ↔ resonance energy (HYP-2155); does the smallest zero track the apex/2q seam?
- The `OCF`/partition-function transfer: which hard-core-gas theorems (cluster expansion, correlation decay) give
  LRC loneliness bounds via `P(z)`?
