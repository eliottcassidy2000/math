---
source: opus-2026-07-11-S269
status: HONEST NEGATIVE (NOT PROVED). Attempting to prove eps_v=O(1) for core>=2 via the cluster/Mayer expansion
  -- the natural rigorous route, built on the PROVEN pairwise bound |Cov(D_v,D_w)|<=1/(3vw) (S262) -- provably
  fails: the leading term does NOT dominate. For core speeds v>=17, |eps_v| exceeds (7/6)Sum_w|Cov(D_v,D_w)| in
  555/658 cases, by up to 47x; the fully rigorous bound (7/6)Sum 1/(3vw) is 113x too small. So eps_v is
  higher-order MULTI-LINEAR dominated, not pairwise. A sharper statement than S266: for v=1 (deep well) the higher
  orders CANCEL a large pairwise; for v>=17 the pairwise is NEGLIGIBLE and the multi-linear resonances ARE eps_v.
  eps_v=O(1) for core>=2 remains the genuine anti-concentration -- NOT reachable by cluster/pairwise/magnitude
  tools, which diverge or fail.
tags:
  - lrc14
  - covering-min
  - anti-concentration
  - cluster-expansion
  - multi-linear
  - honest-negative
  - not-proved
---

# The cluster route to ε_v = O(1) fails: ε_v is higher-order-dominated, not pairwise

**opus-2026-07-11-S269.** Owner: prove `ε_v = O(1)` uniformly for core ≥ 2. This is the cleanest form of the
hard core (S268). I attempted the one rigorous route that rests on a *proven* input — the cluster/Mayer
expansion — and it provably fails. Honest verdict: not proved; the object is irreducibly multi-linear.

## The route and why it was the right one to try

Write `μ_{G'} = ∏_w(1−β_w)/|G'|` (the normalized good-set measure, `β_w = 1_{D_w}`). Then
`ε_v = E_{μ}[g(v·)]`, and expanding the product,

> `ε_v·|G'| = −Σ_w Cov(D_v,D_w) + Σ_{|S|≥2}(−1)^{|S|}⟨g(v·), ∏_{w∈S}β_w⟩`.

The **leading term** is pairwise, `−Σ_w Cov(D_v,D_w)`, and each covariance obeys the **proven** bound
`|Cov(D_v,D_w)| ≤ 1/(3vw)` (S262). If the leading term dominated — i.e. if the higher (connected) correlations
were controlled by the pairwise activity, the standard Kotecký–Preiss / Mayer convergence regime — then
`|ε_v| ≤ (7/6)Σ_w|Cov(D_v,D_w)| ≤ (7/6)Σ_w 1/(3vw) = O(1/v)`, a *fully rigorous* `ε_v = O(1)`. This is the only
route to the O(1) bound built entirely on already-proven inputs.

## The decisive failure

Tested on 658 core speeds `v ≥ 17` across core ≥ 2 covering families:

| test | result |
|---|---|
| `\|ε_v\| ≤ (7/6)Σ_w\|Cov(D_v,D_w)\|` (leading dominates) | **FAILS 555/658**, worst ratio **46.9×** |
| `\|ε_v\| ≤ (7/6)Σ_w 1/(3vw)` (fully rigorous) | worst ratio **113×** |

The leading term does **not** dominate — `ε_v` is typically *far larger* than the pairwise sum. So the
higher-order correlations `⟨g(v·), ∏_{w∈S}β_w⟩` with `|S| ≥ 2` **are** `ε_v`; the pairwise part is negligible.
The polymer gas is not in the convergent low-activity regime for this band (indeed the single-band Fourier
`L¹`-mass `Σ_k|b_k|` diverges, so the naive activity is infinite — the S266 divergence, reappearing).

## The sharper picture (versus S266)

S266 found, for `v=1` (deep well), that the low-order truncation *overshoots* (0.13 vs 0.019) and the higher
orders **cancel** a large pairwise. This session shows the *opposite* mechanism for `v ≥ 17`: the pairwise is
**negligible**, and the higher-order **multi-linear resonances dominate**. Concretely the dominant contribution
is the **noncore-pair resonance** sum — pairs `w_1, w_2` of non-core speeds with `w_1 ± w_2 = ±kv` — the same
additive-energy / resonance object as S263 and LEM-015. Its magnitude sum diverges (harmonic), and its O(1)
value is set by signed cancellation. So across the regimes:

- `v = 1` (deep well, core=1): large pairwise, cancelled by higher orders → **S265 runner-1 lemma**.
- `v ≥ 17` (core ≥ 2): negligible pairwise, **dominated** by multi-linear resonances → this session.

Either way, `ε_v` is **irreducibly multi-linear**, and no bound built from the pairwise covariance (the only
proven input) can reach it.

## Net (honest)

`ε_v = O(1)` uniformly for core ≥ 2 is **not proved**, and this session shows *why* the natural rigorous route
cannot prove it: the cluster/Mayer expansion — the one route grounded in a proven bound (`|Cov| ≤ 1/(3vw)`) —
provably fails because `ε_v` is higher-order-dominated, not pairwise (leading term undershoots by up to 47×).
This is the same anti-concentration wall as S266/S268, now seen from the cluster side and sharpened: the
obstruction is not merely cancellation of a computable leading term, it is that the **dominant** term is the
multi-linear resonance sum itself. Closing `ε_v = O(1)` genuinely requires a bound on that multi-linear
resonance object — the additive-energy / inverse-theorem core the fleet's E3 and Minkowski-tail threads
(LEM-015, tasks #42–#43) target. I did not prove it and I will not claim a route that the computation refutes.

→ opus-S268 (`ε_v=O(1)` = the anti-concentration), opus-S266 (multi-linear, L¹ divergence, v=1 cancellation),
opus-S262 (`|Cov| ≤ 1/(3vw)` — the proven input the cluster route rests on), opus-S263 / LEM-015 (the resonance
/ additive-energy core), opus-S265 (v=1 runner-1 lemma). Files:
`lrc14_cluster_route_fails_eps_is_higher_order_dominated_opus_S269.py` (+`.out`).
