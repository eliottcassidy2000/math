---
id: THM-502
title: The closed-walk census ladder — tr(A^k) = sum over connected cycle-configs of (k/period)*embeddings; gives the exact decompositions for k<=8 and the spectral-horizon conservation law
status: PROVED (coefficient law: rigorous rooting argument; explicit formulas k<=8: verified exact on 300+ random tournaments at n=8, normal-equation exact-rational fits; k=9 coefficient-structure verified 80/80, full triple-config enumeration OPEN)
source: monad-explorer-2026-06-14
depends_on:
  - THM-118   # c_k = tr(A^k)/k for k <= 5
  - THM-499   # first spectral boundary: H non-spectral from n=6 via alpha_2
  - THM-500   # second spectral boundary: alpha_1 non-spectral from n=7 via c7/TQ
related:
  - HYP-2498  # codex: tr(A^6) = 6 c6 + 3 c3 + 6 p33_meet
  - OPEN-Q-093
---

# THM-502 — the closed-walk census ladder and the k/period coefficient law

THM-499/THM-500 located two "spectral boundaries" — where the eigenvalue spectrum
loses the OCF count `H` (n=6) and then the bare odd-cycle count `alpha_1` (n=7). Both
breaks were traced to a *trace decomposition*: `tr A^6 = 6c6 + 3c3 + 6 p33` (codex
HYP-2498) and `tr A^7 = 7(c7 + TQ)` (THM-500). This theorem gives the **single
principle** behind every such decomposition, the **explicit complete formulas through
k=8**, and the **conservation law** that makes the spectral horizon transparent.

## The coefficient law (proved)

In a tournament, the diagonal of `A^k` counts **rooted closed k-walks** (a closed
walk with a marked starting step), and `tr A^k` sums them. Loop-erasing a closed walk
returns a multiset of simple directed cycles whose lengths **partition k**, every part
`>= 3` (a tournament has no 1- or 2-cycles), and which are pairwise **overlapping**
(a single closed walk is connected — vertex-disjoint cycles cannot lie on one walk).
Call such a connected multiset of cycles a **configuration** `C` of total length k.

> **Coefficient law.** `tr A^k = Σ_{configs C, |C| = k}  rot(C) · emb_T(C)`, where
> `emb_T(C)` is the number of embeddings of `C` into `T` and `rot(C) = k / period(C)`
> is the number of distinct cyclic rotations of the closed-walk edge-word that traces
> `C`. Equivalently `rot(C)` = number of distinct rooted closed walks with reduced
> shape `C`.

*Proof of the law.* A rooted closed k-walk is a cyclic edge-word of length k with a
marked position; `tr A^k` counts these. Group them by reduced shape `C`. Words with
shape `C` form orbits under the `Z/k` cyclic shift; the marked position ranges over
all `k` positions, so each *unmarked* shape gives exactly `k/period` distinct marked
words, where `period` is the least shift fixing the word. A simple m-cycle has an
aperiodic word (distinct vertices) so `period = m` and `rot = 1` per *cycle* but the
m rotations are m distinct rooted walks, giving the standard `tr A^m ⊇ m·c_m`; a
configuration of `r` identical primitive copies has `period = k/r`, hence
`rot = r`. ∎

Two immediate corollaries of the law, used below:
- a **simple k-cycle** has `rot = k` → contributes `k·c_k`;
- a **(k/2)-cycle traversed twice** has `period = k/2`, `rot = k/2` → `(k/2)·c_{k/2}`
  (k even); generally a **d-cycle traversed r times** (k = rd) → `d·c_d`;
- **two distinct overlapping cycles** (lengths a≠b, or two distinct equal-length
  cycles) glue into an asymmetric word, `period = k`, `rot = k` → `k·O(a,b)` where
  `O(a,b)` = number of (a-cycle, b-cycle) pairs sharing ≥1 vertex.

## The explicit ladder (k <= 8): the only configs are single + pair

For `k <= 8` no partition of `k` into parts `>= 3` has **three** parts (the first is
`9 = 3+3+3`). So every configuration is one cycle or one overlapping pair, and the law
gives closed forms. **Verified exact** (300/300 random tournaments at n=8, and exact
rational normal-equation fits; the equal-overlap coefficient is uniform across
intersection sizes |∩| ∈ {1,2,3}):

```
tr A^3 = 3 c3
tr A^4 = 4 c4
tr A^5 = 5 c5
tr A^6 = 6 c6 + 3 c3 + 6 p33                     p33 = overlapping triangle pairs
tr A^7 = 7 c7 + 7 TQ                             TQ  = overlapping (triangle,4-cycle) pairs
tr A^8 = 8 c8 + 4 c4 + 8 Q44 + 8 TF              Q44 = overlapping 4-cycle pairs
                                                 TF  = overlapping (triangle,5-cycle) pairs
```

(`tr A^8` is **new** here; `tr A^6` matches codex HYP-2498, `tr A^7` is THM-500.)
The terms are exactly the partitions of k into parts ≥3 carrying the law's
coefficients: k=8 ↔ {8}: `8c8`, {4,4}-doubled: `4c4`, {4,4}-distinct: `8 Q44`,
{3,5}: `8 TF`.

## The conservation law (the spectral horizon, made transparent)

The spectrum fixes every `tr A^k`. Reading the ladder bottom-up, `c3,c4,c5` are
spectral (THM-118), so in the decomposition of `tr A^k` the only **non-spectral**
freedom is a single conserved quantity: the simple top-length count `c_k` trades
1-for-1 against the total overlap count. **Verified identically** (offset = 0 on every
cospectral class, 300/300):

```
c6 + p33        = (tr A^6 - 3 c3)/6      (spectral)   k=6
c7 + TQ         =  tr A^7 / 7            (spectral)   k=7
c8 + Q44 + TF   = (tr A^8 - 4 c4)/8      (spectral)   k=8
```

So **within a cospectral class** `c_k + (overlap count)` is constant; cospectral mates
differ by moving weight between the *simple cycle* and the *overlapping reducible
configs*. This is the exact mechanism of `c_k` non-spectrality for `k >= 6`:

> The spectrum sees the **sum** (all closed k-walks) but not the **split** (simple
> k-cycle vs. reducible overlap). Whenever the overlap count is itself non-spectral
> (first at n=6 for p33), `c_k` is non-spectral.

This pins the **spectral horizon table** (`spectral_horizon_table_monad.py`,
exhaustive n≤6): `c3,c4,c5` spectral at all n; `c6` and `p33` and `alpha_2` (hence
`H`) **first non-spectral at n=6**; `c7,TQ` first at n=7 (THM-500). In particular
`c6` is non-spectral *from its onset* — it does **not** get a delayed break — which is
why `alpha_1` (delayed to n=7) is the unique OCF-derived invariant with a nontrivial
spectral window (THM-500 / the resolution-ladder reflection, claim 3, now verified).

## The Witt/necklace form — why the conserved quantity is spectral

The conserved combinations above are not ad hoc: they are the **necklace (Witt)
transform of the trace sequence**. For a tournament `tr A^1 = tr A^2 = 0`, so

```
W_k(A) := (1/k) Σ_{d | k} μ(d) · tr A^{k/d}
W_6 = (tr A^6 - tr A^3)/6     W_7 = tr A^7 / 7     W_8 = (tr A^8 - tr A^4)/8
```

and **verified 120/120** (n=8): `W_6 = c6 + p33`, `W_7 = c7 + TQ`,
`W_8 = c8 + Q44 + TF`. The reason is structural: `W_k(A)` counts the **aperiodic
closed k-walks up to cyclic rotation** (the standard meaning of the Möbius/Witt
transform of a walk-count sequence). An aperiodic closed walk's reduced shape is
either a **simple k-cycle** (one cyclic word ↔ `c_k`) or an **asymmetric overlapping
pair** (one cyclic word ↔ an overlap config); the periodic ones — a `(k/2)`-cycle
doubled — are exactly the terms the Möbius transform removes. So:

> **The Witt transform of the eigenvalue power-sums equals (simple cycles) +
> (overlapping reducible configs).** It is manifestly spectral (a Z-combination of
> traces), which is precisely *why* `c_k + (overlap count)` is constant on every
> cospectral class. The non-spectral content of the cycle vector is the *partition*
> of the spectral number `W_k` into its simple vs. overlapping parts.

This is the clean, all-`k` statement (the named-piece decomposition `single + pair`
is the `k <= 8` special case; from `k = 9` the aperiodic configs also include
triples, but `W_k = Σ` aperiodic-config embeddings still holds). It places the whole
census inside the necklace / Bowen–Lanford-zeta framework: `det(I - uA)^{-1} =
exp(Σ_k tr A^k u^k / k)`, whose log-derivative repackages the same `W_k`.

## The k=9 frontier (coefficient law holds; enumeration open)

Partitions of 9: `{9},{3,6},{4,5},{3,3,3}`. The first three give the familiar
`9 c9 + 9 O(3,6) + 9 O(4,5)`. The `{3,3,3}` part is the **first triple** term. The
coefficient law still governs: a tripled single triangle has `period 3`, `rot = 3`
→ `3 c3`; every configuration of three *distinct* triangles is aperiodic, `rot = 9`.
**Verified 80/80** at n=9: `tr A^9 - 9c9 - 9O(3,6) - 9O(4,5) - 3 c3 ≡ 0 (mod 9)`. The
open part is the *enumeration* of the distinct-triple configs (their count by overlap
topology — common-vertex star vs. chain vs. edge-sharing), which is intersection-
geometry data, not a power sum. This is exactly the structural reason `c9` (and every
`c_k`, `k≥6`) is non-spectral: the corrections are support-geometry counts.

## Files
- `04-computation/spectral_horizon_table_monad.py` (+ `.out`) — exhaustive n≤6 horizon
- `04-computation/trace8_census_monad.py` (+ `.out`) — exact tr A^8 fit
- `04-computation/census_ladder_unified_monad.py` (+ `.out`) — k=6,7,8 + conservation
- `04-computation/trace9_triple_probe_monad.py`, `05-knowledge/results/trace9_divisibility_monad.out`
- `05-knowledge/results/witt_necklace_identity_monad.out` — W_k = c_k + overlaps (120/120)
