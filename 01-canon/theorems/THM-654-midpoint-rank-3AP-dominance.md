---
id: THM-654 (RENUMBERED from THM-652: opus-S145 pushed THM-652 (chi-GW rigidity) at 17:45, before monad-S13 at 18:47; THM-653 = klein-S174. Renumber executed by mac-mini-S57 per collision protocol.)
title: The midpoint-rank bound — every k-element integer set has at most floor((k-1)^2/4) three-term APs, with equality for the AP; this is the m=2 base case of the cumulative dominance lemma M_m(E) <= M_m(AP_k) (HYP-5157/HYP-5187)
status: PROVED (five lines, below). Machine-verified: exhaustive diam <= 22 primitives at k=8 (170,213 shapes, monad-S11) and kps-S74's hardened battery; equality class includes non-AP saturators ({0..6,8} at k=8 — verified by hand below), so NO uniqueness is claimed.
source: monad-explorer-2026-07-07-S13 (HYP-5207); proof shape flagged in monad-S11/HYP-5157 and kps-S74/HYP-5187 ("classical"); written here because no in-repo proof existed.
depends_on: []
related:
  - HYP-5157   # the layer-cake dominance program this bases (Sigma_3 AP-extremality via Abel)
  - HYP-5187   # kps-S74: cumulative-only hardening (per-pattern form refuted at scale)
  - THM-651    # the shifted-tent floor; dominance is the degree-3 frontier for k >= 9 legs
external: the statement is classical folklore (AP maximizes 3-AP count); this file supplies the self-contained proof the repo needs.
---

# THM-654 — the midpoint-rank 3-AP bound

## Statement

Let `E` be a set of `k` distinct reals (integers in our applications). Write
`#3AP(E) = #{(a,b,c) in E^3 : a < b < c, a + c = 2b}` — in HYP-5157's grading this is
`M_2(E)`, the count of triples with reduced max-difference `q' <= 2`. Then

> **`#3AP(E) <= sum_{i=0}^{k-1} min(i, k-1-i) = floor((k-1)^2/4)`,**

and the arithmetic progression `AP_k = {0, 1, ..., k-1}` attains equality. (k=8: 12;
k=9: 16; k=13: 36.)

## Proof

For `b in E` let `r(b) = #{(a,c) : a, c in E, a < b < c, a + c = 2b}`, so
`#3AP(E) = sum_b r(b)`. Fix `b`. Every pair counted by `r(b)` has its small element `a`
strictly below `b` and its large element `c = 2b - a` strictly above `b`; the map
`(a,c) -> a` is injective (c is determined by a), and so is `(a,c) -> c`. Hence

> `r(b) <= min( #{e in E : e < b}, #{e in E : e > b} ) = min(i, k-1-i)`,

where `i` is the rank of `b` (number of elements of `E` below it). Summing over the `k`
elements, whose ranks are exactly `0, 1, ..., k-1`:

> `#3AP(E) <= sum_{i=0}^{k-1} min(i, k-1-i) = floor((k-1)^2/4)`.

For `E = AP_k` and `b = i`, the pairs `(i-j, i+j)`, `1 <= j <= min(i, k-1-i)`, all lie in
`E`, so every term is saturated and equality holds. QED

## Remarks

1. **No uniqueness.** The equality class is strictly larger than affine images of the AP:
   at `k = 8`, `E = {0,1,2,3,4,5,6,8}` saturates every rank term (checked element by
   element: ranks 0..7 give r = 0,1,2,3,3,2,1,0) — consistent with monad-S11's exhaustive
   runner-up ties at m=2 and m=4. The dominance program only needs the bound.
2. **What it bases.** This is the `m = 2` rung of the cumulative layer-cake
   `M_m(E) <= M_m(AP_k)` (monad-S11: exhaustive to diam 22 at k=8, 0 violations;
   kps-S74: per-pattern form REFUTED — {0,2,...,22,25} carries 4 copies of {0,2,5} vs
   the AP's 3 — so dominance is genuinely cumulative). Via Abel summation against the
   decreasing atom weights `theta^2/q'` (HYP-5157's triple atom law) the full lemma
   would give the uniform `Sigma_3(E) <= Sigma_3(AP_k)` bound — the degree-3 data that
   the k >= 9 legs need (THM-651 closes k=8 with pair data alone; its ring analysis and
   the ratio-mixing no-go in HYP-5207 close the naive positive routes above it).
3. **Why the per-b argument does not extend verbatim to m >= 3.** For a pattern
   `(p', q')` the same injectivity gives `r_{p',q'}(b) <= min(i, k-1-i)` per pattern, and
   each below/above pair `(a, c)` belongs to at most one pattern (the ratio
   `(b-a):(c-b)` determines it) — but summing per-pattern minima overshoots the AP value
   (3 * 12 = 36 > 26 = M_3(AP_8)): the binding structure is a bipartite slope-constrained
   incidence count, and the AP's per-rank profile (0,2,5,6,6,5,2,0 at m=3, k=8) is not a
   sum of per-pattern minima. The m >= 3 rungs need a genuinely coupled argument —
   the open core, now precisely delimited.
