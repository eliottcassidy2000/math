# THM-592: Radius-derivative structure of the lonely measure (slope formula, Farey breakpoint grid, convexity criterion, ladder bound)

**Status:** PROVED (elementary real analysis + exact-rational verification)
**Author:** mac-mini-2026-07-01-S93
**Verification:** `04-computation/lrc_radius_slope_frame_macmini_S93.py` → `05-knowledge/results/lrc_radius_slope_frame_macmini_S93.out` (exact `Fraction` arithmetic; slope formula matches in every linearity cell for AP {1..13}, deep well {1..12,182}, spread {1..12,40}; kink signs match the classification; convexity criterion matches observation on all three).

---

## Setting

For a finite set `S ⊂ Z_{>0}` and radius `r ∈ (0, 1/2)` define on `R/Z`:

- danger set `D_v(r) = {t : ||vt|| < r}` = `v` open arcs `((j-r)/v, (j+r)/v)`, `j = 0..v-1`, each of length `2r/v` (total measure exactly `2r`);
- lonely set `L_S(r) = (R/Z) \ ∪_{v∈S} D_v(r)` (closed);
- lonely measure `m_S(r) = |L_S(r)|`;
- covering-min `M(S) = max_t min_{v∈S} ||vt||`.

Note `m_S(r) > 0 ⟺ r < M(S)`-or-(r = M(S) and positive-measure argmax); `L_S` is nested decreasing in `r`. LRC(N) for `|S| = N-1` is `M(S) ≥ 1/N`, equivalently `m_S(r) > 0` for all `r < 1/N`.

## Statement

**(i) Piecewise linearity on the Farey-type grid.** `m_S(·)` is continuous and piecewise linear on `(0, 1/2)`; every breakpoint lies in the grid
```
B(S) = { d/(v+w) : v,w ∈ S, d ∈ Z_{>0} } ∪ { d/(w−v) : v<w ∈ S, d ∈ Z_{>0} }.
```

**(ii) Slope formula.** In the interior of a linearity cell, `L_S(r)` is a finite disjoint union of closed arcs whose left endpoint is a right endpoint of some `D_{v_L}` and whose right endpoint is a left endpoint of some `D_{v_R}`, and
```
m_S'(r) = − Σ_{components K of L_S(r)} ( 1/v_L(K) + 1/v_R(K) ).
```
In particular `−m_S'(r) ≤ 2·b_0(L_S(r))·max_v(1/v) ≤ 2·b_0(L_S(r))` and `−m_S'(r) ≤ 2|S|` (the union/tangent-at-0 slope; `m_S'(0^+) = −2|S|` exactly when all `2Σv` endpoints are exposed).

**(iii) Kink classification.**
- **Merge events** (`r = d/(v+w)`): a lonely component shrinks to a point and dies; `m'` jumps UP by `1/v + 1/w` → **convex kink**. Components never appear as `r` grows (nesting), so all death events are convex.
- **Overtaking events** (`r = d/(w−v)`, `v < w`): the exposed right (resp. left) endpoint of a component changes owner from the slow-endpoint runner `w` to the fast-endpoint runner `v` (endpoint speeds are `1/v > 1/w`); `m'` jumps DOWN by `1/v − 1/w` → **concave kink**. An overtaking event produces a kink only if the crossing happens on the exposed boundary of `L_S(r)`; interior (covered) crossings are measure-invisible.
- No other kink types exist below `r = 1/2`.

**(iv) Convexity criterion.** For a pair `v < w` the smallest overtaking radius is `gcd(v,w)/(w−v)`. Hence if
```
w − v ≤ gcd(v,w) / ρ   for all pairs v < w in S,
```
then `m_S` is **convex** on `(0, ρ]`. At `ρ = 1/N`: `w − v ≤ N·gcd(v,w)` for all pairs ⟹ convex below the LRC radius. (Sufficient, not necessary: exposure may fail.) Examples: AP `{1..13}` convex below `1/14` (max difference 12 < 14); deep well `{1..12,182}` has its first concave kink at exactly `1/181 = gcd(1,182)/(182−1)` (observed).

**(v) Ladder (tangent-with-defect) bound.** For `0 < r_0 < r_1 ≤ 1/2`,
```
m_S(r_1) ≥ m_S(r_0) + ( m_S'(r_0^+) − K(r_0, r_1) ) · (r_1 − r_0),
```
where `K(r_0,r_1) = Σ (1/v − 1/w)` summed over the (exposed) overtaking events in `(r_0, r_1)`. *(S94 sharpening: weighting each event by its remaining distance, `m_S(r_1) ≥ m_S(r_0) + m_S'(r_0^+)(r_1−r_0) − Σ_kinks Δ(ρ)·(r_1−ρ)`, is strictly better and equally immediate from the BV integration; the unweighted form treats all kinks as if at the anchor. S94 defect sieve: generic covering 13-sets have K(1/16, 1/14) = 0 — their overtaking resonances are covered, not exposed — and the ladder from anchor 1/16 alone certifies m(1/14) > 0 for 22/25 sampled covering sets; the exposed defect concentrates on the deep-well family, `lrc_pairlaw_defect_arcradius_macmini_S94.out`.)* If there are none (e.g. under (iv)), this is the plain convexity tangent bound: a certified pair (measure, slope) at a sub-critical anchor radius pushes positivity of `m_S` up to the tangent root. Verified instance: for the deep well, the tangent at `r_0 = 1/16` alone certifies `m > 0` past `1/14` (i.e. `M(S) ≥ 1/14` follows from data at `1/16`).

## Proof

(i) Every danger endpoint is `(j ± r)/v`, affine in `r` with slope `±1/v`. The combinatorial type of the arrangement (order of endpoints, containment) changes only when two endpoints coincide: `(j + ε_1 r)/v = (j' + ε_2 r)/w` with `ε_i ∈ {±1}` gives `r = (j'v − jw)/(ε_1 w − ε_2 v)`, a positive rational with denominator `v+w` or `w−v` (same-runner coincidences require `r ≥ 1/2`). Between consecutive coincidence radii the exposed-boundary structure is constant and `m_S` is affine. ∎

(ii) On a cell, each lonely component is `[(j+r)/v_L, (j'−r)/v_R]` with derivative of length `−(1/v_L + 1/v_R)`; sum over components. ∎

(iii) Nesting `L(r') ⊆ L(r)` for `r' > r` forbids component birth. A component death is the collision of its two endpoints, closing speed `1/v_L + 1/v_R`, at a radius with denominator `v_L + v_R`; the slope loses that component's term, so `m'` increases: convex. An exposed-endpoint owner change is the collision of two same-direction endpoints, relative speed `1/v − 1/w`, at denominator `w − v`; after the crossing the faster endpoint (owner `v`, the smaller speed value) is the exposed one, so the component's shrink rate increases by `1/v − 1/w` and `m'` decreases: concave. Below `1/2` every breakpoint is one of these (or measure-invisible). ∎

(iv) Realizable coincidence values `j'v − jw` at fixed `(v,w)` are exactly the multiples of `g = gcd(v,w)` (Bézout, with `j, j'` ranging over full residue systems as arcs wrap), so the least overtaking radius is `g/(w−v)`; if it exceeds `ρ` for every pair there is no concave kink in `(0, ρ]`, and a piecewise-linear function whose kinks are all convex is convex. ∎

(v) `m'` is piecewise constant; on `(r_0, r_1)` it equals `m'(r_0^+)` plus the accumulated jumps, of which the negative ones total `−K`; integrate the minorant `m'(r) ≥ m'(r_0^+) − K`. ∎

## Why this matters (the assumption it corrects)

The union bound `m_S(r) ≥ 1 − 2|S|r` (vanishing at `1/(2|S|)`, far short of `1/N`) is exactly the **tangent line at r = 0** of `m_S`. The first-moment/SOS "stall" (HYP-3791, HYP-3822: global moment LP provably gives `min m_0 = 0`) is the statement that the tangent at 0 cannot see the critical radius — not that measure methods fail. Under convexity, **every later anchor gives a better tangent**, and the method of HYP-3824 (Γ₀-localization) is re-read as moving the anchor. The obstruction to riding the tangent all the way to `1/N` is precisely the set of exposed overtaking resonances `d/(w−v)` — a *sparse, enumerable, arithmetic* defect list (the "atoms"), not an analytic wall. Slope/measure certificates at sub-critical anchors + a bound on the exposed-resonance defect `K` = a new route to covering floors.

Consumers: THM-579 (CV(N_R)² covering floor), THM-580 (2-adic descent), OPEN-Q-108 fattening lemma, HYP-3824 (the S92 linear-slope observation, now exact: see THM-593).

-> HYP-3840, THM-593, HYP-3824, HYP-3822, HYP-3791, THM-523, OPEN-Q-108.
