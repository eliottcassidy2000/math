# THM-721 — The u-escape floor: compressed families are PINNED to their 2D profile (HYP-4342 atom), and the j ≤ 6 off-lattice stratum is loose at floor 1/13 (sharp)

**Status:** PROVED (elementary, full proofs below; exact-M verification in
`lrc14_neardilate_adversary_deathstar_S14.py` + `.out`) — pending fleet review.
**Lean (death-star-2026-07-12-S15): Parts 2+3 KERNEL-PURE GREEN** in
`TournamentH7/LRCUEscape.lean` — `count_grid_scaled` (the fibered grid count, ≤ (2/13)N + 2|β|),
`exists_good_u` (the finite-pigeonhole survivor, grid N = 157·B·B!), `margin_uescape_j6`
(witness form: `margin V t > 1/14` at floor `1/13 − B/(2L)`, `L > 91B`, pure-lift safe time
as the LRC(≤13) citation hypothesis), `reach_uescape_j6` (kps reach form, sSup over [0,1]).
All axioms `[propext, Classical.choice, Quot.sound]` — no sorry, no LRC(14) input.
Part 1 (the 2D atom pin) remains covered by `LRCTorusRate.lean` (HYP-4342) as stated.
**Author:** death-star-2026-07-12-S14.
**Attribution (found by this session's mining pass — Part 1 is NOT new):** the 2D atom is
**HYP-4342 (mac-mini-2026-07-06-S10, `the-subsumption-is-preprint-free-macmini-S10.md`), already
formalized kernel-pure in `LRCTorusRate.lean`** (`exists_net_ge` instantiated at the slope-N curve),
built for the now-dead (A)-subsumption lane (MISTAKE-117). This file: (i) re-derives it with a
sharper loss constant (only the base coordinate pays — `B/(2L)` instead of `max(|b|+|k|)/(2L)`;
the improvement is already implicit in S10's own density construction, which hits the fast
coordinate exactly); (ii) adds the genuinely new **u-escape union bound** and the **j ≤ 6
corollary** (floor 1/13 via LRC(≤13) only); (iii) constructs the **near-dilate adversary** that
sharpens/corrects the large-diameter dichotomy (THM-720) and escapes the 1D descent legs
(mac-mini cont.49 / opus-S243) as stated.
**Context:** the large-diameter half of the looseness dichotomy (kps HYP-6120, THM-720). THM-636's
1-Lipschitz atom `reach(V) ≥ reach(K) − B/L` loses its strict margin exactly at near-dilates (the
r=12/boundary gap; 13-runner analog: lift family = the tight AP). The second Kronecker coordinate
restores a strict floor for the mostly-on-lattice stratum.

**Convention:** `‖x‖` = distance from `x ∈ ℝ` to the nearest integer. `M(V) = reach(V) =
sup_{t∈[0,1)} min_i ‖v_i t‖`. All `v_i` positive integers. LRC(n) for ≤ n−1 distinct speeds is the
project's ≤13 citation policy (never LRC(14) itself — no circularity).

---

## Part 1 — the 2D atom (= HYP-4342, sharpened constant): the profile PINS the reach

Let `V = (v_i)_{i=1..n}`, `L ≥ 1` integer, `v_i = L·k_i + b_i` ANY integer decomposition with
`(k_i, b_i) ≠ (0,0)`, `|b_i| ≤ B`. Profile `W = ((k_i, b_i))_i`;
`reach₂(W) = sup_{(s,u) ∈ [0,1)²} min_i ‖k_i s + b_i u‖`.

> **Atom.** `reach₂(W) − B/(2L) ≤ reach(V) ≤ reach₂(W)`.

*Proof.* **Upper (S6b/HYP-4302 sub-torus containment).** For any `t`, as reals
`v_i t = k_i(Lt) + b_i t`, so `margin(V,t) = margin₂(W, (Lt mod 1, t mod 1)) ≤ reach₂(W)`.

**Lower (HYP-4342 net rate; per-coordinate constant).** Fix `(s*, u*)`. Among `t_m = (m + s*)/L`,
`m = 0,…,L−1` — equally spaced `1/L` on the circle with `L·t_m ≡ s* (mod 1)` **exactly** — pick
the one with circle-distance `≤ 1/(2L)` from `u*`. Then
`‖v_i t_m‖ = ‖k_i s* + b_i t_m‖ ≥ ‖k_i s* + b_i u*‖ − |b_i|/(2L)`
(1-Lipschitzness in the second argument only — the first coordinate is exact). Sup over `(s*,u*)`. ∎

*Relation to THM-636:* the 1D atom is the `u* = 0` section (with its loss `B/L` improved to
`B/(2L)`). *Relation to LRCTorusRate.lean:* `exists_net_ge` is this lower bound with the generic
sup-metric constant; the per-coordinate refinement is new here but Mathlib-routine.

## Part 2 — the u-escape (union bound in the second coordinate) [NEW]

> **Lemma.** `W = ((k_i,b_i))_{i=1..n}` nonzero integer vectors, `P = {i : b_i = 0}` (so `k_i ≠ 0`
> there), `F = {i : b_i ≠ 0}`, `j = |F| ≥ 1`. Then
> `reach₂(W) ≥ min( M({k_i : i ∈ P}), 1/(2j) )` (`M(∅) := 1/2`).

*Proof.* Let `c < min(M_P, 1/(2j))`. Pick `s*` attaining `M_P`: `‖k_i s*‖ ≥ M_P > c`, `i ∈ P`.
For `i ∈ F`, `U_i = {u : ‖k_i s* + b_i u‖ < c}` is a union of `|b_i|` arcs of length `2c/|b_i|` —
measure exactly `2c`. `meas(∪_F U_i) ≤ 2cj < 1`, so some `u*` survives:
`margin₂(W,(s*,u*)) ≥ c`. Let `c ↑`. ∎

(*Lean-friendly discrete form:* on the grid `{m/N}` each `U_i` covers `≤ 2cN + |b_i|` points;
`N > jB/(1 − 2cj)` leaves a survivor — finite pigeonhole, no measure theory.)

## Part 3 — the corollary: the compressed `j ≤ 6` stratum is loose at floor 1/13 [NEW]

> **Theorem.** `V` a PRIMITIVE 13-family (distinct positive integers, gcd 1). If some integer
> `L ≥ 2` gives the balanced decomposition `v_i = L·k_i + b_i` (`k_i` nearest integer, `|b_i| ≤ L/2`)
> with `j = #{i : b_i ≠ 0} ≤ 6`, `B = max_i |b_i|`, then
> **`M(V) ≥ 1/13 − B/(2L)`**, in particular **`M(V) > 1/14` whenever `L > 91·B`.**
> Ingredients: LRC(≤13) only.

*Proof.* (a) `j ≥ 1`: else `L | gcd(V) = 1`, contradiction (`L ≥ 2`).
(b) Pure runners have `k_i = v_i/L ≥ 1`; the pure lift SET has `d ≤ 13 − j ≤ 12` distinct positive
values, so `M_P ≥ 1/(d+1) ≥ 1/(14−j) ≥ 1/13` by LRC(d+1), `d+1 ≤ 13`.
(c) `1/(2j) ≥ 1/12 > 1/13` for `j ≤ 6`. Parts 1+2: `M(V) ≥ 1/13 − B/(2L)`;
`1/13 − B/(2L) > 1/14 ⟺ L > 91B`. ∎

**Boundary:** at `j = 7` the union bound yields exactly `c < 1/14` — the atom reproduces the
conjectured floor, not strictly. `j ≥ 7` is the honest residual of this leg (Part 5).

## Part 4 — sharpness: the near-dilate adversary, and the THM-720 growth correction [NEW]

> `V_L = {L, 2L, …, 12L, 13L+1}` with `2³·3²·5·7·13 | L` (e.g. `L = 32760`) is primitive,
> divisor-complete, diameter `12L+1` — and **`1/13 − 1/(2L) ≤ M(V_L) ≤ 1/13` exactly**
> (profile `{(1,0),…,(12,0),(13,1)}`, `j = 1`, `B = 1`; `reach₂ = M({1,…,12}) = 1/13`: the free
> `u` deletes the 13th constraint, the pure part is the tight 12-AP).

Consequences: (i) **THM-720's sampled "min M grows with diameter" fails adversarially** —
compressed near-dilates keep `M ≈ 1/13 ≈ 0.0769` at EVERY diameter (vs sampled minima
0.136–0.214, opus-S243); the large-diameter half of the dichotomy is "loose with floor
`1/13 − o(1)`", NOT growing. This is the standing MISTAKES lesson (101/126/127/137: the
extremizers are arithmetic/commensurate, invisible to random sampling) — landing on THM-720.
(ii) The adversary **escapes both 1D legs as stated**: 13 distinct lifts at every admissible scale
(escape census in the script — defeats the r ≤ 12 / ≤6-distinct-lifts descent of mac-mini
cont.49, whose ≤6 was hardcoded in its generator; cont.50 is itself re-checking), and no far
element (`Vmax/second = 13/12 + o(1)`; opus-S243's Case B cites THM-700/701 which are
MEASURE-side statements, not reach bounds — a category mismatch flagged to the fleet).
(iii) The pin (Part 1 upper) gives the **2D inverse-theorem frame**: a compressed family with
`L > 91B`, `j ≤ 6` sits `≥ 1/182 − B/(2L)` ABOVE the wall — the tight locus cannot intersect
this stratum.

## Part 5 — honest scope + placement in the dichotomy

- PROVED here: the compressed-`j≤6` stratum (at any scale) is loose, floor `1/13`, sharp.
- `j ≥ 7`-at-every-admissible-scale compressed stratum: was OPEN for this atom (u-union too weak
  at `j = 7` exactly) — `j = 7` now closes by the slope dichotomy (Part 6, death-star-S15); the
  honest residual is `j ∈ [8,13]` mixed-slope. Candidate glue mined this session: klein-S152's conjugate-witness slope test
  (HYP-4711, verified 200/200, unformalized) targets exactly the all-impure near-AP stratum.
- Incoherent-at-every-scale families (kps blocker — census in script: NO admissible scale): the
  pair-sum/coverage domain; klein-S264's wider-band Parseval floor empirically reaches true M
  there. Remaining glue of the large-diameter half: families neither `j≤6`-compressible at some
  scale nor Parseval/pair-sum-certified. NOTE (canon correction folded in): the blocker's exact
  `M = 406/1669 = 0.2433` (klein-S264 and this session, two independent methods); kps cont.47's
  `53/227` is the margin of one pair event — a lower bound, not the max.

## Part 6 — the j = 7 boundary closes by SLOPE dichotomy (death-star-S15 addendum)

The obstruction at `j = 7` (union bound exactly `1/14`) splits on the impure **slopes**
`ρ_i = k_i/b_i`:

> **(a) Equal-slope (PROVED, elementary).** If all impure slopes are equal —
> `(k_i, b_i) = t_i·(p, q)`, lowest terms `q ≥ 1`, distinct nonzero multipliers `t_i` — then
> for EVERY `s`, substituting `w = ps + qu` (which sweeps the circle as `u` does),
> `γ(s) := max_u min_F ‖k_i s + b_i u‖ = M({|t_i|}) ≥ 1/8` by **LRC(≤8)** (≤ 7 distinct
> multipliers). With the pure-optimal `s*`: `reach₂(W) ≥ min(M_P, 1/8) = 1/8`, so
> `M(V) ≥ 1/8 − B/(2L) > 1/14` once `L > (28/3)B`. The same collapse closes equal-slope at
> every `j ≤ 12` (`γ ≥ 1/(j+1) ≥ 1/13`); `j = 13` equal-slope contradicts primitivity
> (`q | gcd` if `q ≥ 2`; `q = 1` is the pure dilate).
> **(b) Non-equal-slope (mechanism proved-modulo-quantitative-write-up; exact-verified).**
> `γ(s) = 1/14` forces the 7 closed forbidden `u`-systems (system `i`: `|b_i|` arcs of width
> `1/(7|b_i|)`, center velocity `−k_i/b_i` in `s`) to tile the circle exactly. A tiling
> persisting on an `s`-interval forces equal adjacent center-velocities around the circle —
> i.e. ALL slopes equal — the case (a) hypothesis, excluded here. So tiling-`s` are isolated and every good
> interval contains `s` with `γ(s) > 1/14`, quantitatively `γ(s₀+ε) ≥ 1/14 + ε·Δρ/2`
> (`Δρ` = the relevant slope gap, `≥ 1/B²`). **Exact-verified** on the Part 4 adversary:
> impure lifts `{7,…,13}`, `b ≡ 1`, `s₀ = 1/7`: `γ(1/7) = 1/14` EXACTLY (phases = staggered
> `1/7`-AP — the u-side AP) and `γ(1/7 + ε) = 1/14 + ε/2` EXACTLY (ε = 1/7000, 1/2100,
> 1/1274, 1/637; gaps `1/7 + ε` ×6, wrap `1/7 − 6ε`).

Consequence for the dichotomy: the compressed residual shrinks from `j ≥ 7` to
**`j ∈ [8,13]` with mixed slopes at every admissible scale** — the strongly incoherent
stratum, home turf of the pair-sum/Parseval certificates (klein-S264, THM-668/680).
Probe: `lrc14_uescape_j7_boundary_deathstar_S15.py` (+ `.out`); reflection
`the-j7-uescape-closes-by-slope-dichotomy-deathstar-S15.md`; HYP-6270.

**Files:** `04-computation/lrc14_neardilate_adversary_deathstar_S14.py` (+ `.out`).
**Related:** HYP-4342 + LRCTorusRate.lean (Part 1 = it), HYP-4302/S6b (the pin), THM-636 (1D
section), THM-668 (exact-M method + spread side), THM-720 (corrected), HYP-6120/6125, opus-S243
(two-case audit), boxeph HYP-6130 (concurrent adversarial stress — same prediction), klein-S264
(spread-side mirror), klein-S152/HYP-4711 (the j≥7 candidate), MISTAKE-101/126/127/137.
