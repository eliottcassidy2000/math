# The j=7 u-escape closes by a SLOPE dichotomy — equal-slope collapses to LRC(≤8), non-equal-slope breaks the tiling at first order

**death-star-2026-07-12-S15.** Companion to THM-721 (Parts 2+3 now kernel-pure in
`LRCUEscape.lean`, this session). THM-721 Part 5 left the compressed `j ≥ 7` stratum open:
the u-escape union bound gives `reach₂(W) ≥ min(M_P, 1/(2j))`, which at `j = 7` is exactly
`1/14` — the wall, not strict. This reflection records the mechanism that makes `j = 7`
STRICT, found by asking *when can the union bound be tight?*

## The objects

Profile `W = ((kᵢ, bᵢ))`, impure set `F` (`bᵢ ≠ 0`, `|F| = j = 7`), pure set `P`
(13 − j = 6 distinct lifts, `M_P ≥ 1/7` by LRC(7)). For `s` in the pure good set
`S_good = {s : min_P ‖k_p s‖ ≥ 1/13}` define

    γ(s) = max_u min_{i∈F} ‖kᵢ s + bᵢ u‖   (the u-escape value at s).

Union bound + endpoint argument: `γ(s) ≥ 1/14` for EVERY s (7 forbidden systems, each of
total measure `2c`; at `c = 1/14` total = 1, and a finite union of open arcs of total
measure 1 always misses a point). Tightness is rigid:

> `γ(s) = 1/14` ⟺ the 7 closed forbidden systems at `c = 1/14` tile the circle exactly
> (pairwise interior-disjoint, endpoints touching).

System `i` consists of `|bᵢ|` arcs of width `1/(7|bᵢ|)` whose centers move, as `s` varies,
at velocity `−kᵢ/bᵢ` — **the slope of runner i**. This is the whole story.

## The dichotomy

**(a) All slopes equal** (`kᵢ/bᵢ = p/q` for all `i ∈ F`, lowest terms, so
`(kᵢ, bᵢ) = tᵢ·(p, q)` with distinct nonzero integers `tᵢ`). Then
`kᵢ s + bᵢ u = tᵢ(ps + qu)`, and as `u` sweeps `[0,1)` the variable `w = ps + qu` sweeps the
whole circle (`q ≥ 1`). Hence for EVERY s

    γ(s) = max_w min_i ‖tᵢ w‖ = M({|tᵢ|}) ≥ 1/8,

by **LRC(≤8)** — at most 7 distinct multipliers (± collisions only shrink the count).
The s-freedom is not needed; take the pure-optimal `s*` (`M_P ≥ 1/7`) and the atom gives
`M(V) ≥ min(1/7, 1/8) − B/(2L) = 1/8 − B/(2L) > 1/14` once `L > (28/3)·B`.
The degenerate direction of the tiling analysis is not a hard case at all — it is the
LOOSEST case. (At `j ∈ [8,12]` the same collapse gives `γ ≥ 1/(j+1) ≥ 1/13` — equal-slope
strata close at every `j ≤ 12`; `j = 13` equal-slope is a global dilate, excluded by
primitivity at `L ≥ 2`.)

**(b) Slopes not all equal.** If `γ(s) = 1/14` on a positive-measure subset of an interval
of `S_good`, then (finitely many combinatorial tiling patterns; each pattern's
endpoint-matching equations are affine in `s`) some pattern persists on a sub-interval, and
differentiating the endpoint matches forces adjacent arcs to move at equal velocity —
around the circle this forces ALL slopes equal. Contradiction. So tiling-`s` are isolated:
**every good interval contains `s` with `γ(s) > 1/14` strictly**, and the margin is
first-order quantitative: moving `ε` off a tiling point grows the max circular gap at rate
`= max adjacent slope difference Δρ ≥ 1/B²` (rates are differences of `−kᵢ/bᵢ` and sum to
zero around the circle, so some gap grows):

    γ(s₀ + ε) ≥ 1/14 + ε·Δρ/2  (verified exactly on the near-dilate adversary:
    lifts {7..13}, b ≡ 1, s₀ = 1/7 — γ(1/7) = 1/14 EXACTLY (the u-side AP: phases
    {7s..13s} = staggered 1/7-AP), and γ(1/7 + ε) = 1/14 + ε/2 exactly, the wrap gap
    donating −6ε to the six growing gaps).

The available motion is `ε ≤ headroom/max_P k_p` (pure headroom `1/7 − 1/13 = 6/91`), so the
per-family strict floor is `1/14 + δ(W)` with `δ(W) ≳ Δρ·headroom/(2·max_P k_p)` — positive
for every fixed profile, and beating the Lipschitz loss `B/(2L)` for `L` above an explicit
per-profile threshold (near-dilates: `L³ ≳ B³·Vmax`-scale, comfortably inside the
large-diameter leg).

## Why this is the same object as the residue law (fleet convergence)

At rational `s = a/b'` with unit frequencies (`b ≡ 1`) the phases `{kᵢ a/b'}` are a residue
orbit and `γ(a/b') = (max circular residue gap)/(2b')` — mac-mini-S66's MAX-GAP RESIDUE LAW
(`maxgap{frac(eᵢx)} = G_E(a,b)/b`), one rung down (profile level instead of family level).
The s-motion strictness above is the DERIVATIVE form of "residue-band dodging": the good
period exists because the band cannot follow the sweep. Same physics, two scales — worth
keeping in view for THM-527-A's large-spread half.

## The honest residual after this

- `j = 7`: closed in mechanism (this note); needs write-up as THM-721 Part 6 + the
  quantitative `δ(W)` lemma. The equal-slope leg (a) is complete as stated (pure LRC(≤8)
  citation + the `w`-substitution — two lines, Lean-friendly, same shape as
  `margin_uescape_j6`).
- `j ∈ [8,12]` equal-slope: closed by (a) at floor `1/13 − B/(2L)`.
- `j ∈ [8,13]` NOT-equal-slope at every admissible scale: still open for this atom (union
  bound below `1/14` before strictness even starts) — but this stratum is strongly
  incoherent (8+ of 13 runners off-lattice at every scale), the home turf of the
  pair-sum/Parseval certificates (klein-S264, THM-668/680). The compressed lane's gap is
  now `[8,13]`-impure-with-mixed-slopes only.

Probe: `04-computation/lrc14_uescape_j7_boundary_deathstar_S15.py` (+ `.out`) — exact
Fraction arithmetic: the adversary rate law, adversarial profile search (does anything pin
`sup_{good s} γ = 1/14`?), the equal-lift corner `inf_c max_u min ‖c + bᵢu‖`, and DC
inhabitation of the `j = 7` stratum.

→ THM-721, HYP-6270, LRCUEscape.lean, LRCDecorrelation13.lean, mac-mini-S66 (residue law),
klein-S264 (Parseval side), klein-S152/HYP-4711 (the old all-impure candidate — subsumed on
the equal-slope part by (a)).
