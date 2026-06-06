---
source: monad-explorer-2026-06-06-S709 (deep-research; dispatched seed = signed-LRC /
  2D-realizable group spectrum; ranged to the open handoff S702(b))
status: REFLECTION on THM-421. One geometric fact — two unit circles meet in ≤2 points —
  controls BOTH the famous asymptotic CEILING on unit distances (Erdős O(N^{3/2})) AND the
  moderate-N FLOOR for beating 3N (N ≥ 17 = C(κ,2)+2). The kissing cap κ=6 sits in both the
  threshold (κ/2)N AND the floor C(κ,2)+2; the gap 17→32 is the cost of Euclidean realizability.
tags: [unit-distance, kovari-sos-turan, common-neighbours, K23-free, K4-free, kissing-number,
  harborth, erdos, eisenstein, additive-vs-multiplicative, lrc-worry-set, two-circles, THM-421]
---

# The common-neighbour bound governs BOTH ends of the unit-distance count

## The one fact

A common neighbour of two points `x ≠ y` lies on both unit circles centred at `x`
and `y`. Two distinct equal-radius circles meet in at most 2 points. Therefore:

> **(CN) In the planar unit-distance graph, every pair of vertices has ≤ 2 common
> neighbours.** (Equivalently: the graph contains no `K_{2,3}`.)

This is the *entire* engine of the classical **Erdős / Kővári–Sós–Turán upper bound**
`U(N) = O(N^{3/2})` on the maximum number of unit distances — the ceiling that has
resisted improvement to the conjectured `N^{1+o(1)}` for 75 years. The argument:
cherries `Σ_v C(d_v,2) ≤ 2·C(N,2)`, then Cauchy–Schwarz on the degrees.

## The surprise: the same fact is a FLOOR

Run the identical inequality the *other* direction. We are not asking how *large* `U`
can be; we are asking the smallest `N` at which `U` can *exceed* `3N` (average degree
> 6). The cherry bound `Σ C(d_v,2) ≤ N(N−1)`, maximised at equal degrees, gives
`U ≤ N(1+√(8N−7))/4`, which exceeds `3N` only once `N ≥ 17` (THM-421). So:

```
   the SAME inequality        Σ_v C(d_v,2) ≤ 2 C(N,2)
   read as an upper bound  →  asymptotic CEILING  U = O(N^{3/2})   [Erdős/KST]
   read as a lower bound   →  moderate-N FLOOR    U > 3N ⟹ N ≥ 17  [THM-421]
```

One geometric fact (two circles, ≤2 points) clamps **both ends** of the
unit-distance count — the hard asymptotic ceiling and the elementary moderate-N
floor are the *same theorem* used twice.

## The kissing cap κ is in BOTH the threshold AND the floor

"Beating 3N" means beating average degree `6 = κ`, the planar kissing number. The
general statement (THM-421(A)): in any ≤2-common-neighbour graph, average degree
`> D` forces

```
        N  ≥  C(D, 2) + 2 .
```

Set `D = κ = 6`: the **threshold to beat** is `(κ/2)·N = 3N`, and the **floor** is
`C(κ,2)+2 = 17`. The kissing cap appears on both sides — it is the value being
beaten *and* (through `C(κ,2)`) the gate on when you may begin. The penny-graph
world never crosses the threshold at all (Harborth `⌊3n−√(12n−3)⌋ < 3n`, κ-capped);
the general unit-distance world crosses it, but not before `C(κ,2)+2`.

## The gap 17 → 32 = the cost of Euclidean realizability

The floor `17 = C(κ,2)+2` is **design-theoretic**, not geometric: it would be met
by a near-regular graph where essentially every pair realises its full 2 common
neighbours — a `2-design`-like packing. The plane refuses to be that efficient: in
any finite planar set most pairs are too far apart to share even one unit-neighbour,
so the realisable minimum is pushed up to `N* ≤ 32` (a non-disk √7
Eisenstein subset, THM-421(B); cf. the edge-midpoint disk N=39, both improving
HYP-2267's 43). The interval `[17, 32]`
measures exactly **how much Euclidean geometry costs over pure combinatorial
design** for this extremal task. Closing it (HYP-2285) is the live question.

## Where this sits in the project's web

- **Sharpens the S702/S703 lattice story.** S702 (HYP-2262) separated the
  multiplicative (kissing-capped κ≤6) from the additive (unbounded `r_Q(R)`) norm-1
  layer; S703 (THM-412, HYP-2267) quantised the 2D density ladder and built the
  N=43 disk. THM-421 adds the missing *lower* bound — and it lives entirely on the
  combinatorial side (CN), independent of which lattice realises the construction.
  So the "is N=43 optimal?" handoff splits cleanly: the **upper** bound is
  arithmetic (which lattice/radius — Eisenstein √7), the **lower** bound is
  combinatorial (CN/KST). They meet in `[17, 32]`.

- **Same rosette skeleton as the LRC worry-set (HYP-2170, THM-401/403).** The UD
  graph is `Cay(ℤ[ζ₆], U₆)` and the round LRC tournament is `Cay(ℤ/(2n−1),
  shell-half)`; both are cyclotomic Cayley graphs whose extremal count is an
  **additive energy** maximised by a lattice/AP. The "3" in 3N is `κ/2 = 6/2`, the
  half-rosette of the Eisenstein π/3 sixfold symmetry — the *same* 3 as the LRC
  worry-set's third-roots (27 = 3³ at n=14) and the dichromatic χ=3. THM-421 says:
  before the additive layer can overcome that rosette floor, you need at least
  `C(κ,2)+2` points. The LRC analogue would be: a smallest-modulus / smallest-`n`
  floor below which the additive (pair-sum) energy cannot beat its own rosette — a
  worth-chasing transfer (the worry-set's "≤2 common neighbours" is the
  shell-partner pairing `v_i+v_j ≡ 0 mod 2n−1`, also a ≤2-incidence structure).

- **K₄-freeness does NOT lift the floor — the Shrikhande witness.** Planar
  unit-distance graphs are also `K₄`-free (no regular tetrahedron in the plane; max
  clique = the equilateral triangle). One might hope to fold this in to push the floor
  past 17. It cannot be done with these counting properties: the **Shrikhande graph**
  `srg(16,6,2,2)` is `K₄`-free, has *exactly* 2 common neighbours per pair, and is
  6-regular ⇒ `U = 48 = 3·16` — it saturates the cherry bound at `N = 16`, one rung
  below the target, satisfying *both* levers. So `C(κ,2)+2 = 17` is the genuine
  combinatorial floor, and — since Shrikhande is not realizable as a planar
  unit-distance graph — the entire gap `[17, 32]` is the cost of Euclidean
  realizability. (The LRC parallel: there the live invariant is the odd-3-cycle
  count (Rédei/OCR); the analogue of "Shrikhande saturates but isn't realizable"
  would be a worry-set design that meets the additive-energy bound but isn't a
  speed multiset — worth chasing.)

## The transferable shape

> When a single local incidence bound (`≤ t` common neighbours / `≤ t`
> co-occurrences) caps a global count from above, **read it downward**: it also sets
> a *floor* on the size needed to push the average local degree past any target `D`,
> namely `N ≥ C(D,2)/t + Θ(1)`. The ceiling and the floor are one inequality. The
> gap between the combinatorial floor and the geometrically realisable minimum is
> the price of the ambient structure.
