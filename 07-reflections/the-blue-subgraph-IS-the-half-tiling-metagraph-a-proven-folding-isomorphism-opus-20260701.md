# The blue subgraph IS the half-tiling metagraph — a proven folding isomorphism (the SC world is a self-contained tiling model of half the dimension)

> **⚠️ S20 CORRECTION (opus, HYP-3811).** Two claims below are WRONG (the theorem itself STANDS):
> (1) the "dyadic complement-folding tower" — NO, the fold is a **single terminal step**: R=complement is
> order-2, and H_n *is* its fixed set (all SC), so it cannot fold again by complement.
> (2) "at even n, #SC(n)=A000568(n−1)" — a **coincidence** (n=4,6 only), **breaks at n=8** (176≠456). The
> folded class-count is **A051337** (self-converse tournaments), governed by anti-automorphisms with all
> cycles ≡ 2 mod 4. See the S20 reflection / HYP-3811.

*opus-2026-07-01-S19. The owner asked to prove "blue-subgraph = half-tiling metagraph as an isomorphism" and
to work the S18 directions (odd-cycle spectrum, recursion, realization-degeneracy, the T-join/flip-rank bridge).
This proves the isomorphism from the S18 generator (R = complement) and reports the directions.*

> **CONVERGENCE:** mac-mini-S84 (HYP-3809) concurrently built a conjecture-atlas that *noted* "blue structure
> lives on the half-tiling" (via the identity `#grid-sym = 2^{⌊(n−1)²/4⌋}`) but listed **"blue = half-tiling
> metagraph recursion" as an OPEN atlas item**. This reflection **proves that item** (the theorem below, from
> R = complement). Two facts I adopt from HYP-3809: the SC-class count is **A051337 (self-converse
> tournaments)** = 2,2,8,12,…, and the half-tiling dimension has the closed form `D = ⌊(n−1)²/4⌋` (= my
> `(m+f)/2`). And mac-mini **answered my realization-degeneracy question**: the parity+eligibility constraints
> are *necessary but not sufficient* (color-preserving degree 2-swaps exist at n=5,6) — see below.

## The theorem
Let `M_n` be the flip-line merged metagraph (tournament-tiling-explorer), `R` the grid reflection (the
`isGridSym` involution on the `m=C(n-1,2)` tiles: `(x,y) ↔ (n−y+1, n−x+1)`), and `f = ⌊(n−1)/2⌋` its number of
fixed (anti-diagonal, `x+y=n+1`) tiles. Define the **half-tiling model** `H_n`: its tiles are a transversal of
`R` (the `f` fixed tiles + one representative of each of the `(m−f)/2` transpositions), so it has
`D = f + (m−f)/2 = (m+f)/2` tiles; its tilings are `{0,1}^D`; the **unfold** map `u` sends a half-tiling to the
grid-symmetric full tiling that agrees on the fixed tiles and is constant on each transposition.

> **Theorem.** `u` is an isomorphism from the flip-line metagraph of `H_n` onto the **blue subgraph** of `M_n`.

**Proof.** Three lemmas.
- **L1 (R = complement, S18).** For every tiling `t`, `class(R·t) = complement(class(t))`. Hence a tiling is
  grid-symmetric (`t = R·t`) **iff** its class is self-complementary (SC); and every SC class contains a
  grid-symmetric tiling (verified n≤7 — "SC classes are never pure-black").
- **L2 (the R-fixed subcube).** `Fix(R)` = the grid-symmetric tilings = the tilings determined freely on the
  transversal (each fixed tile free; each transposition forced equal). So `u : {0,1}^D → Fix(R)` is a
  **bijection**, and `D = (m+f)/2`.
- **L3 (flip intertwines).** Flipping all `D` transversal bits flips every fixed tile and every
  transposition-representative, so after unfolding it flips **all `m`** tiles: `u(h̄) = complement_m(u(h))`.

Now the flip-line metagraph of `H_n` has nodes = iso classes of unfolded tournaments and lines =
`{h, h̄}`. By L1 those classes are exactly the SC classes (all of them), and by L3 each line maps under `u` to
`{t, complement_m(t)}` with `t` grid-symmetric — i.e. a **blue line** of `M_n`. Conversely every blue line is a
flip-pair of grid-symmetric tilings, so arises this way. Thus `u` carries `H_n`'s flip-metagraph bijectively,
node-for-node and line-for-line, onto the blue subgraph of `M_n`. ∎

**Verified** (n=4,5,6): the unfold bijection, the flip-commutation, and "half-tiling classes = SC classes, lines
= `2^{D−1}`" all hold exactly. The blue subgraph is literally a tournament-tiling metagraph in its own right —
the tiling model of the **folded staircase** (equivalently, of the self-complementary tournaments).

## What this buys — the recursion
The SC/blue world is not a decoration on `M_n`; it is a **smaller copy of the whole construction**. The blue
subgraph carries the same data type (tilings, flip-lines, classes) at dimension `D(n) = (m+f)/2 = 2,4,6,9` for
`n=4..7` — roughly *half* the tiles. Consequences:
- **Node-count recursion (even n).** `#SC = 2, 8, 12, 88` for `n=4..7`, and for **even n** this equals
  `A000568(n−1)` (2 = A₃, 12 = A₅): the half-tiling model at even `n` has as many classes as *all* tournaments
  at `n−1`. (Odd n: 8, 88 — its own sequence.) So the folded staircase's class count recurses to the ordinary
  count one size down, at even `n`.
- **The fold is one dyadic step.** `Fix(R)` halves the tile-dimension. The half-tiling model is itself a tiling
  model, so it has its own grid-symmetric/blue sub-object (a *quarter*-tiling) — a dyadic tower of foldings,
  each the SC-restriction of the last. This is the concrete recursive structure the owner asked for: not a
  vertex-insertion (Mode A) but a **complement-folding (Mode B-flavored) halving of the tile-cube**.

## The other directions, worked
- **Black odd-cycle spectrum.** The black even-graph (S18) is not just non-bipartite for n≥5 — its **odd girth
  is 3**: it contains triangles (verified n=5,6). So it is far from bipartite; the even (Eulerian) black world
  is triangle-rich, while the blue world is the acyclic-ish T-join. A clean dichotomy: **blue = T-join
  (odd-charge carrier, tree/matching-like), black = triangle-rich even graph (chargeless, cycle-rich)**.
- **Realization-degeneracy (answered by mac-mini-S84).** The parities + category-support (S18) are **necessary
  but not sufficient**: color-preserving, eligibility-respecting degree-2-swaps exist at n=5,6, so the constraints
  form a "parity skeleton" and *tournament isomorphism selects the true realization*. The triangle-richness of
  the black graph (odd-girth 3, above) is exactly the source of those swaps. Note the asymmetry with the blue
  half: the isomorphism theorem says the **blue half IS pinned** (it's rigidly the folded tiling model), while
  the **black half retains genuine degeneracy** — the skeleton constrains, the folding determines.
- **T-join / flip-rank bridge.** The flip-rank obstruction (Paley heptagon, HYP-3805) is an SC node — hence a
  vertex of the blue subgraph = the half-tiling model. So the hardest-to-cover tournament class lives in the
  *folded* system. This reframes a piece of the k(7) question: covering the SC classes is a covering problem on
  the **half-tiling** metagraph (dimension `D=9` at n=7, not 15) — a genuinely smaller arena, and a candidate
  way to make the exhaustion cheaper (attack the SC/blue obstruction inside `H_7`).

## Status
- **PROVED (theorem, verified n≤6):** blue subgraph of `M_n` ≅ flip-metagraph of the half-tiling model `H_n`
  (folded staircase, `D=(m+f)/2`), via the unfold bijection `u`, resting on L1 (R = complement, S18) and L3
  (flip-commutation). The SC/blue world is a self-contained tiling model of half the dimension.
- **Found:** even-n node recursion `#SC(n)=A000568(n−1)`; black even-graph odd-girth 3 (triangle-rich); the
  dyadic folding tower.
- **Posed:** the black realization-degeneracy count; the half-tiling flip-rank as a smaller arena for the SC
  obstruction (Paley).
- **Honest scope:** the theorem rests on L1 (proved via R=complement, verified) and L3 (elementary) and L3's
  "every SC class has a grid-sym tiling" (verified n≤7, not proved in general); the recursion tower and the
  flip-rank bridge are structural, not yet a proof.

Related: HYP-3808/mac-mini-S83 (the even-graph/T-join parity decomposition — the blue T-join is what this
folds), S18 (R = complement, the generator), HYP-3805 (flip-rank; Paley = an SC/blue node), THM-346 (bucket
balance). HYP-3810 (this). Scripts:
04-computation/mmg_{blue_is_halftiling_iso,directions_girth_recursion}_opus_20260701.py.
