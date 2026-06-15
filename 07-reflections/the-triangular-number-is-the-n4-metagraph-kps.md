# The triangular number's two operations ARE the n=4 tournament metagraph

**Source:** kind-pasteur-2026-06-14-S3. Dispatch (user observation): two ways of
writing `T(x) = x(x+1)/2` via `a(x)=x+1`, `b(x)=x/2`; the four compositional
options match the four tournament iso classes on 4 vertices.

## The observation, verified

Write the triangular number with two "alphabetical" operations
`a(x) = x+1` (additive successor) and `b(x) = x/2` (dyadic halving). Then
`T(x) = x · b(a(x)) = x(x+1)/2` for **all** `x` — and the parity-split form
`T = f·g` with `f = a` (x even) / `id` (x odd), `g = b` (x even) / `b∘a` (x odd)
is the same identity, routing the *forced* even factor: consecutive integers
`x, x+1` always contain exactly one even, so `v₂(x(x+1)) ≥ 1` and the `/2` lands
on whichever factor is even. (All four lines check by hand.)

The four ways to apply `{a, b}` each at most once are the subsets
`∅, {a}, {b}, {a,b}` — the **Boolean square** `2^{a,b}`. The user's claim: this is
the structure of the **four tournament isomorphism classes on 4 vertices**
(A000568(4) = 4). Computed (`04-computation/triangular_composition_metagraph_kps.py`),
it is a **graded-poset-with-involution isomorphism**, exact on every feature:

| Boolean square `2^{a,b}` | rank `|S|` | tournament class on 4 vtx | `H` | complement |
|---|---|---|---|---|
| `∅` | 0 | transitive `(0,1,2,3)` | 1 | self-comp (fixed) |
| `{a}` | 1 | dominator + 3-cycle `(1,1,1,3)` | 3 | ↔ `{b}` |
| `{b}` | 1 | sink + 3-cycle `(0,2,2,2)` | 3 | ↔ `{a}` |
| `{a,b}` | 2 | strong `(1,1,2,2)` | 5 | self-comp (fixed) |

- **Grading.** Subset size `|S| ∈ {0,1,1,2}` ↔ Hamiltonian-path count
  `H ∈ {1,3,3,5}` (the OCF; `H` is the metagraph's height function). The two
  rank-1 elements ↔ the two `H=3` classes.
- **Involution.** Tournament **complementation** `T → T^op` fixes the transitive
  and the strong class and swaps the dominator/sink pair. On subsets this is
  **swap `a ↔ b`** (fixes `∅` and `{a,b}`, swaps `{a}`/`{b}`) — *not* subset
  complement (which would swap the extremes). The match is with the leg-swap.
- **Self-complementary ↔ extremes.** The two self-complementary classes are
  exactly the rank-extremes `∅` (apply neither) and `{a,b}` (apply both); the
  complement-pair is the two singletons (apply one). This is precisely the user's
  "one empty, the other both combined, the remaining two a single of either."

## Why the involution is the leg-swap — and why that's the point

`a = +1` is the **successor / additive** operation; `b = /2` is the **dyadic /
halving** operation. In the repo's foundation ("everything is the triangle") the
staircase `δ_{n-2}` has two legs: the **vertical leg** (source column = scores,
the additive `+1`/Mode-A direction) and the **horizontal leg** (sink row =
complement, the SC distinction). Complementation `T → T^op` is exactly the
reflection that swaps the two legs. So:

> tournament complement = swap(`a`,`b`) = swap the staircase's two legs = swap
> the successor operation with the halving operation.

And the composition that *builds* the triangular number, `b∘a` (apply successor
then halve), is the arc-count itself: `T(n-1) = C(n,2)` = the number of arcs in an
`n`-vertex tournament = the staircase area. The user's `T(x)=x·b(a(x))` is the
arc-count written as "one of each leg," and the Boolean square of `{a,b}` is the
`n=4` shadow of the two-leg structure.

## Why n=4 specifically — the honest boundary

The clean Boolean-square match is **special to n=4**. A000568 = 1, 1, 2, 4, 12,
56, … : the values `1,1,2,4` (n=0..3 classes... n=4 gives 4) are the powers of 2,
and **4 = 2² is the last one** — n=5 has 12 classes, not `2^k`, and its H-graded
poset (H ∈ {1,3,5,9,11,13,15}, with degeneracies, 8 self-comp + 2 comp-pairs) is
not a Boolean lattice of any operation set. So n=4 is the unique place where the
metagraph *is* the Boolean square of the two fundamental operations — the smallest
nontrivial level, where "apply each leg at most once" closes up exactly.

What *does* persist to all n (the general shadow, not the coincidence): the
metagraph is graded by `H`; the rank-extremes (transitive `H=1`, regular/strong
`H=max`) are self-complementary; complementation is the leg-swap with the
comp-pairs off the SC-spine. The n=4 case is these general features collapsing
onto a 2-element Boolean lattice because there are exactly two ranks-worth of room.

## The 2-adic coda

`b = /2` is the halving the repo keeps meeting as the **2-adic seam** (THM-469's
parity dichotomy, the dyadic valuation). `T(x)=x(x+1)/2` always absorbs *exactly
one* factor of 2 (consecutive integers), and the four compositions are the four
parity-routings of "where the `/2` and the `+1` live." That the two
self-complementary tournaments sit at the *parity extremes* (neither/both) and the
complement-pair at the *single-operation* middle is the same `ℤ → ℤ/2` parity
structure, now read on the 4-vertex metagraph. The triangular number is where the
additive `+1` and the dyadic `/2` meet — and at `n=4` that meeting *is* the
tournament metagraph.

## Status / honesty

This is a **structural resonance** (a graded-poset-with-involution isomorphism at
n=4), verified exactly, not a new theorem about tournaments — the `H` values and
complement structure are known; the new content is the identification with the
`{a,b}` composition algebra and the reading of complementation as the leg-swap.
Extension to n≥5 is honestly negative for the *Boolean* form (HYP-2515); the
general H-grading + leg-swap involution is the part that lifts. Cross-links:
[[triangle_foundation]] (the two legs), THM-499 (`H = 1+2(c3+c5)+4D`, the grading),
the merged-metagraph reflections (SC-spine = self-comp at extremes), THM-469 (the
`/2` dyadic seam).
