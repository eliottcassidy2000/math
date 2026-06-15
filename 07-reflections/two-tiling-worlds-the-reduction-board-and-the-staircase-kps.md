# Two tiling worlds: the 1-D reduction board (growth constants) and the 2-D staircase (the OCF) — and the bridge is the user's board model

**Source:** kind-pasteur-2026-06-15-S2. Dispatch (user): think about the tiling
interpretations of these phenomena, leverage the concept, explore. The user's
model: `R_d(n,k)` counts tilings of a length-`(n+d)` board with exactly `k` long
tiles of length `d+1` and the rest unit tiles; `S_d(n)=Σ_k R_d(n,k)` satisfies
`S_d(n)=S_d(n-1)+S_d(n-d-1)` because the last tile is a unit (`S_d(n-1)`) or a long
tile (`S_d(n-d-1)`).

## The repo has been tiling in two dimensions at once

The board model is **1-D tiling**. The repo's foundation ("everything is the
triangle") is **2-D tiling**: a tournament = a binary tiling of the staircase
`δ_{n-2}` (the `m=C(n-1,2)` cells, each `±`), with the OCF `H=I(Ω,2)`. These are
two different tiling worlds, and the user's family is the bridge — it is the *1-D
reduction shadow* of the 2-D staircase.

## The centerpiece: Mode-A/Mode-B descent IS the d=1 board tiling

The repo reduces tournaments two ways: **Mode A** (`n→n-1`, hypotenuse removal /
vertex insertion) and **Mode B** (`n→n-2`, both legs / Cayley–Dickson). A reduction
schedule from `n` down to the base is a sequence of `−1` and `−2` steps — i.e. a
**composition of `n` into parts 1 and 2 = a monomer–dimer tiling of a length-`n`
board**, the user's `d=1` model. Verified: the number of Mode-A/B schedules from
`n` is `1,2,3,5,8,13,21,…` = **Fibonacci**, growth `φ`
(`04-computation/pascal_slope_family_growth_ladder_kps.py` companion run).

> **This gives HYP-614 a combinatorial origin.** HYP-614 says `φ` (the Dedekind
> regulator `log φ` of `ℚ(√5)`) governs `H` growth, the Ising free energy, and the
> Lyapunov exponent. In tiling terms `φ` is simply the **entropy of the
> reduction-board tiling**: the descent has two tile types (Mode A = unit, Mode B =
> dimer), and a 2-tile board grows at the golden ratio. The arithmetic constant and
> the tiling entropy are the same number because the reduction process is a
> Fibonacci tiling.

**The Mode-ladder generalization (HYP-2519).** If the repo defined a "Mode-`d`"
reduction (`n→n-(d+1)`) alongside Mode A, the descent would be the user's `d`-board,
with growth `β_d` = dominant root of `x^{d+1}=x^d+1` — the ladder `2 ⊃ φ ⊃
supergolden ⊃ … ⊃ plastic ⊃ … → 1`. The repo currently uses steps `{1,2}` (d=1,
golden). Higher modes would access the higher constants. Mode A is the hypotenuse
(`H=1+2^d` direction) and Mode B both legs; a Mode-`d` would peel a longer
sub-staircase — the question is whether such a reduction is *natural* on tournaments
(does a `n→n-(d+1)` tournament operation with the board recurrence exist?).

## Fugacity: the OCF is the board at weight 2

Weight each long tile by `x`: the board count becomes the **path-power independence
polynomial** `I(P_n^d, x) = I(P_{n-1}^d, x) + x·I(P_{n-1-d}^d, x)`. So the user's
`S_d(n)` is fugacity `x=1`. At `d=1`: `x=1 →` Fibonacci, `x=2 →` **Jacobsthal**
`1,1,3,5,11,21,43,…` (growth 2) — verified. The repo's OCF lives at **fugacity 2**
(`H=I(Ω,2)`): the `2` in the OCF is "each long tile / each odd cycle weighted 2,"
the same `x=2` that turns the Fibonacci board into Jacobsthal. THM-485's "two
temperatures" (claudebox) is exactly this `x`-axis at `d=1`. The user's family
(`x=1`, all `d`) and THM-485 (`d=1`, all `x`) are the two axes of a `(d,x)` grid of
weighted tilings, with `H` at `(d=?, x=2)` on the conflict graph `Ω` rather than a
path.

## Free vs interacting: the board is the mean-field reference for the OCF

The board count `R_d(n,k)` is **always realizable and log-concave** (verified `d=1`,
`n≤8`) — it's a genuine tiling count, the *free* hard-core lattice gas (the column
sequences are the smooth Pascal diagonals). The tournament's OCF-coefficient vector
`α_k` (disjoint-odd-cycle collections, the THM-466 2-adic digits) is the **same kind
of object on the conflict graph `Ω`** — but `Ω` is *interacting*: its Helly /
intersection structure forbids certain `α`-vectors, producing the **forbidden values
`H∈{7,21}`** (THM-029/075; last session's `0+0=1` work). So:

> The board tiling is the **free / mean-field reference**; the tournament OCF is the
> **interacting version on `Ω`**; the forbidden values measure the *deviation* of the
> real conflict graph from the free gas. "Why is 7 forbidden?" becomes "why can't the
> free board's `α=(3,0)` configuration be realized by an actual `Ω`?" — the
> three-pairwise-intersecting-cycles Helly obstruction, an interaction the free board
> doesn't have.

This reframes the forbidden-value program (last session's verdict) cleanly: the
unconstrained model (board, gap-free) vs the constrained model (`Ω`, gaps); the gaps
are pure interaction.

## The growth constants are tiling entropies (the arithmetic of the rungs)

`β_d` is the **topological entropy** of the `{1, d+1}`-tile subshift of finite type;
`log β_d` is the free energy per site. The named rungs are exactly the famous tiling
constants: golden `φ` (`d=1`), the **plastic number** (`d=4`, the smallest Pisot, the
Mahler measure of `x³−x−1`). So HYP-2518's "arithmetic of the higher rungs" has a
clean dynamical answer: they are **subshift entropies / Mahler measures of the tile
recurrence** — and `log β_d` is the per-step growth of the corresponding Mode-`d`
reduction. HYP-614's `R = log φ` is the `d=1` instance of "regulator = tiling
entropy."

## The synthesis

The user's board model is the **1-D reduction world** the repo has been using
implicitly (Mode A/B = monomer/dimer), with `φ` (HYP-614) as its fundamental entropy
and a whole ladder of higher entropies (supergolden, plastic, …) waiting behind
higher reduction modes. The OCF `H` is the same tiling at **fugacity 2 on the
interacting graph `Ω`**; the forbidden values are where the interaction breaks the
free board's smoothness. Two tiling worlds — 1-D reduction (constants) and 2-D
staircase (the OCF) — meet at the fugacity-2 weighting and the reduction schedule.

Cross-links: HYP-614 (`φ` = the d=1 reduction-tiling entropy), HYP-2518 (the rung
ladder; the constants are subshift entropies), HYP-2519 (the Mode-ladder),
THM-485 (the fugacity axis), THM-466 (the OCF `α_k` = the interacting board columns),
THM-029/075 (the forbidden values = interaction deviation),
[[the-slope-ladder-of-pascal-and-the-growth-constants-kps]] (the constant ladder),
[[triangle_foundation]] (the 2-D staircase world).
