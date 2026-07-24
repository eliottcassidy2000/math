# The isoperimetric-dimension axis: the honeycomb's √ is a dimension, and G_n has none

*opus-2026-07-23-S1. Seed: OEIS A263135. See HYP-9022, script
`04-computation/metagraph_isoperimetric_dimension_opus_S1.py`.*

## The resonance

The owner handed me A263135 — max penny contacts on the **honeycomb**,
`a(2n) = 3n - ⌈√(3n)⌉`. It looks like a packing curiosity. It is not. That
`⌈√(3n)⌉` is not a correction term; it is a **dimension**.

Every "max-contact" lattice sequence has the same shape,
`max-edges(N) = (d/2)·N − c·√N`: bulk minus a √-thin boundary (triangular
A047932, square A123663, honeycomb A263135, verified together). The exponent on
`N` in the boundary is `(D−1)/D`. For these planar lattices `D = 2`, so the
boundary is `N^{1/2} = √N`. **A263135's square root literally says "the
honeycomb is two-dimensional."** Same for its diameter, `~√N`.

So the seed hands us a *ruler*: a graph whose isoperimetric dimension is 2, with
which to measure any other graph — in particular the object this whole project
orbits, the tournament metagraph `G_n`.

## What the ruler reveals

`G_n` reads at the **opposite end of the axis** (three independent legs, HYP-9022):

- diameter `~ 0.7·log₂V`, **logarithmic** — not `√V`;
- algebraic connectivity `λ₂(L) ≈ 2`, flat — forcing a **linear** edge-boundary
  at balanced cuts, not `√V`;
- the sparsest cut is a balanced bisection **transverse** to every coordinate we
  privilege — it ignores the H-gradient and ignores the SC/NS split.

`G_n` has **no** isoperimetric √. It is small-world; conjecturally `D = ∞`. The
honeycomb and the metagraph are isoperimetric **antipodes**.

## Two species of square root

The project's canon is saturated with √'s and told to see them as the triangle's
constants (`√2` hypotenuse/leg, and the honeycomb even brings `√3`). This seed
forces a distinction I had been blurring:

1. **Metric √'s** — `√2`, `√3`, the fiber fraction `~1/√(πn)`. These live in the
   *geometry* of the staircase triangle and its Gaussian/Wallis asymptotics. They
   are genuinely present.
2. **Isoperimetric √'s** — the `√N` *boundary exponent*. This encodes the
   *dimension* of a graph. And on the metagraph it is **absent**: `G_n` is not a
   surface, it is an expander.

The mathematics points past the pretty coincidence "√ everywhere" to a sharper
claim: the √'s that matter for *counting* (metric, present) and the √'s that
matter for *cutting* (isoperimetric, absent on `G_n`) are different animals. The
"everything is the triangle" frame captures the first and is silent on the second.

## The one-line consequence for proof strategy

The honeycomb bound is proved by **peeling a thin boundary layer** (perimeter
`~√area`). The antipode result predicts, with a mechanism, that *any* `G_n`
inequality attempted this way must fail: **there is no thin layer to peel.** A
balanced half of `G_n` touches the other half along a front that is a constant
fraction of the whole. This is why the principal line is a *gradient* and never a
*cut* — and it retro-explains MISTAKE-035 (G_n is not an H-DAG) and MISTAKE-037
(H-convexity fails): both are symptoms of a graph with no low-dimensional
boundaries to be convex or monotone across.

## Where it points next

The axis is a *ruler for the whole zoo*. The spine (SC-SC) already reads as the
sparse quasi-1D thread and the sea (NS-NS) as the `D=∞` core — the project's own
anatomy spans the dimension axis. The live question the seed opens: **is any
natural tournament graph genuinely 2-dimensional — a real honeycomb analog?** The
even-graph metagraph `E_n` (the standing opus mandate, and the honeycomb *is* a
tiling by even 6-cycles) is the first door to knock on. If some `E_n`-derived
graph reads `D=2`, the penny sequence stops being an analogy and becomes a
member of the family.
