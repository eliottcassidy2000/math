---
id: HYP-9022
title: "The tournament metagraph G_n is the isoperimetric ANTIPODE of A263135's honeycomb (dim 2 vs dim infinity)"
status: >
  VERIFIED for n=3..6 (exact) and corroborated through n=8 (via the repo's
  corrected diameter). The tournament arc-flip metagraph G_n is small-world /
  expander-like: diam(G_n) ~ 0.7 log2 V, algebraic connectivity lam2(L) ~ 2
  stays bounded away from 0, and the sparsest cut is a balanced bisection
  TRANSVERSE to both the H-gradient and the SC/NS split. This is the opposite
  end of the isoperimetric-dimension axis from the A263135 seed (honeycomb,
  dim 2, boundary ~ sqrt N). OPEN: the true n->infinity expansion constant
  (does conductance stay bounded below, or decay slowly?) and the isoperimetric
  dimension of the SC-spine.
source: opus-2026-07-23-S1
depends_on: []
related:
  - MISTAKE-036   # corrected diam(G_n)=1,2,3,4,7,8 ~ n^2/4 -- reused as the log-V leg
  - HYP-3543      # R-even bulk (A000568+SC)/2 vs R-odd -- the repo's bulk/boundary split
  - HYP-2748      # spine=SC, sea=NS; wiggly metagraph = hypercube Q_m
  - OPEN-Q-108    # A000568 edge-sandwich (the only prior max-edges framing)
seed: OEIS A263135 (max penny contacts on the honeycomb; a(2n)=3n-ceil(sqrt(3n)))
script: 04-computation/metagraph_isoperimetric_dimension_opus_S1.py
output: 05-knowledge/results/metagraph_isoperimetric_dimension_opus_S1.out
prior_art: 04-computation/concrete_cheeger_s92v.py  (exact Cheeger at n=5 only; "Paley=max expander")
---

# HYP-9022 — G_n is the isoperimetric antipode of the honeycomb

## The seed as a dimension statement

OEIS **A263135** = max penny-to-penny contacts of `N` pennies on the **honeycomb**
lattice; owner's form `a(2n)=3n-ceil(sqrt(3n))`. It sits in a family of three
**edge-isoperimetric profiles** of `d`-regular vertex-transitive graphs (all
verified in the script):

| lattice | OEIS | `d` | max-edges(N) |
|---|---|---|---|
| triangular | A047932 | 6 | `floor(3N - sqrt(12N-3))` |
| square | A123663 | 4 | `2N - ceil(2 sqrt N)` |
| honeycomb | A263135 | 3 | `(3/2)N - ceil(sqrt(1.5N))` (even N) |

All are **`(d/2)N - c*sqrt(N)` = BULK minus a sqrt-THIN boundary.** The exponent
`1/2` is `(D-1)/D` with `D=2`: the `ceil(sqrt(...))` in A263135 *is* the
statement "the honeycomb is 2-dimensional." For any vertex-transitive graph of
isoperimetric dimension `D`, the min edge-boundary of an `N`-set scales like
`N^{(D-1)/D}`; and diameter scales like `V^{1/D}`.

## Porting the question to G_n (arc-flip metagraph: iso-classes = vertices, single-arc-flip = edges)

Self-contained rebuild reproduces canon exactly: `V=A000568(n)=2,4,12,56`,
`E=1,5,30,290`, `SC=2,2,8,12`, `Hmax=3,5,15,45`, `diam=1,2,3,4` (n=3..6).

**Three independent legs, all saying "not 2D, but small-world / expander":**

1. **Diameter is LOGARITHMIC, not sqrt.** `diam(G_n)=1,2,3,4,7,8` (n=3..8; the
   corrected sequence, MISTAKE-036) vs `log2 V = 1,2,3.6,5.8,8.8,12.8`. Ratio
   `diam/log2 V = 1.0,1.0,0.84,0.69,0.79,0.63` — bounded. So `diam ~ 0.7 log2 V`.
   A 2D lattice would have `diam ~ sqrt V` (`=1.4,2,3.5,7.5,21,83`) — off by an
   order of magnitude already at n=8. **Logarithmic diameter is the small-world
   signature.**

2. **Algebraic connectivity stays bounded away from 0.** `lam2(L) = 2.0, 2.0,
   1.60, 1.96` (n=3..6). The Fiedler bound `|dS| >= lam2 * |S|(V-|S|)/V` then
   forces min-edge-boundary `>= ~ V/2` at balanced cuts — **linear, not sqrt V.**

3. **The sparsest cut is TRANSVERSE to the project's coordinates.** Exact
   min-conductance cut (n=5, `phi=0.379`): a balanced 6|6 bisection whose two
   sides each span the *full* H-range `[1,15]`/`[3,15]` and each mix SC and NS.
   It is **not** an H-interval (down-set) and **not** the SC-vs-NS split. So the
   **principal line / H-gradient is a gradient, not a bottleneck**; there is no
   thin neck along the repo's natural axes.

## Verdict and the dimension axis

- **honeycomb A263135:** isoperimetric dim **D=2**, boundary `~ sqrt N`, diam `~ sqrt N` — amenable.
- **tournament G_n:** boundary `~ Theta(V)` at balanced cuts, diam `~ log V`,
  transverse sparsest cut — **small-world, conjecturally D = infinity (expander family).**

They are **isoperimetric antipodes.** Structurally expected: the wiggly metagraph
is the S_n-quotient of the hypercube `Q_{C(n,2)}` (HYP-2748), and hypercubes have
linear (not sqrt-thin) edge-boundary at balanced sets — the honeycomb's opposite.

Internal anatomy **spans the dimension axis**: the **SEA (NS-NS)** is itself the
expander core (n=6: 44 nodes, `lam2(L)=1.98`, matching the whole graph), while the
**SPINE (SC-SC)** is the sparse low-dimensional skeleton (n=6: 12 nodes / 13 edges,
maxdeg 3, cyclomatic 2 — near-tree, quasi-1D; though at n=5 it is more branched, so
its true dimension is itself open).

## What is OPEN (do not overclaim)

- **Asymptotic expansion.** The *normalized* gap decays (`2,1,0.53,0.32`) purely
  because avg-degree grows (`~C(n,2)`) while `lam2(L)~2` is flat. Whether the
  scale-free conductance `phi` stays bounded below as `n->inf` (true expander) or
  decays slowly is **not settled** by n<=6. Exact `phi = 0.6 (n=4), 0.379 (n=5)`;
  n=6 bracketed `phi in [0.16, 0.46]`.
- **Isoperimetric dimension of the SC-spine** (candidate finite-D substructure).
- **A 2D-like tournament graph?** Is there *any* natural tournament construction
  whose metagraph has finite isoperimetric dimension (a genuine honeycomb analog)?
  The even-graph metagraph E_n (opus mandate) is the first place to look.

## Why this matters (the structural bridge, per the analogy-bridge pattern)

The seed and the project's central object answer the *same* extremal question
("how thin can the boundary of a k-blob be?") at **opposite extremes**. The
honeycomb's `sqrt` is a dimension; G_n has no such dimension. This reframes the
"principal line" as a gradient direction inside a small-world graph, and predicts
that any G_n bound proved by "peeling a thin boundary layer" (the honeycomb proof
strategy) must FAIL — there is no thin layer to peel.
