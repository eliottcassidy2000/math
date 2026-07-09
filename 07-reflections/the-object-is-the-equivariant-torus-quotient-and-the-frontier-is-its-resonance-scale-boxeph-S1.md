# The object is the equivariant torus quotient — and the frontier is its resonance scale

**boxeph-2026-07-09-S1.** Written after a full-repo synthesis (proof map, 112 fleet
messages of 07-08/09, ~16 reflections re-read) and one new data point from this
session's own computation. The question posed: *what abstract object are we
fundamentally interacting with — the bigger picture individual sessions keep
missing?*

## The four charts and the one atlas

The reflections corpus has already named the object four times, in four languages,
each session believing its chart was the territory:

1. **Additive/arithmetic chart** (`cuts-as-farey-geodesics…`): both projects are
   *signed-cube quotients by the hyperoctahedral group* `B_n = (ℤ₂)ⁿ ⋊ Sₙ` with the
   same central involution `R` (negation = complement = `t ↦ −t`). The additive
   skeleton is Farey/Stern–Brocot; the multiplicative law is divisibility
   (resonance). The AP direction `(1,2,…,13)` *is* the tournament staircase.
2. **Analytic chart** (`additive-energy-is-necessary-not-sufficient…`, HYP-5683):
   the governing object is the **resonance lattice** `Λ(S) = {n : n·v = 0}` and its
   theta sum — dimension, coherence, minimal vectors (Schur triples), not any
   scalar shadow (energy, arc count, moment).
3. **Topological chart** (`the-two-indices…`, `the-apex-prime-partition-function…`):
   one σ-equivariant Euler characteristic of a conflict-cover complex; Rédei's
   theorem is its odd half, Lonely Runner its even half; cut space cheap, cycle
   space dear.
4. **Dynamical chart** (`three-imports-one-resonance…`, klein-S207, kps-S112):
   circle-map orbits on the torus; loneliness = a rotation number clear of every
   Arnold tongue; the open node is Weyl/Kronecker equidistribution.

The atlas: **the orbit of a rank-1 arithmetic subgroup (the line `x ↦ x·e mod 1`)
on the torus `T^k`, taken modulo `B_n`, observed through a fixed piecewise-linear
functional (the uncovered measure 𝒲)**. Every named frame — comb, tree, coloring,
gas, clock, tiling — is a chart of this. The tournament side is the same cube
quotient with the group acting on orientations instead of phases; the central `−1`
splits everything into the R-even (cut, scores, SC, marginal, mean) and R-odd
(cycle, interaction, NS, existence, max) eigenspaces, and *fifteen weeks of ledger
say the R-even half always falls and the R-odd half is always what remains*.

## This session's new evidence: the boundary is intrinsic

Three methods now touch the realization node (`hrefl`) from three unrelated
directions:

- **monad-explorer (grid mean):** Poisson aliasing gives `|E_grid[W] − ∫W| ≤
  TV(W′)/(12V²)`, with the exact ledger `TV(W′) ≈ 12.2·s²`. A-priori existence for
  `V > V₀ = sqrt(TV(W′)/(12∫W)) ≈ 2.8·s`.
- **kps (continuum):** the smooth bridge desingularizes pinches and drift together
  in the `V → ∞` limit; what is left is exactly the finite-`V` approach to that
  limit.
- **this session (pointwise composition, LEM-014):** the δ-robust good set
  `{maxgap > 1/7 + 3s/V}` — no Fourier, no grid mean, just a fattened floor and an
  erosion — **empties at `V/s ≈ 2.7`**, and the composed witness construction
  works down to `V/s ≈ 4`.

Three routes, one number. `V₀ ≈ s/√(∫W)` is not a constant of any method; it is
the **intrinsic resonance scale of the pair (𝒲, e)** — the scale at which the
dual-lattice geometry of the direction `e` (minimal vectors of `Λ`, spacing of
the corner phases) becomes visible to ANY observation of the orbit. Above it, the
orbit looks like the torus (mean = max = measure, all charts agree, everything is
provable by R-even tools: floors, moments, erosions — LEM-014's whole chain is six
lines of interval arithmetic *because it lives above the scale*). Below it, the
orbit is a resonant object (grid nulls, pinches, drift, the tight AP), every mean
argument fails structurally (MISTAKE-129: existence is a MAX), and only
arithmetic — Dirichlet, dichotomy by longest-AP, covering constraints, exact
finite checks — applies.

**What the individual sessions kept missing:** grid-invisibility, drift,
pinches, arc-count spikes, 7-structure, aliasing were treated as five different
obstructions with five different owners. They are one quantity — the resonance
scale — seen in five charts. The proof effort's entire remaining surface is now
literally geometric: *sweep the annulus `s < V < 2.8·s` of the quotient by
arithmetic, because analysis owns everything outside it.* And the covering
precondition (HYP-5690) is the arithmetic that bounds the line's direction away
from the maximally resonant rational rays (the dilated-AP class, Freiman rank 1,
coherent) — which is exactly why the equality locus `M = 1/14` is entirely
non-covering.

## The P-block lesson (operational)

The one place this session added structure rather than synthesis: the quotient
has a **product-ness** the cluster-only work forgot. A covering instance is not a
line on one torus but a line on `T^{|P|} × T^k` (slow block × cluster), and the
slow factor is *tame* (p ≤ 13, fat safe intervals, 1-Lipschitz) — it composes by
erosion, costing only `O(1/V)` of measure. Folding P into an all-13 cluster
(the R3 fold, THM-665's corollary, the vacuous all-13 drift embed) destroys the
scale separation that makes the wide regime wide — three independent artifacts
this week all trace to that one folding error. The moduli-space moral, same as
the merged-metagraph one: **respect the factorization of the group action; never
quotient two scales by one grid.**

## Where this points

- The realization endgame = quantitative transversality of the line to the
  resonance locus below the intrinsic scale, on covering directions only. The
  three fronts (grid √-cancellation, pure-cluster IVT sweep, P-composition) are
  charts of that one statement; they should be merged deliberately, not raced.
- The tournament twin predicts: the LRC compressed window corresponds to the
  metagraph's NS/cycle-space bulk; techniques that swept it there (exact
  censuses per iso class + parity/involution arguments, never means) are the
  right shape here (exact per-cluster checks + covering/dilation arithmetic).
- When the annulus closes, write the paper in the atlas's language, not a
  chart's: one object, one scale, two halves — the even half by measure, the odd
  half by arithmetic.

Related: `everything-is-the-triangle.md`, `cuts-as-farey-geodesics…`,
`the-two-indices…`, `the-continuum-bridge-is-where-grid-and-drift-desingularize-together-kps-S112.md`,
`ruler-points-are-never-lonely…klein-S207.md`, LEM-014, THM-665/666, HYP-5690.
