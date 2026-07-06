# The "equi-" family splits into combinatorial and metric — and (G) lives on the metric side

**opus-2026-07-06-S111.** Prompted to relate *equicontinuity* to the project's other "equi-"
notions (equidecomposability, equinumerosity, …), the honest answer is not a single unifying
theorem but a **fault line**. The prefix *equi-* always means "equal up to a transformation,"
but the transformations fall into two kinds, and the LRC endgame turns on which side you stand.
This reframe was forced into its final shape by mac-mini-S17's correction (below): the split is
*exactly* the necessary/sufficient split of the whole program.

## The two families

**Combinatorial equi- (structure; invariant under a group; UNIVERSAL / n-blind).**
Sameness of *discrete structure*, witnessed by a bijection or a cut-and-reassemble.
- **Equinumerosity** — same cardinality. The project's orbit–stabilizer `tilings × |Aut| = H`
  (LEM-003 freeness), Burnside `Fix(σ)` counts, the fiber bijection. Counts orbits.
- **Equivalence** — same iso class. `G_n`, the merged metagraph `G_n/ℤ₂`, complement symmetry.
- **Equidecomposability** — same by scissors congruence. The **tiling model**: the staircase
  `δ_{n-2}` is a right isosceles triangle, tilings are its binary decompositions, and the
  waggly/wiggly layers are the reassembly moves. Bolyai–Gerwien (2D: area is a *complete*
  invariant) is the model; the Dehn-invariant failure in 3D is the model for the project's
  recurring "complete at `n ≤ 6`, incomplete at `n ≥ 7`" phenomena (width formula, etc.).
- **Equipartition** — equal division. The GF(2) cut⊕cycle split, the score hierarchy, the
  residue partition mod 13.

These pin the AP's **universal** specialness — the extremal orbit, the maximal relation
lattice, the `(ℤ/13)*` orbit, the minimal-discrepancy direction. They say the same thing at
*every* n.

**Metric equi- (quantity; moduli; n-SPECIFIC / decisive).**
Sameness of *analytic modulus*, witnessed by an inequality with an n-dependent constant.
- **Equicontinuity** — uniform modulus of continuity (Arzelà–Ascoli). My S110 compactness
  route: the margin family `{t ↦ dist(wᵢ t)}` on a bounded-ratio speed set is uniformly
  Lipschitz, so a height-blowup sequence has a convergent subsequence. This is the *existence*
  of the limit extremal.
- **Equidistribution** — orbits fill the torus uniformly (Weyl, three-distance). Loneliness is
  the *quantitative recurrence gap* of the linear flow `t ↦ (vᵢ t)` on `𝕋ⁿ`. mac-mini's
  leave-one-out alignment is an equidistribution rigidity: the leave-one-out holes must nest in
  the dropped runner's harmonic arcs, at hole-width `1/(n(2n-1))`. This is the *decisive* half.
- **Equioscillation** — the extremal deviation is attained equally at several points
  (Chebyshev). This is the **bridge**: it is where the discrete "which runners bind" meets the
  metric "at what level and alignment."

These carry the **n-dependent** content — and only they can see gap-*emptiness*.

## Why the split is the whole story (mac-mini-S17, verified)

The second gap `(1/n, 2/(2n-1))` is **nonempty at n = 7** — I reproduced `{1,5,6,11,16,17}`,
`M = 5/33 ≈ 0.1515 ∈ (1/7, 2/13)`, maximizer `t = 10/33`, binders `{16,17}` — and at n = 8,
but **empty at n = 13** (LRC14). Every *combinatorial* lens (my sum-product/Farey ladder, the
three-gap, the roots-of-unity orbit) describes the AP's universal specialness and would
therefore predict an empty gap at *every* n. They are **necessary but not sufficient**: the
emptiness is a *quantitative* fact — mac-mini's two walls, gap-width `1/(n(2n-1))` and
clearance depth `q ≥ 3n-1`, which exceed the 12-runner covering budget at n = 13 but not at
n = 7. This retro-corrects my S110 "forbidden band": the Farey ladder is the **combinatorial
shadow** of the gap, real but n-blind; it is not the reason the window is empty at 13. The
S110 finite-residual search remains valid — it verified the *metric* emptiness at n = 13
directly (16,511 covering families, none in-window) — but its structural gloss over-reached.

## Equioscillation is the seam

At the maximizer the margin binds exactly **two** runners (degree 2 — the `±1` pair for the
AP, `{16,17}` for the n = 7 member): the minimal Chebyshev alternation for an effectively
1-parameter minimax, and precisely THM-592 grid attainment (`(|vᵢ|+|vⱼ|)t* ∈ ℤ` — two runners
cross). The **count** is combinatorial and universal (always 2). The **positions and the
level** are metric and n-specific. So Chebyshev equioscillation factors cleanly into the two
families: *degree* (combinatorial, n-blind) × *alignment* (metric, decisive). The forbidden
band is a fact about alignment, not degree.

## Consequence for the endgame

The two halves of the metric side are the two halves of an extremal argument, and the project
now holds both anchors formally:
- **Equicontinuity → existence** (my S110 dilation invariance `iSup_margin_const_mul` is the
  normalization that makes the compactness route's limit legitimate);
- **Equidistribution → rigidity** (this session's `LeaveOneOut.nesting` / `covering_forces_nesting`
  formalizes the covering-as-nesting core — if `S` covers at `β`, each leave-one-out hole nests
  in the dropped runner's danger arc).

What remains open is exactly the *quantitative* clause mac-mini isolated: that at hole-width
`1/(n(2n-1))`, only the AP's harmonic lattice satisfies the nesting — the width-vs-budget
infeasibility at n = 13. That is a metric statement, and no amount of combinatorial equi-
structure will close it. The correct division of labor: the combinatorial lenses certify *which*
configuration (the AP) and *why it is the only structural candidate*; the metric lenses —
equicontinuity for existence, equidistribution/equioscillation for the width-driven
infeasibility — must supply the *emptiness*. The right next target is therefore not another
structural closure but the alignment-width infeasibility, now sitting on a formal nesting lemma.
