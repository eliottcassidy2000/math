# The converse-parity of the arc cube IS the FKN dictator hierarchy — scores are the level-1 ranking signal, cycles are the even content, and the tiling model breaks the symmetry that separates them

**Source:** kind-pasteur-2026-06-15-S6. Dispatch: consider the Friedgut–Kalai–Naor
(FKN) theorem — degree-1 Fourier concentration ⟹ near-dictator — and apply it to the
tiling model (transitive default + single-arc / sub-tiling perturbations); map tilings to
iso classes; find meaningful metrics; develop the signed sub-tournament inclusion–exclusion
`T = A+B+C−D−E−F+G`. Builds on THM-510 (the n=4 Boolean square), THM-502/507 (the
converse-shift `A↦A−J`), THM-468 (the determinant lens), and last session's
tournament-as-zero-sum-game / bipartisan-Nash work.

## The one-line discovery

> **On the full arc cube, the converse involution `T↦T^op` (reverse every arc) is the
> global sign flip `x↦−x`. So a tournament invariant's Fourier *level-parity* equals its
> *converse-eigenvalue*: converse-EVEN invariants (`H`, the cycle counts `c_k`, the whole
> OCF/conflict graph) live on EVEN levels; converse-ODD invariants (the scores) live on
> ODD levels. FKN's "level-1 vs higher" hierarchy is exactly the project's
> "ranking/transitive vs cyclic" split — and the tiling model's fixed base path is the
> chosen symmetry-breaking that turns the score into a genuine level-1 signal.**

This is THM-511, verified exactly (`fkn_converse_parity_kps.py`, `fkn_tiling_fourier_kps.py`):
on the full cube the odd-level mass of `c_3` and `H` is **zero** (`n=4,5,6`), the scores sit
**purely on level 1**, and `H` is supported on `L0,L2,L4`.

## Why FKN lands here, and what it says

FKN: a Boolean function with `≥1−ε` of its Fourier mass on level 1 is `O(ε)`-close to a
**dictator** `±x_i`. Its historical home is **quantitative Arrow** (Kalai 2002): a
pairwise-majority social-welfare function that *rarely* produces a Condorcet cycle must be
close to a dictatorship. A tournament *is* the pairwise-majority outcome of a profile; a
**3-cycle is a Condorcet cycle**. So the project's central dichotomy —

| transitive / score-ordered | ⟷ | Condorcet-cyclic |
| level-1, converse-odd, "the ranking" | ⟷ | level-`≥2`, converse-even, "the cycles" |
| `H=1` (transitive) | ⟷ | `H=max` (regular/Paley) |
| pure Nash on the dominator (last session) | ⟷ | uniform bipartisan Nash |

— is **literally FKN/Arrow's dictator-vs-spread axis**, read on the arc cube. Near-transitive
= level-1-concentrated = near-dictatorial; genuinely cyclic = even-level mass. The
**FKN-defect** I define, `(Fourier mass at level ≥2)/(total variance)`, is a clean new
**metric on tournaments** = "how much of the structure is irreducibly cyclic rather than
score-explained." For `H` on the tiling cube it grows `0.25 → 0.51 → 0.68` (`n=4,5,6`):
**score-determination of `H` fails progressively and quantitatively** — the precise sense in
which large tournaments stop being rankings.

## The default state and the perturbations (the dispatch's framing, resolved)

The transitive tournament is the **zero tiling** (`H=1`), the converse-symmetric reference
where *all* cyclic (even-level) mass vanishes and only the score line survives. The two
perturbation scales the user names map onto the two cube structures:

- **Flip a single arc** (a wiggly line, Hamming `d=1`): moves one coordinate. On the *tiling*
  cube this is a degree-1 (level-1-touching) move; it shifts two scores by `±1` and may create
  one 3-cycle. This is the smallest "ranking edit."
- **Flip a whole sub-tiling** (e.g. the complement tiling, `d=m`, or a triangle of tiles):
  a coordinated multi-arc move. The full-reversal of a sub-triangle is exactly the move under
  which a 3-cycle indicator is **invariant** — which is *why* cyclicity has no odd-level
  component. Sub-tiling flips are the even-level (cyclic) generators.

So "single arc" lives at the score/level-1 (ranking) scale, and "sub-tiling worth of arcs"
lives at the cycle/even-level scale — the two perturbation sizes the user intuited are the
two Fourier parities.

## The signed sub-tournament sum is the inclusion–exclusion that vanishes on cycles

`T = A+B+C − D−E−F + G` (three size-`(n−1)`, three size-`(n−2)`, one size-`(n−3)`
sub-tournaments) is the **degree-3 finite difference / Möbius derivative** over a distinguished
3-vertex set — exactly the inclusion–exclusion `Δ_{ijk}` whose value is the **level-3 Fourier
coefficient**. The discovery makes its fate precise:

> For any **converse-even** invariant (every OCF/cycle quantity, including `H`), the degree-3
> inclusion–exclusion `A+B+C−D−E−F+G` **vanishes** — there is no irreducible 3-arc
> interaction. The irreducible interaction is **degree-2**: the *pairwise* conflict between
> odd cycles, i.e. an **edge of the conflict graph `Ω`**. And `H = I(Ω,2)` (the OCF, THM-002)
> is the independent-set generating function of exactly this level-2 interaction graph.

So the user's "naive `2^{C(n,2)}` count ignores the structure" gets a sharp reading: the
structure it ignores is the **Möbius/Fourier stratification by converse parity**, and the
correct organizing object is not the count of tilings but the **even-level conflict hierarchy**
`L2 (Ω edges) ⊂ L4 (α_2 disjoint pairs) ⊂ …` — conjecturally the OCF `α_k` strata themselves
(HYP-2532; qualitatively confirmed: `α_k` first appears at level `2k`, since `k` disjoint
3-cycles span `3k` arcs and each contributes up to level 2).

## A clean by-product: the variance of the 3-cycle count

`Var(c_3) = 3·C(n,3)/16` over uniform labelled tournaments, sitting **entirely at level 2**
(verified `n=4,5,6`). Reason: a single triangle's 3-cycle indicator is
`¼ − ¼Σ_{pairs} x_a x_b` (level 0 + level 2, *no* level-1 or level-3 term, by reversal
symmetry), and distinct triples' cyclicities are **uncorrelated** — a shared single arc
carries no covariance because each indicator is even in its own three arcs. The triangle
count behaves like `C(n,3)` independent `Bernoulli(¼)` trials at the level of variance.

## How it threads the repo

- **THM-510 (n=4 Boolean square).** There, tournament complement = the leg-swap `a↔b`, fixing
  `∅` and `{a,b}`. Complement on iso classes is the converse; the **self-complementary classes
  are the converse-fixed points** = the converse-even extremes of the `H`-spine. The n=4
  square is the converse-parity grading collapsed onto two ranks.
- **THM-507 (walk counts are spectral).** Its engine is the **converse-shift**
  `A ↦ A−J = −(Aᵀ+I)` — the same `T↦T^op` seam. The spectral (det-side, converse-even) world
  vs the score/permanent side is the same even/odd split seen here.
- **THM-468/472/507 (determinant lens).** `S = A−Aᵀ` is converse-**odd** (`S(T^op)=−S(T)`),
  but the **skew spectrum** `{μ_j}` and `det(I+S)=∏(1+μ_j²)` are converse-**even** (eigenvalues
  in `±` pairs) — so the determinant lens is an even-level object, consistent with `H`'s even
  support. The OMWU **frequencies = skew spectrum** (last session) are therefore the
  even/cyclic content's oscillation; the **score** is the odd/ranking content; the bipartisan
  Nash interpolates (pure on the score-dominated transitive end, uniform on the cyclic Paley
  end) — last session's dichotomy is this session's parity split.
- **The `28 = C(8,2)` / E₇ antipodal pairs (HYP-2530).** The converse pairing `λ↔−λ` of the
  E₇ 56-rep symplectic form is exactly the global sign flip — E₇'s "tournament" is by
  construction a converse-odd (symplectic) object, so its invariants split by the same parity.

## Status / honesty

PROVED & VERIFIED: the level-parity = converse-eigenvalue theorem (`χ_S(−x)=(−1)^{|S|}χ_S(x)`),
the scores-are-level-1 and cycles-are-even facts (`n=4,5,6` full cube, `n≤7` tiling), the
`Var(c_3)` value, and the vanishing of the degree-3 inclusion–exclusion on converse-even
invariants. CONJECTURED: HYP-2532 (`H`'s even-level weights are the OCF `α_k` strata exactly —
qualitatively confirmed, exact weights intricate: `c_3` alone gives `7.5` of `H`'s `L2=16.875`
at `n=5`, the rest from `c_5` and cross-correlation) and HYP-2533 (the FKN-defect of `H` is
monotone in `n`, with a quantitative-Arrow lower bound on level-`≥2` mass — a real
tournament-flavoured Arrow inequality). Cross-links: THM-511, THM-510, THM-507, THM-468,
THM-002, HYP-2530, [[the-triangular-number-is-the-n4-metagraph-kps]],
[[two-lyapunovs-the-tournament-game-and-the-e7-symplectic-tournament-kps]].
