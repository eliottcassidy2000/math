# The RG flow's local + global structure: the cusp is a RAZOR-SHARP repeller (only 1 of 182 single-element AP perturbations keeps M=1/14 — the razor manifold is ~0-dimensional, steepest perturbation already at 2/27>1/14) AND a rare small basin (CUSP ~2–8%, measure-0); the DOUBLET basin (floor 0.198) DOMINATES (~50%), MULTIPLICITY ~10% (residue collisions), OTHER ~30–40%. The cusp is triply non-generic (measure-0 · sharp-repeller · small-basin), so the conjecture binds OFF-cusp; the basin is a STATIC 2-adic shell decomposition

*opus-2026-06-30. Owner: work both — linearize at the cusp, map the basins. Did both. The cusp's
non-genericity is now quantified three ways; the doublet is the generic attractor.*

## Task 1 — linearize at the cusp (the razor is ~0-dimensional)
The AP `{1..13}` is the global M-minimum (`M=1/14`). Over all **182 single-element replacements** `e→e'`
(`e∈{1..13}`, `e'∈{14..27}`):
> **Only 1 perturbation is FLAT** (keeps `M=1/14`); the other **181 are RELEVANT** (`M` increases). The
> flat one removes `e=12`; the steepest relevant perturbation still lands at `min M = 2/27 = 0.0741 > 1/14`.
So the AP extremal is a **razor-sharp, nearly-isolated minimum** — the "razor manifold" (perturbations
preserving `M=1/14`) is essentially 0-dimensional. In RG terms, **almost every direction at the cusp is
RELEVANT** (drives `M` up, away from the extremal); the cusp has no stable manifold to speak of. This is the
"razor's edge" made quantitative: the extremal is not a basin, it is a point. (And every perturbation keeps
`M ≥ 1/14` — the flow never dips below the conjecture's floor, consistent with LRC.)

## Task 2 — the basins (a static 2-adic shell decomposition)
The basin of a 13-set is fixed by its **2-adic shell profile** (sizes of valuation classes) + the residue
structure mod 7. Classifying `2000` random 13-sets per range:
| range | DOUBLET (0.198) | OTHER (>0.198) | MULTIPLICITY (<0.198) | CUSP (meas-0) |
|---|---|---|---|---|
| dense `{1..20}` | **52.9%** | 28.1% | 10.6% | 8.5% |
| mid `{1..60}` | **53.0%** | 33.1% | 11.2% | 2.8% |
| sparse `{1..300}` | **45.9%** | 42.0% | 10.0% | 2.0% |
- **DOUBLET basin (floor `0.198`) DOMINATES (~50%)** — the doublet is the *generic* attractor/floor; the
  typical 13-set realizes THM-590's `4cos²(3π/7)` exactly. Determinant: a 2-element shell with distinct
  residues mod 7 (common).
- **OTHER (~30–40%)**: no doublet shell, higher floor (`0.308`/`1`) — *easier* sets.
- **MULTIPLICITY (~10%, stable)**: a shell with a repeated residue mod 7 (a birthday-collision rate); `g`
  dips below `0.198` — the looser-bound regime, density-independent.
- **CUSP (~2–8%, density-dependent, RARE)**: needs 7 odds covering `Z₇`; shrinks as sets get sparser.

## The synthesis: the cusp is TRIPLY non-generic
The cusp (the `Z₇` / AP extremal) is non-generic in three independent senses, now all quantified:
1. **measure-0** (apex gap `g=0`, `ρ=0`) — the descent product vanishes there;
2. **sharp repeller** (Task 1: ~0-dimensional flat manifold, 1/182) — the flow leaves instantly and no
   perturbation stays;
3. **small basin** (Task 2: ~2–8%, shrinking with sparsity) — almost no 13-set is a cusp.
> **The conjecture binds OFF-cusp.** The cusp is the global extremal (AP, `M=1/14`) but a measure-0, sharp,
> rare special point — handled by the comb-witness. The *generic* set sits in the DOUBLET basin (floor
> `0.198`, positive measure), and the binding M-constraint is the **off-cusp covering-min `n/Φ₆(n)`** (the
> descent refinement). The RG geometry confirms it: the hard case is not the isolated cusp but the dense
> off-cusp interior where `M` is tightest.

## What the local+global picture buys
- **The razor is a point, not a ridge.** The AP's `M=1/14` is an isolated minimum (1 flat direction). This
  sharpens "razor-thin": the extremal has essentially no neighbors at the floor — perturbing it *anywhere*
  lifts `M`. The LRC bound is tight at exactly one configuration (up to the lone marginal direction).
- **The doublet is the generic attractor.** ~50% of sets flow to floor `0.198` — THM-590's value is not
  exotic, it is what a typical 13-set does. The obstruction atom is the modal floor.
- **The basin is static and computable** — the 2-adic shell profile + residues determine everything (cusp /
  doublet / multiplicity / other) without running the flow. The "dynamics" is a decomposition.

## Status
- **Computed (opus):** Task 1 — AP perturbation landscape (`1/182` flat, steepest relevant `2/27`); the
  razor is ~0-dimensional, cusp = sharp repeller. Task 2 — basin distribution (DOUBLET ~50%, OTHER ~30–40%,
  MULTIPLICITY ~10%, CUSP ~2–8%); basin = static 2-adic shell profile + residues.
- **Synthesis:** the cusp is triply non-generic (measure-0 · sharp-repeller · small-basin) → conjecture
  binds off-cusp (the dense interior, covering-min `n/Φ₆(n)`); the doublet (`0.198`) is the generic floor.
- **Honest:** the multiplicity basin (~10%) has `g<0.198` (looser bound); not the binding case but the one
  regime where the simple doublet floor doesn't directly apply.

Related: the-renormalization-flow-cusp-to-doublet (the flow this localizes/globalizes),
the-descent-product-is-renormalization (the off-cusp binding refinement), per-level-vs-total-doublet (the
0.198 floor), the razor-thin reflections (now: the razor is a POINT); klein THM-590, THM-523/580;
OPEN-Q-108.
