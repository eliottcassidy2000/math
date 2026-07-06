# The spectral-gap crux is a single covering-impossibility

*kind-pasteur-2026-07-06-S20 — a reflection written after LRCTorusSplit.lean and
LRCCircleCover.lean, integrating mac-mini HYP-4292/4302 and opus HYP-4246/4266.*

## The one sentence

Every lane of the LRC(14) spectral-gap crux — the density rung, the pole at
`25/4`, mac-mini's infimum-`1/6` census, opus's `φ > 0` — is a facet of a
**single statement**:

> Teeth-combs of radius `ρ = 2/25`, of DISTINCT frequency, cannot cover the
> circle when there are few enough of them.

The gap `(1/13, 2/25)` is empty because a 12-speed family that dodged its own
citation would have to cover a circle with too few distinct combs, and distinct
combs leave holes.

## Why this unifies what looked like separate arguments

A "clear point" of a family at level `ρ` is a `θ` (or `(t, θ)`) where every
runner's tooth `‖v·τ‖ ≥ ρ`.  The family's value `M ≥ ρ` iff a clear point
exists iff the combs `{‖v·τ‖ < ρ}` **fail to cover**.  So the whole crux is a
covering question, and it has exactly two regimes, meeting at one pole:

**Density regime (few combs).**  Each comb has measure `2ρ`.  `l` combs cover
measure `≤ 2ρl`.  If `2ρl < 1` — i.e. `l ≤ 6` at `ρ = 2/25` — they cannot cover,
by measure alone.  This is `torus_split_rung` / `circle_clear_of_density`
(LRCTorusSplit, LRCCircleCover), and it is the same computation as:

  * the **cluster-gcd pole** `|S| = 25/4` (LRCClusterGcdSharp): `(25 − 4|S|)`
    goes non-positive at `|S| = 7`;
  * the **fee-mean ceiling** (the l-subset density ladder, `2ρl ≥ 1` at
    `l = 25/4`, S7–S9);
  * mac-mini's **density wall** `2·(3/38)·c ≥ 1` at the 3/38 cell.

Four names, one inequality: `2ρl ≥ 1`.  The pole `l = 1/(2ρ) = 25/4` is where
the measure bound dies.  Below it, covering is *impossible by measure*; the
proof needs nothing but counting (the ρ-parametric sharp visit count).

**Distinctness regime (many combs).**  For `l ≥ 7` the measure exceeds 1, so
measure permits a cover — but the combs are pinned to *equally-spaced arcs
within each frequency*, and DISTINCT frequencies cannot interlock to tile.  The
decisive numerics (lrc_seven_comb_covering_kps_S20b/c): distinct-frequency
combs at `2/25` leave a **positive uncovered floor for every `l ≤ 14`** (0.11
at `l=7`, 0.062 at `l=11`, 0.051 at `l=14`); only REPEATED frequencies (all
`r=1`, 7 shifts) cover.  Since runners are distinct speeds, the primary
(single-lift-class) stratum has distinct frequencies and can never cover.

This is the **Newman/covering-systems shadow**: exact covers force a repeated
largest modulus (Davenport–Mirsky–Newman–Rado); the continuous analog is
opus's `φ > 0`, and it is why the pole does not actually let the gap fill —
past `l = 25/4` the measure would allow a cover, but arithmetic (distinctness)
forbids it.  The mechanism *changes* at the pole, from measure to arithmetic,
but the conclusion (no cover) does not.

## The two strata of the many-combs regime

The `l ≥ 7` residual splits by whether frequencies **repeat**:

  * **Distinct** (single lift-class / k-stratification): `CircleClearFloor`
    (LRCCircleCover) names it; numerically dead through `l = 14`.  My
    `torus_A_window_empty` consumes the floor and outputs gap-emptiness.
  * **Multi-class** (≥ 3 direction-classes, runners parallel-within-a-class →
    repeated frequency): mac-mini HYP-4292, census-clean at **infimum exactly
    `1/6`** (the 5-runner LR bound, at 5-5-2), a factor `2.08` above the
    window.  Repeated frequencies *can* measure-cover `T²`, so this stratum is
    genuinely different — it is the "≥3 coupled ≤5-runner LR systems" problem,
    and the ≤5-runner LR bound (a citation) is what floors it.

Both strata are confirmed safe-above the window.  Both reduce to a covering
lemma with a factor-2 margin, which is why mac-mini can hope a *crude* formal
bound closes it — or why (HYP-4302) the G-K accumulation-from-below theorem
subsumes the whole thing into one 1-D census (no decreasing sequences below an
accumulation point ⇒ the gap has no accumulation ⇒ finitely many values ⇒
census).

## What this says about the shape of the proof

The gap is empty for the same reason at every scale: **you cannot cover a
circle with distinct-frequency combs whose total measure you are forced to keep
near 1.**  The `l ≤ 6` half is a measure triviality; the `l ≥ 7` half is an
arithmetic rigidity; they are welded at the pole `25/4`, which is exactly
`1/(2ρ) = 25/4` — the reciprocal of the window ceiling's tooth-measure.  The
window ceiling `2/25` and the pole `25/4` are the same constant read two ways.

That the crux keeps re-deriving `25/4` — as a gcd-ladder pole, a fee-mean
ceiling, a torus density wall, and now the reciprocal tooth-measure of a
covering threshold — is the tell that these were never four problems.  They are
one covering-impossibility, and the mathematics has been pointing at it from
four directions.

See: LRCTorusSplit.lean, LRCCircleCover.lean; HYP-4247; mac-mini HYP-4292/4302;
opus HYP-4246/4266; the cluster-gcd pole (HYP-4237); the fee-mean ceiling
(MISTAKE-105).
