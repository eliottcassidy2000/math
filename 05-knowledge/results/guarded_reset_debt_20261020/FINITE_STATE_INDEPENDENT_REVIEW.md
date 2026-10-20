# Independent review: finite-state debt barriers

Date: 2026-10-20. Reviewed `FINITE_STATE_BARRIERS_ADDENDUM.md` against the guard and ROOT conventions in `RESET_DEBT_COMPOSITION.md`. No base files or repository/remote state were changed.

## Verdict

**Approve the mathematical results.** Theorem 3 now explicitly requires an integer `L>=1`; this correction was reread and verified after the initial audit. The earlier wording had permitted `L=0` when the edge set was empty, which made its witness progression include `n=1` and invalidated only the separate claim that ROOT was not reached. The applied clarification closes that boundary gap without weakening the result. No unresolved mathematical correction remains from this review.

## Proof audit

- **Theorem 1:** Both directions are valid for a finite directed multigraph with positive multipliers and nonnegative debts. The chronological weighting in the telescoping identity is correct. A neutral cycle forces every debt to vanish because all weights in its debt sum are strictly positive. Closed-walk decomposition handles repeated vertices, parallel edges, and loops correctly.
- **Constructive sufficiency:** The finite simple-cycle minimum supplies a strictly positive epsilon wherever needed. The shortest-path inequality has the correct orientation: `r_u-r_v<=w_e`. Exponentiation gives nonnegative edge slack, strictly positive on every positive-debt edge, including such edges outside all cycles. Common scaling therefore pays every debt. Isolated vertices, acyclic graphs, all-zero debts, and neutral zero-debt components cause no exception.
- **Exact rational feasibility:** Scaling a positive barrier to have every coordinate at least one preserves feasibility specifically because debts are nonnegative. The stated rational-polyhedron argument is valid. Floating-point cycle equality must not replace exact rational comparisons.
- **Strict cycles:** Cycle extraction leaves a simple residual path and yields the asserted exponential lower bound on multiplier growth. The acyclic case is vacuous. Integer height consequently enters the finite head at graph/block endpoints; this is not by itself a ROOT theorem or a bound on every intermediate value inside a block.
- **Theorem 2:** The witness sources have exactly the first `L` valuations equal to one, not merely an upper bound on them. Their odd endpoints strictly increase from a source at least three. Thus neither a reset-containing guard nor an early ROOT stop can evade the obstruction. Unit-budgeted physical-height words and genuinely descending words must contain a reset, so both corollaries are correctly scoped.
- **Theorem 3:** With its now-explicit `L>=1` hypothesis, every purported lawful `V`-edge path lies entirely inside the initial valuation-one prefix. Its repeated state creates a closed walk with product below one, contradicting Theorem 1. This works for every initial state, including source-dependent choices. Equivalently, the subgraph of all-ones-labelled edges must be acyclic: each such edge strictly increases its state barrier. Guard restrictions can only remove paths, not rescue coverage.

## Independent controls and scope

Separate exact-arithmetic checks passed 2,000 randomized instances of the path telescoping identity, 384 valuation-one cylinders with prefix lengths from 1 through 64, and all explicitly claimed feasible hostile edge pairs. These are implementation controls, not substitutes for the proofs. The separate LP/cycle checker was not relied upon for this review.

No universal Collatz, arbitrary-source entry, or infinite-bank coverage theorem follows. Unbounded word lengths, additional unbounded state, and a different signed or matrix resource remain outside the finite fixed-word physical-height no-go theorem. The addendum preserves these restrictions correctly.

The subsequently added **Fully rational constructive barriers** paragraph is also approved: the binomial bound validates the integer choice of N, the modified cycle products are at least one, multiplicative shortest paths give positive rational reciprocals with strict slack exactly where needed, and rational common scaling supplies the claimed explicit certificate without logarithmic numerics.
