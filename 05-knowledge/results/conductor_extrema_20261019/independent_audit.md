# Independent audit: exact conductor-defect endpoints

Date: 2026-10-19. **Verdict: PASS in the stated straight-line, characteristic-zero, U(3,r) scope.** No correction to the theorem is required. The intermediate attainable defects and general branch-conductor vectors remain unclassified, exactly as the note says.

## Proof checks

1. **Universal bounds.** Products of at most q pair/singleton lines separate each point within every selected subset of size at most 2q+1. Uniformity ensures a pair's joining line avoids the selected point. Hence min(r,2q+1)≤h_q≤min(r,binom(q+2,2)). Summing r−h_q gives precisely the two claimed endpoint formulas. At r=3 the empty lower sum is zero, and the printed closed formula also gives zero.

2. **Equality bookkeeping.** The two differences from the extremal defects are sums of nonnegative rank discrepancies, so equality is termwise. The adjacent-degree minimum test is sufficient: multiplication by a linear form nonzero at all points propagates injectivity downward and surjectivity upward. There is no cancellation argument hidden in the equality characterization.

3. **Rational minimum witness.** The sparse-positive-root lemma is valid by induction and Rolle's theorem. It proves nonsingularity of every square generalized-Vandermonde minor on positive distinct nodes and distinct exponents. For G_r, the degree-q exponents a+rb are distinct when q<r; s≤r−1, so this covers every degree needed to reach surjectivity. Multiplication by x supplies all subsequent degrees, including those where exponent collisions do occur. Nonzero rational minors remain nonzero over every characteristic-zero extension, so the real-root proof does not impose an ordering on the ground field.

   At degree one the exponent set {0,1,r} proves U(3,r). Before degree s, deleting any branch leaves at least as many nodes as monomial columns, still with full column rank; therefore no premature branch separator exists. Thus every conductor exponent of this specific witness is s.

4. **Maximum and geometric equality.** Restriction to a smooth conic gives all binary forms of degree 2q, hence h_q=min(r,2q+1) and every separator degree is ceil((r−1)/2). Conversely, for r≥6 upper equality forces h_2=5, so a nonzero quadric vanishes on all points. A reducible plane conic contains at most four no-three-collinear points. Thus the conic is smooth, and the one-dimensional quadratic kernel also gives uniqueness. This proves the claimed equivalence, not merely a sufficient condition.

5. **Full Hilbert function versus branch conductors.** D7 and G7 both have ranks 1,3,6,7,7,… and minimum defect five. In D7 the quadratic xz−y² separates the displaced branch, while a quadratic cannot vanish at the other five conic points without vanishing identically on the conic. Thus its conductors are (3,3,3,3,3,3,2), whereas those of G7 are all three. This is a valid counterexample to recovering the branchwise conductor from the full Hilbert function. It does not contradict the common surjectivity threshold rho=s.

6. **Scope and attribution.** The conductor identities concern cones over projective points, not curved branches with the same tangent data. Seminormal and normalization defects differ by r−1. The asymptotics and small-value table agree with the formulas. The note appropriately identifies the interpolation/Hilbert-function machinery as classical and makes no literature-wide novelty claim.

## Verifier audit and independent controls

The producer's script has explicit exception gates, exact integer/Fraction arithmetic, no randomness, and two genuinely different rank implementations. Its row-deletion test is correct: e_i belongs to the evaluation image exactly when deleting row i lowers the rank by one. Its threshold checks suffice because the normalized evaluation spaces are nested. The advertised 224 full matrices are the distinct degree-loop matrices; the additionally repeated threshold matrices do not enlarge that distinct-matrix count.

A fresh producer replay is retained as `independent_producer_replay.json` and compared with the supplied output. The separate `independent_audit_verifier.py` imports no producer code. It checks:

- Both closed formulas against direct sums for r=3,…,10,000.
- Both witness families for r=3,…,30: 560 full-rank and 924 row-deletion modular certificates.
- Every modular certificate reaches the proved rational-rank upper bound, so it is an exact rank certificate, not a numerical approximation or probabilistic test. The prime moduli are independently tested for primality. The first prime suffices in all 1,484 cases.
- Exact rational ranks and every first separator degree for D7 and G7.
- A third family, (1,t,t³), at r=6,7,10,13, including genuine interior defects 12 at r=10 and 22 at r=13. These are finite controls, not a classification of intermediate attainable values.

The independent program is run normally and under Python optimization, with byte-identical JSON outputs. `independent_audit_receipts.json` records the source/output hashes and replay comparisons; `INDEPENDENT_SHA256SUMS` hashes the audit artifacts. All audit files are confined to this staging directory. No checkout or Git writes were made.
