# Independent audit: Fibonacci-colored addition

PASS on the declared canonical-word and occurrence-dependent branch carriers. No source packet changed.

## Analytic proof review

1. The canonical A and C indicators are consecutive increments of the exact prefix-count floors. Their phase intervals are respectively [1−alpha,1) and [1−2alpha,1−alpha); B is the remaining interval. Positive orbit positions meet none of the endpoints.
2. Conditioning the target phase retains addition exactly: {a alpha} must lie in I_c intersect ({N alpha}−I_d). This set has at most two interval components. The count is not generally the product of densities at fixed N.
3. The discrepancy constant is safe. A q-block for a convergent is within circular distance less than1/q of a shifted q-grid, creating at most six changed memberships at the two endpoints; rational-grid discrepancy is at most two. The Fibonacci denominator ratios are at most two, so greedy removal leaves less than half the preceding length. At most bitlength(L) blocks and two interval components give16·bitlength(L). The proof is uniform in starting point and target phase.
4. The AA convolution is the stated triangular tent with support (1−2alpha,1), area alpha², and height alpha. Outside support no AA representation is possible. Removing at most one diagonal and dividing by two gives exactly unordered distinct pairs.
5. The quadratic equidistribution argument is valid: fixed-shift differencing gives an irrational linear phase with vanishing normalized geometric sum; van der Corput followed by H tending to infinity proves every nonconstant quadratic Fourier average vanishes. The cited author-hosted Tao notes were checked and contain the specified Weyl framework. This is not a prime-weighted equidistribution statement.
6. Trimmed support intervals have a positive tent minimum, so the O(log N) error is eventually dominated by N times that minimum. Equidistribution and passage of the trimmed endpoints to the true endpoints give density exactly2alpha for canonical distinct-A square targets. This argument correctly avoids a universal finite arrival horizon near the support boundary.
7. The full branch family contains canonical A pointwise. Thus all branches have lower density at least2alpha for represented square roots, and liminf normalized AA averages at least the canonical limit. An arbitrary branch need not have a limiting density or average. Equality is established for the canonical branch, not every branch.
8. For Amax (all canonical B promoted), inclusion-exclusion gives H−2c(H)+R_CC(N). The stated floor lower bound exceeds1 for H>=11, and equals2 at H=10. Explicit N7..10 pairs complete the proof. Direct checking confirms last distinct-summand failure6.
9. At one fixed target, reflection a↔N−a partitions indices into a matching plus an optional singleton. Hence the selector polynomial factors exactly: AA contributes z, AB contributes1+z, BB contributes3+z, and unconstrained B choices contribute2. Optional C choices multiply by a scalar but cannot create A. The note explicitly rejects multiplying different target polynomials without retaining their shared selectors.

## Independent exact computation

A separate checker generates colors using Beatty address formulas, not the source interval classifier. It matches the substitution on10,000 positions; verifies exact quadratic-field AA tent discrepancies for targets2..2000; checks all Amax failures through2000; and directly enumerates compatible color assignments for every target2..13 to reproduce the AA-count polynomial. Amax failures found:2,3,4,6. Among square roots2..100, the sole positive-support target without a distinct canonical AA pair is k=2 (the diagonal-only small case).

The original source checker was replayed normally and with `-O`; both results match its frozen receipts, with84,812 checks each. All five source-manifest size/SHA pins were checked. These are finite controls supporting the analytic proofs, not evidence of a Goldbach or Collatz proof.
