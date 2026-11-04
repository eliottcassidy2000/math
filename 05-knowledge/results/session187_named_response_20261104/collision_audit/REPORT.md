# Independent audit: transitive-module order separation

**PASS after a scope-wording correction.** Reviewed the complete proof and author checker without importing its implementation. No theorem or numerical correction was needed. The finite ambient universe is 75 labelled quotients with root fixed at 0, covering arbitrary rooted quotients through order four up to relabelling, not every labelled choice of root separately. The author has corrected that wording.

## All-order proof

The literal named query makes the old transitive module and both ports one enlarged module in the unchanged rooted quotient. The color-class argument gives the stated polynomial: when the exceptional old vertex has y's color, a class with k preceding vertices has 2^k+1_[k>0] Hamilton words. Later vertices contribute no extra choices. Multiplying by the source singleton x is legitimate in the color polynomial.

The consecutive-rank difference is 2 times a sum of nonnegative monomials. Conversion through nonnegative Stirling numbers gives coefficientwise monotonicity of every ordered-cover coordinate. Its first coordinate is strictly positive, with the stated values 2 at i=0 and 2^i at i>=1. The occupied-run coefficients of every rooted quotient are nonnegative and the single-run coefficient is H(Q)>0. This proves strict rank separation in arbitrary fixed quotient contexts, independently of any bound on run number. Sorting the r fixed physical queries therefore recovers the transitive module order.

The three-state gate evaluation formula follows by converting q1,q2,q3 to P(1),P(2),P(3). Its potentially negative evaluation coefficient does not undermine the preceding positive-cover proof. The proof correctly limits the theorem to the same exterior quotient and root incidences. It does not assert universal reconstruction, arbitrary-module separation, or a full-response collision.

## Executed verification

- Reran author verify.py normally and under -O. Both outputs agree byte for byte with each other and the frozen author receipt. This reproduces 21 complete profile controls, 750 actual ambient queries, 1,452 interpolation-query checks, and 45 actual gate-query checks.
- Independent check.py imports no author implementation. A predecessor-based Hamilton DP evaluates all 4,096 actual attachment assignments for each of T, U and the literal whole converse of T: 12,288 counts, not merely the interpolation battery.
- Boolean inversion of those complete truth tables exactly reproduces every serialized named contrast coefficient. The arbitrary comparison has degree four and 290 nonzero coefficients; the whole-converse comparison has degree three and 200 nonzero coefficients.
- Independent directed-reachability and pair tests establish that both six-vertex controls are strong and have no common homogeneous pair.
- Literal color-assignment enumeration, using Hamilton permutation counts for every induced color class, checks 45 controls: every rank at r=1,...,5 and t=1,2,3. All agree with the color-polynomial formula.
- The independent checker also reproduces byte-identical results under -O. Its gates use explicit exceptions. Source SHA-256 digests are retained in normal.json and optimized.json.

These finite computations audit the controls and formulas. The arbitrary-order theorem rests on the reviewed color and occupied-run counting arguments, not extrapolation from the finite data.
