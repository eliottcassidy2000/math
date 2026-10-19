# Independent audit: one-sided three-slave theorem

Date: 2026-10-19. Result: approved.

The main proof was checked independently using a different shuffle derivation. Instead of deleting y, retain all special letters x,a1,…,am,y. Choose one of the m+2 placements of y in the fixed x,a1,…,am order. Assign the special letters, in that order, to the m+1 chain-gap blocks, enforcing that aj occupies block index at least j. There are 8,30,112 admissible order-and-block assignments for m=1,2,3.

For each assignment, the product of weak-gap binomial polynomials counts its extensions. If i_x and i_y denote common-chain counts before the roots, the wordwise aggregate deficit is c+3(i_x−i_y). Replacing these ranks by the exact weak-gap order-statistic means gives a direct rational polynomial. This construction is independent of the manuscript's 2,5,14-cut formulas obtained by first deleting y.

For m=1,2,3, the three resulting deficit polynomials were expanded and matched coefficient by coefficient against the main verifier's receipt. All matched exactly. The positive monomial counts for m=1,2 are 13 and47. For m=3, the only negative monomials lie on p=s=0; that boundary factors as

    (t+1)(q+t−1)[6q³+6q²t+18q²+4qt²+14qt+12q+t³+5t²+6t]/24.

It is nonnegative on the required integer domain. The remaining145 coefficients are positive.

The conditioning argument was checked with cumulative thresholds, including propagation of predecessor relations along the imposed upper-slave chain. The exterior tree-deletion identity preserves the family: no lower slave can be created, and deleting upper slaves decreases their number. All zero-length gap endpoints are covered by the binomial-polynomial convention.

The c=0 suffix certificate was independently discovered by direct core-word enumeration. Every projected maximum-span order for m≤3 has a single forward incomparable pair whose 3N−2E is nonpositive in every suffix generator. The main verifier's enumeration has exactly 1,2,6,48 orders for m=0,1,2,3. Projection preserves the root span because no exterior point is between the roots, so this finite test covers actual exterior realizations. This certifies the c=0 exclusion without relying on the aggregate's empty sum.

The four-upper-slave counterexample was checked independently by whole-poset subset DP. For x<a1<a2<a3<a4,z<a4,y isolated, E=35, N(x<z)=28 and N(z<y)=20. Thus the aggregate extension really fails, even though (z,y) is balanced. The manuscript correctly distinguishes this obstruction from a counterexample to the 1/3–2/3 conjecture.

No mathematical statement or proof mismatch was found.
