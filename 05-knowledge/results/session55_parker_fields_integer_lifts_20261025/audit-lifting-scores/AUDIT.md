# Independent audit: integral lifting and residual scores

Status: symbolic proof audit PASS; finite controls as recorded in normal.json and optimized.json. No remote edits, novelty assertion, or resolution of the integer square-magic problem.

## Integral lifting

The parametrization is valid over arbitrary commutative unital rings: adding the two diagonals and middle column and subtracting the top and bottom rows gives T=3c without division. All opposite entries then total 2c, and outside columns recover the other coordinates. The six nonpivot unit vectors give a right inverse to S, proving surjectivity and an integral direct sum, not merely rational rank. The exact carry condition follows by integer cancellation of the nonzero modulus; positivity and distinctness are additional predicates. The p59 grid has ordinary sum 48, with only 4 and 16 integer squares. This is a genuine hostile to residue/no-wrap sufficiency, not an integer solution.

The interval claim requires strict diameter less than m, or one excluded endpoint at diameter m. Distinct nonnegative entries force positive line sums, justifying (0,m]. Zero in an integer solution forces a positive square equal to twice another positive square, impossible by 2-adic valuation. Opposite entries therefore lie strictly between zero and twice the center square. The least quotient modulus is exactly the least common sum under the stated at-most convention. Prime equality would require e=1 and p=3; nine distinct entries below 2 cannot occur. This proves the conditional minimum-height equivalence, without assuming existence. A field extension's cardinality is not the characteristic and its prime subfield has only p possible integer images; characteristic two cannot preserve nine distinct integer entries. The finite-bank hostile expressly does not preserve local distinctness.

## Residual scores

R2 is the exact tail-assignment matching criterion. For any left edge collection A its endpoint set S has all r(S) slots as its neighborhood; the stated cut inequality supplies at least |A|. Total equality saturates all slots. Complementary inequalities give upper bounds and precisely the stated absolute imbalance formulation.

For even n>=6 and k<=n/2, n-k-2>=1. Deleting a cycle removes at most 2k crossing edges, so cut capacity dominates k, and every balanced prescribed +/-1 set works. For a path, its endpoints reduce the removed crossing bound by epsilon and increase the imbalance bound to k+epsilon; both changes have the correct signs. The n4 cycle matching requires one high vertex per edge. The n4 path complement is a tree, whose orientation toward its unique low vertex realizes all four patterns. The odd path excess-unit proof has |b(S)|<=2 and cut capacity >=2 for n>=5; n3 needs an endpoint.

Restoring the fixed directed path adds one outgoing arc except at its end. Balance consequently requires the end high. Native complementation transports endpoint full degrees by +1 at the start and -1 at the end relative to reversal. Both source and partner balanced is exactly both endpoints high within this slice, leaving p-1 high interiors, hence binomial(n-2,p-1). Conversely internal full degrees yield precisely this slice. No counts of individual orientations or unlabelled objects follow. The regularizing parent has order n+1, not n; its cycle residual is Eulerian, distinct from the even carrier.

## Independent controls

check.py does not import either producer checker. It enumerates every orientation of the path and cycle complements for n4 and n6, compares prescribed score realizability including the sharp exception, tests all subset inequalities for every claimed prescribed pattern through n10, and enumerates every 3x3 array modulo 2 and 3 to test the characteristic-independent parametrization against line sums directly. Source SHA256 pins are embedded in both output certificates. Run python3 check.py and python3 -O check.py and compare JSON byte-for-byte. Producer carry verifier is separately run in both modes. General proofs, not finite tests, justify all-size claims.
