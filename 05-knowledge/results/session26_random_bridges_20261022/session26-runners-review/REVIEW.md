# Independent review: runner clearance height curvature

Verdict: PASS, scoped to the precise statements in RESULT.md. No proof-breaking defect found. This is not a proof of LRC(14), a literature-novelty determination, or an obstruction to all reconstruction algorithms.

## Independent mathematical verification

For each affine segment, pushforward density is 1/|s| on its height interval (l,h). Its contribution to signed height curvature is (delta_h-delta_l)/|s|. The superlevel length is therefore sum c_h(h-delta)_+ for delta>=0. The complete signed atomic measure has sum c_h=0 and sum c_h h=1. On positive heights it is the second distributional derivative of superlevel length. To interpret its zero atom literally as a second derivative, extend superlevel length as 1 for negative thresholds; this yields c_0=-density(0+). The text's sidecar interpretation is correct without requiring that extension.

Direct integration gives A_k=sum c_h h^k=k integral f^(k-1), k>=1. There is no 0^0 issue: A_0 is excluded. The recurrence and rational generating function concern A_k, not raw moments. For example raw single-runner moments are 2^(-j)/(j+1), whose generating function is -2 log(1-z/2)/z, not rational.

At M every incident segment has M as its upper endpoint, so all contributions are positive. Equality level sets have measure zero because every segment slope is nonzero. The full distribution detects the support endpoint despite zero mass exactly at M. This is not inference of a maximum from one mean.

Equal-height aggregation is essential and performed by the source. After aggregation, each nonzero height gives a simple pole at z=1/h with analytic residue -c_h. The largest height corresponds to the nearest positive pole, not the largest pole location. Distinct poles cannot cancel. A_k for k>=1 has exactly r minimal constant-coefficient recurrence order by the Vandermonde argument; equivalently its r-by-r Hankel determinant is a nonzero signed coefficient product times a squared Vandermonde determinant. A_0 including the zero atom is not part of this minimal-order assertion.

For N=2m, independent reflection/tooth analysis gives the stated positive atoms at k/(N+1), k=1..m, weight 2+2/N; negative atoms at k/(N-1), k=1..m-1, weight -2+2/N; and c_0=-4+2/N. Coprimality of odd N-1,N+1 makes the two height sets disjoint. Thus r=N-1 exactly, M=m/(N+1). The result rules out a runner-count-only bound for this exponential recurrence representation, not arbitrary finite nonlinear reconstruction.

For (1,4), atoms are (0,-7/2),(1/5,5/2),(1/3,-3/2),(2/5,5/2). Mean clearance is 1/6 while maximum is 2/5. The negative knot reflects a monotone excursion-boundary speed change, not disappearance of a connected superlevel component.

Integer dilation preserves the whole pushforward distribution via the circle covering map. Duplicate speeds are harmless. The stated finite nonempty positive-integer domain excludes empty rows (undefined minimum) and zero speed (identically zero f, requiring separate handling); negative speeds can be reduced to absolute values but are not needed. Arbitrary noninteger dilation on the fixed unit interval is not claimed. There are no plateau degeneracies in the stated domain.

## Source review

probe.py constructs rational tent-wall intervals, clips each affine runner branch against all others, resolves duplicate affine branches by index, merges matching adjacent segments, and checks partition length. Its Fraction arithmetic is exact. family_check.py checks all even N through 80 and independently sweeps threshold walls of the original distance inequalities. Explicit raises remain effective under python -O. The scripts do not directly compute recurrence rank; the all-N proof, rather than the finite output labels, establishes minimal order.

Normal and optimized source reruns and artifact comparisons are recorded in adjacent files. independent.js uses a different envelope construction: enumerate every pairwise intersection within each tent-wall interval, select the minimum at each subinterval midpoint, and integrate unmerged pieces with exact BigInt rational arithmetic. It checks moments of orders 0..10, both atomic normalizations, terminal positivity, and superlevel lengths at every midpoint between consecutive atom heights. It also checks the even-family atom formula through N=40, 320 ordered-with-repetitions triple controls and their doubles, the original 162-row census, and both large controls. Cross-comparison uses full atom maps, not just maxima or counts.

## Optional editorial refinements

- Opening phrase "rational moment generating function" could explicitly say "normalized-moment ordinary generating function" to exclude confusion with raw moment series or exponential moment-generating functions; Section 3 already resolves this.
- JSON field terminal_residue stores c_M, not the complex analytic residue, which is -c_M. Rename to terminal_atom_weight for precision; this does not change a result.
- "Positive-height poles" informally labels poles by their heights. Pole locations themselves are reciprocals. Section 3 correctly states the relationship.

No repository or Git writes were performed by this audit.
