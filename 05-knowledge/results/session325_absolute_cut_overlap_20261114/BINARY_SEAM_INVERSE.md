# A target-first binary decoder in a reusable strong seam

The [ordered seam](READY_OVERLAPS.md) also admits a genuinely branching fixed right bank, with the same updated hidden coordinates for both letters. For any strong rooted left source put (H,h,L) as in that note. Let B be a strong right block whose shared root has indegree one, with its unique in-neighbor y in the terminal SCC of B0=B−root, and with B0 not strong. This is the [left-zero seam condition](../session324_busch_inverse_and_gap_structure_20261113/BUSCH_INVERSE.md). Then its actual start count at the root is h_B=H(B0), because every old Hamilton path starts before that terminal SCC, and the seam formula gives

    H'=h_B H+(H(B)−h_B)L,    h'=h_B h.

Let Q_B0=Q_(B0)(2), and let L_B be the outgoing relaxation at the root on the right block. The two optional slots on the ordinal sum of the two root-deleted blocks give

    L'=Q_B0(L−h)+h L_B.

Indeed an out-cut start in the left block permits an arbitrary two-cover on the right; a start in the right block forces the entire left block into the first slot. An empty second slot contributes h h_B.

## Two literal right letters

Take C=(10,28,25,16,1), with H(C)=9 and Q_C(2)=84. The optional-start response at y is 11 for both y=0 and y=3. The optional-end responses at these two vertices are 11 and 23 respectively. Let B0 be an ordinal singleton followed by C. Add root with sole in-neighbor y and all other old vertices as out-neighbors. This right block is strong, with left defect zero, and has

    h_B=9,    Q_B0=168,    L_B=157,    H(B)=31 or 55.

To see L_B, the relaxed second slot can start anywhere except y, hence L_B=Q_B0−S_y(B0)=168−11. The full count is h_B+E_y(B0), with E_y doubled from E_y(C) by the initial singleton and no old path starting at y. The replay checks both right blocks and both seams with a literal left source.

Thus one may independently choose either physical step

    H'=9H+δ L,    h'=9h,    L'=168L−11h,   δ=22 or 46.

The choice changes H but not the subsequent (h,L). At each step all old arcs are retained, and the strongness of both factors plus their common mark gives strongness.

## Exact inverse including its capacity limit

Fix an initial strong source with L_0>h_0 and H_0≤L_0. The two hidden sequences are independent of choices:

    h_i=9^i h_0,
    L_i=168^i(L_0−11h_0/159)+(11h_0/159)9^i.

For depth k≥1 the attainable counts are exactly

    9^k H_0 + Σ_(i=0)^(k−1) (22+24 b_i) 9^(k−1−i) L_i,   b_i∈{0,1}.

The size shell is unique: 22L_(k−1)≤H_k≤55L_(k−1), whereas L_k>157L_(k−1). These ranges are strictly disjoint for successive k. For a prescribed target locate the only possible k, subtract the all-zero baseline, divide its residual by 24, and use descending greedy on the displayed positive weights 9^(k−1−i)L_i. Success requires an integral nonnegative residual and that the final residual is exactly zero. This test is exact because consecutive weights grow by factor greater than 157/9>2, so every weight exceeds the sum of all preceding weights. Successful digits reconstruct a literal sequence of fixed right letters; failure certifies no hit in this two-letter grammar.

There are exactly 2^k values at depth k, but L_(k−1) grows as 168^(k−1), so this grammar gives only O(X^(log_168 2)) counts through X. The genuine binary freedom carries the endpoint state and admits a target-first inverse, yet does not make the integer intervals overlap.

A general obstacle to independent interval tiling is already visible in H'=h_B H+(H(B)−h_B)L. Suppose a left palette varies only H with [h,L] fixed, while the right palette has fixed positive odd h_B and fixed output-update coordinates. Compare choices of right H(B) whose offsets lie in the same residue class modulo 2h_B after removing the common left parity. For gcd(L,h_B)=1, their H(B) values differ by at least 2h_B, giving shift at least 2h_B L. Even if all left counts between 1 and L were physically realized, their output spread would be less than h_B L. Thus these same-residue intervals are forced disjoint. If gcd(L,h_B)>1 even the full odd residue grid is already unavailable in that one bank. To fill gaps by composing seams one must change the active left L response, or use a different coupling, not only increase the number of independent offsets.
