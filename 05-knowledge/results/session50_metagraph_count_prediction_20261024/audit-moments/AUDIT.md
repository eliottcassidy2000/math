# Independent audit: finite odd-fibre moment reconstruction

## Verdict

**PASS**, with the statistical/computational scope below. The proof and exact computation recover the order-eight ordinary tournament-isomorphism count **6880** from eight supplied word moments. This is not an order-only formula, a previously unknown A000568 term, or a demonstrated acceleration of enumeration.

The audit independently reconstructs and checks every one of the 3,528 source atlas rows. Normal and optimized Python runs and an actual fresh-directory manifest-only replay are provided. The mathematical checks do not use `assert` and remain active under `python -O`.

## Normalization and inputs

For positive order n, fix a directed path on all n owners. There are N = 2^m free words, m = (n−1)(n−2)/2. Each ordinary isomorphism class C has word fibre h_C = H(C)/|Aut(C)|. The automorphism action on directed Hamiltonian lists is free, which proves the quotient formula. The classical odd-Hamiltonian-path theorem and the odd order of a tournament automorphism group make every h_C a positive odd integer. Odd group order follows because an involution would swap two owners and reverse the arc between them. The inherited fibre theorem explicitly records these parity facts; the all-order Hamiltonian parity theorem is inherited mathematical input, not inferred from the finite tests.

A converse-merged node has mass μ = h for a self-converse class and μ = 2h for a distinct converse pair. In particular μ modulo 2 identifies the multiplicity, and h is μ if μ is odd and μ/2 otherwise. Our eight-moment extraction uses the numeric merged masses alone, not the class labels or their number:

K_r = sum_words h(word)^r = sum_classes h_C^(r+1)
    = sum_merged_nodes μ (μ / (1 if μ odd else 2))^r.

The shift is important: K_0 = N, whereas the desired class count is the negative word moment sum_words 1/h. Neither μ nor a transition degree can be substituted for h in the higher moments.

Our independent `reconstruct.py` accepts only the two keys `n` and `moments`. It rejects side-channel fields such as target, hmax and atlas. The producer's narrower list-only moment input is inspected separately. An actual isolated process is run with only our decoder and its moment JSON present, with no atlas or count oracle.

## Congruence proof

Let Q_k(x) = product_{j=0}^{k−1}(x−(2j+1)) = sum_{j=0}^k a_j x^j. For positive odd h,

Q_k(h) = 2^k k! binom((h−1)/2,k).

When (h−1)/2 < k the binomial and the product are both zero. Therefore Q_k(h) is divisible by 2^e, where e = k + v_2(k!) = 2k − popcount(k). Summing over ordinary classes gives

a_0 A + sum_{j=1}^k a_j K_(j−1) = 0 mod 2^e.

The constant a_0 is odd and hence invertible modulo 2^e. The first k word moments determine the displayed residue for A. If e ≥ m, the interval 1 ≤ A ≤ N contains exactly one integer in the residue class. This includes modulus N: the interval has N consecutive integers, and residue zero is represented by N, not by zero. The theorem is stated for n ≥ 1; the empty tournament needs a separate base convention.

For n = 8, m = 21. Twelve moments suffice for the unconditional interval [1,2^21], since e(11)=19 and e(12)=22. Eight moments suffice here only because the next rigorous moment bounds shrink the interval.

## Global exact-rational bounds, no support cap

Set q_L(x) to the degree-three polynomial with coefficients C_L^−1 b_L, where (b_L)_i = K_i and (C_L)_ij = K_(i+j+1), 0 ≤ i,j < 4. Define P_L(x) = 2q_L(x) − xq_L(x)^2. Its pointwise error is

1/x − P_L(x) = (1−xq_L(x))^2/x ≥ 0 for every x > 0.

For the upper bound, let q_U have degree two, with (b_U)_i = K_(i+1)−K_i and (C_U)_ij = K_(i+j+2)−K_(i+j+1), 0 ≤ i,j < 3. Define

P_U(x) = 1 − 2(x−1)q_U(x) + x(x−1)q_U(x)^2.

Then

P_U(x) − 1/x = (x−1)(1−xq_U(x))^2/x ≥ 0 for every x ≥ 1.

Thus L = sum_words P_L(h) ≤ A ≤ sum_words P_U(h) = U. The respective polynomial degrees are seven and six, using K_0…K_7 and K_0…K_6. Exact normal equations give L=b_L^T C_L^−1 b_L and U=K_0−b_U^T C_U^−1 b_U. The matrices are Gram matrices for the actual nonnegative word measure, and they are nonsingular for this input. These pointwise bounds are valid for any q; inversion is only the optimization that chooses q. A general implementation must handle singular matrices rather than claim they are always invertible. The producer does so by stopping at the largest invertible leading block.

The independent audit checks both square identities coefficient by coefficient over rational numbers. Spot tests are only extra controls, not the reason the inequalities hold. No upper support endpoint, observed maximum h, mass-one multiplicity, or histogram enters this bound certificate.

The exact endpoints are:

L = 27035513046341570927644802943075201007181614224844706720 / 4035942584549039167511930263776554850885060915100907

U = 547789015950209817357891493316406828364640 / 20793022244349980456910650230381182077.

These imply the integer interval [6699,26344]. The eight-moment congruence is A = 6880 mod 32768, so the only candidate is 6880. The wider n-only word bound and all decimal approximations are unnecessary for this final uniqueness comparison.

## Shorter order-eight certificate

A separate self-contained decoder avoids the optimized rational coefficients entirely. Fix

q(x) = (x² − 1250x + 465000) / 50000000,

and use P(x) = 2q(x) − xq(x)² + (1−xq(x))². Exactly as above,

P(x) − 1/x = (x−1)(1−xq(x))²/x ≥ 0 for every x ≥ 1.

From K_0,…,K_6, exact evaluation gives

A ≤ 496641366973194539 / 15625000000000 < 32768.

The integer range is therefore [1,31785]. The already verified eight-moment congruence selects 6880 as its only representative in this range. No lower polynomial, matrix inverse, fitted support maximum or small-fibre count is needed. `simple_certificate.py` contains this fixed polynomial and no target constant; it is independently replayed with only the eight integer moments. Its coefficientwise square-identity check ensures that choosing convenient coefficients cannot smuggle in an unproved bound.

## Limited sharpness claim

The exponent e(k) is optimal among degree-at-most-k integer-coefficient polynomials F with odd F(0), when requiring every positive odd value F(1), F(3), … to vanish modulo the same power of two. Write f(t)=F(2t+1). Newton interpolation gives

F(0) = sum_{j=0}^k Δ^j f(0) binom(−1/2,j).

If every positive odd value of F is divisible by 2^b, every finite difference here is also divisible by 2^b. The reduced denominator of binom(−1/2,j) is exactly 2^(2j−popcount(j)). These exponents increase with j. If b>e(k), each term on the right is an even integer, contradicting odd F(0). Q_k achieves e(k), with equality of valuation already at h=2k+1.

This does not prove that eight moments are globally minimal for the actual tournament count. It does not rule out methods using other statistics, discrete-support feasibility, nonlinear identities, or additional tournament structure.

## Independent finite controls and leakage boundary

- Enumerated every fixed-path word for n=1,…,6, independently canonicalized ordinary classes, and checked H/Aut using both subset dynamic programming and direct Hamiltonian-permutation counting. Predictions obtained before comparison are 1,1,2,4,12,56.
- Recomputed H, automorphisms and actual converse isomorphism for all 3,528 order-eight atlas representatives. All source masses pass the normalization; 176 merged nodes are self-converse.
- Checked explicit frozen adjacency matrices to fix the atlas's MSB-first lexicographic-pair convention.
- Tested the falling-factorial identity and Newton denominator valuation through k=40, integer-polynomial sharpness controls through degree four, and congruence reconstruction on synthetic odd-mass multisets.
- Generated the eight input integers independently from μ values; they agree with the producer's frozen moment list. The unique result is compared with the separately computed atlas class count only after inference.
- Replayed a decoder in a fresh directory containing only code and the two-key moment input. Class counts, IDs, number of atlas rows, support caps, and the atlas itself are absent.

The preparation stage nevertheless uses a complete, previously materialized atlas. That source already encodes the class count, and this audit does not claim that summarizing it predicts an unseen A000568 value from n alone. Its contribution is a valid finite-aggregate reconstruction theorem and a much smaller sufficient statistical certificate. Whether the needed moments can be obtained without building the atlas, or faster than direct class enumeration, remains open here. Likewise, the audit's full order-eight quotient completeness is inherited from the source atlas; it rechecks every row and the total mass but does not independently enumerate all 2^21 words at order eight. The n≤6 quotient checks are completely fresh.
