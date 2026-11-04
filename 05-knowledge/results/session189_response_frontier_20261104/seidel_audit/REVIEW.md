# Independent review of positive Seidel energy and carousel separation

**Status: APPROVED for the exact profile-restricted all-order theorems in the frozen source.** The core resolvent, fixed-point, energy, dyadic approximation, ambient, substitution-closure, domination and carousel statements are sound. This is not approval of unrestricted arbitrary-profile reconstruction or an uninflated-source consequence of a thickening operation.

Reviewed source `PROOF.md` SHA-256: `204d7cd868d3d83e019a20ea5bb6254ee455e5f0ecaa01a6635c55c7afecd3d0`.

## Mathematical review

1. For K=(S+J)^(-1), the identity K+K^T=2(K1)(K1)^T implies h=h^2. Reduction modulo two correctly selects h=0 for even induced sets and h=1 for odd ones. Every defining inverse has odd determinant. The even Sherman–Morrison formula and the odd null-vector formula retain their signs and the factor v=1+2t correctly.
2. I independently expanded the two port orders, all prefix/suffix factors and both direct transitions. Their factorization yields W_even=v(1+2A)^2(1+2C)^2 and W_odd=−v^(-3)(1−2A)^2(1−2C)^2. The empty subset has W=v, not one away from t=0; the source wording was repaired before this freeze and the formulas already had the correct factor.
3. Inclusion-exclusion contributes (−1)^|I|2^(−m), cancelling the odd-case negative sign and giving the positive Gram sum. The no-port computation gives exactly1/(1+2t), proving every scalar ordered-cover coordinate is fixed. This is a formal rational/2-adic fixed point, not a realizable finite tournament profile or a real limiting Hamilton count.
4. Boolean extraction of [a_i a_j c_i], with i!=j, gives the factor2^(6−m) and the joint term r_i^2 r_j. Contracting with S_ij gives the positive sum of squared inverse-vector entries. No independently chosen endpoint marginals are multiplied. At least one cubic is nonzero because the energy is strictly positive.
5. Integrality is justified by the convergent integer run functional, despite the displayed power-of-two denominator in the Gram expression. The dyadic TT_N congruence follows from exact surjection counts; odd bases have the required unit exponent and even bases have sufficient valuation. If a normalized profile agrees with g modulo2^L, its k-th coordinate has valuation at least min(k−1,L). Thus only finitely many total-excess shapes contribute modulo the modulus, uniformly across cover coordinates.
6. The contrast is exactly linear in the cut variable. Ambient substitution therefore multiplies its Hamilton contrast by c1+2c2=H(R[TT2]), an odd number. This is a full-cover argument; it does not use an invalid generic lift of a scalar modular gap.
7. Substitution closure divides by the common odd Hamilton factor correctly. The profile-dominating-target argument uses a negative baseline member of a converse-complement pair and nonnegative run coefficients. The transitive-profile minimum follows from the exact partition formula. The prime-root original-module qualification correctly excludes crossing quotient modules at composite roots and all-pairs output reconstruction.
8. For the carousel, the switching order and alternating diagonal signs are correct. The even transitive skew inverse gives the stated coordinates. The independent-sign encoding is a bijection on the even subsets containing the two distinguished endpoints. Its second moment yields exactly−8(n−4), including the n=3 empty-sum boundary. The modulo-sixteen claim for transitive orders divisible by four follows directly from odd fourth powers and even-base vanishing, not from numerical extrapolation.

## Independent computational checks

All four normal/optimized pairs below passed byte-identically; the producer replay also matched its original frozen receipt.

- **Complete small geometric response test:** `audit.py` uses Cramer's-rule determinant expansions rather than the producer's Gaussian inverse solver. A separate C++ cut-aware multiset engine sums geometric run profiles modulo256. Their results agree for all seven quotient isomorphism types at orders2,3,4, both orientations, every Boolean query, and the first five ordered-cover coordinates:2,336 complete query records. A full labeled enumeration independently verifies the small quotient coverage. The test also checks the scalar fixed point, literal eight-corner shared-cubic extraction, positive energy, and every higher contrast coordinate.
- **Independent carousel subset check:** `carousel_check.py` uses fraction-free Bareiss determinants and Cramer's rule on the actual original-label Seidel principal matrices, without the source's switching or sign-average reduction. It checks all5,461 relevant even subsets at orders3,5,7,9,11,13,15 and recovers the exact shared cubic−8(n−4). Complete subset terms and fixed names are retained.
- **Actual finite physical tournaments:** `physical_check.py` counts3,994 literal augmented tournaments using vertex-subset ordered-cover DP, through augmented order14. It retains complete named Hamilton banks and all named coefficient coordinates for C3 with pair children and a two-name quotient with TT4 children, totaling139,264 full physical response values. It separately tests every coarse query for C3 with TT4 children. All288 source/target coarse physical queries agree with the rational Gram functional modulo4 or16 as appropriate, in all first five cover coordinates. These are full named banks only for the stated old orders6 and8; the old-order12 control is a complete coarse bank, not a claimed full physical bank.
- **Producer replay:** the frozen `verify.py` was rerun normally and under optimization, reproducing its exact receipt and rational banks. The producer's finite calculations are distinguished from the unbounded algebraic proof.

All full coefficient banks retain zero, cubic and quartic entries, physical names, aligned source/target rows, and exact integer or rational values. The formal geometric banks are not described as physical banks of large inflated tournaments.

## Reproduction and files

Run `python3 audit.py`, `python3 carousel_check.py`, and `python3 physical_check.py`, and repeat each with `python3 -O`. The programs are self-contained within this audit packet, use only the standard library and an already installed C++17 compiler, and perform no installation. Source/data SHA-256 pins and byte lengths are in `MANIFEST.json`.

No further mathematical repair was required. The separately drafted global thickening integration is not a dependency of this review and requires its own scope check.
