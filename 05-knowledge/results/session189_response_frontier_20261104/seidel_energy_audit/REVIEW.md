# Independent review of the Seidel-energy and dyadic-inflation theorem

**PASS, including finalized producer Sections 1–8.** The frozen producer proof has SHA-256 `204d7cd868d3d83e019a20ea5bb6254ee455e5f0ecaa01a6635c55c7afecd3d0`; its verifier has SHA-256 `bc19002e71f92505237a74633ba01361ab4420ee1c72be4e8575416b85f4e9a6`. This review is independent of the producer's implementation. The only defect found so far was a sentence assigning value one to the empty-subset two-port cut functional. Its value is 1+2t; the author corrected the sentence, and the displayed formulas already used the correct value.

## 1. Main all-order argument

The theorem is valid for every finite quotient of order at least two, under the stated geometric-profile congruence hypothesis. It is not a theorem for arbitrary child profiles.

The run series converges in the 2-adic topology coefficientwise and uniformly modulo each fixed power of two over all cover coordinates. A shape of excess e has a factor 2^e. Modulo 2^L only e<L is relevant; every retained word has at most m+L+1 symbols including the two ports, so even the possible number of cuts is bounded. This justifies the inverse expansion, inclusion-exclusion, and the complete-profile ambient transfer without discarding an uncontrolled tail.

For an even induced subset, the skew Seidel matrix is invertible modulo two: it is I+J, whose square is I. For an odd subset, I+2B is invertible modulo two. Thus every displayed r-vector has odd denominators. The formulas K1=r, K^T1=-r in the even case and K1=K^T1=r/(1+2t) in the odd case have the stated signs. They give the two port orders and cut contributions exactly. The empty subset contributes 1+2t.

The inclusion-exclusion sign and the factor 2^(-m) are correct. The unported functional is exactly 1/(1+2t), so the formal profile g is a fixed point for every nonempty quotient. This is a rational identity interpreted in the 2-adics, not a limit in the real topology or a claim that a finite tournament has profile g.

The source-minus-converse shared cubic is

    K_ij = 2^(6-m) sum_(I even, i,j in I) r_(I,i)^2 r_(I,j).

In particular the factor is not 2^(5-m) or 2^(7-m). Multiplying by S_ij and summing yields a strictly positive rational sum of squares. Two-name subsets already supply strict positivity. Although the displayed prefactor can have a power-of-two denominator, the complete coefficient is a 2-adic integer by the convergent run identity. The same holds for the energy, so L=1+v2(E) is a positive finite integer.

The odd-subset contributions cancel under converse, leaving the exact all-cover contrast (1+2t) times its Hamilton coordinate. The Boolean contrast has degree at most three. Its vanishing quartics refer to the full child-uniform formal bank, not to every individual physical-vertex coefficient of a finite inflated tournament.

The transitive approximation follows from the surjection formula with the stated sufficient bound. Even bases vanish; odd bases have exponent dividing 2^(L-1); k>N is harmless because k-1>=N>=L. Profile matching gives the same odd normalization scalar on both sides. Nonzero cubic finite difference implies an actual corner query, without assuming a coefficient itself is a realizable query.

For the ambient, the augmented core remains a common module precisely because both ports imitate its original exterior incidences. The full profile contrast is (delta,2 delta,0,...), modulo the chosen power of two. Therefore the ambient multiplier is c1+2c2=H(R[TT2]), which is odd. There is no replacement of the complete profile by the Hamilton coordinate alone. The same selected core query survives every fixed ambient.

## 2. Independent finite checks of the proof interfaces

`audit.py` imports no producer code. It uses exact rational Gaussian elimination and independently expands both orientations' complete bidegree-(2,2) coefficient banks. Its controls cover every labeled quotient of orders two through four and selected quotients through eight, 81 quotients in all. Both source and converse banks retain all zero, cubic, and quartic coordinates. There are 13,584 named positions in each orientation.

The same program checks 17,199 inverse/segment identities at rational t values, every individual shared cubic constant, the energy identity, exact cancellation of quartics and odd-subset banks, 648 finite-difference corners, and 3,192 transitive-profile congruences including coordinates beyond the child order.

`run_walk.cpp` is an independent finite-word automaton, not a matrix inverse or fixed-multiplicity kernel. It tracks the set of children already used, the two singly occurring ports, repeat excess, last symbol, and cut count. A repeated old child costs -2, while a new child costs one. Uncut repetitions of the same child are prohibited; cuts permit them. It matches all 53,984 complete profile coordinates in 4,704 source/converse query evaluations modulo the sufficient modulus. Full input and output are retained.

`physical.cpp` counts actual labeled permutations with path cuts by subset dynamic programming. `physical_audit.py` retains six complete physical named response banks, totaling 4,218 coefficients, for three finite inflated source/converse pairs. These are actual physical banks, separately distinguished from the formal child-uniform banks.

The literal nonzero control is C3[TT4,TT4,TT4], with twelve old names and the two actual ports. The query cube on block variables {a0,a1,c0} contains corner (3,1), with contrast eight modulo sixteen. The entire first-five cover contrast is (8,0,0,0,0) modulo sixteen there. Higher coordinates vanish because h! is divisible by sixteen for h>=6. A cyclic exterior ambient is retained on the same physical core and ports, reaching sixteen augmented vertices.

This TT4 control uses the sharper individual-cubic modulus sixteen, not the energy's sufficient modulus thirty-two for C3. It supports the explicitly stated individual-cubic refinement and does not claim to literally enumerate the much larger E-based transitive inflation.

`ambient_audit.py` and `ambient_walk.cpp` independently count old words in six rooted ambient controls, allowing repeated geometric-profile core children but permitting each exterior physical name only once. All 96 source/converse ambient queries agree with the complete-profile transfer. The controls have nonzero c2 equal to one or thirteen, so the second-cover term is exercised rather than silently discarded. The cyclic control has c1=3, c2=1, and multiplier five.

## 3. Substitution closure and profile-dominating targets

The finalized producer Section 6 is also approved mathematically. Applying the unported fixed-point identity to a quotient of children in the same normalized class gives q_h=Lambda g_h and H=Lambda modulo the modulus. Division by the actual odd H restores the normalized class. This is closed under recursive substitution without asserting that all actual tournaments belong to it.

For target profiles larger over the integers, pair each baseline query with its core-only converse-complement query, leaving exterior attachments fixed. The entire baseline augmented profile contrast changes sign, so the ambient contrasts are opposite too. Choose the strictly negative source-minus-baseline value. Nonnegative run coefficients preserve that strict inequality after increasing any target child profile. This is a source-relative conclusion, not arbitrary pairwise injectivity.

Every tournament on N names has profile coordinate at least h! S(N,h), since every block of every set partition has at least one Hamilton path. Thus transitive source children permit arbitrary target interiors on the same physical sets. The prime-quotient original-module consequence is valid because every proper original module lies within a child. The producer correctly excludes arbitrary composite-quotient root-family changes.

`extensions_audit.py` retains a nontransitive actual source C3[TT4,TT4,TT4] whose normalized first-five profile is (1,14,4,8,0) modulo sixteen, together with a changed target child of strictly larger profile. It checks complete first-five baseline duality, both fixed exterior dual queries, and persistence of the negative separator. All fifteen physical queries and exact counts are saved. Higher closure coordinates vanish as before by factorial divisibility.

## 4. Independently derived odd-carousel cubic

For n=2d+1, order the physical carousel names as

    0,d+1,1,d+2,...,d-1,2d,d.

Switch by alternating signs. The resulting Seidel matrix is the transitive skew matrix. Its even-order inverse has entries (-1)^(i+j) above the diagonal and their negatives below it.

For an even selected subset containing the first and last displayed names, let u be the alternating switched-sign sum of the selected interior positions. The two corresponding coordinates of r are u-1 and u+1. Hence the cubic summand is u^3-u^2-u+1. Reflection of the odd-length interior reverses u, cancelling the first and third moments.

There are 2^(n-3) even interior subsets. Their second moment is (n-3)2^(n-3). To check this without an extrapolation, expand the square. For interior length ell>=3, each diagonal term occurs 2^(ell-2) times. For a pair, parity-preserving toggling cancels its sum unless the positions are adjacent or are the two endpoints of the interior. The adjacent pairs contribute +2^(ell-3) each before the factor two in the square; the endpoint pair contributes -2^(ell-3). The result is (ell-1)2^(ell-1). The n=3 empty-selection boundary is separate and gives zero second moment.

Therefore the coefficient on the original physical quotient labels is

    [a0 a_d c0] Delta_Q^g = -8(n-4).

This has valuation three for every odd n>=3. Modulus sixteen consequently suffices for this coefficient. Positive transitive child orders divisible by four already realize g modulo sixteen, by the even/odd base argument with odd fourth powers. This supports unequal such child orders and the same ambient transfer. It does not include every even child tournament.

`carousel_audit.py` checks the label order, switching and explicit inverse identities, exact subset moments through n=19, and independent full rational coefficient crosschecks at n=3,5,7,9. The structural derivation above is the proof; those checks are controls, not the reason for the all-order conclusion.

## 5. Reproduction and scope

Run each of `audit.py`, `physical_audit.py`, `ambient_audit.py`, `carousel_audit.py`, `extensions_audit.py`, and `producer_audit.py` with normal and optimized Python. Each script compiles only its own standard C++17 source using the existing compiler, and uses explicit checks active under optimization. No software is installed. Receipts pin sources and generated artifacts. Normal and optimized receipts must agree byte-for-byte.

The literal subset-DP maximum order is sixteen; its first-five ordered-cover counts are bounded by 16! binomial(15,4), safely below 2^64. The finite-word automata reduce modulo powers of two throughout. All complete source and converse banks preserve the same named child sets and fixed port roles; there is no independent relabeling between sides.

Approval does not extend to arbitrary child profiles, full unrestricted response injectivity, exact collisions, all prime orders without profile restrictions, or any complexity lower bound. The complete producer bank was independently compared: all 10,098 named positions in each orientation, including 6,133 quartics per side and every zero, agree with this auditor's separately recomputed exact subset inverses. Both clean producer replays reproduce the frozen receipt and all data byte-for-byte. The producer receipt SHA-256 is `2c66bc0a74e8245ee6a5ef1fe62f93a6912e9214c511cd84af3923d3ecf06b4d`. All six independent audit programs also have byte-identical normal/optimized receipts. Source and artifact hashes were rechecked after an accidental overwrite of one audit source was restored byte-for-byte; no data were changed or lost. The unlisted empty normal.json from the other worker is not part of this packet. No proposed thickening appendix or other future addition is approved by this review.
