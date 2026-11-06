# Independent audit: digital gcd inverse and targeted certificates

**PASS for the stated numerical domain and original indices v>=2.** The all-order inverse, search-node bound, and fixed-finite-template membership complexity are sound. Exact targeted certificates pass. A small public-helper domain issue is recorded below: the frozen helper assumes, but does not enforce, its original-index lower bound. This does not affect the theorem, default calls, or complete template caller.

Universal twelve-template occupancy and the original spectrum remain **OPEN**. A membership algorithm is not an occupancy theorem.

## 1. Frozen inputs and independent implementation

`INPUTS.json` records all twenty frozen source artifacts. Principal hashes:

- `PROOF.md`: `23a4d730689a643daf46b75356b4e8c7dc39af904015f44502db9d57b017221b`.
- `verify.py`: `dff630f52771985da4bf5eb6d9a3b347c5f68af671ea0ab2632851c493012740`.
- `check_targets.py`: `1761949b4f9ffc1a58462fc19077f7b6294f8f96c971eabedb06aa223b8abcc9`.

The core proof, matrix header, core code, menu, and factorization proof match their previously independently audited dependency hashes. No source file was edited.

`independent.py` imports no author code. It computes g through literal binary source operations, reconstructs the full central-word matrices, uses a uniform-depth inverse instead of the author's adaptive stopping, recomputes local modular inverses, independently builds target adjacencies, and checks the recursive Lucas chains. `compare_author.py` separately treats the author helper as the implementation under test and compares it with direct physical-digit enumeration.

Author, target-certificate, independent, and API-comparison normal/optimized runs all pass with byte-identical paired receipts. Author receipts also exactly match the frozen originals. The exhaustive author run at 2134117 was separately replayed and exactly matches its frozen receipt.

## 2. Numerical versus physical domains

The coefficient assumptions are sufficient numerical conditions, not a characterization of all physical sources. The proof keeps that distinction explicit. For an actual source, the row is exactly

    (alpha,beta,gamma,delta)=(d,U-d,E-d,D-E-U+d).

The reviewed nonnegative coordinates give the claimed inequalities. Odd D,d and U=E modulo two are preserved by each source operation, so alpha is odd. The adjoint row for every central word comes from the actual source obtained by applying its dagger to tau(u), not an independently selected cone vector.

The independent checker verifies this row identity, cone membership, and beta>=u-1 on all 35 full-menu matrices for u=2 through 128. The all-order reason for the latter lower bound is correct: the empty word gives beta=u-1; otherwise the final J in the dagger makes beta equal to its predecessor's U, and physical operations never decrease U.

## 3. Peeling, fixed gcd, and original indices

Substituting h=2j and h=2j-1 gives exactly the two stated coefficient updates. In the odd case the new lower delta bound is -2 beta, which equals alpha-(alpha+2 beta), so the entire numerical domain is preserved. Alpha stays fixed. Because alpha is odd, multiplication of beta by two is invertible modulo alpha; both branches preserve gcd(alpha,beta). Without odd alpha this reasoning would not hold.

The index relation v=s h-r and both offset updates have the correct orientation. At depth k the scale is 2^k and 0<=r<s. These maps partition original indices by binary residues. The independent checker tests the exact original-row identity at internal h=1,2,3,7 for all visited nodes in its fixed-row controls, not just at returned witnesses.

An internal scalar h=1 is legitimate. The odd recurrence at that point is the equality 1=1; no physical singleton-Z construction is inferred. The independent controls exercise 218 scalar singleton candidates. A particularly small positive boundary is the row (1,1,0,0) at target three: internal h=1 reconstructs the legal original index v=2.

### Public-helper precondition

The original theorem requires v>=2, and all default/complete-template calls meet that condition. However, the frozen public helper permits the out-of-domain call

    fixed_inverse((1,1,0,0), 2, vmin=1)

and returns [1]. This is scalar arithmetic, not a physical-source certificate. A source-level API precondition vmin>=2 should be made explicit or enforced in the release. This is not a counterexample to the theorem or to the published target searches. The author was notified; any later guard must be separately archived and reviewed rather than silently changing frozen inputs.

## 4. Interval, progression, and leaf bounds

The lower inequality follows from g>=1, gamma>=-alpha, and delta>=alpha-beta; the upper inequality follows from g<=h and nonpositive gamma,delta. Their rearrangement gives the interval endpoints with the stated ceiling and floor. Intersecting with h>=1 and the exact transformed requested-index window preserves completeness. A target below alpha is impossible.

The width is strictly bounded: before integer rounding, its upper estimate is

    1 + n alpha/[beta(alpha+beta)] - alpha/beta,

which is strictly less than 1+n alpha/beta^2. The strict sign matters at the smallest modulus.

For h=2t+p, reducing modulo alpha gives exactly the displayed progression. The divisibility check is necessary. After division by d0, the coefficient 2 beta/d0 is invertible modulo m=alpha/d0 because m is odd and the reduced gcd is one. At a leaf beta^2>=n d0, the width is less than 1+m<=2m. Thus each parity progression has at most one candidate, including m=1. Every candidate still requires an exact g evaluation and equality check.

The author code uses the correct parity-specific RHS, transformed bounds, modular inverse, and candidate retest. Independent fixed-depth tests also compare every locally recomputed inverse with the prepared depth formula, verifying that the factor of two in the progression coefficient is retained.

## 5. Node bound and bit complexity

Every branch has beta_k>=2^k beta_0, hence stops by the stated least K. The complete binary tree therefore has at most 2^(K+1)-1 nodes and at most 2^(K+1) candidate evaluations. This proves the one-row bound O(1+sqrt(n d0)/beta_0); pruning is not needed for that asymptotic claim.

For a fixed finite word menu, restoring the full transpose-closed family justifies u<=v. The reduced twelve representatives alone would not justify that ordered-index restriction. The independent reconstruction finds exactly 35 matrices at length five and verifies transpose closure and complete positive-ray product bounds.

Using d0<=alpha<=beta and beta>=u-1 gives the per-row bound O(1+sqrt(n/(u-1))). Summing for u<=sqrt(2n) gives O(n^(3/4)) arithmetic nodes. With a fixed word-length cutoff, starting coefficients are O(u), and the numerical cone bounds all later coefficients by the growing beta. Depth, index sizes, and coefficient bit lengths are O(log n). Exact arithmetic, modular inverse preparation, and g evaluation therefore contribute only polylogarithmic overhead. This supports the claimed O(n^(3/4) polylog n) bit complexity for complete membership in a fixed finite family.

The constants are allowed to depend on that fixed family. The statement does not apply with a word cutoff growing unrestrictedly with n, claim average-gcd savings, or imply that every odd target has a witness.

## 6. Complete finite controls

The independent uniform-depth implementation checks 225 valid numerical rows, including interior gamma/delta values, with 94,725 complete target/window comparisons. Its direct expected answers are generated from actual binary-source g values. It checks strict intervals, invariant gcd, depth inverses, reconstructed indices, duplicate exclusion, and complete tree/candidate budgets.

A separate test of the frozen author API supplies 58,725 complete/window comparisons and 13,500 first-result comparisons on those 225 rows. First-result checks establish that a returned value is exact and that an empty result agrees with the complete finite control; they do not assert smallest-index ordering.

For finite-template membership, independent direct product-bounded enumeration and the independent inverse agree for every target 1 through 501. The independently reconstructed full length-five inverse at 2134117 exhausts all 25,400 rows and finds no witness. The uniform-depth implementation visits 1,103,398 nodes; this larger number is expected because it does not use the author's adaptive early leaf rule and parity pruning.

The author's separately replayed exhaustive run visits 97,131 nodes, evaluates 38,217 scalar candidates, and also exhausts 25,400 rows with no witness. These results certify absence in the selected finite menu, consistent with the older complete finite prefix. They do not certify absence among arbitrary tournaments or longer menus.

## 7. Exact prime chains and positive physical certificates

All 26 recursive Lucas certificates pass independently, with 117 exact modular-power/gcd checks. Every p-1 factorization is complete, all factors have already been certified, and two is the explicit base case. The full-order argument is sound: for every prime divisor r of p, the supplied conditions force p-1 to divide the multiplicative order of the base modulo r, hence p-1 divides r-1 and r>=p. Thus p itself is prime. No probabilistic primality result or SymPy assumption is used.

The nine prime targets/companions exclude the two factorized templates at every source height by the separately audited divisor theorem. The following remaining-menu witnesses are exact:

- 100000837 = Psi_JOZJJ(3,1890994), with 29 vertices.
- 10000001533 = Psi_empty(2,7840792990), with 35 vertices.
- 1000000001911 = Psi_JOZJ(8,22798770134), with 43 vertices.

The independent checker verifies both g coordinates, matrices, full source profiles, canonical chronological graph words, complete adjacency orientations, and actual graph orders. It verifies every recorded normalization trace and additionally proves the nontrivial boundary absorption used in the last witness by an exact row identity, valid for every right index.

The Hamilton-path counts follow from the independently reviewed physical compiler and the exact profiles. This audit does not claim literal enumeration of the large graphs' paths. Their one-cycle property likewise uses that compiler. The long first-witness search statistics remain author-only: the positive certificates were checked without replaying those searches, and no minimality or complete representation enumeration is inferred.

## Reproduction and limits

From this audit directory:

    python3 independent.py > check-independent.json
    python3 -O independent.py > check-independent-optimized.json
    cmp check-independent.json INDEPENDENT_NORMAL.json
    cmp check-independent.json check-independent-optimized.json
    python3 compare_author.py > check-comparison.json
    cmp check-comparison.json COMPARISON_NORMAL.json
    sha256sum -c MANIFEST.sha256

The independent checker uses copied, frozen input data. The API comparison locates verify.py in a containing release packet or in the original neighboring source packet. The original unguarded API boundary is deliberately recorded in its receipt, so a later reviewed guard should change only that boundary status when replayed.

This audit proves no occupancy theorem, extends no interval census, asserts no optimal graph order, and does not turn a two-template prime obstruction into a full-menu counterexample.
