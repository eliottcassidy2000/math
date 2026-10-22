# Independent audit: reciprocal powers and bounded bridge claims

Date: 2026-10-22 UTC. Scope and input hashes: INPUT_SHA256.txt. No repository or Git writes were made. Reviewed source files were not changed.

## Verdict

PASS: RECIPROCAL_POWER_THEOREM.md is mathematically sound as stated. No mandatory theorem correction found. The author's exact verification script reran successfully with all 229 assertions; replay.txt is byte-identical to the supplied verification.txt.

PASS with separate scope: the recurrence, K3 attribution, null quotient, and E8 facts checked in PRIOR_1729_AND_E8_HOLES.md. The inherited tournament computations, remote provenance/history, and complete inherited 27-dimensional script were not independently replayed in this audit. This report does not certify the entire modular directory.

Recommended strengthening: add the independent E8 stabilizer justification below. The conclusion is correct, but an intermediate sentence in the cited paper confuses the full symmetry group of the cross-polytope with its index-two subgroup. No claim of novelty or all-degree literature completeness is certified.

## 1. Universal theorem review

1. The companion matrix correctly advances (u_n,u_(n+1),u_(n+2)); its characteristic polynomial is (x+1)(x²−Tx+1), and its determinant is −1.
2. Over Q(T), the three roots are distinct, including separation from −1. Its tensor powers are diagonalizable. Tensor eigenvalues are (−1)^a λ^(b−c), with a+b+c=d. Since a≡d−(b−c) modulo 2, they are exactly μ_r=(−1)^(d−r)λ^r for −d≤r≤d. Every such r occurs by taking a=d−|r|. No tensor dimension is being mistaken for a count of distinct eigenvalues.
3. Pairing r and −r gives product 1 and trace (−1)^(d−r)L_r(T). The displayed polynomial therefore kills the full 3^d-dimensional tensor matrix, not merely its symmetric quotient or pure powers.
4. Every entry of F(M) lies in Z[T]. The embedding into the splitting field is injective. Proving zero there proves a polynomial identity before specialization; arbitrary ring maps Z[T]→R preserve it. Thus repeated roots, positive characteristic, zero divisors, nilpotent parameters, and the zero ring cause no denominator or diagonalization problem. This is the crucial universal argument.
5. Applying coordinate functionals to tensor products of independently chosen state vectors proves mixed monomials, including repeated factors. Linear combinations prove every stipulated homogeneous polynomial. No division by d! or polarization identity is used; that distinction matters over small characteristic.
6. The d=0 case is F=x−1 and a constant observable. Fixed shifts remain legitimate solutions. Nonhomogeneous observables are correctly excluded from the stated single-polynomial result.
7. F is monic with constant coefficient (−1)^(d+1). This proves two-sided propagation and unique backwards extension. A block of 2d+1 zeros is sufficient, not necessarily minimal for the chosen residual.
8. At t=−2 all modes coalesce and F=(x−(−1)^d)^(2d+1). At t=2, counting positive/negative modes gives (x−1)^(d+1)(x+1)^d. The manuscript retains these necessary multiplicities. Dropping multiplicities fails in our direct tensor controls.

## 2. Sharpness and observable-space interpretation

For real |t|>2 choose |λ|>1. The |λ|^r are pairwise distinct, so the 2d+1 modes are distinct and independent by Vandermonde. All occur among monomials of the three eigenmode solutions, while all degree-d monomials lie in their span. Therefore the common annihilator has exactly degree 2d+1. This proves a common-space statement, not a claim about an individual sequence, residual, or chosen input family.

The rational claim follows by using rational initial-state basis solutions, then extending scalars to the splitting field. A polynomial annihilating all rational-basis degree-d monomials would annihilate their extended scalar span. The relation XY=Z² reduces the homogeneous coordinate dimension to 2d+1; Vandermonde prevents additional homogeneous relations in this range. No modular, K3, or tournament map follows from this fact.

## 3. Exact application and independent controls

The independent script uses Chebyshev polynomials 2T_r(t/2), rather than importing the reviewed Lucas implementation. It tests full symbolic tensor matrices for d=0,1,2,3; compares the squarefree characteristic polynomial at t=83 for d=1,2,3; and tests dual-number rings (Z/mZ)[epsilon]/(epsilon²) for m=4,6, t=epsilon,2+epsilon,−2+epsilon, d=2,3. The latter explicitly exercise nilpotents combined with zero divisors. These controls supplement, not replace, the universal proof.

The t=83,d=3 descending coefficients independently agree:
[1, −564733, −3889880904, 322812677460, 322812677460, −3889880904, −564733, 1].

Direct rational-function series expansion reproduces all seven displayed triples and zero residuals. The backwards triple is (−9,12,10), its cubic difference is −1, and the rearrangement is 9³+10³=12³+1=1729. The script includes hostile checks that dropping repeated-root multiplicities fails, and that a nonhomogeneous observable need not satisfy the degree-one candidate. See independent_controls.txt for exact results.

Harness development note: its initial squarefreeness check compared the integer polynomial gcd with 1, overlooking the harmless constant content 63. The corrected check asks whether that gcd has degree zero. This was an audit-harness error, not an input theorem failure.

## 4. Primary-source attribution review

[Hirschhorn 1996](https://web.maths.unsw.edu.au/~mikeh/webpapers/paper40.pdf), author PDF pages 1–3, explicitly reduces the same three cubes to seven spectral modes, establishes seven initial checks as sufficient, and computes those seven checks. The manuscript correctly calls this recovered prior mathematics. The bibliographic reference is independently repeated in the 2006 paper.

[Han–Hirschhorn 2006](https://web.maths.unsw.edu.au/~mikeh/webpapers/paper128.pdf), author PDF pages 2–3, provides quadratic forms in consecutive solutions of h_(n+2)=9h_(n+1)+h_n, the stated 3×3 evolution, and explicit negative-index validity. Thus neither the negative extension nor the 83=9²+2 mechanism is newly established priority here. The exact recurrence theorem is a synthesis; the identical all-degree formula's first appearance was not investigated exhaustively.

[Ono–Trebat-Leder v5](https://arxiv.org/html/1510.00735v5), Theorems 1.1 and 1.3 and Section 3.2, support k(T), rational function-field rank 2, rank 4 over Q(sqrt(−3))(T), geometric Picard number 18, six type-IV fibers, and rho=r+14. The curve has j=0. Our script separately verifies both displayed sections and squarefreeness of the degree-six k. This is an attribution audit plus polynomial checks, not an independent Magma rank computation or full reproof of K3 classification.

[2017 erratum](https://link.springer.com/article/10.1007/s40993-017-0076-8) adds omitted historical acknowledgement of Ramanujan and Hirschhorn. It does not retract the rank/Picard theorems. The addendum correctly separates this elliptic K3 construction from the cubic equation in four homogeneous coordinates and from unsupported tournament-to-period rhetoric.

## 5. Null quotient review, separate from recurrence

The II_(9,1) argument is correct. Unimodularity identifies L with its integral dual, and primitivity of v makes the pairing functional primitive, hence surjective onto Z. Choose (v,w)=1. Evenness makes w'=w−(w,w)v/2 integral; its norm is zero and (v,w')=1. The span U of v,w' is unimodular, so it splits integrally. The complement has positive signature (8,0), is even unimodular, and is E8 by the explicitly imported classification. Quotienting v-perp by Zv gives that complement. Mapping hyperbolic bases and complements proves the full-orthogonal-group orbit assertion. Orientation/cone conventions must remain specified if the acting group changes.

The same argument in II_(25,1) leads to positive rank-24 complements; their 24 classes produce 24 cusp types. That does not say there are 24 Leech deep-hole types. [Borcherds 1985](https://math.berkeley.edu/~reb/papers/leech/leech.pdf), introduction and preliminaries, supports the indefinite classification, the Leech covering radius sqrt(2), and 23 deep-hole types corresponding to the non-Leech Niemeier lattices. The rank-eight classification is an imported classical theorem, not proved by the introductory historical sentence alone.

## 6. E8 nearest shell and stabilizer

The E8 normalization and elementary nearest-point argument in the addendum are correct: at e1 the nearest set is e1±e_i, 16 points at distance 1. Integer candidates are checked by parity; half-integer candidates are at squared distance at least 2. The cited [Cambridge paper](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/polarization-of-lattices-stable-cold-spots-and-spherical-designs/16B869EA8216480737BB2BD00B115F9D), Lemma 3.6(3) and Appendix A.3.1, supports the deep-hole orbit, radius, shallow/deep distinction, and affine stabilizer W(D8).

Source caveat and independent repair: the proof near Lemma 3.6 says the full symmetry group of a standard cross-polytope is W(D8). In fact it is W(B8), all signed permutations. To recover the correct E8 conclusion, write an affine stabilizer at c=e1 as phi(x)=c+g(x−c). Since g permutes the centered nearest set {±e_i}, g is a signed permutation. The translation c−g(c) belongs to D8, regardless of its signs. Thus phi preserves E8 if and only if g does. A signed permutation with k sign flips sends h=(1/2,...,1/2) to h plus an integer vector with coordinate sum −k. Consequently it preserves h+D8 exactly for even k. Every such g works, and no odd one does. This yields precisely (Z/2Z)^7⋊S8=W(D8), independently of the erroneous intermediate source sentence.

No hole census was recomputed, and nearest-shell geometry is not identified with the Lorentzian null quotient. The addendum correctly keeps these constructions distinct.
