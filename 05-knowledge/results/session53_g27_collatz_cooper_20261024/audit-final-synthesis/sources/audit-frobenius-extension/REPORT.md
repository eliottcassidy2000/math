# Independent audit: two-step Frobenius and modulo-nine digits

Verdict: PASS for the new deductions, explicitly conditional on the inherited modular inputs listed below. This is not a fresh independent certification of inherited Lemma T, its integral model, or the global modularity/pole calculations. Finite replay is supplemental evidence, not the proof of an infinite quantifier. No originality claim is made.

## Audited mathematical statements

1. For every prime p and integer n>=0, c(p²n)=c(n) modulo p² follows from inherited Theorems B/E/D and the new constant-term argument. Here G=sum c(n)t^n is the normalized formal square root of the s7 generating series, G(0)=s(0)=1, and s²=(1+t)(1-27t).
2. Modulo9, G=P9 G(t^9), P9=1+2t+4t²+3t³+3t^4+6t^5, follows from inherited B at3 and the exceptional excellent-lift identity. It implies the stated base9 product and exact one-defect support criterion.
3. C3G=(1+3t)G(t³), H=(1+6t)G(t³), and C27G=C3G modulo9 are correct. C3G=H is false modulo9.

## Independent proof audit of the inert extension

The coefficient change of variables is exact: dt/t=theta²s dq/q and G(t(q))=theta imply c(n)=CT_q(theta³s t^(-n)). For f*=theta³a with a=1-3t/4, f* dq/q=G a/s dt/t=H a dt/t. Thus the correction term is a coefficient of H a epsilon, not of G a epsilon. This distinction is essential.

The exact divisor-sum identities in the producer proof are correct, including the V_p correction when p divides the coefficient index. Modulo p², U_p acts on A1 with eigenvalue1 and A2 with eigenvalue chi(p), not the reverse. Constants are preserved. Therefore inert U_pf=f* and U_p²f=f modulo p². The companion f*=1-7A1/8+49A2/8=theta³(1-3t/4) has matching coefficients (1,21/4,105/4) through q². Given the inherited standard modularity of the weight3 Eisenstein forms and theta³t, both sides have the same weight, level7, and character. The Sturm bound2 proves the identity (multiply by8 first). Every denominator is a unit for inert primes, which are odd.

Writing T2=t(q^(p²)), the identity T2=t^(p²)(1+p epsilon(t(q^p))) modulo p² has the correct sign and precision: the other factor is (1+p epsilon(t))^p=1 modulo p². Inverting and raising to the integer exponent n gives t^(-p²n)=T2^(-n)(1+pn epsilon(t(q^p))) modulo p². No division by n is used, including when p divides n, and n=0 gives zero correction.

The adjoint relation CT(A(q)B(q^p))=CT((U_p A)(q)B(q)) holds on formal Laurent series with finite negative tails. U selects coefficients at all integer multiples of p. Products have well-defined constant terms since the negative tails are finite. Applying it twice gives the first term c(n); applying once to the correction and reducing modulo p gives [t^(pn)]H a epsilon.

Inherited B gives H=H_p G(t^p) modulo p. The proof of inherited E supplies epsilon H_p=P_p with P_p(0)=0 and degree<=d_p-1 for inert p>=5. Multiplication by a raises this bound only to d_p<p. A polynomial of that degree and zero constant has no exponents divisible by p, even after multiplication by G(t^p). This proves vanishing at every pn, not merely a finite set. At3, epsilon=0 modulo3 proves vanishing separately, without an invalid weight2 Eisenstein argument. Split primes (including2) and7 use iteration of inherited D, not the inert proof.

The result is same-precision two-step return. It is not the open split-prime modulus p^4 refinement and is not a scalar sign rule or naive G/H swap.

## Independent proof audit of modulo-nine digits

Modulo9 the normalized s satisfies s²=1+t, and H_3=1+6t. Inherited B and t^sigma=t³ give G=s(1+6t)G(t³). Twice substituting this exact modulo9 formula gives s(t)s(t³)(1+6t+6t³)G(t^9). The candidate normalized square root 1+5t+t² squares to (1+t)(1+t³) modulo9. Uniqueness follows because2 is invertible; this supplies an elementary alternative to the producer's valid binomial perturbation. Multiplying gives exactly P9 above.

Since deg P9<9, coefficient extraction has a unique base9 address with no carries. Digits0,1,2 contribute units; digits3,4,5 contribute valuation1; digits6,7,8 contribute zero. Consequently modulo9 nonzero means no forbidden digit and at most one valuation1 digit. For 0<=n<9^r, this also gives exactly 3^r units and r*3^r nonzero multiples of3, hence (r+1)3^r nonzero coefficients. This count is an elementary consequence, not an assertion of exact integer 3-adic valuations when the residue is zero. Leading zeros do not matter.

The hostile c9=26337699740, c1=2 gives difference18 modulo27. The first companion mismatch is c3-d1=336-15=6 modulo9.

## Independent finite path and guards

independent.py constructs s7 directly by its binomial sum through240, then constructs its normalized formal square root by convolution, independently of both producers' three-term recurrence. It tests every p²n<=240, the modulo9 digit law through240, digit-budget counts through five base9 places, and the exact finite polynomial square and product. Controls reject modulus27 strengthening, naive companion swap, deletion of negative Laurent indices, a corrupted digit weight, and coefficientwise numerical square roots. Normal and -O outputs agree. These tests do not prove infinite statements.

Startup read AGENTS, navigation, guardrails, research protocol, integer-sequence frontier, and core-paper status. The shared math tree was already dirty and remains untouched. Agent message processor was unavailable because no .machine-id is registered; no identity was invented. Anchor: conditional inert proof. Niche: base9 support budget. Wildcard: Laurent operator adversary. Closest inherited mechanism B/E; hostile naive swap; corrected near miss naive Frobenius replacement; underused sidecar weight3 companion.

No new external-literature priority claim was attempted. The producer's bounded search and inherited source remain provenance, not a proof of novelty. Outstanding frontier: independently audit inherited Lemma T/integral modular geometry and stronger p-adic precision, rather than representing this consequence audit as that missing work.

## Appendix: exact residue distribution (independently audited supplement)

Let r>=1, j in {0,...,5}, and h=(2,1,-1,-2,-1,1), with cyclic indices. Among 0<=n<9^r the number with c(n)=2^j modulo9 is

U_j(r)=(3^r+(-1)^j+2^r h[(r-j) mod6])/6.

The number with residue3 is (r*3^r+r)/2; with residue6 it is (r*3^r-r)/2; with residue0 it is 9^r-(r+1)3^r. At r=0 there is just n=0, residue1. These are exact consequences of the inherited-dependent digit theorem, not extrapolated counts.

Proof: unit digits0,1,2 give weights2^0,2^1,2^2. Their residue counts are the coefficients of (1+x+x²)^r modulo x^6-1. Sixth-root Fourier evaluation has eigenvalues 3, 2exp(i*pi/3), 0, 1, 0, 2exp(-i*pi/3). Inversion gives U_j(r). This uses characteristic-zero complex numbers only to evaluate exact integer coefficients; no reduction of roots of unity modulo9 is invoked. For the single defect digit, normalized defect weights3,3,6 divided by3 have signs +,+,- modulo3, sum1. The unit-digit signs are +,-,+, also sum1. Choosing the defect position gives signed difference r between residues3 and6, while their total is r*3^r. This proves all formulas.

For explicit rational generating functions in the place-count variable z, set h_k=h[k mod6]. Then

sum_{r>=0} U_j(r)z^r = delta_{j,0} + (1/6)(3z/(1-3z) + (-1)^j z/(1-z) + (2h_(1-j)z-4h_(-j)z²)/(1-2z+4z²)).

The residue3/6 generating functions are (1/2)(3z/(1-3z)² +/- z/(1-z)²). The residue0 generating function is 1/(1-9z)-1/(1-3z)². The r=0 correction in the unit formula is essential because the zero Fourier eigenvalues contribute only there.

Consequently integers n with 9 dividing c(n) have natural density1. Indeed, if 9^r<=X<9^(r+1), the number of exceptions below X is at most (r+2)3^(r+1), whose ratio to X is at most 3(r+2)/3^r, tending to0. This concerns coefficient divisibility, not Collatz or orbit coverage.

## Final packaging acceptance

Both producer packets were hardened and refrozen before final acceptance. Their authoritative manifest hashes and independently rerun /tmp replay receipts are in producer-replay-receipts.json. The two-step packet now carries its inherited sources and recursively closed manifest; the modulo-nine packet now rejects nested undeclared files. Independent optimized-mode nested-file attacks are rejected for both. The two-step producer hostile suite also rejects corruption, missing files, and nested extras normally and under -O. No mathematical blocker remains for the explicitly conditional statements. This audit packet itself has a recursive closed manifest, copied proof snapshots, and independent normal/-O replay.
