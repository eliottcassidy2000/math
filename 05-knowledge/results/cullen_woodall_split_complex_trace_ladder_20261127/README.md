# The split-complex Cullen–Woodall trace ladder

This extends the [two prime-certificate machines](../cullen_woodall_prime_certificate_and_sieve_barrier_20261127/README.md) and the [shared Collatz countdown](../collatz_cullen_woodall_post_exit_inverse_20261127/README.md). The point is an exact common recurrence, not a claim of infinitely many prime terms.

## 1. The split-complex unit exists in the twin's residue ring

Write X=n·2^n, C=X+1, W=X−1, and M=CW=X²−1, with n≥2. In R=ℤ/Mℤ the element L=X satisfies L²=1. The factors 1+L and 1−L are both nonzero, but their product vanishes. Since 2 is invertible, set e_C=(1−L)/2 and e_W=(1+L)/2. They are orthogonal idempotents summing to 1. Under R≅ℤ/Cℤ×ℤ/Wℤ they are (1,0) and (0,1), respectively. Thus one can run a common integer recurrence over R and read different rungs on the two channels. This splitting is distinct from the prime-by-prime splitting below.

## 2. One trace, two consecutive stopping rungs

Write X=k·2^e with k odd. For any integer t let α satisfy α²−tα−1=0 in its rank-two companion algebra; conjugation sends it to t−α=−α⁻¹. Let B_j(t)=Tr(α^(k·2^j)). Since k is odd, B_1=B_0²+2. Thereafter B_(j+1)=B_j²−2. This is the usual Lucas trace mechanism (compare the [Lucas–Lehmer–Riesel framework](https://www.mi.uni-koeln.de/~mmertens/thabit2.pdf)), applied here with the same norm-minus-one seed on both candidates. The offset by a rung is essential.

**Theorem.** For Woodall n≥2, W is prime iff some t gives B_(e−1)=0 modulo W. For Cullen n≥3, C is prime iff some t gives B_(e−2)=0 modulo C. More precisely, for prime C exactly (C−1)/4 seed residues work, and for prime W exactly (W+1)/2 seed residues work. At the excluded n=2, C=9 is a genuine false positive: t=4 gives B_1=18. The trivial Cullen prime n=1 is outside the two-channel formulation.

**Soundness.** If q is an odd prime factor of the candidate and a positive-index B_j vanishes modulo q, let β=α^(k·2^j). Its norm is one, so its characteristic identity gives β²=−1. In particular the 2-part of the order of α is exactly 2^(j+2). Put D=t²+4. If D=0 modulo q, reducing α to its repeated root gives α²=−1 and cannot supply the required order ≥8. If D is a nonzero square, the companion algebra splits as two copies of 𝔽_q and the order divides q−1. Otherwise conjugation is Frobenius in 𝔽_(q²); α^(q+1)=−1 and the order divides 2(q+1). Consequently:
- For W, j=e−1 and q≥2^e−1>√W. For the strict inequality, k is odd with k<2^e−1, so (2^e−1)²−W=2^e(2^e−2−k)+2>0.
- For C, j=e−2 and q≥2^(e−1)−1>√C when n≥4. Indeed the squared difference is 2^e(2^(e−2)−1−k), positive directly for n=4 and from 2^(n−2)−1>n for n≥5. At n=3, C=25, q=5: neither q−1 nor 2(q+1) is divisible by 8. Thus every passing candidate in the stated ranges is prime.

**Completeness and counts.** For prime C, C≡1 mod 8. Inert or ramified seeds cannot provide the required order because v₂(2(C+1))=2 and the ramified root has order at most four. In the split case the roots are a and −1/a; the positive-index trace vanishes at j=e−2 exactly when a is a nonsquare. The other root is then also a nonsquare. No nonsquare is fixed by a↦−1/a: a²=−1 forces a to be a square when C≡1 mod 8. Pairing the (C−1)/2 nonsquares gives (C−1)/4 distinct traces t=a−1/a. This also converts any usual Cullen Proth witness a to a seed t with no discrete logarithm.

For prime W, W≡7 mod 8. No split seed can work since v₂(W−1)=1, and there is no ramified seed because −1 is not a square. Every nonsplit seed works: α^(W+1)=−1, which is precisely β²=−1 at j=e−1. The norm-minus-one fiber has W+1 elements, all nonscalar (a scalar would square to −1). Conjugation pairs them, with different pairs giving different traces. This proves the seed count. In particular the golden seed t=1, with α=(1+√5)/2, certifies exactly those Woodall primes congruent to 2 or 3 modulo 5, by quadratic reciprocity. For example it certifies n=2,3,6,123,249,362,384, but not n=30,75,81,115. Failing a particular seed is not evidence of compositeness.

## 3. A simultaneous certificate and its limitation

Let r=B_(e−2) evaluated modulo M. The split-complex element e_C·r+e_W·(r²−2) vanishes exactly when both tagged certificates pass. Equivalently, without ever dividing by 2, test
\[
(X-1)r+(X+1)(r^2-2)\equiv0\pmod{X^2-1}.
\]
For n≥3 the existence of a seed passing this test is equivalent to both terms being prime. The forward direction is soundness and the reverse direction is CRT; on a hypothetical prime pair the fraction of successful t modulo M is exactly (C−1)(W+1)/(8CW). This is a compression of two local certificates, not a way to make the index prime-producing. The previous [factor-planting theorem](../cullen_woodall_prime_certificate_and_sieve_barrier_20261127/README.md) still forces arbitrarily long composite blocks along every periodic phase. The common 2-adic clock also appears in the signed Collatz runs, but in the trace ladder the two stopping rungs measure multiplicative order, not integer descent.

Run **node replay.js** to check the two-channel identity and small candidates, reproduce the exceptional false positive at 9, and check two Cullen and eleven Woodall prime certificates.
