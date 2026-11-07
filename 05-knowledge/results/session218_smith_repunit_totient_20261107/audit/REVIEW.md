# Independent mathematical review

Status: PASS. Bridges A and B independently proved below; full Smith and dyadic producer proofs reviewed.

## A. Simultaneous Mersenne/Wagstaff primes force a nontotient exponent

Let q>5 be an integer and suppose M=2^q−1 and W=(2^q+1)/3 are primes. Then q is prime: if q=ab with a,b>1, 2^a−1 is a proper divisor of M. In particular q is odd.

For any prime q>3, the complete inverse-totient fiber of 2q is {2q+1,4q+2} if 2q+1 is prime, and empty otherwise. Indeed each prime r|x must satisfy r−1|φ(x)=2q, so r is among 2,3,q+1,2q+1. Here q+1 is even and greater than2. If 2q+1 is absent, φ(x) has only prime factors2 and3, impossible. If present, it contributes the entire factor2q; no factor3, no additional power of2 beyond the first, and no higher power of2q+1 is allowed. Thus x is exactly2q+1 or2(2q+1).

Suppose ℓ=2q+1 were prime. Fermat gives ℓ | 2^(2q)−1=3MW. Since ℓ>3, primality gives ℓ=M or ℓ=W. But 2^q>6q+2 for every integer q>5, by equality at5 and strict induction thereafter. Hence W>ℓ and M>ℓ, contradiction. This proof is simpler than Euler's criterion.

The original Euler proof also works: q≡3 mod4 gives ℓ≡7 mod8 and ℓ|M; q≡1 mod4 gives ℓ≡3 mod8 and ℓ|W. The exceptional equalities occur at q=3 and q=5 respectively.

Boundaries: q=3 has φ⁻¹(6)={7,9,14,18}; q=5 has φ⁻¹(10)={11,22}. Both M,W are prime there, so q>5 must not be dropped. A nontotient has zero preimages; this is not evidence against the Carmichael totient conjecture, which excludes a fiber of size one.

## B. Base-dependent Smith bridges

Use the composite-only definition and prime factors with multiplicity. Write s_b(n) for the base-b digit sum.

For prime p≥5, set b=p−1≥4. Then p=b+1 has digit sum2, p²=b²+2b+1 has digit sum4, and2p=2b+2 has digit sum4. The prime-factor digit sums are respectively2s_b(p)=4 and s_b(2)+s_b(p)=4. Both numbers are composite. The restriction p≥5 matters: p=3 gives neither9 nor6 Smith in base2.

If R_b(k) is prime, k≥2, and2≤c<b, then s_b(cR_b(k))=ck, since this product has k digits all equal to c. Its prime-factor digit sum is k+sopfr(c), where all prime factors of c are single digits. Since sopfr(c)≤c, the Smith equation (c−1)k=sopfr(c) forces2c−2≤c, so c=2 and then k=2. Conversely these parameters give equality. This theorem says nothing about multi-digit c; omitting c<b is invalid.

For b≥8, s_b(b²−1)=2b−2. Prime-factor sums are additive across (b−1)(b+1), even when the factors share primes. The contribution of b−1 is at most b−1. If b+1 is prime, its digit sum is2. Otherwise every prime divisor is less than b, and sopfr(b+1)≤(b+1)/2+2<b−1. The bound sopfr(n)≤n/2+2 for composite n follows by choosing a factorization n=uv, 2≤u≤v, and bounding sopfr(n)≤u+v≤n/2+2. Consequently equality is impossible. At b=2 the number3 is prime and excluded; for b=3,4,5,6,7 the digit sums are4,6,8,10,12 and the factor sums6,5,9,7,11, respectively. Thus no base b≥2 makes b²−1 Smith. A blanket assertion that its factor sum is always smaller would be false at bases3 and5.

## Finite independent verification

Run `python3 session218-audit/verify_independent.py` from the common workspace, or run the script directly from any location. It imports only the Python standard library and no producer code or data. Exact ranges and counts are recorded in `verification.json`. Finite checks support, but do not replace, the proofs above; the inverse-totient scan alone does not exclude preimages beyond its bound. No infinitude claim, novelty claim, or claim resolving Carmichael's conjecture is made.

## Full producer artifact review

Reviewed `session218-smith/SMITH_REPUNIT_BRIDGES.md` in full. All six sections pass, including the finite all-base sieve with unique unbounded exception4; carry-free repunit powers; the general separated-block seed; decimal multiplier1540; binary complement obstruction; carry mass and grouped-base transport; and prime-length/order restrictions. Definitions must retain integer bases b≥2, lengths k≥2, positive x,y, and grouping exponent h≥1.

For the1540 seed, factor multiplicity gives A_10(1540)=18, t=171, and R_10(4)−t=940, so9+13−4=18 exactly. The identity is unconditional for k≥4; the Smith conclusion requires a prime repunit or independently equal repunit factor weight. The composite k=6 hostile is correct. The binary complement proof is valid including the endpoint c=2^k. The small c bound is essential:238 is a valid binary Smith counterexample to an unrestricted obstruction.

Reviewed `session218-scaling/DYADIC.md`. The dyadic existence criterion and complete source-shape parameterization pass: exactly one factor supplies q, either through p−1 for p=2^j q+1 or via q² when q−1 is a power of2. All remaining odd factors are distinct Fermat primes. The residual two-adic contribution r has precisely the stated two choices when r=0 and one otherwise. No existence of a repair height for every q follows. Suggested wording clarification: q=17,a=4 demonstrates that q² must be retained to enumerate all preimages; it does not show that criterion(ii) is needed for membership there, since j=3 gives137 prime and272=phi(548) also.

The independent extended verifier directly propagates raw convolution carries (rather than defining carry mass from the claimed identity), checks grouped-base identities and factor-weight transport, the binary complement range, the finite all-base sieve, and decimal seed lengths4..200. Its separate finite dyadic membership check uses explicit sufficient witnesses, so it does not silently claim its scan cutoff exhausts all inverse fibers. The frozen core audit and receipt are separate files and unchanged.
