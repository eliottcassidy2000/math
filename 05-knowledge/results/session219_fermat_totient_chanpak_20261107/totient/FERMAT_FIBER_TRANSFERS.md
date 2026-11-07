# Fermat holes, totient multiplicities and dyadic supplier convolution

**Status: PROVED elementary identities and window theorems; FINITE-EXACT regression.** No claim of historical priority. Chan–Pak completeness and Carmichael's totient conjecture remain OPEN. This extends Session218's support classification from membership to exact fiber counts.

## Sources, inheritance and the typed map

Write A(m)=#{n>=1: phi(n)=m}. [OEIS A007374](https://oeis.org/A007374) is a(n)=min{m:A(m)=n}, indexed by n>=2. It is neither the nontotient sequence nor the Camion exponent sequence. Its Ford citation concerns existence of every positive multiplicity >=2; no such theorem for tournament fibers is imported. The independent small computation below verifies its first ten terms, not the complete published table.

Inherited result: Session218 `totient/DYADIC.md` classifies all preimages of 2^a q via a unique q-supplier. The closest mechanism is support-aware Euler-product scaling. Hostile: odd/even preimages can scale differently, so neither equal fibers nor multiplicity is preserved by multiplying every source by two. The least-used sidecar is the complete prime-valuation vector of the OUTPUT phi(n).

Our concept board: zero fibers; exact multiplicity; Fermat binary subset support; prime q-suppliers; mixed odd-core sectors. Source is prime-power factorization of n, map is its totient valuation vector, preserved predicate is exact phi(n)=m, lost data after scalar projection is prime support and supplier identity. Restore the latter to count without collisions. No map from these fibers to tournament isomorphism classes is asserted.

## 1. Pure dyadic generating function

Let F_i=2^(2^i)+1, let I={i>=0:F_i is prime}, and set

    P(x)=product_{i in I}(1+x^(2^i)),
    B(x)=(2-x)/(1-x)=2+x+x²+... .

Then, as a formal power series,

    sum_{a>=0} A(2^a)x^a = B(x) P(x).                 (1)

Proof. If phi(n) is a power of two, every odd prime p|n has p−1 a power of two, and occurs to exponent one. If 2^h+1 is prime, h must itself be a power of two: an odd factor of h would give the ordinary odd-power factorization. Thus odd support is a subset of Fermat primes, contributing distinct weights 2^i to v_2(phi(n)). The 2-part 2^e contributes max(e−1,0): two choices e=0,1 at weight zero, one at every positive weight. This is B. Prime factorization gives a bijection, not merely a count estimate. All coefficients are finite; only i with 2^i<=a matter.

Equivalently, if c(h) is the indicator that the binary support of h lies in I, then

    A(2^a) = 2c(a) + sum_{0<=h<a} c(h).             (2)

Binary subset sums are unique, so c(h) is an indicator, not an arbitrary partition number.

## 2. A missing Fermat prime causes a precise drop and plateau

Suppose F_0,...,F_(r−1) are prime and F_r is composite; put L=2^r. For the entire window 0<=a<2L, regardless of later Fermat numbers,

    A(2^a)=a+2                 if 0<=a<L,
    A(2^a)=L                   if L<=a<2L.         (3)

Proof. Below weight 2L, permitted binary supports give precisely h=0,...,L−1. Insert this into (2). Later possible prime generators have weight at least 2L and cannot affect the window.

In particular F_0,...,F_4 are prime and

    F_5 = 2^32+1 = 641 * 6700417,

so A(2^a)=a+2 for 0<=a<=31 and exactly 32 for 32<=a<=63. The exact multiplicity drops from 33 to 32 under doubling 2^31 to 2^32. This is not a numerical artifact: the lost final binary subset had two 2-part choices, and once it is below the target it would have only one; the new endpoint subset is absent. Membership remains true everywhere in this pure dyadic family, while multiplicity decreases. The boundary a=0 gives A(1)=2.

Formula (3) supplies the unconditional nonminimal bound a(n)<=2^(n−2) for 2<=n<=33 in OEIS notation. It makes no claim that these are minimum targets. Already n=6 has a(6)=12 rather than 16. It supplies no proof of all-multiplicity existence beyond its window.

## 3. Exact prime-odd-core convolution

For an odd prime q define

    S_q(x)=sum_{j>=1 : 2^j q+1 prime} x^j.

Then

    sum_{a>=0} A(2^a q)x^a
      = B(x) P(x) S_q(x)
        + 1_{q=F_r for some r} B(x) x^(2^r) P(x)/(1+x^(2^r)).  (4)

Proof. Exactly one prime-power factor of n supplies the single q in phi(n). One possibility is a prime p=2^j q+1 appearing to exponent one; it contributes x^j, and all other odd support consists of distinct Fermat primes. Such p is never a Fermat prime since q is odd. The other possibility is q², which is allowed precisely when q=F_r, contributes x^(2^r), and excludes q from remaining Fermat support. The cases are disjoint and exhaustive. The 2-part is counted by B. This is the full Session218 source-shape classification with each case counted once.

Set w(t)=0 for t<0; w(t)=t+2 for 0<=t<=31; w(t)=32 for 32<=t<=63. For non-Fermat q and a<64, (4) becomes the explicit finite convolution

    A(2^a q) = sum_{1<=j<=a : 2^j q+1 prime} w(a−j).            (5)

For Fermat q, add the exact excluded-support coefficient in (4), not the full w term. At a=1, q=3 yields A(6)=4, whereas a naive nonexceptional supplier count yields only 2. At q=17,a=4 the q² supplier produces additional preimages even though an earlier prime supplier already established membership. Thus a least repair height is not enough to reconstruct multiplicity.

All missingness statements along a fixed odd-prime ray reduce to absence of all suppliers up to a (and absence of a permitted q² supplier). Once a supplier exists, membership persists at all greater a, although the number of preimages need not be monotone. No eventual supplier existence for every q is asserted.

### Exact deconvolution: the count sequence recovers its prime suppliers

For non-Fermat q, write C_a=A(2^a q), s_a=1 if 2^a q+1 is prime for a>=1 and s_a=0 otherwise (also set C_a=s_a=0 at negative a, and s_0=0). For 0<=a<64,

    C_a−2C_(a−1)+C_(a−2)
      = 2s_a−s_(a−1)−2s_(a−32)+s_(a−33).         (7)

Indeed in this window B(x)P(x)=(2−x)(1−x^32)/(1−x)^2. Multiply (4) by (1−x)^2 and compare coefficients. Thus the multiplicity sequence recursively determines every supplier indicator in the window. Below32 the formula simplifies to second_difference(C)_a=2s_a−s_(a−1); the F5 hole introduces a precisely delayed correction, not arbitrary noise. For Fermat q the same statement holds after subtracting its explicit q² contribution from C.

This gives a stronger boundary than a heuristic analogy: exact multiplicities contain the primality-supplier information, whereas a zero/nonzero membership sequence records only the first supplier height. A method which replaces counts by membership cannot reconstruct later births without additional data.

## 4. Where 641 and 6700417 actually go

These are NOT Fermat primes. Their arithmetic retains the lost odd coordinates:

    phi(641)=640=2^7*5,
    phi(6700417)=6700416=2^7*3*17449,
    phi(F_5)=2^14*3*5*17449=4288266240.

The verifier proves primality of 641,6700417,17449 by integer trial division, and checks the factorization exactly. Both factors have ord_p(2)=64, since 2^32=-1 modulo p and the order divides64 but not32. This common multiplicative order is compatible with DIFFERENT totient odd sectors; it cannot identify those sectors.

641 is exactly the j=7 supplier for q=5 in (4). For 7<=a<=38, its disjoint contribution to A(2^a*5) is a−5; for 39<=a<=70, its contribution is exactly32. These ranges follow by applying the complete pure window w(a−7), not by claiming the whole q=5 fiber has that size. The contribution beyond70 requires further Fermat data. The other F5 factor lies in odd core52347=3*17449, outside the single-prime sector; it cannot substitute for 641 or for missing pure generator weight32.

More generally, use variables X_l for each prime l and write W(m)=product_l X_l^v_l(m), with W(1)=1. Then the coefficient identity

    sum_{m>=1} A(m) W(m)
      = product_{p prime} [1 + W(p−1)/(1−X_p)]                (6)

is valid coefficientwise. The p-factor chooses absence or exponent e>=1, with phi(p^e)=(p−1)p^(e−1). The p=2 factor is B(X_2), including its doubled constant term. For each fixed m, only p−1|m may contribute, so the coefficient is finite. Taking all odd-variable exponents zero gives (1); selecting q-exponent one and every other odd exponent zero gives (4). Replacing a composite Fermat number by its prime factors changes the monomial, so no false generator substitution survives this multigraded model.

## 5. Constructive inspiration and its boundary

The useful operation is not “totients and tournaments share gaps.” It is: choose a primitive supplier, retain its support/valuation interface, freely attach independent compatible generators, and prove source uniqueness before counting. In the totient setting, prime factorization proves uniqueness and binary weights make an explicit coefficient kernel. In a tournament grammar, ordered-sum multiplicativity alone does not provide unique factorization of the chosen source category or ensure every target has a supplier. Any transplant must establish those two properties for its own endpoint/profile/isomorphism-class states.

The Fermat-hole theorem is a hostile control against a proposed “more scaling always gives more realizers” proof. It also suggests separating three obligations: supplier existence, extension-kernel coverage, and collision control. Equation (4) solves all three here only because the complete inverse-fiber classification is already proved. It is a template for a precise constructive lemma, not a Chan–Pak theorem.

## 6. Reproduction, universes and controls

Run `python3 verify.py` from this directory. `receipt.json` is generated beside it.

* All primality checks use exact integer trial division, not probable-prime tests.
* Pure coefficients through a=63 use complete allowable Fermat support; no claim about F6 or later is needed. An independent candidate-divisor/prime-power recursion checks a=0..18.
* The same recursion checks (4) at a=1..18 for q=3,5,7,13,17,19,31:126 comparisons. It derives possible primes from divisors d|m with p=d+1, rather than importing the supplier formula.
* A complete phi sieve to20000 verifies that recursion for every target1..100, including the even nontotient prefix through76. Completeness follows from phi(n)^2>=n/2, hence a preimage of m<=100 has n<=20000.
* This complete window certifies A007374 minima for multiplicities2..11:1,2,4,8,12,32,36,40,24,48.
* Exceptional q=3 and q=17, F5 mixed sectors, and the multiplicity33→32 drop are explicit hostile controls. Finite checks support the proofs but do not replace them.
