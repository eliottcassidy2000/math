# Independent audit: exact dyadic totient fibers

**PASS.** The statements below are elementary proved identities, with independent exact computation. They do not establish a new value of A007374, minimality of any target, the existence of only five Fermat primes, or Chan–Pak completeness.

## Definitions and proof audit

Write A(m)=#{n>=1: phi(n)=m}. An odd prime dividing a preimage of 2^a must have p−1 a power of two; primality then forces p to be a Fermat prime. Its exponent in the preimage is one. Conversely every product of distinct Fermat primes, times an admissible power of two, is a preimage. If n has 2-adic exponent e, the cost in v2(phi(n)) is zero for e=0,1 and e−1 for e>=2. This gives the exact formal power series

    sum_(a>=0) A(2^a)x^a = (2−x)/(1−x) product_(i: F_i prime)(1+x^(2^i)).

Formal infinite products are legitimate coefficientwise: only finitely many weights contribute to each exponent. F_i denotes 2^(2^i)+1, not its individual prime factors.

The five primes F0,...,F4 have weights 1,2,4,8,16. Their subset sums are exactly 0,...,31, once each. F5=641·6700417 is composite; no next candidate can have weight below64. Therefore, unconditionally,

    A(2^a)=a+2 for 0<=a<=31;
    A(2^a)=32 for 32<=a<=63.

The boundary a=0 is correctly included: phi^{-1}(1)={1,2}. The drop A(2^31)=33 to A(2^32)=32 is real. The factorization does not license an all-exponent plateau without further Fermat compositeness information.

## Unique odd-prime supplier

For an odd prime q, a preimage of 2^a q contains exactly one source of the odd factor q. That source is either a prime p=2^j q+1, j>=1, appearing to the first power, or q itself appearing to the second power, with q−1 a power of two. In the second case q is a Fermat prime. No other odd prime may have exponent greater than one; no two suppliers can occur together. The remaining primes are distinct Fermat primes. In the q² case the residual Fermat subset must exclude q.

With B=(2−x)/(1−x), P=product_(i:F_i prime)(1+x^(2^i)), S_q=sum_(j>=1:2^j q+1 prime)x^j, the exact identity is

    sum_(a>=0) A(2^a q)x^a = B P S_q
      + 1_(q=F_r) B x^(2^r) P/(1+x^(2^r)).

Both sides count integers, not supplier descriptions with multiplicity: the unique-supplier proof ensures disjoint cases, and unique factorization ensures injective descriptions within each case. The a=0 coefficient is zero, as required for odd targets q>1.

## Independent checks

- `python3 verify_dyadic_independent.py`: enumerates every admissible support and 2-part for all 64 exponents a=0,...,63, checks exact phi values and uniqueness, certifies all displayed primes by trial division and F5 by multiplication. A complete totient sieve to 2,097,152 independently recovers every preimage for a=0,...,10. Completeness uses phi(n)>=sqrt(n/2).
- `python3 verify_supplier_independent.py`: builds complete inverse fibers through a general divisor/candidate-prime recursion, without assuming unique suppliers, then compares the formula on 130 targets: q in {3,5,7,11,13,17,19,23,29,31}, a=0,...,12.
- Both scripts pass. JSON receipts preserve their exact universe and results. Finite checks support the proofs, not a conjectural extension.

## Scope and provenance

This audit is independent of the producer's source code. The supplier statement was supplied in a worker message after the independent verifier had been written; its formula agrees exactly. The prior Session218 support-aware classification already supplies the existence criterion. The present new refinement is exact fiber multiplicity and its Fermat-hole convolution, not reinvention of that existence criterion. Sources and literature claims belong to the producer's source note and are not independently certified by these arithmetic checks.

## Producer-artifact addendum

Read `session219-totient/FERMAT_FIBER_TRANSFERS.md` in full. Equation (6), the multigraded Euler product, is valid coefficientwise: for a fixed target m only primes with p−1|m can occur, and the sole nontrivial constant multiplicity comes from p=2. This is a well-defined coefficient model, not an assertion of analytic convergence. The term 641 contributes a−5 for7<=a<=38 and32 for39<=a<=70 to A(2^a·5), since its cost is7 and the *remaining* pure-dyadic exponent stays in0,...,63. These ranges do not claim that the whole q=5 fiber has that size. All displayed identities and status boundaries pass.

The later deconvolution identity (7) also passes symbolically. Modulo x^64, P=(1−x^32)/(1−x), so (1−x)^2 BPS_q=(2−x)(1−x^32)S_q. With all negative-index sequences zero, comparing coefficients gives exactly the stated second-difference recurrence. Fermat-q exceptional counts must be subtracted first, as the producer states.
