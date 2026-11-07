# Session219: Fermat factors, totient fibers and Chan–Pak constructions

**CHECKPOINT1: audited totient results. Constructive tournament extensions are in progress. Carmichael’s totient conjecture and Chan–Pak completeness remain OPEN.**

## The three sequences have different types

- Even nontotients begin14,26,34,38,50,62,68,74,76. These are missing OUTPUTS of Euler’s phi.
- The corrected Camion exponent set through64 is4,6,10,12,14,26,34,38,62. These k make composite2^k−1 pass the strong-tournament congruence. The shared fragment is not equality of sequences;50 separates them.
- [OEIS A007374](https://oeis.org/A007374) gives the least target having exactly n totient preimages, indexed by n>=2. It measures multiplicity, not missingness.

[Session218](../session218_smith_repunit_totient_20261107/README.md) proves the genuine common subfamily: simultaneous Mersenne/Wagstaff primality at q>5 makes2q both an even nontotient and a Camion exponent. It also gives the exact first dyadic repair heights. The present session extends membership to full fiber counts.

## A Fermat factorization causes an exact multiplicity drop

Write A(m)=#{x:phi(x)=m}. The full inverse-fiber generating function is

    sum A(2^a)x^a = (2−x)/(1−x) product_(F_i prime)(1+x^(2^i)).

Since F5=2^32+1=641·6700417, this proves

    A(2^a)=a+2 for0<=a<=31,
    A(2^a)=32 for32<=a<=63.

Thus doubling2^31 to2^32 decreases the exact multiplicity33→32. No unknown Fermat primality is needed. A complete prime-supplier convolution counts every preimage along2^a q, including the exceptional q² supplier for Fermat q. Second differences recover supplier births. The multigraded Euler product explains why641 and6700417 cannot replace the missing Fermat generator: their totients carry odd parts5 and3·17449.

[Full proof and exact formulas](totient/FERMAT_FIBER_TRANSFERS.md) · [Independent audit](audit/DYADIC_REVIEW.md).

## Reproduce checkpoint1

From this directory:

    python3 totient/verify.py
    python3 audit/verify_dyadic_independent.py
    python3 audit/verify_supplier_independent.py

Producer:245 complete recursion comparisons and a complete small-target sieve. Independent audit: explicit pure-dyadic inverse fibers through exponent63, a complete bounded totient sieve, and130 prime-core inverse fibers. JSON receipts state the universes. The first ten A007374 minima are recovered by complete small-target computation; no new OEIS term or historical priority is claimed.

## Research direction

The useful Ford–Konyagin inspiration is controlled source construction plus exclusion of every unintended preimage, rather than copying a multiplicity slogan. The next checkpoint investigates exact repunit-factorization control, gap absorption under tournament products, and the order64 signatures of641 and6700417 in physical tournament grammars. These extensions are not premises of this checkpoint.
