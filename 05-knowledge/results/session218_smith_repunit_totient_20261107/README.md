# Session 218: Smith numbers, repunits and totient scaling

**ACTIVE research. Carmichael's totient conjecture and Chan–Pak completeness remain OPEN.**

## First audited result

If both M=2^q−1 and W=(2^q+1)/3 are prime and q>5, then **2q is an even nontotient**. This connects the repo's Mersenne/Wagstaff Camion family to the even-nontotient sequence without identifying their full sets.

Proof in one paragraph: if ell=2q+1 were prime, Fermat would make it divide 2^(2q)−1=3MW. Thus ell would equal M or W. Both exceed ell for q>5, a contradiction. For an odd prime q>3, the complete preimage set of 2q under Euler's phi is {2q+1,2(2q+1)} when 2q+1 is prime, and empty otherwise. The cases q=3 and5 are genuine exceptions.

[Full proof and scaling law](totient/PROOF.md) · [Independent audit](audit-core/CORE_REVIEW.md).

## The prefix needs a definition

14,26,34,38,... starts the even nontotients, followed by50,62,68,... . The same four appear in the repo's Camion Mersenne exponent table, but that full corrected table through64 is **4,6,10,12,14,26,34,38,62**. These are different objects. The theorem above proves a common subfamily, rather than choosing an interpretation from four numbers. Historical Carmichael, repunit and tournament sources are being reconciled for the second checkpoint.

A nontotient is a missing OUTPUT of phi. Carmichael's totient conjecture says every nonempty output fiber has at least two INPUTS. Empty fibers cannot refute it. Carmichael numbers and the Camion strong-tournament congruence are separate predicates.

## What scaling preserves

The exact law is phi(tn)=t phi(n) times the product of (1−1/p) over primes p dividing t but not n. A fiber can be scaled uniformly only with the relevant support sidecar. The examples phi(3)=phi(4)=2 but phi(6)!=phi(8), and14 missing while28=phi(29), prevent unconditional scaling claims.

## Reproduce

Source-path mentions inside frozen producer and audit reports record their original provenance. Use the following portable commands from this directory:

    python3 totient/verify.py
    python3 audit-core/verify_core.py

The producer uses SymPy; the independent verifier uses only the standard library. The complete totient sieve cutoff is justified in the proof, not guessed. Receipts are included beside each verifier. Historical Camion censuses are inherited, not newly re-enumerated. No infinitude of simultaneous Mersenne/Wagstaff primes is assumed. Further Smith, base-change and dyadic results are under final integration.
