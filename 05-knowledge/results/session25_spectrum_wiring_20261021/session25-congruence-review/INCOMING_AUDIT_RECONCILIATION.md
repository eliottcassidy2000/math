# Reconciliation with incoming independent audit

## Subsequent additive closure (Session25)

A separate argument, now fully written and independently audited in ADDITIVE_PRIME_POWER_PROOF.md, proves s(n+phi(p^r))=(2/p)s(n) modulo p^r for every odd p, r>=1, and n>=2r-1, with the same law for D. Thus the incoming additive conjecture is now proved by this new filtered-Hurwitz/Newton argument. The earlier caution that the multiplicative theorem alone does not settle it remains correct. The new proof also gives arbitrary-unit-base closure and the sufficient connected-graph bound n>=2r. No minimal-threshold or literature-priority claim is made.

**Current source:** commit 3596c5d0b1990b05ae960e980a863956b47cc6b5, MISTAKE-651, revised THM-4746 and THM-4747, and `05-knowledge/results/tournament_pseudoprimes_audit_20261021.md`.

The incoming audit already corrects the 262-class denominator, the Z_3 representative boundary, the Mersenne table, and the false connected-graph negative control. It also proves the square-modulus additive period and deduces the prime-square equivalence. These corrections and their provenance must be preserved. Our pre-audit patch is obsolete and must not be applied.

## Connected graphs: the initial term is the entire obstruction

Write b(n) for the count of connected labelled simple graphs, with b(1)=1. Let an odd-prime cycle act on a selected p-element block among m+p-1 labels. A fixed graph has a circulant simple graph inside that block and uniform adjacency between the block and each outside vertex. The internal block has a=2^((p-1)/2) choices.

If m>=2, collapse the block. The original graph is connected if and only if the quotient graph on m labels is connected, for every internal choice:
- Every nonempty circulant graph of prime order is connected, since any nonzero step generates the group.
- If the internal graph is empty, a connected quotient gives the contracted vertex an outside neighbor. Its uniform adjacency connects all p internal vertices through that neighbor. Quotient paths then lift.
- Conversely, an original connected graph has connected quotient.

Thus there are exactly a*b(m) fixed connected graphs, yielding b(m+p-1)=(2/p)b(m) modulo p for m>=2.

At m=1, the empty internal graph is disconnected and all other a-1 choices are connected. Hence b(p)=a-1 modulo p, differing from (2/p)b(1) by -1. Testing only from the initial index falsely suggests no tail law.

This elementary check agrees with incoming Lemma C, which gives a Hurwitz-ring proof and credits Mani–Stones 2016. The literature attribution is preserved from that source, not independently bibliographically verified in this audit. Our independent recurrence tests cover 990 tail cases and all ten corresponding initial obstructions for odd primes through 31.

## What remains genuinely additional

Incoming Lemma B proves additive quasi-periodicity modulo p^2 with shift p(p-1) and threshold n>=3. Its higher additive law modulo p^k remains a conjecture. Our PRIME_POWER_PROOF.md proves the distinct multiplicative law

s(mp^r)=(2/p)^m s(mp^(r-1)) modulo p^r,

and the same for D, for all m,r>=1. It supplies an alternative square proof and an all-prime-power dilation result. It does not prove the additive higher-power conjecture.

The exact closed-presentation congruence and cycle-stabilizer obstruction remain independent useful consequences. Aggregate count congruences still do not determine individual H(T), S2/S3 spectra, or wiring.

## Remaining editorial boundaries

The revised canon still benefits from the usual explicit singleton Hamiltonicity convention and from distinguishing multiplier classes from full isomorphism classes at composite order. The p=3 endpoint of the 3p family is justified by the now-proved square law, even though its displayed representative cannot literally equal 3 modulo p-1. These are residual wording issues, not grounds to replace the incoming audit wholesale.

The incoming Mersenne table is source-qualified as corrected: exponents 4, 6, 10, 12, 14, 26, 34, 38, 62 pass through 64, and exponent 59 is decided. Our independent checks confirm the small boundary 63; they do not independently reproduce the entire corrected table or the expanded search through prime 10000.
