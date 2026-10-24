# The two-step prime-power tower: exact evidence and an open boundary

Status: FINITE-EXACT certificate; the general higher-depth statement is OPEN. This packet does not alter the separately frozen conditional base theorem.

## Precisely stated hypothesis

For every prime p, integer r≥2, and n≥0,

c(p^r n) ≡ c(p^(r−2)n) (mod p^(2r−2)).

Here c(0)=1,c(1)=2 and (k+1)²c(k+1)=(26k²+13k+2)c(k)+3(3k−1)(3k−2)c(k−1). The base r=2 is the conditional theorem in sources/base-two-step-proof.md. At odd depth r=3 the target is c(pn), not c(n). The n=0 case is trivial.

Equivalently it suffices to consider n prime to p: writing n=p^a m yields depth r+a. The full hypothesis would then imply the stronger valuation bound 2r+2a−2 for the displayed difference. The verifier checks this strengthened bound whenever n is divisible by p, rather than omitting such n.

## Exact bounded evidence

verify.py builds all c(k), 0≤k≤20000, using arbitrary-precision integers and checks exact division at every recurrence step. It tests every prime p≤floor(sqrt(20000))=141, every depth r≥2 with p^r≤20000, and every positive n with p^r n≤20000. It records minimum valuations, the first minimizing n, and a nonzero normalized residue modulo p certifying exactness. All checks survive Python optimization; normal.json and optimized.json are identical.

Selected sharpness and exceptional extra precision:

- p=3: every depth r=2,...,9 has valuation exactly 2r−2 at n=1. Thus both parities support the claimed exponent sharply; a uniform additional digit is already false at each of these tested depths.
- p=2: r=2 has minimum valuation2. Every tested depth r=3,...,14 instead has minimum valuation2r−1, attained at n=1. This is an observed extra digit, not a proved dyadic refinement.
- p=5: depths2,4,6 have minimum valuations2,6,10; depths3,5 have minima5,9, one above the proposed bound.
- p=7: depths2,4 have minima2,6; depths3,5 have minima5,9.
- p=11: depths2,3,4 have minima2,4,6. Inert p=13,17,19 at depth3 have minimum4.

All available n, including p-divisible n, are included. These are overlapping tests: the same numerical difference can be represented at more than one (r,n), and the comparison count is not a count of independent samples.

## Relationship to existing conjectures: no relabeling

The inherited root note, section8, already records the split/ramified one-step conjecture

c(p^j n) ≡ c(p^(j−1)n) (mod p^(2j)).

On those primes, its cases j=r and j=r−1 sum to precisely the two-step bound p^(2r−2). For r=2 the lower step is the inherited proved mod-p² law. Thus the tower above is not a new stronger conjecture at split primes or p=7; it is a weaker consequence of that existing one-step conjecture. The substantive extension is a common two-step formulation that also covers inert primes, where the one-step law fails even modulo p.

The repository search in repository-search.txt found the prior one-step higher-power conjecture and its explicit unproved status in the root and Beukers-extension notes. A bounded public-source search on 2026-10-24 found no exact two-step formulation, but is not a novelty proof:

- OEIS A279619, https://oeis.org/A279619 , identifies the root sequence and links O'Brien's thesis and Chan–Cooper–Sica Conjecture5.4; the retrieved entry discusses Lucas congruences, not this tower.
- Alinquant–Osburn, On sporadic sequences, https://www.ucd.ie/mathstat/t4media/AlinquantOsburn.pdf , discusses the original sporadic sequences, including s7. Those are not automatically statements for the square root G. No transfer from F=G² to G's two-step tower is asserted.
- Searches used “A279619 congruence prime powers”, “A279619 supercongruences”, and “Cooper square root supercongruence”, plus an O'Brien conjecture query. This is a narrow search, not a completed literature review.

## Why the base proof does not prove the tower

The base proof has two precision-limited inputs: U_p²f≡f mod p² and an excellent-lift correction ε_p H_p congruent modulo p to a small-degree polynomial. Applying that proof at n replaced by p^(r−2)n proves only modulus p²; it does not produce p^(2r−2).

At precision p^4 and beyond, the first-order expansion of t(q^p)/t(q)^p is insufficient. Higher powers and products of ε, higher derivative terms, and the nonzero p² corrections in the exact Eisenstein U_p relations enter. The existing degree bound kills the first-order Cartier contribution only. It provides no degree/Cartier bound for the higher error terms, nor their required cancellations with U_p²f−f.

A plausible route is a depth-dependent Frobenius or crystalline congruence compatible with the two-dimensional companion system f↔f*, retaining all higher corrections and controlling their Cartier image. The inherited note itself flags the missing higher-precision truncation-ratio congruence. No such theorem has been derived in this packet. Ordinary iteration, a formal squaring of the root series, and the numerical valuation pattern are not proofs of the proposed precision gain.

## Reproduction

Run python replay.py from any current directory using its path. It validates exact recursive membership and SHA256 hashes, then reproduces normal and optimized fixtures. Run python hostile.py to test rejection of a corrupted proof, an omitted fixture, and an undeclared nested file. No assertions are used for critical checks.
