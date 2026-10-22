# Independent audit: odd-exponent finite-descent densities

2026-10-22. **Verdict: ACCEPTED within the stated scope.**

Audited the phase2 `DENSITY_ONE_DESCENT.md` and `verify_density.py`. Exact frozen copies are in `frozen_inputs/`; their SHA-256 values are recorded by `results.json`. No producer functions are imported by the independent checker.

## Proof review

1. **The distributional input is valid.** There is an independent short proof of the exponent/source residue bijection: for distinct positive odd k,l,
   
   v2(u_k−u_l) = v2(k−l), where u_k=(3^k+1)/4.
   
   Indeed, take k>l and k−l=2d. Then v2(9^d−1)=3+v2(d), so division by4 gives1+v2(d)=v2(k−l). The identity for9^d follows by factoring off the odd part of d, then repeatedly using v2(9^(2^j)+1)=1. Thus the induced map on odd residue classes modulo2^L is injective and hence bijective. Combining this with the exact valuation-word source cylinder modulo2^(A+1) gives density2^(−A) among odd exponents. This is finite-parameter counting, not independence along an individual integer orbit.

2. **The carry threshold is uniform.** The recurrence beta'=(3 beta+1)/2^a, a>=1, implies0<beta_h<=(3/2)^h−1. For s_h=min{s:2^s>3^h}, every source above T_h=((3/2)^h−1)/(1−3^h/2^s_h) has endpoint descent exactly when A>=s_h. All actual valuations, however large, obey the same threshold. Since u_k is increasing, the exceptional sources induce finitely many exponents for each fixed h.

3. **The exact density is a finite-cylinder calculation.** The bad endpoint set differs by finitely many exponents from the finite union of words whose total A is in[h,s_h−1]. The composition count is binom(A−1,h−1). Its density is q_h=sum binom(A−1,h−1)/2^A. Taking a complement gives1−q_h. No infinite union of natural-density sets is used.

4. **The exponential estimate is rigorous.** Compare this finite bad-word sum with the convergent formal geometric generating series at t=3/4. Since3^5<2^8, bad words satisfy5A<8h. Hence q_h <= (3/5)^h(3/4)^(−8h/5), whose fifth power is(65536/84375)^h. The formal probability language is optional; this is also a direct inequality of positive series.

5. **The density-one passage is correct.** For every fixed h, endpoint descent is contained in the event of some finite strict descent. Its lower density is therefore at least1−q_h. Letting h tend to infinity gives natural density1, without interchanging an infinite union and density.

6. **Each of the six registers is retained.** A word fixes one odd exponent class modulo2^(A+1), including its mod4 value. A compatible fixed odd class modulo12 adds only the independent mod3 constraint. The word has conditional density2^(1−A); an incompatible word has density0. Consequently each register's bad density is at most2q_h. These are source-register constraints, not an ownership or arrival-time theorem.

## Independent tests

Run `python3 audit_density.py` and `python3 -O audit_density.py`. The outputs are byte-identical.

- 20,622 modular exponent/source cases: complete exponent periods modulo lcm(12,2^s_h), for h1..8, checked against the full bad density and all six conditional register densities. This uses modular powers and literal valuation calculations, not the producer's composition enumeration.
- 12,040 literal family endpoint checks for odd k1..601, h1..40: affine identity, carry bound, and descent classification above the threshold.
- Exact tail identity, rational Chernoff inequality, and sharp all-ones carry bound through depth300.
- Replayed both normal and optimized producer checkers. Both match the supplied producer verification exactly.

`results.json` records the independent results. `producer_replay.json` and `producer_optimized_replay.json` record the separate replay, not a replacement for the independent tests.

## Interpretation and limits

No correction is required. In this packet, A is the sum of valuations; the ordinary unaccelerated Collatz clock is A+h, as specified in phase1. The phrase 'clock threshold' must not be interpreted as A ordinary Collatz steps.

The reported endpoint densities need not increase with h and are not a first-descent CDF. Finite descent does not prove convergence: the descended integer need not remain in the family. Density one does not prove every exponent descends. Exponent density cannot be transferred to source-height density or to the153N ledger. The report explicitly maintains all these boundaries.

The all-order audit rests on the proofs above, not on the finite checks. No new literature-priority claim or external theorem beyond the frozen report's stated scope is endorsed.
