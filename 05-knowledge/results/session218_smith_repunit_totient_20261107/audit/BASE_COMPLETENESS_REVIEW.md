# Independent audit: complete moving-base Smith existence

PASS for the proposed theorem: a composite positive integer n is Smith in at least one integer base b≥2 if and only if n is not6,8,9,12. This is an elementary proof review, not a historical novelty claim.

Write F(n) for the sum of prime factors with multiplicity. For every composite n other than4,6,8,9,12, the specific base b=n−F(n)+1 works. The case4 instead works precisely in bases b≥5.

If n=2p for prime p≥5, then F(n)=p+2, b=p−1≥4, and n has digits22. Both prime factors are single-digit except p, whose digits11 give weight2; the total factor weight is4. Thus this branch needs its own proof and must not be included in the later assertion that every factor is below b.

For all remaining n it suffices to consider n>12: the other composites below13 have already been handled or are the four proposed exceptions. If n is odd, choose a prime divisor p with3≤p≤sqrt(n). Subadditivity under factorization gives F(n)≤p+n/p≤3+n/3. If n is even but not twice a prime, n/2 is composite and F(n)=2+F(n/2)≤n/4+4. Both upper bounds are strictly less than(n+2)/2 for n>12. Therefore b>n/2, every prime divisor is at most n/2<b, and0<F(n)−1<b. The representation n=b+(F(n)−1) has digit sum F(n), equal to the prime-factor digit sum. This proves existence.

Nonexistence uses the necessary congruence b−1 | n−F(n). Candidate bases and (number digit sum, factor digit sum) are:
- n=6: b=2 gives(2,3).
- n=8: b=2 gives(1,3); b=3 gives(4,6).
- n=9: b=2 gives(2,4); b=4 gives(3,6).
- n=12: b=2 gives(2,4); b=6 gives(2,7).
These exhaust the candidates, so the exclusion list is exact. The strict inequality at n>12 cannot be weakened blindly: n=9 and12 hit the relevant failure boundary.

Independent finite check: `python3 check_base_existence.py`, receipt `base_existence.out`. A smallest-prime-factor sieve checks all composites through200000 using direct digit expansions, the specific constructive base, branch conditions, and complete exceptional candidate sets. This supports the algebraic proof but is not the reason the global claim holds.

Full producer artifact `session218-smith/BASE_COMPLETENESS.md` subsequently reviewed: PASS with no mathematical changes required. Its proof correctly isolates the twice-prime branch and all exceptional candidate bases. Producer verifier `check_base_completeness.py` was also rerun successfully; receipt is `producer_base_completeness_rerun.out`.
