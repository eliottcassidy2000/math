# Independent reconciliation: nine open orders, inherited degree 35

2026-10-21 UTC. Verdict: PASS. This review intersects the incoming degree-35 result with the separately reviewed Session 22 two-shore insertion exclusion. It does not claim a new degree-35 theorem or a completed center census. All work here is read-only with respect to the repository; review files and isolated replay inputs are under this review directory.

## Exact theorem and attribution

For a strong self-converse tournament with black Hamilton-path-orbit deficit six, the inherited small-order coefficients are one at orders 5, 6, 7 and zero at 8, 9. Every still-undetermined coefficient is supported on

    {12, 15, 18, 20, 21, 24, 25, 28, 35}.

Every coefficient in this nine-element list remains OPEN. In particular deg R6 <= 35; the paired-black trivalent series has minimal eventual recurrence order 18 and a guaranteed homogeneous recurrence onset n >= 48. This is not an earliest-onset assertion.

The incoming COMBINED_CUTOFF.md established degree 35 and ten residual orders by composing the older closing-cycle theorem, the Session 22 odd-shore/score-run proof, and the frozen order-ten certificate. Its attribution must remain intact. The new intersection here removes order 14 using the independently reviewed q=2 insertion theorem; it does not lower the degree or recurrence onset. The standalone degree-47 and later degree-45 arguments are valid weaker results, not the sharpest current bound.

## Why the older closing theorem is genuinely general

The cycle-cut note starts with deficit two, but its section 5 explicitly leaves that assumption and proves: every strong self-converse tournament of odd order n >= 7 admits a Hamilton cycle with nonautomorphic successor. Section 3 and the explicit prefix source support this independently of deficit two.

Take any Hamilton cycle. If its successor is nonautomorphic, stop. Otherwise it gives circulant coordinates with 1 in the skew connection set S. Every positional reflection in these coordinates is anti. The prefix classification supplies a nonblue closing Hamilton word, and hence another actual Hamilton cycle. Its successor cannot be an automorphism, since that would again make all its positional reflections anti.

Closure is not an implicit bypass-to-cycle conversion: F1 starts at 2 and ends at n-1, so its return increment is 3, allowed by F1. F2 starts at k-1 and ends at n-1, with return increment k, allowed by F2. Both F3 variants start at 1 and end at n-1, with return increment 2, allowed since the first negative sign has index at least 3. Unit reindexing is multiplication by an invertible allowed step; it preserves every arc, including the final return arc, and preserves failure of anti-reversal. The four order-seven paths have respective return increments 4, 4, 3, 2, each in its displayed S. Thus all branches yield cycles, not merely bypass paths. The three residual sign patterns are eliminated by these same reindexed forbidden-prefix families, including the m=4 case B two-step reindexing.

The theorem is existential in the cycle. It does not assert q>1 for every cycle. The incoming Paley7 hostile correctly retains a q=1 standard cycle and a q=7 nonblue closing cycle on the same source.

A further novelty boundary is important: the older section 5 already proves Delta >= p(n)-1 for odd n >= 7. At prime order this reads Delta >= n-1 and already dominates the Session 22 auxiliary Cauchy–Davenport prime-center exclusion. That auxiliary connectivity lemma remains valid, but its prime-center exclusion is not a newly improved center theorem.

## Exact intersection proof

For the chosen cycle let m=|Aut(T) intersect <R>| and q=n/m. Tournament automorphism groups have odd order, so m is odd. The n cycle cuts meet exactly q path orbits; at most gcd(2,q) <= 2 are blue. Hence q <= 8. Even n forces q >= 2 by oddness of m. Odd n >= 7 permits the chosen nonautomorphic-successor cycle above, also forcing q >= 2.

The inherited odd-shore/score-run result applies to that actual chosen cycle: m lies in {1,3,5,7}; q=6,8 with m>=5 is excluded; and (q,m)=(7,7) is excluded. At orders >=10 the complete surviving parameter list is

    (n,q,m) = (10,2,5), (12,4,3), (14,2,7),
              (15,3,5), (15,5,3), (18,6,3), (20,4,5),
              (21,3,7), (21,7,3), (24,8,3), (25,5,5),
              (28,4,7), (35,5,7), (35,7,5).

The reviewed insertion theorem proves Delta >= 2m-2 for q=2 and odd m>=5. Therefore both q=2 rows are impossible at deficit six. Removing them leaves exactly the nine stated orders. The frozen 128-frame order-ten census independently agrees with the analytic order-ten exclusion, but is no longer necessary to remove ten once insertion is used. Neither parameter rows nor this arithmetic establish existence at any surviving order.

## Recurrence boundary

With p=1-z-z^3, u=p(z^2), the inherited numerator over p^2 u^2 is

    N = z^5(u^2-p^2u)/2 + p^2u R6 + z^10 p^2(1+z+z^3).

The degree bound 35 gives deg N <= 47, hence recurrence from n=48. The coprime squarefree factors p,u do not cancel: at a root of p the first term is nonzero, and at a root of u the final term is nonzero, independently of R6. Thus the minimal eventual order remains 18. The formal control R6=z^5+z^6+z^7+z^35 still respects the nine-order support and has numerator degree 47 with leading coefficient -1. This certifies the sharp uniform algebraic bill for the stated support, not existence of an order-35 center or the actual earliest onset.

## Consumer guards, pins, and reproduction

The incoming consumer uses ck raising ValueError, not assert, so its guards survive Python optimization. It pins eight inherited proof/receipt inputs after LF normalization; verifies the complete 128-mask domain; independently reconstructs every source; checks all stored automorphisms and anti-involutions; verifies free divisions, Hamilton words, source-bound targets and the inherited deficit histogram; checks the Paley7 cycle-choice hostile; and checks exact support and formal recurrence arithmetic. It does not recompute exhaustive H or group completeness, and correctly declares itself a consumer of the pinned producer and receipt, not a new census. Some properties of black ports and completeness rely on those exact pinned producer bytes rather than independent new tests, consistent with its scope.

Copied unchanged incoming consumer inputs were replayed in an isolated tree. Normal and optimized --verify both passed 16,889 gates and matched combined receipt LF SHA256 9a1bdb64fdc717cb3af1aa9c2c65eef5b59978786f161969fbe1147411d8e048. The older cycle-cut source was initially sparse-absent and read successfully using git show HEAD:path; no proof blocker remains. Source hashes, including the q=2 proof and independent review, are in source_pins.json. The nine-order intersection checker is separate, leaves the incoming ten-order receipt unchanged, and runs normally and under -O.
