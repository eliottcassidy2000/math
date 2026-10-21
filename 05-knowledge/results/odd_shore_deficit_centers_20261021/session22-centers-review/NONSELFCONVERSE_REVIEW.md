# Independent review: non-self-converse quadratic corollary

Date: 2026-10-21 UTC. Verdict: PASS as a necessary order bound, not an atom classification or census.

Reviewed session22-root/NONSELFCONVERSE_QUADRATIC.md separately from the self-converse center theorem. The proof uses strong-tournament Hamiltonicity/pancyclicity, the explicitly demonstrated rotation-shore construction, and odd order of tournament automorphism groups.

1. For a fixed marked Hamilton cycle, equivalence of two cut words requires the unique positional correspondence, namely a rotation. Thus the number of cut orbits is exactly q and q<=h. The automorphic rotation subgroup has odd order m and n=qm. No blue word exists for a non-self-converse tournament.
2. If q=1 then successor rotation is an automorphism, so the tournament is cyclically invariant. Negation of the cyclic coordinates exchanges each directed difference with its negative and is an anti-automorphism. Hence non-self-converse T necessarily has q>=2.
3. For odd m>=r+1, the regular strong rotation shore contains an r-cycle, which is proper. Literal marked-cycle matchings show its complement is strong. Another shore gives equal row sums k, so exact half-balance would imply rk=(r/2)m and hence 2k=m, impossible.
4. At least r/2+1 qualifying exits give that many Hamilton words with a common suffix. Any two differ by a nonidentity power of an r-cycle, whose order is even because r is a power of two. Such a permutation cannot be an automorphism, so these are distinct orbits. No blue/mate analysis is required, and no n>=2r half-support hypothesis is used.
5. Entries are handled in the opposite tournament. Reversal of words gives a bijection of path orbits back to the original tournament; Aut(T)=Aut(T^op), h is unchanged, and no ambient anti-automorphism is assumed. Thus r/2+1>=d+1 contradicts h<=d in both cases.
6. Odd m cannot equal even r; consequently m<=r-1 and n<=d(r-1). For d>=2, least-power-of-two minimality at threshold 2d gives r<4d, implying the stated strict bound below 4d^2.

The final reviewed version settles d=1 internally: q>=2 and q<=h<=1 directly contradict one another. This incorporates the simplification suggested during review and no longer invokes an inherited rigid classification at that endpoint. The source status is now PROVED; this is an editorial/dependency simplification, not a change to the bound.

For combining with the self-converse bounded-defect theorem, the needed link is Delta=h-b<=h, since b>=0. Hence bounded h is covered by that theorem on the self-converse side and the present theorem on the non-self-converse side. This gives a quadratic search cutoff, not the unknown atom polynomials or a classification of extremizers.

No finite census or new computational authentication is claimed. No core proof was edited. Source hash is recorded separately in nonselfconverse_source_pins.json. Only the review directory was written; no Git or remote writes were performed.
