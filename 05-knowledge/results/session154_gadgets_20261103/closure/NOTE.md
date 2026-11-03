# Prime-lane obstruction to the transported family and adaptive two-branch closure

Status: rigorous consequences of the verified Session 153 transfer identities; not a proof of the global tournament-spectrum conjecture. Source: `golden_inverse_defects_20261103/TRANSPORT.md`, especially T1--T3, T5 and T7, at commit f6d5304849e50b29df6e368239e92fffa92e57e3.

## 1. Explicit matrix proof of growth

The actual seed U has state (Y,X,g)=(25,34,6), Hamilton count h=Y+g=31. The source-prefix symbols q2,b act on this state by

    q2: (Y,X,g) -> (11Y+4X,18Y+7X,g),
    b : (Y,X,g) -> ( 5Y+2X, 8Y+3X,-g).

For every prefix, Y>=25, X>=34 and |g|=6. These hold by entrywise positivity and induction. The exact differences are

    h(q2 w)-5h(w)=6Y+4X-4g > 0,
    h(b w)-5h(w) =2X-6g > 0.

Therefore h(wU)>=31*5^|w| for every prefix w in {q2,b}*, regardless of parity. This is an actual count estimate, with no eigenvalue approximation, norm cancellation, or parity assumption. For even b-parity, T3 identifies actual strong-tournament outputs 3h+2 and 3h+4, arising from literal suffix replacements U2 and U4. Their source words are wU2 and wU4 respectively.

For X>=93, if either positive output is <=X then length |w|<=floor(log_5(X/93)). At length l there are at most 2^l prefixes. Including both suffixes gives at most

    2 sum_{l=0}^{floor(log_5(X/93))} 2^l = O(X^alpha),
    alpha=log_5 2<1.

The estimate counts words with multiplicity: collisions and the even-parity restriction only reduce the number of distinct H values. The same bound holds for any fixed finite set of replacement outputs satisfying a uniform lower bound c*h(wU)-C with constants c>0,C; it does not extend automatically to arbitrary suffixes, an unbounded bank, or transformations which shrink H.

## 2. Both prime lanes are missed infinitely often

Chebyshev's classical lower bound gives pi(x)>=c x/log x for all sufficiently large x. Since X^alpha=o(X/log X), the number of primes p<=X exceeds the number of values in the family up to X for all sufficiently large X. Thus infinitely many primes are absent from this family. Likewise, there are pi(X/7) integers of the form 7p with p prime and 7p<=X, which also eventually outnumber the family's O(X^alpha) values. Thus infinitely many 7p are absent too. Removing the finitely many exceptional p specified in the two-prime-lane reduction changes nothing.

This strengthens the warning that the family has zero integer density: zero density alone would not exclude a set containing all primes. It says nothing about which missing lane values other strong tournaments realize.

## 3. Adaptive two-branch no-go, independent of finite-state assumptions

Let a rooted tree of candidate constructions have finitely many roots, each labeled by a positive integer H, and at most r children per node. Assume every child has H_child>=a H_parent, with a>r>1. The rules may depend on the full history and need not be regular or finite-state. If H<=X, depth is at most log_a(X/H_min), so the number of constructed values at most X is O(X^(log_a r)). By the same Chebyshev argument, neither all but finitely many primes nor all but finitely many 7p can occur.

In particular, start from any finite set of seeds and at each state allow at most two realized updates selected adaptively from H->3H, 3H+2, 3H+4. Every child is >=3H. There are O(X^(log_3 2)) reachable values <=X and the two lanes are both missed infinitely often. This applies even if multiple distinct guarded realizations have the same H, because the count is of construction paths, provided the branching bound is per full state. It does not rule out three branches, an infinite seed bank, reverse moves, or operations with insufficient uniform growth. It also does not claim that any of those three updates is universally realizable.

## 4. Independent concrete limitation of a different real tournament family

For T(w) on x and a transitively oriented chain 1<...<n, with x->i exactly for bit b_i=1, strongness is equivalent to first bit 1 and last bit 0. A path-counting formula from Session 25 gives a positive term 2^(n-2) from the first 1 and last 0, plus two endpoint paths. Thus H>=2^(n-2)+2. If H=13 then n<=5. Direct exact recursion from (H,A,p)=(1,0,0), with append-1 (H,2A+1,1) and append-0 (H+A+p,2A,0), yields these attainable counts for n=2,3,4,5 respectively:

    {3}; {5}; {9,11}; {15,17,19,23}.

Consequently this entire strong feedback-vertex-one family does not realize prime 13. It is an explicit guard failure, despite its proved complete coverage of parity-admissible residue classes modulo every fixed modulus. No claim is made that 13 is globally absent; it is known realized by other tournaments.
