# Two moving-guard choices form a positive three-state algebra

**Status.** Exact physical transfer, a strong recurrence family, and an all-depth power-of-two modular realization theorem. These do not supply a complete ternary bank or exact odd-spectrum classification. No literature-priority claim is made.

## Sources

This extends [S163's moving-guard lemma](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/strong_returns/MOVING_GUARD.md) and its [unconditional first-step refinement](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session163_ternary_triples_20261103/strong_returns/FIRST_STEP.md). The latter already supplies the scalar correction 3H+2d, propagation to the a guard, and the alternative b-guard defect. The result here packages **both choices** into a closed, positive state and iterates arbitrary words. The characteristic polynomial that appears below is already present in [S159's overlapping-cycle transfer](https://github.com/eliottcassidy2000/math/blob/d9d9f3e483c54d83ebd569b91a85b7d015923cfe/05-knowledge/results/session159_constructions_20261103/transfer_knots/NOTE.md), with different initial data.

## 1. A positive three-state transfer

Let q have unique in-neighbor r in X. Define

    h = H(X),
    u = # Hamilton paths of X not starting at q,
    d = # Hamilton paths of X−q starting at r.

Add a,b with q→a→b→q and a,b→X−{q}. Both a and b have indegree one in the enlarged graph, with respective unique in-neighbors q and a. Use A to mean this attachment with new guard a, and B to mean the same attachment with new guard b. Then

    A(h,u,d) = (3h+2d, 2h+2d, 0),
    B(h,u,d) = (3h+2d, 2h+d, u).                              (1)

Equivalently, on column vectors,

        [3 0 2]             [3 0 2]
    A = [2 0 2],        B = [2 0 1].                         (2)
        [0 0 0]             [0 1 0]

**Proof.** S163 gives total H=3h+2d. The paths starting at a are precisely a,b followed by an X-path, so there are h. The paths starting at b are counted by h+d in S163's path decomposition. Subtract these counts from the total to get the two u coordinates. After deleting a, b is a source, so no Hamilton path starts at a's unique in-neighbor q; the new d is zero. After deleting b, a must be followed by an X-path whose initial vertex is not q, giving precisely u for b's defect. ∎

All entries in (2) are nonnegative. The subtraction in the intermediate description d_b=H−H_start disappears once u is made a state variable. Every finite word in A,B therefore has a literal guarded-tournament realization. If X is strong, all descendants remain strong, because the old strong graph and q→a→b→q form one strong component. Each letter costs exactly two vertices.

The state is exact for these two choices; it is not asserted to determine all future tournament attachments. Arbitrary marked-arc weights are also not automatically preserved by this bijection, because inserting after q can remove an old outgoing arc of q.

## 2. A rational strong family with an existing denominator

Take X_0=C3 with any vertex as guard. Its state is

    v_0=(3,2,0).

Apply B repeatedly and let a_n=H(B^n X_0). The matrix B has determinant 1 and characteristic polynomial

    λ³−3λ²−λ−1.

Consequently

    a_0=3, a_1=9, a_2=31,
    a_(n+3)=3a_(n+2)+a_(n+1)+a_n,
    sum_(n>=0) a_n z^n = (3+z²)/(1−3z−z²−z³).                 (3)

These are strong tournaments on 2n+3 vertices. The beginning is

    3, 9, 31, 105, 355, 1201, 4063, 13745, 46499, 157305, …

This denominator matches S159's overlapping reversed-interval family, whose numerator is (1+z)² and whose initial counts are 1,5,17. Equality of scalar denominators alone would not identify graphs or full interfaces. Section 5 goes further by factoring the guard step into one-vertex moves and proving an explicit graph isomorphism for the intervening even orders.

The real dominant root is approximately 3.3829757679. The polynomial has one real root (its discriminant is −176), and the other two roots have modulus 1/sqrt(λ)<1, so its dominant root is a cubic Pisot number. This observation concerns this recurrence's growth, not exact integer coverage.

The first 31 is 7 mod 12, and 46499 is 11 mod 12. Thus the family escapes both congruence classes forbidden to the entire ordinal-sum/triangle profile grammar of [PROFILE_ALGEBRA.md](PROFILE_ALGEBRA.md).

## 3. Exact multiplicative checkpoints and all odd 2-adic residues

Call a state canonical when it equals t(3,2,0). For every positive integer t,

    A[t(3,2,0)] = 3t(3,2,0),
    A(A(B[t(3,2,0)])) = 31t(3,2,0).                           (4)

In temporal notation the second word is BAA. These are literal source returns with the current physical guard, not just arithmetic multipliers. Starting from C3 and composing A and BAA therefore gives strong tournaments with

    H = 3^(a+1) 31^b,       order = 3+2a+6b,    a,b>=0.       (5)

**Theorem.** For every k>=1 and every odd residue c modulo 2^k, one of the strong sources (5) has H congruent to c modulo 2^k.

**Proof.** For k>=3, 93=3·31 satisfies v_2(93−1)=2. Hence 93 has order 2^(k−2) modulo 2^k and generates all units congruent to 1 mod 4. Multiplication by 3 supplies the other coset. Therefore 3 and 31 generate the full unit group modulo 2^k. Since it is finite, every inverse is also a nonnegative power, so the monoid generated by 3 and 31 equals that group. Multiplication by the initial factor 3 is a bijection on the group. Cases k=1,2 are immediate. ∎

For completeness, the order claim follows by induction from

    v_2(93^(2^j)−1)=j+2.

The base is 93−1=92. At every doubling, the other factor 93^(2^j)+1 has 2-adic valuation 1, so the valuation increases by exactly one.

This is modular universality within a very restricted smooth-number family. It does not say that every odd integer is one of its Hamilton counts, or that exact +2 and +4 branches have become available.

## 4. Fixed observer queries and finite modular searches

Let p≠q be an observer port, and let X_0=X, X_1=X^a_p, X_2=X^aa_p be S163's successive queries. Form the three raw states (h_i,u_i,d_i), with the same q and its same unique in-neighbor r. Both attachment choices commute with either query: every attachment vertex beats every query vertex, while the query vertices beat only their designated ports. Therefore the same A or B acts simultaneously on all three raw states.

This is an exact nine-coordinate nonnegative physical interface for arbitrary guard words. If the usual observer is h=h_0, x=(h_1−h_0)/2, y=(h_2−h_1)/2, its value can be recovered after each word. This does not reduce the nine coordinates to the original three, nor prove the required affine ternary offsets are reachable.

Reducing (1) modulo any fixed M gives a complete finite search with at most M³ states for this guarded language and seed. The following saturation results are finite controls only:

| M | Reachable states | Hamilton residues |
|---|---:|---|
| 3 | 27 | all |
| 5 | 125 | all |
| 7 | 217 | all |
| 9 | 573 | all |
| 12 | 216 | all odd |
| 15 | 3375 | all |
| 16 | 512 | all odd |
| 21 | 5859 | all |
| 25 | 15625 | all |
| 27 | 13911 | all |
| 32 | 4096 | all odd |

A useful next target is an odd-modulus reachability theorem or a certified exceptional modulus. The displayed finite controls do not prove either. For ternary-bank work, the stronger target is a physical return for the nine-coordinate observer state with the desired affine correction; scalar modular coverage is insufficient.

## Reproduction

Run `python3 check.py` and `python3 -O check.py`. Independent subset/endpoint DP verifies (1) on all 110 qualifying guarded labeled tournaments of orders 1 through 4, including non-strong and boundary cases. It verifies strength and actual guard states of B^n(C3) for n=0,…,6 (through 15 vertices), the multiplier identities, the stated finite modular saturations, and unit-generation controls through modulus 4096. All checks remain active under `-O`. The general results follow from the path bijection, matrix identities, and elementary unit-group proof above. Sections 5 and 6 add further checks described at the end.

## 5. One-vertex factorization, a tribonacci root, and the literal S159 bridge

The parallel one-vertex guard calculation in this session refines the two-vertex language. Add only a with q→a and a→X−q. Either retain q (call this D) or mark a as new guard (call this C). Both choices are valid unique-in-neighbor guards. For a strong source the output remains strong: the cycle r→q→a→r joins a to the old strong graph.

Use the nonnegative coordinates (s,u,d), where s counts Hamilton paths starting at the guard, u counts the others, and d is as above. A direct decomposition gives

    C(s,u,d)=(u,s+u+d,s),
    D(s,u,d)=(s+d,2u,0).                                     (6)

Indeed, paths starting at a correspond to the u old paths not starting at q. Paths in which a is not first either insert it immediately after q into an old Hamilton path, or use the d exceptional q,r joins. Thus total H=s+2u+d. For C, the new start count is u, and deleting a gives defect s. For D, paths starting at q count H(X−q)=s+d; deleting q makes a a source, so the defect is zero. This proves (6).

The fixed-port query commutation in Section 4 also holds for each one-vertex move. Consequently C and D act simultaneously on the three raw endpoint triples, giving the same exact nine-coordinate observer lift at one-vertex resolution.

Processing letters from left to right, CC is exactly the two-vertex B, and CD is exactly A, including their physical guards. The matrix of C in (s,u,d) coordinates is

    [0 1 0]
    [1 1 1],
    [1 0 0]

with characteristic polynomial λ³−λ²−λ−1. Thus a_n in (3) is the even subsequence of a tribonacci family. If t_j=H(C^j(C3)), then

    t_0=3, t_1=5, t_2=9,
    t_(j+3)=t_(j+2)+t_(j+1)+t_j,
    sum_(j>=0) t_j z^j=(3+2z+z²)/(1−z−z²−z³).                (7)

More strongly, the odd subsequence is **literally isomorphic** to S159's overlapping reversed-interval family. This upgrades the scalar comparison in Section 2.

Write the initial triangle as g_0→s→r→g_0, and let g_1,…,g_(2m−1) be successive new guards under C. Relabel the tournament by the following vertex order:

    (g_(2m−2),g_(2m−1), g_(2m−4),g_(2m−3), …, g_0,g_1, s,r).

Every edge is forward in this order except precisely

    (2i+3)→(2i),       i=0,…,m−1.                             (8)

To check this, each displayed guard pair has its internal edge forward. A later-created guard beats every older vertex except its immediate predecessor in the guard chain. Between consecutive displayed pairs, that sole exceptional predecessor edge runs from position 2i+3 to position 2i; nonconsecutive pairs have only forward edges. The final pair (s,r) has s→r, and its only exceptional edge to the preceding pair is r→g_0, again of the form (8). All other edges involving s,r point forward. This verifies every edge, proving the isomorphism.

Consequently S159's H_m equals t_(2m−1) for m>=1. To include m=0, start instead from a directed edge with its sink as guard, state (s,u,d)=(0,1,1). One C move gives the canonical C3, so S159's entire family is the even subsequence C^(2m) of this edge seed. Its cycle-packing transfer and the moving-guard endpoint transfer describe the same underlying graphs on even orders. The odd-order cap t_(2m) gives the B family of Section 2. No independent-new-family claim is appropriate for the even-order subsequence.

## 6. A smaller canonical multiplicative bank

In the one-vertex coordinates the canonical state is t(1,2,0). For n>=0,

    D^n[t(1,2,0)] = t(1,2^(n+1),0),
    D(C(s,u,0)) = (s+u)(1,2,0).

Therefore the temporal word D^n C D returns to a canonical source with multiplier 1+2^(n+1), using exactly n+2 new vertices. In particular CD multiplies by 3 and DCD multiplies by 5. This yields the more economical strong sources

    H=3^(a+1)5^b,      order=3+2a+3b,      a,b>=0.             (9)

The same power-of-two modular theorem holds with this smaller bank. For k>=3, 5 generates the subgroup of units congruent to 1 mod 4, because v_2(5^(2^j)−1)=j+2, and 3 supplies the other coset. Hence 3 and 5 generate all odd units modulo 2^k. As before, finite-group inverses are nonnegative powers, and the leading factor 3 does not reduce coverage. This improves the vertex bound from Section 3; the earlier 3/31 return remains valid inside the stricter two-vertex A/B language.

The checker additionally verifies both one-vertex transfers on the same 110 qualifying sources, checks the literal CC factorization, verifies the explicit isomorphism permutation in (8) through m=10 (22 vertices), and controls the general multiplier identity for t=1,…,29 and n=0,…,7. Formulae (6), (8), and the canonical identity prove the corresponding all-depth assertions.

## Full-profile extension

[POINTED_PROFILE_ALGEBRA.md](POINTED_PROFILE_ALGEBRA.md) lifts this scalar calculus to every ordered slot count and to positive unordered-cover polynomials. It gives a closed endpoint module combining C/D, ordinal grafting, and triangle resets, with exact finite modular quotients.
