# Ordinal-pair insertion and a two-parameter arithmetic route

Status: the formulas and strong-construction statements below are proved; the prime-coverage condition is open, and no infinite prime claim is made. This builds on the rooted block insertion formula of Session 150 and does not establish priority over the literature.

For a tournament A set h(A)=H(A), with H(empty)=1, and define the ordered two-path partition count

    P(A) = sum over S subset V(A) of H(A[S]) H(A[V(A)−S]).

The two paths are distinguished and may be empty. If C,D are tournaments and C→D denotes their ordinal sum, then

    H(C→D)=H(C)H(D),    P(C→D)=P(C)P(D).               (1)

Proof: a directed path in C→D has its C vertices before its D vertices. This gives a bijection, independently in each of the two distinguished paths, with a path pair of C and a path pair of D. The first identity is the analogous one-path bijection. For nonempty A, P(A) is even, by exchanging the two paths, since no pair is fixed when their disjoint union is nonempty. Furthermore P(A)≥2^|A|, with equality exactly for transitive A: each induced tournament has a path, and a nontransitive tournament contains a directed triangle, whose three cyclic Hamiltonian paths make its subset summand strictly larger than 1.

Let B be any tournament, v a vertex of B, q=H(B−v), and F the Session 150 number of ordered path pairs in B−v whose first path is empty or ends at an in-neighbor of v. Take nonempty A, form T=A→B, and adjoin a transitive k-vertex block X, k≥1, with v→X→T−v. Then

    H(T_X) = h(A)H(B) + (2−2^k)h(A)q + (2^k−1)P(A)F.  (2)

In particular, for k=1,

    H(T_x) = h(A)H(B) + P(A)F.                           (3)

Proof: H(T)=h(A)H(B), H(T−v)=h(A)q, and no Hamiltonian path of T starts at v, since A is nonempty and dominates B. We claim F(T,v)=P(A)F(B,v). Split each path of a pair on T−v into its A and B−v segments. There is no cross-boundary constraint other than their inherent order because A dominates B. If the first B segment is nonempty, the first whole path ends at its endpoint, so it must end at an in-neighbor of v in B. If the first B segment is empty, then any first A segment is admissible: if nonempty it ends at an A vertex, which beats v, and if empty the whole first path is empty. Thus admissible pairs are exactly arbitrary A path pairs times admissible B path pairs. Substitute these identities into the Session 150 rooted block formula h−s+q+(2^k−1)(F−q), yielding (2) and (3).

If B is strong, T_X is strong: every old vertex reaches v, v reaches X, and X dominates all old vertices except v. X is nonempty, so these routes join all vertices. More generally it is enough that every vertex of B reaches v.

## Concrete seed and arithmetic family

For C on vertices 0,1,2,3 use adjacency rows 0001, 1000, 1100, 0110 (equivalently lexicographic upper-triangle mask 4 with bit 1 meaning lower-label→higher-label). Direct path counting yields H(C)=5, P(C)=32, and at root v=1, F(C,1)=9. A short check for P is that the subset summands by subset size are 5,8,6,8,5. The eight F contributions, with first path's vertex set ordered as empty, {0}, {2}, {0,2}, {3}, {0,3}, {2,3}, {0,2,3}, are 3,0,1,0,1,1,1,2. C is strong (it contains the directed spanning cycle 0→3→2→1→0).

Put A=C→...→C (m copies), B=C, and use root 1 with k=1. Equations (1) and (3) construct a strong tournament with

    H = 5·5^m + 9·32^m,    m≥1.                          (4)

For m=1,2, these are 313 and 9341, respectively, both prime by trial division through their square roots. No assertion is made that (4) contains infinitely many primes or covers either required prime lane.

There is a serious seed barrier: if B is the directed triangle, then H(B)=F(B,v)=3 for every root v. Formula (3) gives H=3(h(A)+P(A)); hence for every nonempty A the output is composite and cannot realize any prime. In particular, arbitrary variability of the initial ordinal-sum seed A does not rescue this one-vertex triangle-root template for the prime lane. This statement does not cover larger k, other roots, or other gadgets.

## Exact conditional reduction

Suppose that for every odd prime p>199, there are nonempty tournaments A_p and D_p satisfying

    p = 5H(A_p)+9P(A_p),
    7p = 5H(D_p)+9P(D_p).                                (5)

Then the full H-spectrum conjecture follows. Indeed, (3) and the strong seed C provide strong witnesses for both numbers. For the remaining odd primes p≤199, the spectrum through 100001 certified in Session 150 contains both p (except the known exclusion 7) and 7p (except 21), since 7p≤1393. These attainments also provide strong witnesses: any tournament with count p has a strong component with count p; any tournament with count 7p has a strong component of count 7p because 7 is excluded (including p=7). The Session 150 two-prime-lane reduction then supplies all remaining allowed odd values via ordinal products; Rédei parity and the known exclusions 7 and 21 give the reverse inclusion. Condition (5) is an explicit sufficient arithmetic target, not an established fact or a claimed equivalence.

## Verification controls

Run `python3 check.py`. The script builds labelled tournaments independently, computes directed-path counts by subset/last-vertex dynamic programming, computes P and F from their definitions, explicitly forms all extensions, and checks (2) in 957 cases for all labelled A,B of orders 1–3, all roots, and k=1,2,3. This finite check supports the proof but does not replace it. The output also lists h,P for all such A. Additional finite examples above are checkable from the subset contribution tables; primality is elementary trial division. 

## Varying the root

For roots v=0,1,2,3 of C, direct enumeration gives F(C,v)=7,9,6,5, respectively. Allowing all roots in the construction with A=C^→m gives four families 5·5^m+c·32^m, with c in {5,6,7,9}. The c=5 family is always divisible by 5. The c=7 family is always divisible by 3, since 5≡32≡−1 (mod 3) and 5+7≡0 (mod 3). Both families consist of composite numbers for m≥1. The c=6 family gives 217=7·31 at m=1, 6269 at m=2, and 197233 at m=3; 31, 6269 and 197233 are prime by trial division. Thus root variation supplies both a specific seven-prime witness and additional prime examples, but neither establishes a prime-lane coverage theorem. Allowing all four roots in condition (5) gives a weaker sufficient condition: for each indicated target N, require nonempty A and c∈{6,7,9} with N=5H(A)+cP(A). The c=5 option can be omitted because the targets in (5) are coprime to 5.
