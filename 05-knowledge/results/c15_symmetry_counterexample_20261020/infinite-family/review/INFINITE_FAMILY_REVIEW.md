# Independent proof review: infinitely many mixed-prime symmetry deficits

**Verdict: PASS for the strict inequality stated below.** This review does not claim an exact deficit, a smallest example, or a full classification. It supplements the independently checked n=10 C15 counterexample.

Let a<b<2a be distinct primes, both congruent to 3 modulo 4. Put k=b-a<a and n=2b. On disjoint blocks A,B,F of sizes a,b,k let

    H=H_a times H_b,
    H_p={x -> ux+v on F_p: u a nonzero square, v arbitrary},

and let H fix every vertex of F. Then H has odd order and

    tau_E(H)<tau_T(H).

There are infinitely many such pairs (a,b), hence counterexamples at infinitely many vertex orders.

## 1. Internal actions and the normalizer

For p=3 mod4, H_p has order p(p-1)/2, which is odd. It is transitive on unordered pairs: translating one endpoint to zero leaves a difference up to sign, and exactly one of d,-d is a square. Consequently its invariant internal graphs are empty or complete, and its invariant tournaments are the Paley tournament and its converse. A nonsquare multiplier interchanges these two tournaments, so they form one orbit under AGL(1,p).

The translation subgroup C_p is the unique Sylow p-subgroup of H_p: its number divides (p-1)/2 and is 1 modp. It is characteristic. Thus any permutation normalizing H_p normalizes C_p and is affine. Conversely all affine permutations normalize the square-multiplier subgroup, so

    N_Sp(H_p)=AGL(1,p).

The H-orbits are A, B and the k singleton vertices. Their two nontrivial sizes are distinct, so the full normalizer is exactly

    N=N_Sn(H)=AGL(1,a) times AGL(1,b) times S_F.

No normalizer element exchanges A or B with a union of smaller H-orbits.

## 2. Why the internal prime-order tournaments are prime

A nontrivial module M of a regular tournament induces a regular tournament: all vertices of M have the same outside out-neighbors, so equal global outdegrees imply equal internal outdegrees. Every nontrivial regular tournament is strongly connected. For example, if it had distinct source and sink strongly connected components, regular indegree forces the source and regular outdegree forces the sink each to have more than half the vertices, a contradiction.

A strongly connected tournament module is strong. If two modules properly overlap, their intersection and difference in either module form a nontrivial uniformly directed cut, contradicting strong connectivity. This was proved explicitly in the n=10 review.

Now a vertex-transitive tournament of prime order p has a transitive C_p-action in the present case (translations). Every translate of M is a strong module of the same size. Equal-sized strong modules are equal or disjoint. The translates cover all vertices, so they form a partition into sets of size |M|. Thus |M| divides p, forcing |M|=1 or p. The Paley tournament has no proper nontrivial module.

This proof uses the prime order and the translation action, not an unproved assertion about the full automorphism group of a Paley tournament.

## 3. Unique large modules remove tournament multiplicity

Every H-fixed tournament is a substitution with internal prime strong tournaments on A and B, arbitrary orientation between the two blocks and the k fixed vertices, and an arbitrary tournament on F. Both A and B are strong modules in the full tournament.

Any other prime strong module of size a is nested with or disjoint from A and B. It cannot be a proper nontrivial module inside prime B, cannot contain B because a<b, and cannot fit inside F because k<a. Thus A is unique. A second prime strong module of size b would be disjoint from B and therefore equal to its complement A union F (whose size is b); but that tournament contains the proper nontrivial module A, so is not prime. Hence B is unique as well.

Every full isomorphism between two H-fixed tournaments therefore preserves A, B and F. This **does not require the given isomorphism itself to normalize H**. Replace its restriction on A by identity or a chosen nonsquare multiplier, according as the two internal Paley orientations agree or are opposite; do the same on B, retaining its original permutation of F. These replacements preserve all external uniform orientations and produce an isomorphism in N.

It follows that each actual H-symmetric tournament type meets Fix_T(H) in exactly one N-orbit:

    tau_T(H)=#(Fix_T(H)/N).                              (1)

The last replacement step is necessary; unique modules alone identify blocks but do not identify normalizers of their internal groups.

## 4. The one-sided graph count is sufficient

Apply inherited THM-4721 Theorem R with odd P=H and Gamma=N. It gives

    #(Fix_T(H)/N)=#({N-even graphs in Fix_G(H)}/N).       (2)

A globally sign-even H-fixed graph is N-even. Also every globally even H-symmetric graph type has an H-fixed representative, so there is a surjection from N-orbits of globally even H-fixed graphs onto those actual types. Therefore

    tau_E(H) <= #({globally even H-fixed graphs}/N).     (3)

Consider K_b,b with parts B and A union F. It is H-fixed. Every element of N preserves these two sides separately; orienting every edge from the first side to the second shows that all these normalizer automorphisms have sign +1. Thus it belongs to the right side of (2).

But swapping the two sides pairs the vertices into b adjacent pairs. Its sign is (-1)^b=-1, because b is odd. Hence K_b,b is globally odd and is absent from (3).

The globally even subset is N-invariant and excludes at least this one whole N-orbit. Combining (1)–(3) proves

    tau_E(H) <= tau_T(H)-1.

No injectivity for graph types is needed. Additional graph fusion or additional normalizer-even/global-odd graphs can only increase the deficit. The proof therefore establishes a strict inequality, not deficit exactly one.

## 5. Infinitely many admissible prime pairs

The authoritative formula [NIST DLMF 27.11.9](https://dlmf.nist.gov/27.11.E9), with h=3,k=4, gives divergence of the reciprocal sum over primes congruent to 3 modulo 4. I independently opened the formula and checked the coprimality hypothesis and the coefficient 1/phi(4).

List these primes increasingly as q_1,q_2,.... If only finitely many adjacent pairs satisfied q_(j+1)<2q_j, then the sequence would eventually at least double at each step. Its reciprocal tail would then be bounded by a convergent geometric series, contradicting divergence. Thus infinitely many adjacent pairs satisfy a<b<2a. Their larger members b are unbounded and distinct, giving infinitely many orders n=2b.

This needs more than mere infinitude of primes in the progression. Reciprocal divergence supplies exactly the required multiplicative-gap conclusion. No explicit bound on the first such pair after a given threshold is claimed.

## 6. Finite examples and exact scope

The first pair in the displayed construction is (a,b)=(7,11), with k=4,n=22, followed for example by (11,19) and (19,23). These are examples, not a minimal-counterexample claim; the separate cyclic construction already gives the n=10 counterexample.

For k fixed vertices, normalization of the two internal Paley orientations leaves the count of tournaments on two individually distinguished anchors and k interchangeable vertices. Its elementary Burnside formula is

    tau_T(H)=(1/k!) sum over odd-order sigma in S_k
                  2^(1+2 c(sigma)+c_2(sigma)),

where c counts cycles on F and c_2 counts orbits on its unordered pairs. At k=4 the identity contributes 2^15 and the eight 3-cycles contribute 2^7 each, giving tau_T=1408. The theorem then gives tau_E<=1407 for the (7,11) action; it does not determine equality.

The family proof relies on the inherited restricted twisted Burnside theorem, which has already been independently audited in the session. The internal prime-module lemma and the normalizer-to-full-type argument are supplied here, and the global odd witness is explicit.
