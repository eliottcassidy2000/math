# Bounded merge slack has a finite core: inverse unit-amplifier certificates

Session 197, 2026-11-05. **All-order rooted-profile normal form; exact finite-core classifications; all-order binary-cloning obstructions at H=13 and H=15. The full {7,21} tournament Hamilton-path spectrum remains open.** No literature-priority claim is made.

## 1. Sources, notation, and the sign convention

Named starting sources are the [S195 guide](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/00-navigation/S195.md) and [S196 guide](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/00-navigation/S196.md). In particular:

- [S195 zero-slack theorem](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/zero-slack/PROOF.md) gives the positive two-endpoint extension lemma and the boundary-cover formula for merge slack.
- Its [universal refinement](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/zero-slack/UNIVERSAL_BOUND.md) proves a finite-core/dyadic-padding classification **when slack is zero**. The result below extends this to **every bounded slack**, without a finite-gadget-bank hypothesis.
- [S196 transitive cloning](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/spectral-moves/PROOF.md) gives the exact star synthesis and classifies binary cloning geometrically. [S196 unit determinants](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/clone-determinants/PROOF.md) treats deletion count one and prime amplitudes.

For a rooted tournament (T,r), put

    h=H(T), e=ends_r(T), d=H(T-r),
    R=number of ordered spanning covers (P,Q), Q nonempty starting at r,
    f=R+e-h-d, U=R-d, beta=d-f.

Empty paths are allowed in P, and H(empty)=1. The standard physical inequalities are e<=d<=h<=R, f>=0. In particular

    U=h-e+f >= f.                                    (1)

The stationary-root operation K adds a vertex z with z->r and every other old vertex ->z. It has

    K(h,e,R,d)=(h+R-e,d,2R-d,d).

Consequently, for m>=1,

    H(K^m(T,r))=U*2^m+beta,
    e_m=d, R_m=d+U*2^m, d_m=d, f_m=f.                (2)

The new vertices form a transitive clone class. A **unit offset** means beta=+1 or -1 in the exponential formula. The iterative scalar operation is H -> 2H-beta. Thus beta=-1 means the +1 binary step; beta=+1 means the -1 binary step. This convention matters at H=13.

## 2. Two elementary transport facts

Set S=T-r and O=N_T^+(r). For x in S write R_x(S) for the same rooted cover count in S. Directly remove the initial r from the second path of an R-cover. A singleton second path contributes d; otherwise its next vertex is a uniquely determined x in O. Therefore

    U = sum_(x in O) R_x(S).                         (3)

Also e=d-sum_(x in O) ends_x(S).

For distinct x,y in S let B_yx(S) count ordered nonempty spanning two-path covers whose first path ends at y and whose second starts at x. The S195 boundary decomposition gives

    f=sum_(x,y in O, x!=y) B_yx(S).                  (4)

The S195 positive two-endpoint lemma states B_xy+B_yx >= H(S-{x,y}) >=1. Summing over unordered pairs gives the useful quantitative strengthening

    f >= binom(|O|,2).                              (5)

In particular, with

    k(f)=floor((1+sqrt(1+8f))/2),

we have |O|<=k(f). When f=0 this says |O|<=1, agreeing with the exact inherited zero-slack theorem. The argument would also handle f=1, although physical f is even: the singleton-birth tournament has odd count h+R-e, and h,d are odd, so f=R+e-h-d is even.

Write the strong components of S in their unique order C1 -> ... -> Cq. For x in Cj,

    R_x(S)=(product_(i<j) H(Ci)) R_x(Cj)
             (product_(i>j) P_Ci(2)),                (6)

where P_C(2)=sum_A H(C[A])H(C-A) counts ordered two-bin path covers. A singleton has H=1 and P(2)=2. Formula (6) is inherited from the S195 universal refinement, and follows by cutting the two paths at the strong-component boundaries.

Every nontrivial strong component C has H(C)>=|C| and H(C)>=3: every vertex is the endpoint of some Hamilton path. If the nontrivial component counts are d1,...,ds, their product is d. The elementary inequality sum di<=product di for di>=3 shows that the total number of vertices in all nontrivial components of S is at most d. If there are none, that total is zero.

## 3. Universal finite-core normal form at bounded nonzero slack

Define ell(f)=floor(log_2(max(1,f))).

**Theorem 1 (exact profile normalization).** Every rooted tournament (T,r) has a rooted tournament (C,c), and an integer t>=0, with

    |C| <= 1+d+k(f)+ell(f),
    H(C-c)=d, f(C,c)=f,
    profile(T,r)=profile(K^t(C,c)).                  (7)

The normalization is an equality of the four rooted counts, not an assertion that the original tournament is isomorphic to a padded core. C is an induced subtournament of T; some neutral source-side singleton components are discarded, and certain remaining singleton components are relocated to a terminal K-chain.

**Proof.** If O is empty, r is a universal sink. Delete every singleton strong component of S and keep r and all nontrivial components. The four counts all equal d before and after this deletion, f=U=0, and t=0 works.

Now suppose O is nonempty. Let C_a and C_b be the first and last strong components of S meeting O. A singleton component meeting O is itself an out-neighbor and must be retained. There is no ambiguity about a singleton lying in the same component as an out-neighbor: if a singleton component meets O, its unique vertex belongs to O.

Keep all nontrivial components, every singleton in O, and every singleton strictly between C_a and C_b. Delete all other singleton components. Let C be the induced tournament on these retained vertices and r. Let t be the number of deleted singletons after C_b. Those deleted before C_a are neutral and are not counted in t.

The deleted prefix singletons precede every possible x,y in O. In (3) and (6) they contribute a factor H(K1)=1, not a two-bin factor. In every boundary cover (4) they must lie in its first path: the second path starts at an out-neighbor in a later component. They therefore contribute one and do not alter f. Endpoint counts at out-neighbors are unchanged by this prefix deletion as well.

Every deleted suffix singleton follows every x in O, so by (3) and (6) it contributes the same factor two to U. Its position among later nontrivial components is immaterial because these factors multiply. In every cover (4) it must lie in the second path: the first path ends at an earlier out-neighbor. Hence it contributes one to f, not two. Deletion count d is unchanged throughout. If t>0, the sink component of S is disjoint from O, so e=d. Appending the t singletons at the end as a K-chain has exactly the same U, d, f and e. If t=0, only prefix deletions occurred, and e is unchanged. Since h=U+e-f and R=d+U, all four profile coordinates agree.

It remains to bound the retained singleton vertices. At most |O|<=k(f) are in O. Let u be the number of retained singleton vertices outside O, all strictly between C_a and C_b. If a=b, then u=0. If a<b, choose x in O intersect C_a and y in O intersect C_b. Choose a Hamilton path of C_a starting at x, and a Hamilton path of C_b ending at y. Such paths exist by strongness, including singleton components. Put all components before C_a into the first path, all components after C_b into the second. Put C_a entirely in the second path and C_b entirely in the first. Each of the u indicated singleton vertices can independently choose either path. Put each remaining intermediate component entirely in the first path. Concatenating within each slot in component order gives 2^u distinct covers counted by B_yx(S). Thus f>=2^u. In particular u<=ell(f). This argument is needed only when a<b, which implies |O|>=2 and f>0; it never takes log 0.

The nontrivial component vertices number at most d, the retained O singletons at most k(f), and the other retained singletons at most ell(f). Add the root to obtain (7). QED.

**Corollary 1.1.** For fixed integers D,F, all profiles with d<=D and f<=F belong to finitely many K-rays from actual rooted cores. Their scalar counts at most X number O_(D,F)(1+log X).

Indeed the core-order bound is uniform. There are finitely many labelled cores of that bounded order. Each contributes its t=0 value and, for t>=1, the dyadic family d-f+U*2^t. A core with U=0 contributes a constant only. This recovers the S195 zero-slack universal sparsity result when F=0 and extends it to arbitrary fixed F.

## 4. An inverse-target consequence for unit offsets

For singleton-partner star synthesis, equivalently (2), beta=±1 forces f=d-beta. Therefore fixing d<=D automatically bounds f<=D+1. This implication is **not asserted** for arbitrary two-piece synthesis, where beta=d*H(Y)-f*H(Y-c) and the partner Y may vary.

**Corollary 2 (bounded-deletion unit rays).** All singleton-partner unit-offset amplifiers with d<=D, allowing arbitrary rooted source tournaments at every order, lie on finitely many dyadic rays. Their resulting Hamilton counts at most X number O_D(log X). Infinitely many prime targets, and infinitely many targets 7p with p prime, are absent from this restricted synthesis class.

For the last statement, the primes cannot have counting function O(log X). Otherwise summing reciprocal primes over dyadic intervals would converge, whereas the elementary Euler-product argument from the divergent harmonic series proves their reciprocal sum diverges. The same applies to 7 times those primes. This does not assert a missing tournament count.

There is also a sharper finite inverse certificate for odd prime amplitudes.

**Endpoint insertion bound.** If U>0, the root has at least one out-neighbor, by (3). Every Hamilton path ending at r has the form P,r, where P is a Hamilton path of S ending at an in-neighbor. Remove r and reinsert it immediately before the first out-neighbor along P. This is legal: before that first out-neighbor every vertex beats r. The new path does not end at r, and deleting r recovers P, so this is an injection. Therefore e<=h-e. As h is odd,

    h <= 2(h-e)-1 = 2(U-f)-1.

For p=U an odd prime and beta=±1, f=d-beta and d<=h consequently give

    3d <= 2p+2beta-1.                               (8)

This is stronger than the immediate bound d-beta<=p from (1). It holds for arbitrary positive U with the identical expression replacing p, whether prime or not.

Moreover the K-tail in Theorem 1 must have length zero, because a positive tail makes U even. Thus such a prime amplitude is represented by an actual source of order at most

    1+d+k(d-beta)+ell(d-beta).                       (9)

For a fixed d one obtains a crude explicit size restriction on p. Let n=d+k(d-beta)+ell(d-beta). There are at most n! rooted two-path covers starting at a specified vertex of an n-vertex tournament: choose the length of the first path and an ordering of the other n-1 vertices. Equation (3) then gives

    p <= k(d-beta) * n!.                            (10)

In particular, unbounded odd prime amplitudes require unbounded deletion data. This is a rigorous finite inverse reduction, not an unproved claim that every surviving core target exists.

## 5. Exact all-order classifications at deletion counts 1, 3 and 5

The following table lists all possible positive U=R-d in a unit source. Each entry means the displayed finite set multiplied by {2^j:j>=0}.

| d | beta=+1, f=d-1 | beta=-1, f=d+1 |
|---|---|---|
| 1 | {1} | {3,5} |
| 3 | {6} | {9,15} |
| 5 | {11,12,16} | {15,39} |

The d=1 row is inherited from S196. The d=3 and d=5 rows follow from Theorem 1 and the complete finite enumeration specified below. They are all-order statements, not a vertex-order census standing in for a theorem.

For d=1, S is transitive. For d=3, S has exactly one nontrivial component, a directed triangle. For d=5, S has exactly one nontrivial strong component of count 5. Such a component has at most five vertices. Exhaustive labelled enumeration through order five, independently counted, proves that its unique isomorphism type is the strong four-vertex tournament with rows [8,1,3,6].

Accordingly a normalized deletion core is a transitive sequence of singletons, then this one fixed component (if any), then a transitive sequence of singletons. Enumerate at most k(f)+ell(f) singleton positions, every position of the fixed component, and every root incidence word. Retain only the reduced cores in Theorem 1 and those with exactly the required f. This is an exhaustive finite list by (7), not a guessed bound. Its positive U sets are precisely those in the table. Every retained core is an actual graph, and K realizes all its dyadic multiples, proving both necessity and sufficiency. The verifier saves representative rows and profiles for every displayed entry.

The corresponding actual K-families, m>=1, are therefore

    d=3: 6*2^m+1; 9*2^m-1; 15*2^m-1,
    d=5: 11*2^m+1; 12*2^m+1; 16*2^m+1;
         15*2^m-1; 39*2^m-1,

with optional extra powers of two absorbed into m. These are physical constructions, not primality claims. Some scalar rays overlap earlier constructions; the new assertion is their completeness under the specified deletion data.

## 6. Global obstructions to returning clone buttons at H=13 and H=15

Use the complete physical star synthesis from S196. Every positive binary cloning amplifier, after possibly taking a converse, has

    alpha=U_X L_Y,
    beta=d_X H(Y)-f_X H(Y-c),                       (11)

where all pieces are actual rooted tournaments, U_X=R_X-d_X>0, and L_Y is the root-ending two-cover response. If alpha=p is an odd prime and beta=±1, S196's prime-amplitude reduction gives

    U_X=p, L_Y=H(Y)=H(Y-c)=1,
    beta=d_X-f_X.                                  (12)

For completeness, the excluded factorization U_X=1 would force f_X=0 and h_X-e_X=1. Zero-slack geometry then gives (h_X,e_X,R_X,d_X)=(1,0,2,1). The unit determinant forces H(Y)=1; Y is transitive, and L_Y is a power of two, not an odd prime. Thus L_Y=1, which forces Y transitive with its root at the source and establishes (12).

**Theorem 3.** There is no tournament T and vertex r for which either of the following holds:

    H(T[r:m])=7*2^m-1 for every m>=1;
    H(T[r:m])=7*2^m+1 for every m>=1.                (13)

Here T[r:m] replaces r by a transitive m-chain with the same exterior incidences.

**Proof.** Such an amplifier has prime alpha=7 and beta=-1 or +1. By (12) it has a singleton-partner source with U=7 and f=d-beta. For beta=-1, (8) gives 3d<=11, hence d is 1 or3. For beta=+1, it gives 3d<=15, hence d is 1,3, or5. None of the corresponding exact table rows contains U=7. This is a contradiction in both cases. QED.

Equivalently, **no H=13 tournament has a transitive-cloning root that indefinitely returns H -> 2H+1, and no H=15 tournament has one that indefinitely returns H -> 2H-1.** These are global obstructions within the fully classified binary-cloning architecture, with no graph-order bound and no fixed-template-bank assumption. They do not exclude an isolated clone step in a higher-run-degree family, another graph operation, or H=13 and H=15 themselves. No global exclusion of the count seven is used in this argument.

The sign control is explicit. The known H13 tournament

    E=[2,12,9,17,7], root 0

has alpha=6, beta=+1, with counts 13,25,49,97,... . Its returning operation is H -> 2H-1. This is consistent with Theorem 3, whose forbidden family would be 13,27,55,111,... . The verifier recomputes the former directly from cloned adjacency matrices.

## 7. What this does and does not settle

The exact normal form supplies a finite inverse search for bounded deletion/slack, and the small-deletion table gives new complete physical lanes and two specific all-order failures of a putative binary induction. It strengthens an architecture obstruction, not the full {7,21} conjecture. The broad prime and 7p problem still needs an operation that handles unbounded deletion/slack, a different amplitude factorization with varying partners, or a different construction entirely. Coprimality, residues, a large finite covered interval, and these sparse families do not provide that missing realization theorem.

## 8. Reproduction

Run `python3 verify.py > receipt.json` and `python3 -O verify.py > receipt.optimized.json`; compare stdout. All checks use explicit exceptions. The verifier checks the normalization and quantitative bounds on every labelled tournament through order five and deterministic larger controls, independently recomputes small Hamilton counts by permutations, classifies all possible H3/H5 nontrivial strong cores, enumerates the exact normalized unit-source families, constructs their K extensions, and checks the H13 sign control. A finite check supports the displayed base table; the universal completeness is supplied by Theorem 1 and the strong-core bound.
