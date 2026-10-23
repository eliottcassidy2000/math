# Two exterior vertices buy one additional tournament path-cover coordinate

Status: PROVED elementary continuation theorem; exact finite tests are separately reported in exact.json. No literature-priority claim. This develops, rather than replaces, the existing substitution kernel in THM-3121.

## Inheritance and the genuinely new target

The requested run-word/path-cover algebra is already canon:

- 01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md proves the complete factorial-Hadamard walk kernel for substitution.
- 01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md identifies its sidecar with an endpoint Taylor jet and gives a same-H hostile pair.
- 01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md supplies the order-join product and records loss of SCC order.
- THM-2195-transitive-quotients-exactly-control-universal-substitution-products, THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity, THM-2249-directed-triangle-forced-quotient-frustration, and THM-2256-automorphism-contact-dichotomy-for-quotient-frustration concern a DIFFERENT consumer, reversal distance. Their transport matrices must not be substituted for Hamiltonian run counts.
- 05-knowledge/results/zero_heavy_ports_20261020g/profile_growth.md already proves the exact quadratic endpoint-tensor response for arbitrary one-vertex attachment and a reinforced same-profile six-core hostile. Our four-vertex converse example below is an illustration, not a claim of new endpoint necessity. The corrected THM-1900/MISTAKES warning forbids assuming that deleting an internal vertex preserves a directed path.
- Session32 PORT_TRANSFER.md supplies the hostile boundary: arbitrary graph interfaces need named port pairings, not just component counts.

The new question is necessity, with a resource budget: exactly which factor information can all small external tournament contexts observe? The answer below turns the known sufficient profile into the precise equivalence classes for this continuation family.

Concept board: maximal module runs; quotient walks; endpoint jets; bounded experimental context; symmetry loss; named ports. Anchor: exact continuation equivalence. Niche: sharp context budget. Wildcard: intrinsic observer dimension versus ambient jet dimension.

## 1. Definitions

For a nonempty finite tournament S, let pc_S(c) count spanning covers by c unordered nonempty directed paths, singleton paths allowed. Define F_S(c)=c! pc_S(c), extended by zero for c>|S|. Thus F_S(1)=H(S).

A module context with m exterior vertices is a tournament Q on {*,1,...,m}; C_Q[S] substitutes S for * and a singleton for each other vertex. Every exterior vertex has a uniform orientation to all of S. Exterior vertices need not have any symmetry. Include m=0, for which C_Q[S]=S.

Define S equivalent_m T when H(C_Q[S])=H(C_Q[T]) for every such Q with at most m exterior vertices. No equality of factor orders is required.

## 2. The bounded-context equivalence theorem

For every integer m>=0 and nonempty tournaments S,T,

    S equivalent_m T
    iff F_S(c)=F_T(c) for 1<=c<=floor(m/2)+1.                 (1)

In particular equivalent_(2k) and equivalent_(2k+1) are identical relations. Equality of full path-cover profiles is exactly indistinguishability by ALL finite module contexts. This is a semantic minimality statement: any invariant sufficient for every such response must determine the profile. It does not mean the profile entries are algebraically independent on tournaments of a fixed order.

### Forward-response proof

Compress each maximal S-run of a Hamiltonian path to *. Its exterior vertices appear once each, and * appears c times, with no consecutive stars. Let a_Q(c) count these directed words, where star occurrences are indistinguishable. For any fixed word, decorating its star slots by an ordered path cover of S has F_S(c) choices. Hence

    H(C_Q[S]) = sum_(c>=1) a_Q(c) F_S(c).                   (2)

Between consecutive stars there cannot be just one exterior vertex v: the transitions would require both *->v and v->*, impossible in a tournament. Each of the c-1 gaps therefore contains at least two distinct exterior vertices. Thus

    2(c-1)<=m,
    a_Q(c)=0 for c>floor(m/2)+1.                            (3)

Equations (2)-(3) prove that the profile prefix suffices. This also explains why a single additional exterior vertex never adds information at an odd budget.

### Converse proof by explicit triangular experiments

Let T_j be the transitive tournament of order j, j>=1, and let

    R_j(S)=H(C3[S,T_j,T_j]),
    R_0(S)=H(S).                                          (4)

This is a valid module context with 2j exterior vertices. Put G_j(b)=b! StirlingSecondKind(j,b), zero outside 1<=b<=j. The inherited C3 word kernel gives

    R_j(S)=sum_(a=1)^(j+1) B_(j,a) F_S(a),                (5)

    B_(j,a)=sum_(b,c=1)^j w(a,b,c) G_j(b)G_j(c),

where w=3 if a=b=c, w=1 if max(a,b,c)-min(a,b,c)=1, and w=0 otherwise. In the last column only b=c=j contributes, so

    B_(j,j+1)=(j!)^2 != 0.                                (6)

Consequently R_0,...,R_k recover F_S(1),...,F_S(k+1) recursively:

    F_S(1)=R_0(S),
    F_S(j+1)=[R_j(S)-sum_(a=1)^j B_(j,a)F_S(a)]/(j!)^2.    (7)

Equality in all contexts of budget m includes these tests for j<=floor(m/2), proving the other implication of (1). The division is over integers for actual tournament input, but the invertible triangular map is naturally over Q.

The first experiments are

    R_1=3F(1)+F(2),
    R_2=11F(1)+17F(2)+4F(3),

with higher rows supplied by the exact script.

## 3. Sharpness and the correct meaning of minimality

No module context with fewer than 2j exterior vertices has any coefficient on F(j+1). The context in (4) with 2j vertices has coefficient (j!)^2. Thus 2j is the exact cost to expose that coordinate in the ambient linear response space.

For factors of a fixed order n, however, THM-3134 proves parity identities among the full profile entries: E_S(-u)=(-1)^(n-1)E_S(u), where

    sum_c F_S(c) z^c = z E_S(z+1/2).

Together with the fixed leading coefficient n!, these reduce the number of free coordinates. We therefore do NOT claim that k+1 independent measurements are necessary on the set of realizable n-vertex tournaments, or that every strict prefix refinement has a realizable hostile pair at every n. Formula (1), the sharp ambient coefficient bound, and the explicit reconstruction are the claims proved here.

There is a smaller sufficient observer budget when the common order n is known. For n>=3 put h=floor((n-1)/2). The first h profile coordinates determine the entire profile, so contexts of budget 2h-2 suffice for equality under every module context. To prove this, the difference E_S-E_T has parity n-1 and its leading coefficient cancels. Write it u^epsilon P(u^2), where epsilon is 0 or 1 according to the parity and deg P<=h-1. Equality of F(1),...,F(h) makes E_S-E_T vanish to order at least h at u=1/2, by the endpoint-jet formula. Since u is nonzero there and u->u^2 has nonzero derivative, P has a root of multiplicity h at 1/4, forcing P=0. Thus orders 3 and 4 need only H; orders 5 and 6 need at most the two-vertex probe; orders 7 and 8 need at most four exterior vertices. Orders 1 and 2 have a unique profile. This is a sufficient finite-size budget, not a sharpness claim on every realizable order.

The inherited order-five same-H pair has F=(15,78,198,240,120) and (15,90,210,240,120). All contexts with <=1 exterior vertex agree, but the two-vertex triangle context returns 123 versus 135. This is a realized sharp first-depth witness, and order five is the smallest possible common factor order: the fixed-order parity argument above shows that H already determines every profile through order four. With the script's edge-bit convention (lexicographic unordered pairs, bit 1 meaning lower-labelled vertex beats higher-labelled), the exact engine finds masks 69 and 76. They are invisible to every budget-0 or budget-1 module context, and the directed-triangle context at budget 2 separates them by 12 paths. Thus both minimal order and minimal external separating budget are proved, not just observed.

## 3a. The strongest first-possible observer is a permanent

The budget bound has a sharper, context-sensitive form. Let O be the exterior vertices beaten by the module, and I those beating it. Each return gap starts in O and ends in I, so

    c <= 1+min(|O|,|I|).                                  (9)

Suppose the context has exactly 2r exterior vertices, r>=1. Its highest possibly visible coefficient a_Q(r+1) is zero unless |O|=|I|=r. In the balanced case let B be the r-by-r zero-one matrix B_(o,i)=1 if o->i. Then

    a_Q(r+1)=r! permanent(B).                             (10)

Indeed, attaining r returns uses every exterior vertex and forces the full word to start and end in the module. Every gap has exactly two vertices, hence is o,i with o->i. The gaps form a perfect matching from O to I, and the r matched pairs can appear in any of r! orders. Conversely each such ordered matching makes a valid word. This is a direct bijection, not a determinant approximation.

Consequently

    a_Q(r+1)<=(r!)^2,

with equality precisely when every vertex of O beats every vertex of I. Equality in permanent(B)<=r! means every bijection is allowed; every matrix position belongs to a bijection, hence B is all ones. The orientations inside O and inside I are completely invisible to this top coefficient. Thus C3[S,T_r,T_r] is an extremal first-possible observer. One may replace either transitive side by ANY r-vertex tournament without changing the top coefficient, though lower response coefficients generally change.

The earliest way an environment can see the next path-cover coordinate is therefore governed by a bipartite matching count. This is a literal source-target map: maximal-return module words map bijectively to ordered perfect matchings. No matching statistic is asserted to determine lower coefficients or the entire tournament.

## 4. Composition, higher cover outputs, and reversal

THM-3121 shows that equality of full profiles is a congruence for arbitrary tournament substitution: replacing any factor by another with the same profile preserves the full output profile, and hence every nested substitution expression's H. This substitutes S or T as a FACTOR. It does not assert that equal-profile tournaments can be interchanged as QUOTIENT SKELETONS when their vertices are independently inflated; THM-3235 explicitly distinguishes that different operation. Associativity here is inherited from literal substitution and maximal-run splitting; no new matrix identity is needed.

There is also a useful filtered output bound. In a spanning d-path cover of C_Q[S], let c be the total number of S-runs. If p of the d global paths meet S, the return gaps number c-p. Each consumes at least two exterior vertices, so

    c <= p+floor(m/2) <= d+floor(m/2).                      (8)

Thus computing output pc(d) for fixed context requires no factor profile entry above d+floor(m/2). The argument holds even when some exterior paths avoid S. It is a dependency bound, not a claim that all these coordinates are independent or observed by every context.

Reversing all arcs of S reverses each path in any path cover, so F_(S^op)=F_S. Equation (1) implies that S and its converse are indistinguishable under every uniform module context, even without reversing that context. Likewise, factors reordered in an order join have equal profiles by THM-3166. These are exact information losses for this consumer, not recovery of tournament isomorphism, SCC order, or arc-reversal distance.

## 5. Why uniform modules and named ports are different

The uniformity hypothesis is essential. Let A=K1 joined before C3 and B=A^op, with corresponding vertex labels. Their full profiles agree. Attach one new vertex with a NONuniform pattern to these four named vertices; the pattern in which only vertex 0 (the distinguished singleton) beats the new vertex gives H=3 for A and H=15 for B. The new vertex beats all three triangle vertices. For A the extension is the order join K1, K1, C3 and has 3 paths. For B it is a cyclic substitution with the triangle as one factor and two singletons, giving 3F_C3(1)+F_C3(2)=9+6=15. Thus the discrepancy is proved directly; exact.json also records both adjacency matrices and the attachment pattern. Therefore the same full-profile invariant is insufficient for arbitrary pinned-vertex continuation.

For modules, external orientation cannot tell endpoints within S apart: only how many paths are available matters. For arbitrary boundaries, it can tell endpoints apart, so endpoint allocations and pairings survive. This is the precise bridge to Session32's port theorem: the path-cover profile is the symmetry-compressed response only when the whole block really is an orientation module. It is not a universal graph-contraction invariant.

Source: tournament block with a uniform exterior interface.
Target: scalar H under every bounded-size module continuation.
Map: maximal-run compression followed by contraction with a_Q.
Preserved: exactly the prefix in (1), and full substitution congruence in the unbounded limit.
Destroyed: endpoint identities, converse orientation, SCC ordering, and all other distinctions inside an equivalence class.
Sidecar when uniformity fails: named boundary endpoint/pairing response, not more scalar H values.
Cheapest decisive test: two exterior vertices in a directed triangle; hostile control: one nonuniform pinned exterior vertex.

## 6. Exact audit design and limitations

check.py uses only the Python standard library and explicit exception checks (valid under python -O). It independently enumerates quotient words and compares formula (2) with vertex-level Held-Karp DP, for every context tournament with 0..4 exterior vertices and every labelled factor of order 1..3. It checks the fixed-order parity-based observer sufficiency and reconstructs the complete profiles of every labelled tournament through order five using the triangular experiments; these profiles are independently computed from permutations and mandatory backward cuts. It checks ten exact triangular rows and the converse/nonuniform hostile, as well as the two-side bound and all 1,032 permanent-top formulas for contexts with two or four exterior vertices.

Reproduce from the workspace root:

    python3 tournament/check.py > tournament/exact.json
    python3 -O tournament/check.py > tournament/optimized.json

Neither the existing rational walk kernel nor endpoint jet is claimed as new. The new application is the exact bounded-context equivalence and two-vertex observability law. No {7,21} classification, LRC statement, Hamiltonicity algorithmic collapse, or literature priority follows.

Startup contract, relevant guardrails and core-paper routes were read; the bounded startup tool found the repository clean at e5b42c67e5d2. The coordination processor is unavailable because .machine-id is absent; no registration was invented. No Git mutation was made by this worker. Publication and independent review belong to the parent task.

## 7. A bounded negative probe beyond first-depth sharpness

The parity proof shows that a pair with equal F(1),F(2) but different full profiles must have common order at least seven. A deterministic sample of 20,000 labelled order-seven tournaments (seed 330071) found no such pair. This is a finite sampled negative, NOT exhaustive nonexistence. The profile calculator in find_twins.py uses a subset Hamiltonian DP followed by an anchored set-partition recurrence, independent of check.py's permutation-cut method. twins.json records the sample size, seed, distinct profile/prefix counts, all observed profiles, and a digest of the sampled masks. The actual existence and minimal order of a second-depth hostile remain OPEN here. No claim that every formal jet depth is realized by a sharp tournament pair is justified.

Reproduce this deliberately bounded probe with:

    python3 tournament/find_twins.py > tournament/twins.json

This stopping point is mathematical rather than numerical enthusiasm: the new theorem already proves the exact continuation relation, the first-depth hostile is minimal, and the deeper search exposes a potential realizability constraint rather than warranting an unbounded census.
