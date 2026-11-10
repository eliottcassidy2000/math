# Structural continuation: strong roots of outdegree one

**INHERITED classification, independently recovered; scoped corollaries proved below.** The complete nonsink zero-defect classification, exact factorization, sharp h≥3e inequality, and endpoint-equals-one dyadic slice are already proved in `05-knowledge/results/session231_overnight_completion_20261109/descent/ROUND2.md`, Sections1–2 (independent audit in its `descent_audit/ROUND2_AUDIT.md`). The outdegree obstruction is older still: THM-4873(ii)(D) gives the stronger δ≥r(r−1), and THM-4875(A1) supplies a weighted version. The reachability argument here is an alternative proof, not a new global structural theorem. The e=3 divisibility, small-count exclusions, and prime gcd filters below are immediate scoped applications; no external priority is asserted.

For all nontrivial strongly connected rooted tournaments, The saturation condition U=H forces outdegree one, with the exact SCC decomposition below. The actual output roots of the new five-letter compiler also show this directly: each last J root initially beats only the preceding root, and the final Z beats it, so its outdegree remains one.

## 1. Saturation forces outdegree one

Split each Hamilton path immediately before its root r. This injects the paths into the covers counted by U: the distinguished suffix starts at r and the complementary prefix is allowed to be empty. A cover is in the image iff its complementary path Q is empty or its last vertex points to r. Thus U−H counts precisely the covers with Q nonempty and r pointing to its last vertex.

We need a reachability lemma. Suppose r reaches every vertex of a tournament. Among the outneighbors b of r, at most one has the property that r fails to reach every vertex after b is deleted. Otherwise let b,c be two such outneighbors, and let S_b,S_c be the nonempty sets unreachable from r after deleting b,c respectively.

These two sets are disjoint. If x lay in both, every r-to-x path would contain both b and c. Choose a simple such path. Replace its prefix ending at the later of b,c by the direct arc from r to that vertex. The resulting path avoids the earlier vertex, a contradiction. Also every vertex of S_b points to every vertex reachable from r without b; otherwise an arc in the opposite direction would make it reachable. Since S_b and S_c are disjoint, each vertex of S_c is reachable without b, and conversely each vertex of S_b is reachable without c. This forces both S_b→S_c and S_c→S_b, impossible in a tournament.

A vertex reaches all vertices of a tournament iff there exists a Hamilton path starting there. One proof uses the ordered strong components: the vertex belongs to the first component; cut a Hamilton cycle there at the specified start, then concatenate Hamilton paths of all later components. Singleton components are immediate.

If r has at least two outneighbors, choose an outneighbor b whose deletion preserves full reachability from r. Take a Hamilton path P of T−b starting at r and Q=(b). This cover contributes to U−H, so U>H. In a nontrivial strong tournament the root cannot be a sink. Therefore saturation implies that its outdegree is exactly one.

## 2. The exact decomposition

Let T be strong, r its root, and N+(r)={a}. Write S=T−r. Every Hamilton path starting r consists of r followed by a path of S starting a. Removing r from the distinguished cover path gives

    U_r(T)=H(S)+B_a(S).

The inherited singleton graft law gives

    H(T)=H(S)+B_a(S)−e_a(S),   e_r(T)=H(S)−e_a(S).

Therefore U_r(T)−H(T)=e_a(S). In particular r is saturated iff e_a(S)=0, and then e_r(T)=H(S)=d_r(T).

Hamilton paths respect ordered strong components, and a vertex is the endpoint of at least one Hamilton path exactly when it belongs to the last strong component. Strongness of T forces a to belong to the first strong component of S: otherwise r cannot reach the preceding components. Thus, in the saturated case,

    S=A⇒B,

where A is the first strong component, a belongs to A, and B is nonempty. Conversely, any strong A with chosen a and nonempty tournament B, followed by adding r that beats a and loses to all other vertices, is strong and has saturated outdegree-one root r.

Let Q(B)=sum_X H(B[X])H(B[V(B)\X]), with H(empty)=1. Ordered two-path covers of A⇒B factor by splitting each path into its A and B portions. The distinguished path must have a nonempty A portion starting a, but either of its B portions can be empty. Hence

    B_a(A⇒B)=B_a(A)Q(B).

The full saturated face parameters are consequently

    e=H(A)H(B),
    h=H(A)H(B)+B_a(A)Q(B).                 (1)

Combined with Section1, this is an exact constructive characterization of all saturated roots in nontrivial strong tournaments. The two factors are not freely selectable arithmetic coordinates.

Since B_a(A)≥H(A) and Q(B)≥2H(B), equation(1) proves h≥3e. If |B|≥2, the proper-subset terms of Q are positive, so Q(B)>2H(B). Equality h=3e therefore holds exactly when B is a singleton and the chosen root of A is itself saturated. This is the OJ (equivalently scalar JZ) recursion (h,e)=(3H(A),H(A)).

## Small endpoint rigidity

If e=1, then H(A)=H(B)=1. Strong A must be a singleton and B must be transitive. Thus

    h=1+2^k, k=|B|≥1.

If e=3, one of H(A),H(B) is3 and the other1. A tournament of count3 consists of a single directed-triangle strong component with transitive singleton components before and after it (inherited strong-fibre classification). Q is multiplicative on ordinal sums, Q(C3)=12, Q(singleton)=2, and B_a(C3)=3. Therefore h is divisible by3 in either case. In particular no saturated strong tournament has H=13: h≥3e forces e∈{1,3}; neither case can yield13. This last exclusion is only for this rooted class, not total tournaments (13 is of course attained).

## Consequence for an induction attempt

The pair e,h is not an arbitrary odd lattice point in the half-bound cone. Even the stronger h≥3e cone contains physically forbidden states, for example (13,1) and(13,3). A proof based on inverse word maps must preserve the decomposition(1), or use a stronger correct physical invariant, rather than selecting an integer predecessor solely by parity and inequalities. An exact response realization for A or a controllable pair(H(B),Q(B)) is still required.

## 3. Prime obstruction and arbitrary tournaments

If h=p is prime in(1), any common divisor of H(A) and B_a(A) divides p. It is at most H(A)≤e≤p/3, so it must be1. The same argument on B gives

    gcd(H(A),B_a(A))=1,   gcd(H(B),Q(B))=1.        (2)

Because both H and Q multiply over ordinal sums, each strong component C of B must satisfy gcd(H(C),Q(C))=1. In particular B has no directed-triangle strong component: H(C3)=3 and Q(C3)=12. Also A cannot be a directed triangle, whose every root has response3. These necessary conditions constrain a prime-target search; they do not imply that a candidate sum is prime or realize any chosen prime.

The saturation classification extends cleanly to tournaments not assumed strong. Write T=P⇒C⇒B, where C is the strong component containing the root, and P and B may be empty. The same cover decomposition gives

    H(T)=H(P)H(C)H(B),
    U_r(T)=H(P)U_r(C)Q(B),

where Q(empty)=1. If B is nonempty, Q(B)≥2H(B), so U_r(T)≥2H(T)>H(T). Thus saturation forces C to be the terminal component. If C is a singleton, r is a sink and saturation is automatic. Otherwise C is a nontrivial saturated strong tournament, and both h and e are multiplied by H(P).

Consequently the only way H=11 or H=13 can have a saturated root is at a sink. Indeed, for a nonsink saturated root its terminal component has count at least3 dividing the prime H, so would itself be a saturated strong tournament of count11 or13, excluded above. This is an all-order exclusion for a rooted interface, not for ordinary Hamilton-path counts: both11 and13 are attained.

The prior Section2 proof for13 applies verbatim to11, since in either case e≤h/3 permits only e=1 or3, and neither the dyadic family1+2^k nor the divisible-by3 family contains the target.
