# Balanced black-core geodesics, with distinct named starts

**Session 75, 2026-10-27. PROVED in the stated scopes, independently audited; finite heads explicitly certified.**
The regular theorem is the user-supplied claim, reconstructed here; it is not
presented as a newly discovered statement. The even-order extension and target
separation are new work in this packet. No all-order black connectivity claim.

## 0. Inheritance, source types, and concept board

Pinned repository source: `96e35642110104cc547516a63f37eebcf930f707`.
`SOURCE_LEDGER.json` lists exact paths and file hashes. The closest mechanism
is the named endpoint current in `black_descent_20261020c/scc.md` and
`merged_metagraph_endpoint_transport_20261020/endpoint-transport-and-complement-leakage.md`.
The regular claim itself is trusted transcript history; its exact older file
location was not yet identified when this packet was frozen.

The black-path existence input is the corrected, audited uniform prefix
classification, not the old affine-reflection shortcut:
`pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`
and `blue_quotient_20261020b/structure.md`, Section 2. Their nontrivial strong
pure-blue exceptions are C3 and the regular cyclic R5. The even extension
additionally uses the **strong black descending-path theorem**, specifically
`black_descent_20261020c/cycles.md`, Sections 5–6. Mere existence of a black
path does NOT suffice for that endpoint-score selection.

Canonical hostiles: C3 and R5; C3 -> C3 and C3 -> 1 -> C3; and the new
same-start route below. Corrected near miss: a cut proof alone permits
source + strong middle + sink; named starts rule out the second singleton.
Least-used sidecar: a literal second start chosen from a single retained
Hamilton cycle of the deleted core.

Board: (1) signed cut current, (2) fixed endpoint charge, (3) two distinct
starts, (4) core deletion connectivity, (5) actual R/K edge words,
(6) nonregular distinction lambda versus minimum semidegree delta.
Anchor: reconstruct regular geodesic. Niche: even balanced extension.
Wildcard: compare exact all-colour cut distance with black attachment.

## 1. Carrier and exact gauge conversion

A tournament has fixed named vertices V. A directed Hamilton word
P=(p0,...,p(n-1)) is blue iff its positional reversal is an anti-automorphism:

    T(p_i,p_j) = T(p_(n-1-j),p_(n-1-i)) for every i<j.

It is black otherwise; every certificate retains an actual failed equation.
R_P reverses precisely the n-1 path arcs. K_P complements precisely the
NONpath pairs, so K_P(T)=(R_P(T))^op. These are not single-edge reversals.
The original P remains directed across K; reversed P is directed across R.
Blackness is stable at both ends of either edge.

For an R route U_j -> U_(j+1) along P_j, the literal K route is

    T_j = U_j              if j is even,
          U_j^op           if j is odd;
    Q_j = P_j              if j is even,
          reversed(P_j)    if j is odd.

Thus each R word supplies an exact named K word, without relabeling or
forgetting path starts. Both nonstrongness and the target

    C = {1 -> S, S -> 1 : S strong on n-1 vertices}

are invariant under converse. Every distance statement therefore holds for
literal labelled K routes, and also after ordinary/converse isomorphism
quotients. It concerns C, not the forest or transitive targets of Session 74.

## 2. General cut and endpoint bounds; an exact ALL-COLOUR theorem

For nonempty proper A put c_T(A)=number of arcs A -> V\A. Under R_P,

    c_(R_P T)(A)-c_T(A) = 1_A(last P)-1_A(first P).       (1)

Indeed the path's signed crossings telescope. For singleton A this is the
exact degree current. Define

    lambda(T)=min_(empty != A != V) c_T(A),
    delta(T)=min_v min(d_T(v), n-1-d_T(v)).

**PROVED, all colours.** Distance to the nonstrong locus is exactly lambda.
For a strong T the lower bound follows from (1): a nonstrong tournament has
a zero directed cut. For the upper bound fix a minimizing A. While fewer
than lambda moves have been made every directed cut is positive, so the
current tournament is strong. A Hamilton cycle meets both shores and has
an arc from V\A into A. Cut at that arc: the resulting Hamilton path starts
in A and ends outside A. Equation (1) lowers c(A) by one. Repeat lambda
times. At the last step c(A)=0. For a nonstrong source lambda=distance=0.
No assertion that these chosen paths are black is made.

**PROVED lower bounds, black or all colours.** Distance to C is at least
delta. Any terminal member has a source or sink, and a single named degree
can change by at most one per R move. A fixed named source target v requires
at least d^-(v) moves; a fixed named sink requires at least d^+(v).
The lambda bound alone is not an exact core-attachment formula (Section 5).

## 3. Regular odd order: exact black distance m

**Reconstructed theorem.** Let n=2m+1>=7 and T be regular. Its black distance
to C is exactly m. The construction is explicit and retains labels.

For |A|=a, regularity gives

    c_T(A)=a(2m+1-a)/2.                                (2)

Every proper cut is at least m. When 2<=a<=n-2 it is at least 2m-1>m.
Consequently every state before step m of ANY R route is strong.

Choose a black Hamilton path P1; its existence at these orders follows from
the corrected pure-blue classification. Write a=first(P1), v=last(P1).
We will always end subsequent paths at v. After j such moves,

    d_j(v)=m+j,
    d_j(u)=m-N_j(u) for u!=v,                          (3)

where N_j(u) counts starts at the actual named vertex u.

**The second start exists.** In U=T-v all semidegrees are >=m-1. For a cut
whose smaller shore has size s<=m, either the outdegree sum on that shore
or the indegree sum on the other orientation gives

    c_U(A) >= s(m-1)-s(s-1)/2 >= m-1.                  (4)

The last inequality follows because this concave quadratic on 1<=s<=m
has minimum at an endpoint, with values m-1 and m(m-1)/2. Restricting R_P1
to U reverses a Hamilton path, so its cut counts change by at most one.
Thus U1=T1-v has every cut >=m-2>=1 and is strong.

Fix ONE actual Hamilton cycle C1 of U1. Vertex v has m-1>=2 inneighbors in
U1. For each such inneighbor u, rotate C1 into a Hamilton path beginning
at succ_C1(u) and ending at u, then append v. Distinct u give distinct
named starts. Choose P2 with first(P2)!=a. This is the retained port guard.
At time1 only a has degree m-1; v has degree m+1, while P2's start has
score m. Blue reversal would require its endpoint scores to sum to 2m,
whereas they sum to 2m+1. Hence P2 is black.

After the split, for every j>=2, N_j(u)<=j-1. At time j any path ending at v
has endpoint-score sum at least

    (m+j)+(m-j+1)=2m+1>n-1.                            (5)

So EVERY such path is black. As long as j<m, the whole tournament is
strong by (2), and a Hamilton path ending at v exists. Continue to m moves.

At time m, v is a source, every other vertex has score at least1, and no
nonsingleton cut can be zero by (2). The remaining core is strong: otherwise
an ordered SCC suffix would yield either a forbidden nonsingleton zero cut
or a second singleton sink. Thus the endpoint lies in C. Lower bound m
comes from either delta=m or (2).

**Why the split does two jobs.** It prevents later blue paths AND prevents
a second extreme vertex at the terminal. Excluding nonsingleton cuts alone
would still allow 1 -> strong-middle -> 1.

### Literal same-start hostile

Start with the regular seven-tournament

    rows=(112,37,73,67,14,28,50).

These three directed R paths have common start6 and end0:

    (6,1,5,4,2,3,0)  black,
    (6,4,5,3,1,2,0)  blue,
    (6,5,4,3,2,1,0)  blue.

Their final shape is 1 -> R5 -> 1, not C. Thus a greedy fixed endpoint
route with an unrecorded/reused start does not prove the black theorem.

## 4. Extension: every even balanced tournament

**Theorem (new here; structural tail plus finite heads).** For every even
n=2m>=4, if all outdegrees are m-1 or m, the black distance to C is m-1.
The proof below is structural for n>=8. The n=4,6 heads are complete exact
finite proofs with explicit labelled universes and independent orbit cover.

Assume m>=4 and put k=m-1. For a cut with smaller shore s<=m, both
minimum indegree and minimum outdegree are k, giving

    c_T(A) >= s(m-1)-s(s-1)/2.

Its minimum is k, and for s>=2 its minimum is 2m-3>k. (For m=4 the
endpoint s=m also gives6>=5; in general m(m-1)/2>=2m-3.) Hence all states
before k moves are strong and all nonsingleton cuts stay positive at k.

By the inherited **descending black path theorem**, choose black P1 with
score(last)>=score(first). If its last vertex has score m, keep it. If not,
both endpoints have score m-1; pass to the converse and reverse P1. Its new
last vertex then has score m. This step uses the direction inequality;
an arbitrary black path with scores (m,m-1) would not give the argument.

Now v has indegree k and every other vertex has score >=k. As in (3),
d_j(v)=m+j and d_j(u)=d_0(u)-N_j(u). In the deletion T-v, minimum
semidegree is >=m-2. Its smaller shores have size<=m-1, so its cuts are
>=m-2 by the same concave bound. After the first move all core cuts remain
>=m-3>=1. Vertex v has k-1=m-2>=2 inneighbors. Rotating one fixed core
cycle again gives P2 with a genuinely different start.

At time1 this different start has degree>=k; its sum with d_1(v)=m+1 is
at least2m>n-1. After the split, for j>=2 every start has degree
>=k-j+1=m-j; again its endpoint sum with m+j is >=2m>n-1. Thus every
chosen path is black. At k steps v is a source and every other score>=1;
the nonsingleton-cut margin makes the remaining core strong as before.
The lower bound is delta=k.

### Exact finite heads

All 2^(n choose2) labelled tournaments are filtered ONLY by balance. There
are 24 retained objects at n4 and 2640 at n6. `check_black_core.py` gives a
length1 or length2 black route for every one, with exact K conversion and
all cut-current equations checked. The second starts differ when length2.

The independent audit takes one n4 and five n6 seeds, generates their full
label-permutation orbits, and proves coverage against an independently
computed degree-sequence recurrence. It obtains 24 and2640. Every seed's
literal K arcs and colours are checked using arc sets rather than producer
bit operations. Thus these heads are finite proofs, not random evidence.

## 5. Sharp nonregular separation of the two targets

Take two regular cyclic R5 blocks A={0,...,4}, B={5,...,9}, initially A->B,
and reverse only 0->5. The resulting strong tournament has

    rows=(966,1004,1016,1009,995,193,384,768,544,96),
    scores=(6,7,7,7,7,3,2,2,2,2), lambda=1, delta=2.

Thus its all-colour distance to nonstrong is1 but core attachment requires
at least2. These two BLACK R words achieve that lower bound exactly:

    (2,3,4,6,7,8,9,5,0,1),
    (3,0,2,5,9,8,7,6,4,1).

The source labels and every literal K edge are in the receipt. This
REFUTES any inference `distance to C=lambda` beyond the balanced lane.
A direct complete first-step probe also finds black nonstrong exits, so this
particular example has black nonstrong distance1 and black core distance2.
The two distances are already separated in the black graph itself.

## 6. Controls, receipts, and scope

Commands, all from this directory:

    python3 -B check_black_core.py > normal.json
    python3 -B -O check_black_core.py > optimized.json
    cmp normal.json optimized.json
    python3 -B independent_audit.py > audit-normal.json
    python3 -B -O independent_audit.py > audit-optimized.json
    cmp audit-normal.json audit-optimized.json

Producer checks complete balanced labelled universes n4=24, n6=2640,
n7=2640. It also constructs 42 specified larger examples through n13,
including Paley11's required five-black route. These examples are finite
controls for the proof, not its source. The independent audit checks 210
actual K edges, named labels, directedness at both ends, blackness at both
ends, exact target strongness, and endpoint lower bounds. It verifies
complete orbit coverage at4/6/7 against a memoized score-vector counting
recurrence. Both programs use explicit failures, not Python assertions.

Pure-blue controls independently enumerate ALL permutations: C3 has3
Hamilton paths, R5 has15, C3->C3 has9, C3->1->C3 has9; none is black.
The corrected HYP-9348 filters remain unchanged. The old unfiltered
endpoint-law claims are not reinstated.

Auxiliary noncanonical probes found: full n6 black BFS gives distances to
both nonstrong and C equal lambda for 19,680 strong objects withlambda1
and2,640 withlambda2; inherited complete n7 class-graph probing gives
266/84/3 strong classes at distances1/2/3. These are finite observations,
not an all-order equality and not needed for the proofs above.

**OPEN / next gate, not a theorem:** for strong black-active T of order>=4,
is black distance to nonstrong always lambda(T), and is black distance to C
always delta(T)? These are separate questions. First inspect prior sources
and corrected hypotheses, then test one-back-arc pure-blue block joins,
cyclic substitutions with unequal block sizes, and prescribed-endpoint
bottlenecks. Maintain actual path colours and labels at every step. The
n10 separation above is a positive two-target control, not a counterexample
to either correctly separated question. All-order black-core connectivity
and arbitrary module-substitution lifting remain OPEN.
