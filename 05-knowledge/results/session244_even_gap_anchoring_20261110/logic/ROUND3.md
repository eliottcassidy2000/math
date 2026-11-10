# Round 3: a finite forbidden-core grammar, and the necessary join repair

Status: **PROVED structural grammar equivalence and finite obstruction bound**;
**FINITE-EXACT** six-core classification and scalar witnesses. Arbitrary odd
full-count coverage remains **OPEN**. No literature-priority claim.

## 1. Inheritance and objects

The parent's structural observation is the starting point: unrestricted
free-reroot O/Z additions add a vertex with outdegree zero or one into the
old tournament. We distinguish that graph operation from its **guarded**
version, which permits Z only when its unique old out-neighbor reaches every
old vertex. O adds a sink. A root may be freely chosen before each addition.
The inherited O/Z/J state laws still require their correct marked root;
free rerooting is a graph operation, not a scalar four-coordinate formula.

**Direct forest inheritance.** Session199
`session199/hidden-structure/PROOF.md`, §3, already proves that fixed-root
O/Z/J histories reverse a chronological caterpillar forest, with at most one
nontrivial component and a spine of successive J roots. We do not claim the
forest encoding or chronological reversal observation as new. Free rerooting
allows arbitrary earlier parents and arbitrary recursive forest components,
which is a larger construction class. The new statements here concern its
hereditary degeneracy characterization, the exact six forbidden cores, and
the join repair of the reachability guard. They do not replace Session199's
stricter fixed-root chronological-caterpillar theorem.

THM-4868,
`01-canon/theorems/THM-4868-the-noncritical-run-of-a-uniquely-hamiltonian-tournament-local-roundness-one-run-cut-arcs-plus-noncritical-equals-n-plus-1.md`,
records the known sole-exit theorem for uniquely Hamiltonian tournaments
and its minimum-outdegree consequences. THM-4869,
`01-canon/theorems/THM-4869-sole-exit-free-tournaments-with-at-most-five-hamiltonian-cycles-terminal-cores-fibonacci-counts-and-the-arc-connectivity-route.md`,
studies sole-exit-free tournaments and small-cycle terminal cores. Neither
inspected statement supplies the hereditary graph class classified below.
These are relevant inherited mechanisms, not claims newly proved here.

Call a tournament **outdegree-1-degenerate** if every nonempty induced
subtournament has a vertex of outdegree at most 1. This is a graph property;
no statement about its Hamilton-path scalar image is built into the name.

## 2. Minimal forbidden cores have at most six vertices

More generally, let T be inclusion-minimal among tournaments with minimum
outdegree at least k>=1. Let S consist of its vertices of outdegree exactly
k, and put s=|S|, n=|T|. For each vertex v, minimality implies T-v contains a
vertex u of outdegree at most k-1. Since deleting v changes an old outdegree
by at most one, u has outdegree exactly k in T and u->v. Thus every vertex
is dominated by some vertex of S.

In particular at least n-s arcs go from S to its complement. The total
number of those arcs is exactly

    ks - s(s-1)/2.

Therefore

    n <= (k+1)s - s(s-1)/2 <= (k+1)(k+2)/2.                 (1)

The last integer maximum is attained at s=k+1 or k+2. For k=2, n<=6.
Every minimum-outdegree-2 tournament has at least five vertices, since its
sum of outdegrees is n(n-1)/2>=2n. Thus an outdegree-1-degenerate tournament
is characterized by avoiding the inclusion-minimal minimum-outdegree-2
cores of orders 5 and 6 as induced subtournaments.

This is a genuinely finite forbidden-graph criterion. It does not assert
that avoiding these graphs forces an arbitrary desired Hamilton-path count.
For general k, (1) similarly gives a finite obstruction-size bound for
outdegree-(k-1)-degeneracy. No optimality claim is made for k>2.

## 3. Exact six-core classification

The independent standard-library checker enumerates all labelled tournaments
through order 6. An n=5 or n=6 candidate has minimum outdegree at least 2,
and deleting every vertex must reduce the minimum below 2. At these orders
this is sufficient for inclusion-minimality, because every smaller possible
minimum-outdegree-2 subtournament has order at least 5.

The resulting isomorphism representatives, encoded by outgoing-neighbor
bitmasks, are:

    (6,12,24,17,3)                orbit 24,  H=15;
    (6,12,48,5,11,27)           orbit 240, H=33;
    (6,12,24,33,11,23)          orbit 720, H=37;
    (6,12,24,33,43,7)           orbit 720, H=37;
    (6,24,26,33,41,7)           orbit 240, H=45;
    (6,24,42,33,13,19)          orbit 240, H=45.

All permutations of each representative are generated. The six orbits are
disjoint and equal the exhaustive candidate sets: 24 labelled order-5 cores
and 2160 labelled order-6 cores. Together with (1), this finite check proves
that these are the complete forbidden list at every order.

## 4. Guarded grammar is smaller, already on a strong five-vertex graph

The unrestricted free-reroot O/Z grammar equals outdegree-1-degeneracy:
reverse a low-outdegree elimination ordering. Conversely each new vertex
has outdegree at most 1 in every induced subtournament containing it, while
old induced subtournaments remain valid.

The reachability guard changes the graph class. Starting from a singleton
(or the equivalent edge seed), every guarded Z addition with at least two
old vertices makes the entire graph strong. Indeed its new vertex reaches
the old root and hence all old vertices; each old nonroot reaches the new
vertex directly, and the root reaches the new vertex through an old nonroot.
O only appends a sink. Therefore every graph produced is transitive or a
strong initial component followed by a transitive sink tail. This is a
necessary condition, not a sufficient characterization of the grammar.

Already source=>C3, of order 4, is outdegree-1-degenerate but not guarded:
its nontrivial strong component is terminal, not initial. A strong obstruction
of order 5 is obtained from source a=>C3 by adding r which beats only a.
Its adjacency rows can be

    (14,20,24,18,1).

It is strong and its unique sole-exit vertex is r. Deleting r leaves
source=>C3, which is not guarded. A strong graph cannot end with an O step,
and its only possible final Z vertex is r, so no alternate guarded history
exists. The graph is outdegree-1-degenerate and has H=15. The independent
census confirms there is no strong counterexample through order 4.

The parent's order-7 variant, with T-r=C3=>C3 and r beating only a vertex of
the first triangle, is also valid and has H=45, but is not minimal. The
smaller obstruction matters when choosing a seed bank or claiming equality
of construction grammars.

## 5. Positive bridge: ordinal joins repair the guard exactly

Let G be the smallest graph class containing a singleton and closed under:

- O (append a sink);
- guarded Z at any freely selected reachable root;
- ordinal joins A=>B, with every cross arc from A to B.

Then **G is exactly the class of outdegree-1-degenerate tournaments**, hence
exactly the class avoiding the six cores in §3.

**Closure direction.** O and Z preserve the hereditary property by inspecting
whether a chosen induced subtournament contains the new vertex. Joins also
preserve it: take the last occupied strong/join block of an induced subset;
a low-outdegree vertex there has no outgoing arcs to preceding blocks.

**Construction direction.** Induct on graph order. For a nonstrong tournament,
its proper strong components form an ordinal join. Each is hereditary and
therefore already constructible by induction. For a strong tournament on at
least three vertices, degeneracy gives a vertex r of outdegree exactly 1,
say r->a. The parent T-r is hereditary, hence constructible by induction,
including any joins needed there. Since T is strong, a reaches every vertex
of T-r: a simple path from r to any old vertex starts with a and does not
return through r. Thus adding r to this parent is a permitted guarded Z at
a. The singleton and edge cases are bases. This completes the construction.

This gives a correct finite structural certificate for the **join-augmented**
free-reroot grammar. The join is already a legitimate physical operation
whose path counts multiply. It repairs the missing graph histories, not just
the scalar count set. The theorem needs free rerooting and joins **inside**
the recursion; merely multiplying the final scalar outputs of the original
fixed-seed grammar is not asserted to reproduce the same graph family.

## 6. Forbidden graphs do not forbid their scalar values

Every distinct count of the six excluded cores already has a literal guarded
fixed-root witness from the edge 0->1 rooted at 0:

    H=15: OJOZ;
    H=33: OOOOZ;
    H=37: OZZOZ;
    H=45: OZOZZ.

The checker constructs these graphs, checks the reachability guard at every
step, and independently counts Hamilton paths by subset dynamic programming.
Their orders are 6,7,7,7 and their final graphs are strong. Thus the forbidden
core count set {15,33,37,45} is wholly contained in even the smaller guarded
fixed-root spectrum. Graph-family strictness does not prove scalar-spectrum
strictness.

Nor can an excluded induced core simply be replaced inside a larger tournament
by one of these same-count witnesses: external attachment incidence and
path-cover responses need not agree. A positive all-value bridge would need
a replacement preserving the exact attachment sidecar, or a target constructor
inside this six-forbidden family. Neither follows from the scalar matches.

## 7. Verification

Run `python3 -I -B verify_round3.py` and compare with `ROUND3_CHECK.json`;
repeat under `-O`. The checker uses explicit exceptions, not assertions.
It independently enumerates 1,2,8,64,1024,32768 labelled tournaments at
orders 1,...,6, tests elimination and guarded histories by exhaustive cached
deletion recursion, verifies all core orbits, and computes the literal scalar
witnesses above. For reference the degeneracy-class sizes are

    1,2,8,64,1000,26720,

whereas the guarded-grammar sizes are

    1,2,8,56,760,18960.

These are labelled graph counts, not counts of distinct Hamilton-path values.
The infinite equivalence and size bound are proved above; the finite
classification is exhaustive within the proved bound. Chan–Pak remains open.
