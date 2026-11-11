# The 3–4–5 blocks need named path ports, not three supervertices

**Status:** PROVED elementary transfer and obstruction statements below;
FINITE-EXACT local and global census at n=25. No novelty claim for general
boundary-state dynamic programming. No asymptotic Hamiltonicity claim follows.

## 1. Inheritance, scope, and the exact graph

Let G_n have vertices 1,...,n, with distinct u,v adjacent precisely when
u+v is an integer square. The closest inherited mechanism is the
endpoint-aware reduction in
`05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md`.
Its corrected near miss is forcing both edges at a degree-two vertex without
first knowing that it is not an endpoint. The exact n=24 obstruction and
compulsory n=25 ear are already proved in
`05-knowledge/results/decoder_prime_square_20260925.md`, sections 4–5.
The n=15 onset is already audited in
`THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md`.
The least-used relevant sidecar here is the complete named-boundary
continuation relation, analogous in purpose, but not in its mathematical
state space, to `session25-ports/PORT_SIGNATURES.md`.

The concept board is: edge-colour reflections; the proposed tree blocks;
endpoint-aware linear forests; named terminal pairings; global connectivity.
The anchor is exact Hamiltonian feasibility, the niche is a lossless block
transfer, and the wildcard is the cycle-rank loss ledger for contraction.
No colouring, Collatz, tournament, or analytic conclusion is imported.

The proposed disjoint blocks are

    A = {2,7,9,16,18},  B = {4,5,11,12},  C = {1,3,6,8,10,17}.

Their induced internal edges are respectively

    A: 2–7, 7–9, 9–16, 7–18;
    B: 12–4, 4–5, 5–11;
    C: 10–6, 6–3, 3–1, 1–8, 8–17.

Thus the internal edge counts are 4,3,5, and the vertex counts are 5,4,6.
All three blocks are trees, but A is branched, not a path. Its vertex 7 has
internal degree three. B and C are paths. In particular, no Hamiltonian
path of G_25 can visit all of A in one consecutive run.

The remaining core is

    K = {13,14,15,19,20,21,22,23,24,25},
    E(K) = {13–23, 14–22, 15–21, 24–25}.

All 16 boundary edges are listed here, with the block vertex first:

    A: 2–14, 2–23, 16–20;
    B: 4–21, 5–20, 11–14, 11–25, 12–13, 12–24;
    C: 1–15, 1–24, 3–13, 3–22, 6–19, 10–15, 17–19.

There are no A–B, A–C, or B–C edges. The counts are

    |E(G_25)| = (4+3+5) + 16 + 4 = 32,
    cycle rank = 32 − 25 + 1 = 8.

The graph is connected, as also certified by the explicit path below.

## 2. The actual contraction is a counterexample

Contract every entire block A,B,C to one supervertex, retaining parallel
edges. The resulting multigraph H has 13 vertices, 20 edges, and cycle
rank 8. This is the exact invariant preserved by contracting these trees:
each tree-edge contraction reduces both the edge and vertex counts by one
and preserves connectivity.

But **H has no Hamiltonian path, although G_25 does**. To prove the negative
statement, delete its three block supervertices. The remaining core has
six connected components: the four disjoint edges above and isolated
vertices 19,20. Every graph with a Hamiltonian path satisfies

    c(H−S) ≤ |S|+1

whenever vertices remain: deleting |S| vertices from the spanning path
leaves at most |S|+1 path pieces, and extra graph edges can only join
components. Here 6>3+1, a contradiction. Parallel edges do not change this
argument. This is an elementary certificate, not merely a failed search.

If the quotient is simplified to a simple graph, two parallel pairs are
lost, namely C–15 and C–19. It then has 18 edges and cycle rank 6, not 8.
Therefore even the cycle-rank claim requires preserving multiplicities.

The opposite implications fail in general too: P4 and the three-leaf
star K_(1,3) are connected trees with three edges and cycle rank zero, but
only P4 has a Hamiltonian path. K_(1,3) is the smallest nonpath tree;
trees with at most three vertices are paths. Contracting a whole tree to
one vertex plainly forgets that obstruction.

**Loss ledger.** The numbers 3,4,5, together with the knowledge that each
block is a tree, determine block sizes and their contribution to edge
bookkeeping. They do not determine port locations, compatible endpoint
allocations, path-cover components, or the ordering of those components
along a global path. Nothing here gives a Pythagorean geometric
interpretation of 3²+4²=5²; equal small counts are not such a map.

## 3. Exact finite transfer theorem

The following statement applies to any finite simple graph with at least
two vertices, partitioned into a nonempty core K and pairwise nonadjacent
blocks X_1,...,X_t. The single-vertex path is a separate trivial case.
A port is a **named boundary edge**, not merely its block vertex: a vertex
may have two different outgoing ports, as vertex 2 does in A.

For each block X, enumerate pairs (F,Q), where F is a spanning linear
forest on X using internal edges, and Q is a subset of its boundary ports.
Impose:

1. For each v in X, d_F(v)+d_Q(v) is 1 or 2.
2. At most two vertices have total degree 1; record those vertices as
   prospective global endpoints.
3. Every component of F has at least one selected port.

Singleton components are allowed. Their two terminals may consist of two
different boundary edges incident to the same vertex, or of one boundary
edge and an endpoint at that vertex. A component with no port is forbidden
because K is nonempty and the desired global path must be connected.

Each component has exactly two terminals, counting selected ports and
prospective endpoints. Indeed, if it has m vertices, then F contributes
2(m−1) to its degree sum; if e of its vertices have desired degree 1 and
p ports are selected, then

    2(m−1)+p = 2m−e,  hence p+e=2.

The complete transfer signature is the selected named port set, together
with its pairing into components and any endpoint terminal attached to
each component. Retain a count of realizing forests when signatures
coincide. Endpoint labels can be retained when prescribed endpoints or an
endpoint spectrum matter. For bare existence, labels may be forgotten,
but the port-to-endpoint pairing and number of endpoints cannot.

To glue signatures, choose a core-edge subset J. Replace each two-port
component by one virtual edge between the two corresponding core
vertices. Replace each one-port/one-endpoint component by one virtual
edge from its core vertex to a new marked endpoint vertex. Retain loops
and parallel edges if they arise. Accept precisely when:

    every core vertex has degree 1 or 2;
    the total number of degree-one vertices is 2;
    the reduced multigraph is connected.

**Theorem.** Accepted realizations are in bijection with Hamiltonian
paths of the original graph, considered as unoriented edge sets.

**Proof.** Restrict a Hamiltonian path to each block. These restrictions
are spanning linear forests and satisfy conditions 1–3. Contracting their
components to virtual edges preserves both connectivity and all exposed
vertex degrees, and yields an accepted reduced path. Conversely choose a
realizing forest for each accepted signature and expand each virtual
edge by its recorded component and boundary edges. All block vertices
are covered once, core vertices have the required degrees, and expansion
preserves connectivity. A finite connected graph with exactly two
degree-one vertices and every other vertex of degree two is a path.
Restriction and expansion are inverse operations. For signatures with
multiple realizations, multiply their weights and sum over accepted
gluings. This proves the counting version too.

**Why pairing is necessary.** Expose four ports 1,2,3,4. Two local
forests can pair them as (1,2)(3,4) or as (1,4)(2,3), with the same used
ports, degree data, and component count. Attach the external path 2–5–3.
The first gives the spanning path 1–2–5–3–4. The second gives the separate
edge 1–4 and cycle 2–3–5–2. A state that remembers only degrees or component
counts cannot distinguish this continuation.

This is an exact finite representation, not a claim that the chosen
encoding is the smallest possible state space. Restricted continuation
families may admit additional safe identifications.

## 4. The five exact states of the branched block A

Write a=2–14, b=16–20, c=2–23. The complete endpoint-labelled table is:

| Selected ports | Component terminal pairs | Endpoint set |
|---|---|---|
| a,b | a↔18; b↔9 | {18,9} |
| c,b | c↔18; b↔9 | {18,9} |
| a,b | a↔2; b↔18 | {18,2} |
| c,b | c↔2; b↔18 | {18,2} |
| a,b,c | a↔c; b↔18 | {18} |

These five states can be proved without enumeration. Vertex 18 is a
global leaf, so 18–7 is selected and 7 is not an endpoint. If 7–2 is
selected, 7–9 is not, forcing 9 to be the second endpoint and 9–16–20 to
be selected. Vertex 2 must select exactly one of a,c. If 7–9 is selected,
the path 18–7–9 must continue through 16 to 20: stopping at 9 or 16 would
make a component with no boundary port. Vertex 2 then selects either
one of a,c and is the second endpoint, or both and is interior.

Consequently every global path uses exactly two internal A-components,
always selects 7–18 and 9–16, and selects exactly one of 2–7 and 7–9.
If the other endpoint lies outside A, its state is uniquely the last
row: replace A exactly by the virtual edge 14–23 and the endpoint path
20–16–9–7–18. The first four rows use both global endpoints already.

For unprescribed endpoint feasibility, the first and third rows merge,
and the second and fourth rows merge. A then has three anonymous-endpoint
states with counting weights 2,2,1. It would be incorrect to claim that
endpoint *identity* is essential to bare Hamiltonian existence: the
endpoint allocation and terminal matching are what must survive.

## 5. Complete n=25 verification and an explicit witness

`verify_ports.py` enumerates the complete finite universe, with no
random sampling and no imposed endpoint spectrum. The only inherited
filters are the three proved local conditions and the global two-endpoint
budget. All nonreversal paths are independently enumerated by ordinary
vertex DFS from 18, whose global leaf status makes the starting choice
complete rather than a restriction.

| Block | Endpoint-labelled signatures | Anonymous-endpoint signatures |
|---|---:|---:|
| A | 5 | 3 |
| B | 101 | 89 |
| C | 149 | 131 |

For these tree blocks the labelled signatures have one realization each:
the unique tree path between each pair of named terminals determines the
component, including the possible zero-edge singleton path. The script
also checks uniqueness explicitly. Anonymous counts refer to distinct
port pairings after forgetting local endpoint identities, with their
realization multiplicities retained in the certificate histogram.

The total endpoint budget leaves 1543 signature triples. Testing the 16
subsets of the four core edges gives 19 degree-feasible assemblies. Nine
are disconnected and contain cycles. Exactly 10 are connected paths,
agreeing edge-for-edge with independent vertex DFS. Their second endpoints
are precisely

    {2,3,4,8,9,10,11,13,22,23},

with one unoriented Hamiltonian path for each endpoint. This is a
FINITE-EXACT result for G_25, not a general uniqueness principle.

One complete path is

    18,7,9,16,20,5,11,25,24,1,3,22,14,2,23,13,12,4,21,15,10,6,19,17,8.

Its A-pieces are 18–7–9–16 and singleton 2. Its B-pieces are 5–11 and 12–4.
Its C-pieces are 1–3,10–6, and 17–8. Thus it uses the blocks in multiple
runs, exactly the information erased by the three-supervertex quotient.
Also 8 has graph degree two, with neighbors 1,17, but is a path endpoint
and omits 8–1. This is an explicit hostile to unguarded degree-two forcing.

The complete certificate contains all 255 local states, all 10 global
solutions, profile IDs, exact graph edges, and three rejected
degree-feasible disconnected examples. It also verifies the P4/star
hostile, the four-port pairing hostile, the contraction cut obstruction,
and the distinction between 20 multigraph and 18 simple quotient edges.

## 6. Square colours as partial reflections

For a square q, the q-coloured edges are the matching induced by

    R_q(v)=q−v,

where both endpoints belong to 1,...,n and R_q(v)≠v. Thus each colour is
a partial fixed-point-free involution. On the full integer line,

    R_b R_a(v)=v+(b−a).

For a valid walk v_0,...,v_k with square edge sums q_1,...,q_k,

    v_k = (−1)^k v_0 + Σ_(j=1)^k (−1)^(k−j) q_j.

This follows by iterating v_j=q_j−v_(j−1). A closed even walk must have
zero alternating sum of its square colours. For an odd closed walk,
twice v_0 equals the corresponding alternating sum. These equations
encode return; they do not enforce positive bounded vertices, distinct
vertices, or Hamiltonian coverage. Those are essential additional
conditions. Two consecutive equal colours return immediately to the
previous vertex and cannot occur in a simple path.

The arithmetic maps explain the zigzag appearance of monochromatic
reflections and two-colour translations. They do not make a branch into
a path and do not recover discarded port pairings. The graph's colour
incidence and the transfer signature are complementary coordinates.

## 7. What 15,23,25 do and do not signify

The independently established classifications belong to different
questions: n=15 is the first nontrivial Hamiltonian-path order;
n=23 has three paths after failures at 18–22; n=24 fails again; n=25
repairs that failure through the compulsory ear 11–25–24. The exact
endpoint-aware repair is inherited from the decoder note cited above,
not inferred from the numbers 15,23,25 themselves.

For adjoining a new vertex with d neighbors in an already connected
graph, cycle rank increases by d−1. This simple topological count says
nothing by itself about spanning paths. In particular, adding the
degree-two vertex 24 increases cycle rank but destroys Hamiltonian-path
existence; adding degree-two vertex 25 increases it again and restores
existence. The 3–4–5 tree-block contraction makes the distinction even
sharper: it preserves rank 8 while destroying Hamiltonian feasibility.

## Reproduction and limits

From the workspace root:

    python3 ports/verify_ports.py
    python3 -O ports/verify_ports.py

Both executions agree. The script uses the Python standard library;
checks use explicit exceptions rather than removable assert statements.
`certificate.json` includes the script SHA-256. `verification.out` and
`verification-optimized.out` preserve the identical short receipts.
The representation theorem is proved for all finite graphs under its
stated partition hypotheses; the state counts and path counts are only
the stated finite census. No literature priority claim is made.

Repository instructions, relevant current guardrails, and the targeted
prior notes were read. The coordination status command was unavailable
because `.machine-id` is absent; no identity or repository state was
changed. This worker made no Git mutation. Publication belongs to the
coordinating task.
