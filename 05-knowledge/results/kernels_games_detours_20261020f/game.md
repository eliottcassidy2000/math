# The complement game: all four variants and sharp edge budgets

**PROVED:** all four feasibility classifications and the minimum number of
red edges in each winnable variant. One sharp lower bound uses a **CITED**
classical theorem of Erdős. **FINITE-EXACT:** independent labelled census
through six vertices, with a second path test through five. These are
undirected graphs, not tournaments or tournament-metagraph lines.

The user's clarification is load-bearing: Player 2 finds a path among the
uncoloured edges; the path need not use every uncoloured edge. Its start
can be chosen. No Eulerian condition belongs in this game.

## Definitions and complete answer

Let G be Player 1's red spanning subgraph of K_n, and B its complement.
Write H for existence of a Hamiltonian path and W for existence of a
spanning walk, with vertex and edge repetitions allowed. W is exactly
connectivity. A one-vertex graph satisfies both. A win in variant X/Y
means that G fails X and B fails Y.

For every n>=4 the answer is:

| Forbidden in red / required in blue | Win exists? | Minimum red edges |
|---|---|---:|
| H / H | yes | n-1 |
| H / W | yes | n-1 |
| W / H | yes | b(n), below |
| W / W | no | -- |

For n<=3 no variant has a win. If n=2m or 2m+1 respectively,

    b(2m)   = min(4m-5, binom(m+1,2)),
    b(2m+1) = min(4m-3, binom(m+2,2)-1).                 (1)

For n=4,...,12 these are `3,5,6,9,10,13,15,17,19`.
Equivalently b(n)=2n-5 except b(6)=6 and b(8)=10. For even n=2m,
the difference between the two candidates is `(m-2)(m-5)/2`; for odd
n=2m+1 it is `(m-2)(m-3)/2`, proving this simplification.
In particular the asymmetry between H/W and W/H matters: complementing
exchanges their winning positions, but not their number of red edges.

## Elementary feasibility and the first two budgets

The red star K_(1,n-1) has at least three leaves when n>=4, so it has no
Hamiltonian path. Its complement is a clique plus an isolated vertex.
This wins H/H and H/W with n-1 red edges. Reversing the colours gives a
W/H win. For n<=3, a disconnected graph's complement is traceable, and
a connected graph on at most three vertices is traceable; this excludes
all wins. More generally the complement of a disconnected graph is
connected: vertices in different components are adjacent in the complement,
and vertices in the same component have a common neighbour in another one.
This excludes W/W at every order.

To disconnect K_n costs at least n-1 deletions: a nontrivial cut of size s
contains s(n-s)>=n-1 edges. This proves the H/W lower bound.

The H/H lower bound follows from the elementary path-closure criterion:
if every nonadjacent u,v in an n-vertex graph satisfies
`d(u)+d(v)>=n-1`, the graph has a Hamiltonian path. For completeness,
adjoin a universal vertex, so the degree threshold becomes N=n+1.
Adding a nonedge uv at degree sum at least N cannot create a first
Hamiltonian cycle. If that cycle used uv, its deletion gives an ordering
`u=w1,...,wN=v`. Among the N-1 cuts, the sets
`{i:u~w_(i+1)}` and `{i:w_i~v}` have total size at least N and therefore
intersect. The corresponding cross-edges close a cycle without uv.
Add all nonedges and reverse this argument. Remove the universal vertex.

Now if red has m<=n-2 edges, for each red edge uv,
`d_G(u)+d_G(v)<=m+1<=n-1`; hence blue's degree sum is at least n-1.
Blue is traceable. Thus at least n-1 red edges are necessary even without
the restriction on red's own paths.

## The W/H budget: a separator, not an Eulerian obstruction

**CITED input.** Erdős's sharp minimum-degree edge bound says that a
non-Hamiltonian graph on N vertices with minimum degree at least k,
`1<=k<N/2`, has at most

    max_(k<=s<=floor((N-1)/2)) [binom(N-s,2)+s^2]          (2)

edges. The quadratic is convex, so only the endpoints are needed.
See P. Erdős, [Remarks on a Paper of Pósa (1962), pp.227–228](https://users.renyi.hu/~p_erdos/1962-17.pdf).
The modern restatement in Liu–Ning–Wang,
[Theorem 1.1](https://arxiv.org/html/2604.01068v1), agrees; no new theorem
from that preprint is needed here.

In W/H, G is disconnected, so B is connected and has minimum degree at
least one. Adjoin a universal vertex to B. The new graph J has N=n+1
vertices and minimum degree at least two. A Hamiltonian cycle in J is
equivalent to a Hamiltonian path in B; thus J is non-Hamiltonian.
Apply (2) with k=2. Its endpoints are s=2 and s=floor(n/2), and

    e(G)=binom(n+1,2)-e(J).

Subtracting the two endpoint bounds gives exactly the two quantities in
(1). This proves the lower bound under the weaker assumption that B is
connected, hence certainly for the game's disconnected-red condition.

Both endpoints are attained within the game. For
`1<=t<=floor(n/2)-1`, take

    B_(n,t)=K_t join (I_(t+1) disjoint-union K_(n-2t-1)).  (3)

Deleting the t universal vertices leaves t+2 nonempty components. A
Hamiltonian path with t vertices removed leaves at most t+1 path pieces,
so B_(n,t) is not traceable. Its complement has t isolated vertices and
is disconnected. Its red edge count is

    r(n,t)=binom(t+1,2)+(t+1)(n-2t-1).                   (4)

This is concave in t, so the minimum over the allowed interval occurs
at t=1 or t=floor(n/2)-1. Those values are (1). At the first endpoint
red is one isolate plus two adjacent hubs joined to all n-3 other
vertices; at the balanced endpoint red is a clique on m+1 vertices
(even n), or a clique on m+2 vertices with one edge removed (odd n),
together with isolates. These are constructions, not a uniqueness
classification of all minimizers.

## What transfers to the tournament questions

The operative invariant is the separator deficit
`components(B-S)-|S|-1`, which is positive in (3). It obstructs a spanning
path while permitting connectivity. It is a necessary-condition witness,
not a characterization of all nontraceable graphs. A directed tournament
always has a Hamiltonian path; neither this complement game nor its
separator deficit directly supplies a black metagraph component invariant.

The relevant operation to borrow is **retain the boundary incidences before
contracting a piece**. Component count alone loses which ports can be
entered and exited. The same port loss is visible in the undirected
line-graph hostile in [walks_knots.md](walks_knots.md), and in the
endpoint-frozen strong tournaments of [detours.md](detours.md).

Reproduce from repository root:

    python 05-knowledge/results/kernels_games_detours_20261020f/game_check.py

The universe is every labelled simple graph on n=1,...,6; there are no
inherited filters. Subset dynamic programming tests Hamiltonian paths,
literal permutations independently agree through n=5, and breadth-first
reachability tests walks. Complements exchange H/W and W/H counts.
Constructions (3) are checked through n=100, with direct path tests through
n=12. Positive controls include complete graphs and the winning stars;
the hostile W/W control is tested on the entire census. Assertions do not
disappear under optimized Python.

## Incoming independent work: connected colours and exact path-cover joins

Commit `9852737522`, integrated during close-out, independently obtained
the same game budgets. Its [reviewed packet](../complement_path_games_20261020/AUDIT_AND_PROVENANCE.md)
adds prescribed-start quantifiers, a self-contained Pósa-type degree
obstruction, and an exact join theorem using minimum path-cover numbers.
The present packet retains its independently derived free-start proof.

It also answers the next hostile question: **both colours can be connected
and nontraceable at every n>=7**. Take a blue clique C on n-2 vertices,
delete one edge hv, and attach two blue leaves a,b to h. Deleting h leaves
three blue components. In red, h is adjacent only to v; a,b are adjacent
and connect to every other clique vertex. Red is connected, but deleting
{a,b} leaves the edge hv and n-4 isolates, hence n-3>3 components. In each
colour this exceeds the number of pieces a Hamiltonian path could leave
after those deletions. The construction uses 2n-4 red edges, without any
claim here that this is the minimum under the connected-both restriction.
Our census excludes connected-both wins through six: its H/H win count
is exactly the sum of the disjoint H/W and W/H win counts. Thus seven is
the first order, with the exclusion explicitly finite-exact.

The new family is independently constructed in `game_check.py` through
n=100; direct path tests cover n=7,...,12. The more general incoming join
criterion retains the path-cover number of each block, and its counting
version retains the entire oriented path-cover profile. This supplies an
actual map to the maximal-run substitution mechanism in THM-3121, cited
by its full slug in that packet. A scalar path count loses this interface,
just as a gap-only defect kernel loses external matching partners.
