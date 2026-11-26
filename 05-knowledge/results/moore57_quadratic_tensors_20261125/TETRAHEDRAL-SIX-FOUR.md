# Six and four: Wilson collision flags, the tetrahedron, and the first ten primes

2026-11-26. **PROVED** finite structural connections and a sharp tetrahedral
action bound. **FINITE-EXACT** replays. Extends the
[Wilson chart action](WILSON-ACTION-AND-RIGID-MOORE-CHARTS.md).

## 1. The numbers count different things

The Wilson action is 6T+4S for a complete rooted sheet chart, where T and S
count triangles and squares on leaves. Compare a direct and a two-step route:
a triangle collision is based at one of its three unordered vertex pairs.
Compare two two-step routes: a square collision is based at one of its two
opposite pairs, *not* one of the four sides. Expanding the squared counting
defect doubles each collision. Equivalently there are six **ordered**
triangle endpoint flags and four **ordered opposite** square endpoint flags.

A four-vertex tournament has six edge slots and four three-vertex faces
before choosing orientations. These are not the same objects as the two
kinds of Wilson flags. One localized identification is natural: on a
chosen square on four vertices, the ordered opposite pair (x,x') maps to
the tetrahedral face omitting x, equivariantly under the square symmetries.
But a chosen triangle's six ordered pairs cannot be identified
equivariantly with the six tetrahedral edges: under permutations of the
triangle's three sites the ordered pairs form one free orbit, while the
tetrahedral edges split into two orbits of three (in the triangle and
incident with its omitted vertex). This explains exactly where the
equal cardinalities do and do not give the same structure.

The [tetrahedral Haar work](../../../07-reflections/tetrahedral-k4-haar-xor-cycle-bridge-codex-20260816.md)
already pairs the six edge slots into three opposite matchings and takes
their three complementary squares. The [simplicial Wilson note](../session34_gauge_representations_20261023/gauge/SIMPLICIAL_WILSON_FANO.md)
records the distinct relation among four signed triangular plaquettes.
For actual directed tournaments the four score types give:

| Sorted outdegrees | Directed triangles | Hamiltonian paths | Labelled orientations |
|---|---:|---:|---:|
| 0,1,2,3 | 0 | 1 | 24 |
| 0,2,2,2 | 1 | 3 | 8 |
| 1,1,1,3 | 1 | 3 | 8 |
| 1,1,2,2 | 2 | 5 | 24 |

In particular H=1+2 times the number of directed triangular faces here.
The factor two resembles the Wilson collision doubling, but here the
tournament count uses only directed triangles; Wilson defects of the sheet
graph count undirected loops. The replay checks all 64 orientations.

## 2. A sharp finite Wilson result on the six-link tetrahedron

Set k=4 in the rooted chart, so there are exactly four branches of three
sheets. There are six base links, four base triangle faces and three base
squares. Gauge-fix the three links out of one branch to the identity,
and call the other three permutations a,b,c on links 12,13,23, respectively.
Write f(g) for the number of fixed sheets; write products in traversal order.
The rooted count is

    T = f(a)+f(b)+f(c)+f(a c b^-1),
    S = f(a c)+f(b c^-1)+f(a^-1 b).

Each term counts lifts of the corresponding simple base loop, and all three-
and four-cycles upstairs project to simple base loops by the matching
property. For these actual integer counts one has the sharp inequality

    2T + S >= 6.

**Proof.** A transposition in S3 fixes one sheet, a three-cycle none,
and identity three. If one or two of a,b,c are odd, two of the four
triangle holonomies and two of the three square holonomies are odd;
T>=2 and S>=2. If all three are odd, all four triangles are odd.
If a,b,c are all even, T is a multiple of three, so only T=0 needs
checking. In that case write their rotation exponents x,y,z as ±1
modulo three. If x=z and y=-x then x+z-y=0 gives a triangle fixed point,
contradicting T=0. If x=y=z then both y-z and y-x vanish, giving six
square lifts. If x=-z then x+z vanishes and y equals one of x,z, so
again at least two square holonomies are identity. This proves the
inequality in every case.

The exact undominated pairs (T,S) are (0,6), (2,2), and (3,0).
They are realized, for example, by (a,b,c) equal to a common nontrivial
three-cycle; to ((12),(01),(021)); and to (identity,(012),(021)),
respectively. The parity cases above exclude T=1, and the inequality
excludes any point below the three stated ones. Thus the positive weighted
minimum for α T+β S is min(6β,3α), and all three listed points are
optimal at α=2β. The Wilson weights 6 and 4 give **exact minimum 18**.
If instead you insist on no triangles, you must pay at least six squares
and action 24. In particular this four-branch, three-sheet chart cannot
have zero action; no larger-degree chart is ruled out by this calculation.

All 6³=216 rooted-gauge triples were replayed. The eight minimizers all
have exactly one identity face holonomy and no square fixed sheets.
There are two choices at each bad face; global conjugation interchanges
them by inverting the three-cycles. Thus there are **four** minimum
classes with branch labels fixed, literally indexed by the bad face.
This is a structural appearance of four within the same Wilson model.

## 3. A literal six-edge plus four-face Moore graph

Distinguish one point ∞ in a five-element set. Use the ten two-subsets as
vertices, with disjointness adjacency. The six pairs not containing ∞ are
the tetrahedron's edges and the four pairs {∞,i} are its faces, labelled
by the omitted vertex i. The graph is Petersen, with the following exact
6+4 structure:

- Two faces never share a graph edge.
- An edge-vertex is adjacent to each of its two containing faces.
- An edge-vertex is adjacent to its one opposite edge-vertex.

So deleting the three opposite-edge links leaves the subdivision of the
tetrahedron's dual K4. The graph is Moore: disjoint pairs have no common
neighbor, and overlapping distinct pairs have exactly one, the remaining
two-subset. This is a genuine way to turn the ten slots into a graph
without relying on the number ten alone.

## 4. Why the prime 6+4 split cannot be identified with the Petersen cut

The [prior prime selector](../session322_target_inverse_prime_channels_20261113/PRIME_SELECTION_4_PLUS_6.md)
uses the square-sum graph on the first ten primes. It is exactly

    (3—13—23—2—7—29)  ⊔  (5—11)  ⊔  (17—19).

The marked component through 2 has size six. Its complement has size
four but is not connected. This is not the first-six cutoff, the residue
bipartition (5+5), or the Fi22/Baby/Monster layer split (6+3+1).
Under any bijection that sends the marked prime component to Petersen's
edge layer and its complement to the face layer, at most **three of
the seven square-sum edges** can survive: only the three opposite edges
are available inside the six-layer and none inside the four-layer.
Three is attained by taking the unique perfect matching of the six-site
path onto the three opposites. There are exactly 3!·2³·4!=1,152 such
maximizing assignments. If the layer condition is dropped, all seven
edges *can* survive as a spanning subgraph (the replay supplies a witness).
Thus it is the marked 6+4 identification, not ordinary graph embedding,
that causes the failure. The usable comparison is that one must carry
the incidences as well as the counts.

Replay in this directory: node tetrahedral_six_four.js.
