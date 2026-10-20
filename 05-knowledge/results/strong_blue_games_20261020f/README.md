# Strong blue repairs, complementary traversal games, and boundary forcing

**Current scope (2026-10-20):** all-order black-core connectivity remains
**OPEN**. This packet proves a general two-block black detour and closes an
infinite strong-to-strong family that defeats the arc-disjoint square method.
It also settles the four stated undirected games, with Player 2 choosing
the start, and proves boundary-colouring extension/stress criteria.

## Closest inheritance and the remaining anchor

The [preceding packet](../spectrum_basins_20261020e/README.md) proves black
detours at every nonstrong mixed endpoint, connects the active {1,C3,R4}
SCC words across feedback inventories, and joins the neighbours of pure-blue
vertices. Its exact remaining obligation is the strong-to-strong mixed blue
case. The [concept board](BOARD.md) keeps that anchor separate from the game
and geometric forcing lanes. No mixed-to-mixed blue bridge is found here.

## Proved mechanisms

| Object | Result and scope | Proof / checks |
|---|---|---|
|Blue path P=(L,a,R)|If (L,R,a) is directed, both reversed blocks have HPs with the original endpoints, and d(first L)+d(a) is neither n-1 nor n, three black moves reach the exact labelled partner|[Two-block repair](detours.md), detours_family.py/json|
|Arc-disjoint P,Q|Transported P remains blue iff the two positional reversals coincide; gives black-square and partner-joining criteria|[Square law](detours.md), detours_square.py/json|
|T_k, two return arcs in TT_(2k+1)|Every marked standard blue edge has a black detour of length <=3; its unique arc-disjoint HP is blue, so squares alone fail for every k|[Detours](detours.md), [hostile proof](double_return_counts.md)|
|T_(a,b), two transitive blocks|H=(2^(a-1)+1)(2^(b-1)+1); for a=b=k, H=q^2, blue count q, deficit q(q-1)|[Bijection and primitive triples](double_return_counts.md), double_return_counts.py/json|
|Open path partners|Unordered boundary/cut pairs give common black hubs, allowing the roles of the two edges to exchange|[Boundary transport](transport.md), boundary_transport_probe.py/json|
|Interior HP relocation|Two directed triangle reversals sharing one vertex; preserves the entire score vector|[Cycle identity](transport.md), transport_identity.py/json|
|Undirected complement game|Three variants involving H admit a winning selection iff n>=4; W/W never does|[Game proof](game.md), game.py/json|
|Complete joins|Traceability iff each block's path-cover number is at most the number of vertices outside it plus one|[Game, section 3](game.md)|
|Pinned geometric exterior|Independent pinned rigidity rows imply every proper four-colour boundary extends; a minimal six-colour failure requires at least h+ceil(a/2) pinned stresses|[Six-colour extension](six_colour.md), six_colour_check.py/json|
|Minimum palette blocker banks|For k>=3 only concentrated and common-star banks attain the union-bound threshold; distinct unit-distance ports exclude the common-star case|[Palette classification](six_colour.md)|
|Marked directed p-by-q torus grid|Euler length 2pq, winding(q,p); the straight parametrization with this winding is embedded iff gcd(p,q)=1|[Typed torus map](torus.md), torus.py/json|

H means a Hamiltonian path; W means a vertex-spanning walk allowing repeated
vertices and edges. A spanning walk is not an Euler trail. The torus grid
adds an orientation, embedding, and marking; its Euler circuit is not itself
an embedded knot. The marked winding, not a scalar path count, carries p,q.

## Cross-lane conclusions and boundaries

1. **Sidecars become constructive.** The black detour needs internal endpoint
   paths; a join needs internal path covers; geometric forcing needs the
   actual port relation; a torus construction needs marked winding. Counts
   alone erase the data consumed by the next operation.
2. **A failed method need not expose an obstruction.** The two-return family
   has no black arc-disjoint replacement, yet a different three-step repair
   works at every order. Likewise a six-colour precolouring can fail on a
   bipartite star: failure for one boundary assignment is not a six-chromatic
   graph. Both examples are explicit and checked.
3. **The square has a combinatorial origin.** The T_k path set is a q-by-q
   array with a blue diagonal. Its primitive triples are
   (q,(q^2-1)/2,(q^2+1)/2), beginning (3,4,5). This arithmetic map does not
   preserve Petersen-minor or three-edge-colouring obstructions.

The earlier D_2(z)=z^4 theorem, fixed-count rational series, and Fibonacci
matching constructions retain their prior scope. The exact raw HP forbidden
set {7,21} remains **OPEN**; these new explicit path-count families do not
prove that spectrum conjecture. The torus-grid correspondence does not
establish a game-winning predicate or a tournament black-detour criterion.

## Audit and reproduction

Each proof was reviewed against its necessary endpoint, symmetry, and
geometry hypotheses. The root reviewed the game, torus, detours and
six-colour proofs; independent agents reviewed the square counts, cycle
transport, pinned rank proof and palette equality cases. Checkers use
explicit failures, so validation also runs under Python -O.

The all-order statements rest on the proofs. Bounded controls include:
all standard-blue presentations through order seven for the square law;
literal family certificates through order201; the 1,252 nonempty graph
atlas classes through order seven for the game; 737 proper boundary
assignments for the rigidity extension; and all 32,768 six-colour pair
families for the extremal cross-intersection inequality. Their universes
and narrower independent checks are stated beside each script.

See [AUDIT.md](AUDIT.md) for final commands, source hashes, correction
mechanisms, and exact remaining obligations.
