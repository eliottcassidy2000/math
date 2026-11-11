# From the original three colors to exact permutohedral and arithmetic maps

**Session257, 11 November 2026. COMPLETE research deliverable.** New proofs
in notes and finite exact experiments, with inherited results at their recorded
scope. **Chan–Pak remains OPEN.** This is a creative synthesis of actual
carriers, not a broad verification pass or a literature-priority claim.

## 1. The original objects, recovered rather than guessed

The requested references have three distinct meanings.

- **The Q25 three-color decomposition** groups edges by the role of the larger
  endpoint when it was born: traceable births15,16,17,23,25; failed
  double-neighbor births13,14,19–22,24; and the remaining scaffolding. It is
  neither the three old Q12 components nor the square labels4,9,16.
- **4→2→5→3→7→13=4** displays the small-prime valuations of |M|/|B|, with
  the explicit identification4↔13 closing the chain. It is not a valid
  polynomial recurrence modulo9. We retain the construction and find an
  actual geometric operation nearby, rather than infer one from the labels.
- **{17,19,23,31,47} versus {29,41,59,71}** are B minus Fi22 and M minus B.
  Their inherited meanings include modular levels1|2|6 and fixed-point
  regimes of the Monster action on2A involutions. The stabilizer is2.B,
  not B itself.

The [source atlas](history/README.md) pins the historical definitions,
original color list, later corrections and eight operation maps. It recovers
Sessions32/40/210/211 as well as Session256's triple/tree constructors.
Inherited color-polynomial and port-transfer results are credited, not
presented as rediscoveries.

## 2. One weighted polytope genuinely connects path blocks and group orbits

For path pieces of sizes m_1,...,m_k, concatenate the oriented pieces in a
chosen order. If L_i vertices precede piece i, its centered position sum is

    y_i = m_i L_i.

The convex hull of all block orders is exactly the weighted permutohedron

    Z(m) = sum_(i<j) m_i m_j [e_i,e_j].

Named first and last vertices impose the real square-sum seam test.
Geometry alone forgets those endpoints and the orientation of each piece.

Now take any finite group acting on X, with orbit sizes m_i. In the fixed
section of the ordinary permutohedron, let r_i be the common coordinate on
orbit i. The explicit affine change

    y_i = m_i (r_i − (m_i+1)/2)

identifies that fixed section with the SAME Z(m). This is equality of the
geometric carriers under a written map, not just a matching vertex count.

Coarsening the orbit partition is Reynolds averaging. Cutting path pieces
is refinement. Their coordinate maps agree exactly: sum the fine y values
in each coarse block and subtract its internal pair-weight constant. Thus
group averaging and path-piece refinement fit one commuting geometric
diagram, while retaining different graph/port interpretations.

For Q12, its actual reversal group C2^3 has seven orbits of sizes
1,1,2,2,2,2,2, hence a six-dimensional fixed section. Averaging further over
S3×S4×S5 gives the weighted(3,4,5) hexagon. This larger group is NOT a graph
automorphism group: it loses endpoint incidence. The two fixed vertices3,7
are midpoints, not endpoints. Exactly this marking matters for gluing.

[Block geometry and refinement](permutohedron/PROOF.md),
[fixed-section identification and orbit averaging](filtration/PROOF.md).

## 3. Q25 reveals which geometric moves are actually useful

All ten unoriented Q25 Hamiltonian paths share14 edges. Contracting their
common linear forest gives11 path blocks and a ten-dimensional weighted
permutohedron, with orientation sheets. A direct vertex search and an
independent contracted seam search give the same ten paths.

The new finite result is striking: **all ten feasible chart states are
isolated under adjacent block swaps and individual block reversals.** Yet
one-edge surgery has two components of sizes6 and4, and two-edge surgery
connects all ten paths.

The proof-level explanation works for every graph. If two paths replace
r edges, their common forest has exactly r+1 pieces; in that pair-adapted
chart the surgery is a reorder/reorientation of those pieces. A small graph
move can be a long signed interval move in a fixed fine permutohedron.
This suggests choosing the chart dynamically rather than searching only
geometric neighbors in one canonical chart.

The original color specialization is retained exactly:

    7 x^7 y^11 z^6 + 2 x^7 y^10 z^7 + x^6 y^12 z^6.

It is a specialization of the inherited edge-marked polynomial. Seven paths
share its first monomial, so the complete paths, endpoints and orientations
remain in the certificate. Original birth colors, square-target colors and
component colors remain separate throughout.

A further finite construction retains the original colors in a genuine
navigation process. Choose another one of the ten paths uniformly and accept
it only if surgery distance is at most d, otherwise stay. For d=1 and d=2,
each kernel has a minimal color-respecting **nine-state** Markov quotient,
but the merged pair is different. A controller allowed to choose either
operation needs all **ten** states. The complete positive quotient matrix and
stationary weights are recorded. Thus even the smallest sufficient positive
state depends on the operations one intends to perform.

[Proof, full ten-path certificate and surgery data](permutohedron/PROOF.md).

## 4. The global square-sum permutohedron differs from the tournament one

An ordering of the vertices of an undirected graph G induces an acyclic
orientation of G. The fiber of that orientation is its set of linear
extensions, connected by adjacent swaps of NONedges.

**New explicit reformulation:** Hamiltonian words of G are exactly singleton
fibers, hence isolated vertices of the nonedge-swap subgraph of the
permutohedron. In the graphical zonotope they are precisely the normal
cones containing one braid chamber, each with probability1/n! under an
exchangeable continuous ordering. Q25 has20 such cones, paired by reversal.

This complements Session210's tournament construction, where directed
Hamiltonian words are SINKS of a coherent orientation on the permutohedron.
The underlying permutation carrier is shared; the operation is different.
One must not orient missing square edges artificially and claim the count
survives.

The necessary geometric information is also sharp. P4 and the four-vertex
star have the same chromatic polynomial and both graphical zonotopes are
combinatorial cubes. Their braid-chamber fiber sizes are respectively

    P4: two1, four3, two5;
    star: six2, two6.

Only P4 has singleton fibers. The abstract polytope loses the alignment
with the labelled braid fan that detects Hamiltonicity.

[Singleton-chamber theorem, controls and coherent-sink comparison](chambers/PROOF.md).

## 5. Primitive square chains now decode into actual tournament sources

Take a primitive triple s²+t²=u² with s<t. The inherited three-reflection
construction makes a path on1,...,t²−1 using targets s²,t²,u². Let a,b,c
count edges with those three ORDERED target colors. Then the new decoder is

    epsilon = b−a−c in {0,1},
    s² = 2a+1+epsilon,       t² = 2b+2−epsilon.

So the abstract target-color census recovers the primitive triple, without
knowing the original arithmetic vertex labels. Its exact membership test
constructs a chain with that census; it does not certify an arbitrary
submitted color itinerary, and an explicit equal-census hostile records
that loss.

Recover the Euclidean spinor(m,n), put x=m−n,y=n, and form the actual
Session256 tournament C3[TT_x,TT_y,K1]. Its count is

    H = 2+(2^x−1)(2^y−1).

This gives a commuting primitive-triple map from square reflection chains
to literal tournament sources and their Berggren branches. Rebuilding a
child carrier is exact but is not a bounded local graph edit. The primitive
triple is preserved; Hamiltonian count is not. The sorted family still
misses class1 modulo3.

That distinction is necessary, not excessive caution: Session211 already
proves no map can satisfy H(T)=2h(Q_n)+1 for all square-sum graphs, since
h(Q23)=3 and h(Q25)=10 would force the forbidden tournament counts7 and21.
The present construction takes a different, fully specified invariant.

[Decoder, reconstruction boundary and actual graph compiler](gram/PROOF.md).

## 6. A positive seam interface, and its exact first failure

For a port-labelled graph G let h_G({a,b}) count spanning paths with those
named endpoints. Glue two copies of graphs along matching port names and
count Hamiltonian cycles using exactly two seam edges. The answer is

    C2(G,H) = sum_{a<b} h_G({a,b}) h_H({a,b}).

This is a positive Gram kernel. With k ports its rank is sharply binom(k,2),
and every basis probe is realized by the SAME primitive square-reflection
chain with different port labellings. A finite actual graph probe bank
therefore recovers the full response, paralleling Session256's fixed-observer
construction with a different interface.

At four seams, connectedness changes the form. The three perfect matchings
on four ports give the Hamiltonian-cycle matrix J−I, with eigenvalues2,−1,−1.
All three sources are port-renamings of the same induced square block
{1,3,7,9}. This is an actual negative kernel inside the square-block category,
not a formal vector counterexample. The assembled seam edges need not
remain square sums, which is stated explicitly.

Thus Session256's tournament Gram positivity transfers to a precise
restricted operation, not to arbitrary path gluing. Endpoint pairs,
component pairings and seam count are the missing distinctions.

[Sharp probe bank and four-seam boundary](gram/PROOF.md).

## 7. The Monster references acquire two precise geometric operations

For the action M/(2.B), a subgroup P fixes a port iff P lies in a conjugate
of2.B. A fixed point becomes a marked singleton orbit in the fixed
permutohedron. This turns the inherited prime shells into port-depth regimes:

- {17,19,23,31,47}: full Sylow depth admits a fixed port;
- {2,3,5,7,11,13}: some proper p-subgroups do, but full Sylow groups do not;
- {29,41,59,71}: no nontrivial p-subgroup does.

Actual marks, not just valuations, determine the fixed geometry. An explicit
S4 example has identical group/stabilizer valuations but2 versus0 fixed
ports, proving that limitation.

The original valuation vector on(2,3,5,7,11,13) is

    v=(5,7,3,4,1,2).

Correcting B to its actual double-cover stabilizer changes only one coordinate:

    w=(4,7,3,4,1,2).

Along this segment the valuation permutohedron hits exactly the equal-weight
wall for prime coordinates{2,7}. **Exactly360 edges contract, reducing720
vertices to360.** Under the original4↔13 identification, the five-state
cycle simultaneously becomes a two-cycle with a three-step tail, changing
its operator polynomial from z^5−1 to z^3(z²−1). These are two exact effects
of one normalization, not a claimed equivalence of their invariants.

The inherited Monster-local AGL(1,7) quotient also acts on the genuine
Paley7 Hamilton subset of the permutohedron: nonsquare affine maps require
position reversal. All42 signed maps preserve its189 Hamilton words.
The orientation sign, like the endpoint mark, cannot be discarded.

[Fixed faces, normalization contraction and signed Paley action](filtration/PROOF.md).

## 8. Incoming work changes the constructive frontier

Concurrent main `f322658908c42e677cfc19826b51805fbf34dbf9` proves an extension
content theorem: if C_k(T) is the gcd of counts of ALL k-vertex tournament
extensions, then C_(k+1) divides C_k and C_(n−1)=1 for every n-vertex T.
The top coefficient of an alternating old/new path is±2, while all tournament
counts are odd. Fixed k requires only polynomially many attachment samples
for this gcd, although counting each sample retains its own cost.

This is inherited concurrent work, not a new theorem of this packet. It
sharpens how our port/mark obstructions should be used: finite extension
symmetry can obstruct a given budget, but cannot imply a permanent common
prime divisor across unrestricted additions. Gcd1 still supplies neither
an exact target nor a return to a selected source grammar.

Session211's apex flip likewise preserves old sinks and creates exactly its
ordered two-path-cover count. That is a genuine construction; the square-sum
birth polynomial uses NAMED endpoints in DIFFERENT components. Their cover
interfaces are related but not interchangeable, and no H-preserving map has
been claimed. These distinctions point toward the next repair problem:
control actual source, boundary response and target simultaneously.

[Incoming proof pin](history/sources/INCOMING_EXTENSION_CONTENT.md),
[operation comparison and next questions](OPEN_QUESTIONS.md).

## 9. What is delivered and what remains open

Four new-claim engines replay in ordinary and optimized mode with explicit
checks. They cover all stated small universes, including all ten Q25 paths,
all5040 orbit orders in the averaging experiment, all720 valuation orders,
31 primitive spinors,64 independent two-seam cycle counts, and complete
small-graph chamber fibers. Historical large-group and global square-sum
inputs are inherited, not re-proved by those runs.

The session supplies a common exact geometric map, concrete graph-operation
charts, a finite sharp response bank, a primitive-triple compiler, and a
normalization contraction. It does not solve Chan–Pak, make the original
valuation chain a canonical arithmetic dynamics, or identify Monster
symmetries with square-sum adjacency.

[Open constructive gates](OPEN_QUESTIONS.md), [concept board](BOARD.md),
[replay](REPLAY.json), [optimized replay](REPLAY_OPTIMIZED.json),
[manifest](MANIFEST.json).
