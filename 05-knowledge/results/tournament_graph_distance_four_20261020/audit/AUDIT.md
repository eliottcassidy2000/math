# Independent audit of the order-five distance-four identification

2026-10-20. **PASS. No required mathematical correction or scope narrowing.**

This audit covers [the distance-four addendum](../DISTANCE_FOUR_SIX_CUBE.md),
SHA-256 `143f7009085bd69e7dfe594236f4f291b51e53d00fa5831727df332a98f0ce61`.
Its inherited definition is Theorem 3 of the
[trijections and tree-gauges packet](../../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md),
SHA-256 `4a19202298e373dff64cb9ca4b887642c415b560102ff40d34d01735ce0660e7`.
Both producer notes remain unchanged. The complete structural identification is
restricted to order five; it is not an all-order halved-cube claim.

## Verdict

Each 32-vertex component of the order-five cut aggregate is exactly the graph on
even six-bit vectors with adjacency at Hamming distance four. Its positive
aggregate edge weight is 8. It is not the standard halved six-cube, whose
adjacency uses Hamming distance two.

The explicit inverse, six complementary-cycle coordinates, exceptional six-point
S5 action, actual adjacency, spectrum, local graph, and clique obstructions all
pass. In particular the structural argument proves the graph identification;
the note does not merely infer isomorphism from a matching spectrum.

## 1. Coordinate domain and invertibility

At order five every cut has even cardinality. After star normalization, edge
parity is exactly the sum of the six free bits, so the selected component really
is the five-dimensional even subspace displayed in the note.

The six row masks `(13,14,19,22,26,28)` sum to zero over F2, forcing even output.
Direct substitution verifies the displayed inverse on both even spaces. Although
the formula does not use input bit `f`, this is not rank loss on its stated
domain: `f=a+b+c+d+e` there. On the full six-dimensional input space the formula
would have a nontrivial kernel, so the even-domain restriction must remain.
The note retains it consistently.

For a cycle C, cut-cycle orthogonality makes `|x intersect C| mod 2` independent
of the cut representative. Replacing C by its complement changes the value by
`|x| mod 2`; therefore independence of which member of a complementary pair is
chosen holds exactly on the even domain. Both well-definedness statements are
used separately and correctly. The other graph component is identified after a
translation, not by silently extending this unmarked pair-coordinate definition
to odd classes.

The six written cyclic words give exactly the stated row masks and seed
intersection counts `(2,1,1,0,2,2)`. Complementing the seed path therefore gives
output `(1,0,0,1,1,1)`, of weight four.

## 2. Six complementary pairs and six Sylow subgroups

There are twelve undirected Hamiltonian five-cycles. The complement of each is
another five-cycle, giving six complementary pairs. A regular cyclic group of
order five has precisely two edge orbits, the two step classes of the pentagon;
these are a complementary cycle pair. Conversely the rotation subgroup of either
cycle recovers that same cyclic group. This proves a bijection, rather than just
an equality of counts, with the six Sylow-5 subgroups of S5.

The pair action is transitive. Its stabilizer has order 20 and has a unique
subgroup of order five: its number of Sylow-5 subgroups divides 4 and is 1 modulo
5, so must be 1. The order-five subgroup stabilizing one pair cannot stabilize a
different pair, because the two stabilizers would then have the same unique
Sylow subgroup, contradicting the established pair/subgroup bijection. Its
orbits on the other five points have size one or five, and size one has just
been ruled out. Thus it is transitive on those five points and S5 is doubly
transitive on all six.

This is not the intransitive action obtained by fixing a sixth point. The latter
has two orbits on coordinate pairs, of sizes 5 and 10, and cannot justify the
required transitivity on all fifteen weight-four supports.

The equivariance claim concerns the even displacement subspace modulo cuts.
Relabelling cycles may exchange the two members of a complementary pair, but
that has no effect on even displacements. Thus it permutes the six output
coordinates exactly. It does not remove the affine observer correction for
actual tournament-state relabellings. The note explicitly preserves this boundary.

All Hamiltonian paths are relabellings of the seed. Equivariance makes every
image have weight four; double transitivity realizes every possible pair of
zero positions. There are exactly fifteen such outputs. Invertibility identifies
exactly fifteen generators, and transitivity of the 120 ordered Hamiltonian words
makes each generator fibre have size 8. The generator and adjacency statements
therefore follow without spectral classification or an extrapolation from data.

## 3. Actual adjacency and the halved-cube obstruction

Linearity sends the difference of two component vertices to the difference of
their outputs. Since the generator images are exactly the weight-four vectors,
this proves both directions of the claimed adjacency equivalence. It is a graph
isomorphism of the cut-aggregate support, with the scalar edge weight 8 retained.

On the even six-bit space the antipodal operator R satisfies `A4=A2 R=R A2`.
Multiplying an adjacency matrix on one side by a permutation matrix is not
conjugation; this identity therefore does not assert an isomorphism between A2
and A4. The distinction is correctly made in the note.

The spectra are:

- actual support: `15^1, 3^10, (-1)^15, (-5)^6`;
- halved six-cube: `15^1, 5^6, (-1)^15, (-3)^10`.

Multiplication by 8 gives the stated aggregate spectrum. Fourier character
indexing modulo the all-one mask and the multiplicities `1,6,15,10` are correct,
including the pairing of the twenty weight-three masks into ten character
classes.

At a vertex of A4, identifying each neighbour with its missing coordinate pair
makes local adjacency exactly disjointness. The local graph is `KG(6,2)`;
its clique number is 3, and hence the graph clique number is 4. At a vertex of
A2, the local relation is intersection in one coordinate, giving `T(6)` and
local clique number 5, hence graph clique number 6. Translation justifies moving
an arbitrary clique to one containing zero, so these are global maxima.

A further elementary separation from the independent audit is the triangle
count: the actual component has 480 triangles, while the halved six-cube has
640. Their adjacency-cube traces are 2,880 and 3,840, respectively. This is
consistent with both the local graphs and the two spectra.

## 4. Exact replay and independent controls

The producer checker was replayed normally and under Python optimization.
Both outputs match the saved receipts byte-for-byte:

- producer checks: **5,977**;
- checker SHA-256:
  `25d79b7968ee098412822e71f7e1d53fc32f87ebebd0cdcf8f26ec0df1517304`;
- receipt SHA-256:
  `425d0a4a2f2c38d80c757b97fb5842f384a568efc2b79cc129c03eb2b7c3cea6`.

The separate [independent checker](independent_controls.py) imports no producer
code and passes **35,980** live checks in normal and optimized modes, with
byte-identical [receipts](independent_controls.json). Its universe and methods are:

- all 1,024 labelled tournaments and all 120 Hamiltonian words;
- all 7,680 ordered path-complement incidences, with witness uniqueness checked;
- all 64 switching classes and all 4,096 ordered aggregate-vertex pairs;
- the actual weighted aggregate compared entrywise with weight-four adjacency;
- all twelve Hamiltonian cycles generated independently, then paired by
  complement without starting from the producer's coordinate masks;
- explicit order-five rotation subgroups and conjugation under all 120 elements
  of S5, yielding six subgroups with normalizers of order 20;
- the faithful 120-element, doubly transitive coordinate action and its
  equivariance on all 32 even classes;
- every formula/inverse input and exact local adjacency and clique maxima;
- both spectra checked independently using integer adjacency-matrix products,
  an annihilating polynomial with the four asserted distinct roots, and exact
  trace moments, rather than repeating the producer's Fourier calculation.

The actual trace moments in orders 0 through 4 are
`32,0,480,2880,55200`; the halved-cube values are
`32,0,480,3840,55200`. The annihilating polynomials and moments jointly certify
both the eigenvalue supports and multiplicities.

All six producer `SHA256SUMS` entries pass, and the inherited main-note hash is
unchanged. There were no repository or remote writes. These are exact finite
controls supporting an independently audited structural proof, not a formal
proof-assistant certification or a literature-priority claim.

## 5. Portable reproduction and publication scope

Run from the published addendum directory:

    python3 check_distance_four.py > /tmp/distance-four-replay.json
    cmp verification.json /tmp/distance-four-replay.json
    python3 -O check_distance_four.py > /tmp/distance-four-replay-optimized.json
    cmp verification.optimized.json /tmp/distance-four-replay-optimized.json
    python3 audit/independent_controls.py > /tmp/distance-four-independent.json
    cmp audit/independent_controls.json /tmp/distance-four-independent.json
    python3 -O audit/independent_controls.py > /tmp/distance-four-independent-optimized.json
    cmp audit/independent_controls.optimized.json /tmp/distance-four-independent-optimized.json

The publication manifest preserves the producer files byte-for-byte and adds
this audit plus its reproducible evidence. The exploratory `probe.py` is not
required by either checker or proof and is deliberately excluded.

Keep the scope at the order-five aggregate quotient. No assertion here concerns
regular-return icosahedral adjacency, lifting quotient paths to prescribed
tournament representatives, or original/black-metagraph connectivity. No
required edits were found.
