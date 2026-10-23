# Arbitrary regular-five attachments: one owner per connected support

23 October 2026. **PROVED:** the prime-rotation attachment lemma, exact
shared-chart potential criterion, and one-owner propagation rule below.
**FINITE-EXACT:** all 4,320 regular nonuniform five-by-five blocks, all
97,336 declared triangle relation systems, and the actual tournament/path
controls. These are named-cell statements. They do not reconstruct missing
blocks from anonymous motif counts or classify arbitrary metagraph wiring.

**Stronger final result:** an explicit intrinsic observer gives the exact
adaptive external-owner budget **zero or one per connected five-cell
support**. A rigid block is detected by owned row/column pair flags or the
inherited displacement word. Otherwise a nontrivial cycle gain selects an
actual intrinsic fixed owner. Only a genuinely balanced C5 component needs
one external actual owner. See section4A; this is more constructive than
knowing abstractly that the full chart group is trivial.

The [prime-refinement extension](PRIME-REFINEMENT.md) proves that every
prime-size tournament whose full automorphism group is a regular C_p has
only indiscrete/discrete equitable partitions. Thus the one-owner
ordinary-refinement propagation rule extends to all such prime cells.

## Inheritance and the new boundary

The closest mechanism is the partial-bijection chart gluing theorem in
[odd-cell gluing](../../odd_cell_gluing_20261023/README.md). Its sufficient
binary-block test required distinct rows and columns. The
[five-cell repair](../../session40_marked_ports_square_sum_20261023/FIVE_CELL_MOD8_REPAIR.md)
classifies matching/complement-matching attachments, while the
[prime-cell necklace](../../session40_marked_ports_square_sum_20261023/PRIME_CELL_BURNSIDE_FORMULAS.md)
retains their full rotation-double-orbit. The new step handles **all regular
nonuniform blocks**, including degree two/three and repeated rows/columns.

The canonical hostile is still independent choices of charts around a
cycle: every pair can be isomorphic while the whole system is not. The
corrected near miss is to use rank or distinct rows as a necessary condition
for chart transport. The least-used sidecar is the action of the actual
internal automorphism group, rather than all permutations of a cell.

Concept board: actual internal rotations; arbitrary regular cross blocks;
shared affine chart potentials; ordinary equitable refinement; literal
Hamiltonian-path owners. The connection map is stated explicitly: an
actual binary block maps to its chart-compatibility relation; this preserves
all cell-preserving isomorphisms between the two supplied systems. A bare
relation does not recover the block, and a family of unaligned relation
witnesses does not recover a shared isomorphism.

## 1. Prime-rotation kernel lemma

Suppose a named source cell and corresponding target cell each have internal
automorphism group exactly C_p, acting regularly, with p prime. Choose an
internal isomorphism as reference, so its possible charts form a C_p torsor.
Let the block to another named cell be equitable and nonuniform. If two
compatible chart pairs agree on the other cell, their quotient is an
automorphism of the C_p cell preserving each individual block row/column.
A nonidentity element of C_p is transitive. Hence every such row/column
would be constant. Equitability then forces the entire block to be uniform,
contrary to the hypothesis. Thus the projection to the other chart is
injective.

If **both** cells have regular automorphism group C_p, both projections are
injective. The complete compatibility relation is a partial bijection even
when the binary block has repeated rows or columns. This removes the
distinct-row/column prerequisite in this carrier. Equitability, nonuniformity,
prime order, and the stated internal automorphism group are load-bearing.
For example a nonuniform block made of all-zero and all-one rows has a
vertical C_p kernel if row equitability is dropped. Uniform blocks permit
independent charts and must not be treated as functional edges.

Let B be a p-by-p block in these cyclic coordinates. Its rotation
stabilizer is a subgroup of C_p x C_p with injective projections. It is
therefore either trivial or a line

    H_B = {(a,s a): a in F_p},  s in F_p^*.

The latter case occurs exactly when, for a nonempty proper subset S of F_p,

    B(i,j) = 1 if and only if j-s i belongs to S.                 (1)

Indeed the line's action preserves j-s i; each of its p orbits has that
form. Conversely (1) has the displayed stabilizer, and injectivity rules
out a larger one. Its row/column degree is |S|. The slope and subset are
properties of the actual named block, not a claim that every permutation
of p symbols is affine.

For two supplied blocks B and B', the chart relation is empty if they are
not equivalent. Otherwise it is a coset of H_B: either one pair of chart
rotations, or exactly p pairs satisfying

    a_j - s_ij a_i = t_ij.                                     (2)

All offsets retain the same chosen reference chart at a shared cell.

## 2. An all-size compatible attachment rule

Consider any number of named C_p cells, with corresponding internal
isotypes supplied on the target. Join cells whenever their source block
is nonuniform; require the same pattern of uniform/nonuniform target
blocks, and check all uniform block values as well. In a connected support
component, choose one root chart. Every nonuniform edge forces at most one
neighbor chart, so there are at most p global candidate isomorphisms.
Checking all edges, rather than only a spanning tree, is necessary and
sufficient for a candidate to be a genuine cell-preserving isomorphism.

The complete answer is **zero, one, or p** isomorphisms. To see this without
a count heuristic, a nonempty isomorphism set is a torsor for the source's
cell-preserving automorphism group. Restriction to one root C_p is faithful;
its image is a subgroup of C_p and therefore has order one or p.

There is an exact potential implementation. On edges of type (2), propagate
along a tree and express every chart as a_i=alpha_i x+beta_i, alpha_i nonzero.
Each closing edge supplies one linear equation in the SAME root x. A
singleton edge pins both of its incident chart values; this also becomes
a root equation after propagation, or can initialize a new propagation
segment directly. Contradictory equations reject the comparison. A
nonzero equation coefficient fixes x uniquely. If no equation fixes x and
all constant equations vanish, all p root choices survive.

For self-isomorphisms, offsets vanish. A singleton edge rigidifies the
whole component. If all edges have line stabilizers, a cycle whose gain
product is not one also rigidifies it. Otherwise there is exactly one
shared C_p symmetry; adding **one actual owner anywhere in the component**
kills it. This is sharp: a balanced component genuinely has p automorphisms
before that owner is added. The cost is not one mark per cell. Independent
components joined only by uniform blocks can have independent C_p factors;
their marking costs add unless another retained constraint couples them.

This describes an automorphism base / actual chart decoder, not a claim
that ordinary refinement discovers every rigid component without a mark.
An unbalanced three-cell example below is rigid while its named five-cells
remain equitable. Parent-level residual-cycle transport potentials must
still be conditioned on the same compatible chart choice; unioning
allowed cells across incompatible choices is invalid.

## 3. Exact regular-five census and a repeated-row repair

Every regular tournament on five vertices is R5, with i->j when j-i is
1 or 2 modulo5; its internal automorphism group is exactly the rotations.
These inherited facts are checked independently for the concrete R5 in
the script. The complete block census is:

| Cross degree | Labelled blocks | C5-stabilized | Rigid | Rotation-double-orbits | Repeated-row blocks |
|---:|---:|---:|---:|---:|---:|
| 1 | 120 | 20 | 100 | 8 | 0 |
| 2 | 2040 | 40 | 2000 | 88 | 600 |
| 3 | 2040 | 40 | 2000 | 88 | 600 |
| 4 | 120 | 20 | 100 | 8 | 0 |
| total | 4320 | 120 | 4200 | 192 | 1200 |

The 192 orbits consist of 24 with C5 stabilizer and 168 rigid ones. Each
nontrivial-stabilizer orbit has five labelled members; each rigid orbit
has25. The degree-two symmetric blocks are exactly the 4*binom(5,2)=40
patterns (1), not every degree-two regular block.

Here is a literal degree-two rigid block, given by the selected target
owners in its five rows:

    {0,1}, {0,1}, {2,3}, {3,4}, {2,4}.                         (3)

It has repeated rows and columns. Its bipartite support is K_(2,2) disjoint
union a six-cycle. Every degree-two repeated-row block has this structure:
the duplicated two-entry rows saturate their two columns, leaving a
three-by-three degree-two block. There are binom(5,2)^2*3! = 600 of them.
All600 have trivial rotation stabilizer. One way to see this immediately
is that the unique repeated-row pair cannot be preserved by a nonidentity
five-cycle. Complements give the other600 degree-three examples.

Thus repeated binary rows do not imply multiple compatible chart choices.
The internal tournament forbids arbitrary row permutations. In the
two-cell tournament, cross degree2 gives distinct total outdegrees4 and5,
so the two cells in this example are intrinsic even without externally
naming them. This supplies a concrete order-ten rigid attachment that the
old distinct-row test did not certify.

## 4. One mark also suffices for ordinary refinement

R5 has no proper nondiscrete equitable partition. Every equitable cell in
a tournament has odd size: its internal outdegree is constant, hence
(size-1)/2 is an integer. The only possible proper nondiscrete partition
of five vertices is3+1+1. A singleton in that partition has all or no arcs
toward the three-cell, and at most one arc toward the other singleton, so
its internal outdegree is0,1,3,or4, never the required2. Contradiction.

Consequently individualizing any one actual vertex splits its R5 cell
completely under ordinary directed equitable refinement. Once a neighboring
cell is discrete, the nonuniform block gives at least one distinguishing
row on the target R5. Its stable restriction cannot stay a whole cell and
therefore also becomes discrete. Induction along connected support proves:

**For any connected network of named regular five-cells with nonuniform
equitable blocks along its support, one actual mark in any cell makes
ordinary directed equitable refinement discrete on the entire network.**

The script checks all five possible root marks in every one of the4320
two-cell blocks, and every actual mark in four three-cell controls:21660
refinement runs. It does not claim a constant-round bound independent of
the support diameter. An isolated R5 needs its own one owner.

## 4A. Intrinsic flags pay for every rigid five-cell component

The one-owner upper bound can be made adaptive without asking for that
owner when the attachment data already name one intrinsically.

### Degree two/three: two actual pair flags

For degree two use the selected target pair in each row; for degree three
use its two-element complement. Give the row owner one bit according to
whether this pair is an edge of the target's intrinsic short pentagon.
Equivalently the pair has one common outneighbor, rather than zero, within
its R5. This is a cell-relative rooted four-vertex incidence predicate.
Give each column owner the corresponding bit for its selected source pair.
These are actual-owner colours, not two independent aggregate histograms.

**Claim.** Both row and column flags are constant if and only if the block
has a nontrivial C5 rotation stabilizer.

The reverse direction follows from transitivity. For the forward direction,
suppose every row pair is a short-pentagon edge, or every row pair is a
long-pentagon edge. Write m_i for the multiplicity of each of that pentagon's
five edges among the selected row sets. Column degree two says the two
incident multiplicities at every vertex sum to2. Alternation around an
odd cycle forces every m_i=1. Thus the row sets enumerate a pentagon's
edges exactly once, with no duplicates.

The two rows meeting at a target vertex are the corresponding column
pair. Constant column type says their adjacency graph is also one of the
source's two pentagons. Hence the row assignment is an isomorphism between
two five-cycles. Such an isomorphism sends cyclic edge indices to
s i+c, with s nonzero modulo5. The block therefore has form(1), and its
C5 stabilizer is real. Complements preserve the argument.

For each of degrees two and three the exact flag table is:

| Stabilizer | Row flag | Column flag | Blocks |
|---|---|---|---:|
| trivial | nonconstant | nonconstant | 1600 |
| trivial | nonconstant | constant | 200 |
| trivial | constant | nonconstant | 200 |
| C5 | constant | constant | 40 |

Every asymmetric degree-two/three block thus produces a proper intrinsic
colouring in at least one cell. Section4 makes that cell discrete, and
then propagates through its whole connected support with **zero added
external owners**. The checker confirms this on all4000 such labelled
blocks. Looking in only one direction misses200 blocks at each degree;
both endpoints of the shared incidence matter.

### Degree one/four: root the inherited displacement carrier

For a matching (or complement matching) let q be its actual selected
owner bijection. The source R5 has an intrinsic directed short successor
c. At owner i retain the directed short-cycle displacement from q(i) to
q(c(i)) in the target, an element of {1,2,3,4}. This is the rooted form of
the incoming complete displacement necklace. Its colours are constant
exactly when q is affine. Every rigid matching block therefore gives a
proper intrinsic source colouring and the same zero-owner refinement
propagation. All200 rigid labelled degree-one/four blocks are checked.

### When all blocks are symmetric: an actual cycle selects the owner

Every remaining block has form(1). Select its target singleton if its
degree is one/four, or select the unique vertex beating the other member
of its selected two-set if its degree is two/three. This gives an actual
bijection f_ij between the cells. It is affine, with exactly the same
nonzero slope s_ij as the block's rotation-stabilizer line.

Orient each support edge once and retain this choice. On backward traversal
use **f_ij inverse**, not a newly selected winner in the transpose block.
A minimal hostile is B(i,j)=1[j-i in {0,1}]. Its forward winner map is
the identity, but independently selecting the reverse winner gives j->j-1.
The resulting out-and-back walk falsely acquires a nonzero translation.
The exact inverse gives the required identity backtrack.
A closed actual transport then has form x->A x+B. If A differs from1, its
unique fixed point is an intrinsic actual owner of the root cell. Marking
that derived owner is not an additional externally supplied owner; it
is computed from retained block incidences and the chosen named cycle.
Section4 again makes the entire component discrete.

If every cycle has gainA=1, the chart group really has order5, so one
external actual owner is both sufficient and necessary for pointwise
decoding. Crucially, a cycle with A=1 and B nonzero has **no parallel owner
section**, but still retains all five chart automorphisms. It must not be
declared rigid or invalid. The inherited order31 systems with empty owner
sections are precisely why this predicate distinction is retained.

Combining the cases gives an effective exact criterion: first try the
owned local flags on every nonuniform block; if all blocks are symmetric,
retain their actual transports and test cycle gains; use an external owner
only if every gain is balanced. This realizes the abstract0/1 group-base
budget using concrete owned records and ordinary refinement.

All8000 actual affine-matching triangles are checked independently of the
abstract compatibility alphabet. There are6000 with a nontrivial gain and
a unique intrinsic cycle-fixed owner,400 balanced flat triangles with five
owner sections, and1600 balanced nonflat triangles with **no** owner section
but five genuine chart automorphisms.

The selector is orientation-sensitive. Under converse the two-set winner
and loser exchange. Retain both alternatives, or carry the selection bit
to the converse record; do not silently recompute a winner and claim that
its fixed owner must be the original owner. The set of possible selected
transports is converse-closed, the slope test is unchanged, and already
chosen actual ports/owners are transported literally. This is an explicit
sidecar, not an extra external owner or a replacement for native mod8 guards.

## 5. A local success can fail globally

Take three R5 cells A,B,C. On AB and BC use matching j=i. On AC compare
the following cases:

- Matching j=i: the whole tournament retains the common five rotations.
- Matching j=2i: the equations a_B=a_A, a_C=a_B, a_C=2a_A force identity.
- Matching j=i+1: compared to the first system, every pairwise attachment
  has five valid local charts, but the three chart equations are globally
  inconsistent. The two full tournaments are not isomorphic.

Their cell outdegrees are4,7,10, so any ordinary isomorphism must preserve
the cells. The third claim is not merely an artifact of imposed names.
The first two examples show the exact point where the cycle gain pays
for an otherwise necessary actual owner. The first versus third shows why
one must retain the translation part as well as the gain.

All97,336 triangles of the relation alphabet consisting of20 affine
five-point relations,25 singleton pairs, and the empty relation are
checked against an independent direct relational join. They give:

| Global chart choices | Triangle systems |
|---:|---:|
| 0 | 83311 |
| 1 | 13625 |
| 5 | 400 |

This abstract alphabet is used only to test the gluing rule. It is not a
claim that independently selected relations always come from one fixed
source-target tournament pair.

## 6. Converse, ports, and scope

Negating coordinates restores R5 after taking converse. The cross block
becomes B^c(i,j)=1-B(-i,-j), and chart rotations become their negatives.
The checker retains this correspondence in all actual network fixtures.
It does not fold the two charts of a block independently or discard the
offset in (2). Existing mod8 incidence guards remain valid additional
coordinates; their residue alone is not substituted for the full relation.

For every checked network, an actual directed Hamiltonian path is built.
Every surviving chart map is checked entrywise on the tournament, and its
transported path satisfies f(K_P(T))=K_(f(P))(f(T)). This is a path transport
statement, not preservation of an unmarked path's original endpoints. Any
already retained actual endpoints restrict the chart domains further and
must be honored. Adding a global source and sink gives fixed external
endpoints without altering the internal rotation analysis.

The prime-rotation lemma excludes P7 with automorphism group C7 semidirect
C3 and all other internal groups not equal to regular C_p. It does apply
to the cyclic R_p family in the incoming prime-cell note. Arbitrary
nonfunctional attachments remain the parent packet's next boundary.

## Reproduction

    python -B 05-knowledge/results/odd_cell_completion_20261023/five/check.py
    python -O -B 05-knowledge/results/odd_cell_completion_20261023/five/check.py

`--write` intentionally regenerates receipt.json. The finite universe,
actual blocks, paths, and hostile chart relations are retained there.
All checks use explicit exceptions, so optimization cannot disable them.
