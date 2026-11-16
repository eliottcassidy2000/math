# The Collatz tree under changes of order, label, scale, and geometry

**Scope.** These results concern the actual known first-hit basin of 1,
without assuming that it contains all positive integers. They neither prove
Collatz convergence nor close the Chan–Pak coverage conjecture.

The best invariant depends on what is being changed. The modulo-nine
hexagon sees inverse exits. The intrinsic metric sees a branching tree,
not a flat lattice. The full uncolored rooted shape remembers arithmetic,
even when the printed labels are changed. Choosing how to complete the
tree to a tournament is a separate degree of freedom; that choice can
destroy every triangle, or encode the whole finite tournament spectrum.

## 1. Hexagon versus flat geometry

[HEX.md](HEX.md) recovers the archived mod-nine hexagon and lifts it:
for ordinary Collatz, among six consecutive doublings of a vertex
prime to three, exactly three give odd inverse exits, of which two stay
prime to three and one starts a unary doubling ray. Away from the cut
at 1 these exits have no exceptions. Iterating the two continuing exits
from 8 gives 2^m distinct descendants within tree distance 7m+3 from 1.
This is a lower bound, not an asserted growth rate. Compare
[Applegate–Lagarias on actual preimage counting](https://websites.umich.edu/~lagarias/doc/applegateI.pdf).
No injective map of this tree into the square, triangular, or honeycomb
lattice can have uniformly bounded edge stretch: radius-R balls in each
lattice contain O(R²) vertices. The argument also works against every
uniformly discrete subset of fixed-dimensional Euclidean space. The
hexagon is a local arithmetic clock, not the global metric of the tree.

The geometric contrast is constructive. Label children of any rooted
binary tree by left and right. Starting at [0,1], assign them the outer
thirds of the parent's interval. Embed vertex v at its interval midpoint
x_v and height y_v=3^(−depth(v)) in the hyperbolic upper half-plane.
For an edge, cosh(distance)=11/6. More strongly, for any two vertices,
with tree distance D and hyperbolic distance d,
```
D ln 3 − ln 18 ≤ d ≤ D ln 3 + ln 3.
```
If their last common ancestor has depth l, both midpoints lie in an
interval of length 3^−l. For incomparable vertices their gap is at
least 3^(−l−1), from the deleted middle third; for ancestor pairs let
the height ratio supply the lower bound. Insert these in
cosh d=(y_u²+y_v²+(x_u−x_v)²)/(2 y_u y_v), with
D=depth(u)+depth(v)−2l. This yields 3^D/18≤cosh d≤3·3^D/2 and hence
the inequality. Using halves rather than thirds can put distant
branches arbitrarily close across a shared interval boundary; the gap
is doing actual work.

The unlabeled end space is a Cantor set after deleting a countable
dense set of isolated ends, exactly the rays that eventually enter a
multiple of three. Continued live branching and sterile exits in every
cylinder are both supplied by the same six-clock. [HEX.md](HEX.md)
also gives a literal 3-to-5 generalization: when 2 generates units
mod p², a p(p−1)-clock offers p odd exits, p−1 live and one sterile.
For p=5 this is four live plus one sterile in twenty phases, not a
claim of convergence for 5n+1.

## 2. The other three-to-five: change the equation and retain the flags

[LATTICE.md](LATTICE.md) recovers the three flat regular maps and five
spherical ones, including the old local-balance counterexample. Setting
D=(p−2)(q−2), the flat choices have D=4, hence τ(4)=3 ordered
factorizations. The spherical choices have D=1,2,3 and give
τ(1)+τ(2)+τ(3)=5. Quotienting by duality swaps p and q, so there are
two flat and three spherical types after forgetting the distinction.
There is no canonical bijection replacing three lattices with five solids.

There is also a literal 3→5 transition inside the Collatz inverse
macro-tree, not identified with the Platonic one. Retain two live
exits and one marked sterile exit at each expansion. Replacing one
live frontier vertex by these three changes the total frontier by
two. After m expansions it contains m+1 live and m sterile vertices,
regardless of the choices or batching, so the totals are 1,3,5,7,...
If the live child slots are marked in exponent order, the visible
m-expansion shapes have Catalan count C_m while the sequential
histories number m!: at step j there are j live leaves to expand.
The hook formula below gives the multiplicity of each visible shape.
This distinguishes an exact frontier growth law from the static
curvature classification.

A second literal three-to-five split is in the local *states*. Modulo
three the shortcut inverse types are sterile unary (0), live unary (1),
and binary (2). Keep the first two, but refine the binary state by the
type of its odd inverse child. Residues 2,5,8 modulo nine give
respectively a live unary, sterile unary, or binary odd child.
Every second doubling is the next fork; multiplication by four cycles
2→8→5→2, so every block of three fork sites sees each type once.
Delete the sterile rays and then suppress unbranched stretches:
the resulting full binary tree is highly symmetric, even though the
original rooted basin is rigid. At a live fork x≡2 mod9 the compressed
outgoing lengths are {2,2} or {2,4}; for x≡8 mod9 they are {1,4}
or {3,4}. Doubling twice takes x to the next fork, skipping another
pair of doublings only when that fork is sterile. The odd-child path
starts at a live unary or at a fork respectively, and likewise may
have to bypass one sterile fork. Thus forgetting lengths and sterile
branches is a real loss of information, not just a prettier drawing.

In fact **the lengths alone tell us where to put the sterile rays
back**. A fork with lengths (2,2) needs none. In (2,4) and (1,4) put
an infinite unary ray two steps along the length-four edge. In (3,4)
put that ray on the length-four edge and another one step along the
length-three edge. Expand the edge lengths and reattach the initial
three-edge stem 1→2→4→8 (inverse direction). This recovers the entire
uncolored rooted tree. Thus four possible edge lengths on a full binary
skeleton store enough information to recover the rigid arithmetic
object; dropping those weights leaves the ordinary symmetric tree.

The weighted skeleton has a literal six-phase hexagon too. Always
take the branch containing the double. At a live fork the residues
modulo 27 cycle
```
2 → 8 → 20 → 26 → 11 → 17 → 2
```
with respective edge lengths 2,4,2,4,2,4. Thus one compressed hexagon
traverses eighteen uncompressed shortcut edges, exactly the order of
2 modulo 27. The four fork types along it are (2,2), (3,4), (2,2),
(1,4), (2,4), (1,4). It is a hexagon of live phases after the sterile
sites are skipped, not a flat embedding of the binary tree.

Taking the dual of the flag triangulation gives cubic figures
(4,2p,2q). If all three face sizes at each cubic vertex are distinct,
the two neighbor face sizes must alternate around every face, so every
length is even. Nonnegative curvature leaves exactly (4,6,8),
(4,6,10), (4,6,12) in that case: write the lengths as 2a<2b<2c and
solve 1/a+1/b+1/c≥1. The (4,6,10/12/14) sequence crosses from
positive to flat to negative, in keeping with the earlier
[5–6–7 flag carrier](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/geometry_567_flag_carriers_20261004.md).
The same way to find another pattern: change a rule, specify whether
flags or duality labels survive, then solve the new constraints.
For example every square shell D has exactly one self-dual type, even
arbitrarily far into the hyperbolic regime.

## 3. Erasing labels does not erase the infinite inverse arithmetic

Here use the shortcut F: n/2 on evens and (3n+1)/2 on odds. It has the
same positive first-hit basin. Cut the outgoing edge of 1. Apart from
that cut, inverse children are D(n)=2n and O(n)=(2n−1)/3 when n≡2 mod3.
The archive already proves rigidity for the ordinary functional graph
[THM-4523's source](https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/collatz_functional_uniqueness_20261001.md).
Here is a direct finite-resolution version for the **uncolored shortcut
inverse tree**, for unit n≥3. Cutting 1's outgoing edge affects none of
these subtrees.

A depth-(4r−2) view determines n modulo 3^r. A unary vertex followed
by another unary is type 0 mod3, unary then binary is type 1, and a
binary root is type 2. This is the base. To lift from r to r+1, first
take the unique double if the vertex has type 1. At a binary vertex x,
the double child is live and has type 1. If the other child is sterile,
x≡5 mod9 and doubling twice through the live child gives
4x≡2 mod9, now binary with both children live. Thus at most three
doublings reach a binary vertex x* whose two children are live.
The view identifies the corresponding moves without labels: in the
sterile case the live child is unique and then has one child.
Now read the two unordered child residues modulo 3^r by induction.
The double's residue 2x* is known from the parent already modulo 3^r.
Cancel it from the multiset, leaving the odd inverse's residue.
Since 3O(x*)+1=2x*, that lifts x* to precision 3^(r+1).
Undo the doublings modulo 3^(r+1). Each lift costs at most four
more levels. The first-hit roots 1 and 2 are separately distinguished
by their three and two steps to the first fork; ordinary unit types
above them branch within one step.

It follows that no two different positive units in this tree have
isomorphic uncolored rooted subtrees. At every fork the two child
trees differ, either by sterility or by the two different unit
labels. So the full rooted basin is rigid. Multiples of three have
identical unary futures viewed in isolation, but their labels are
recovered by entering the ray at the uniquely decoded unit parent
and continuing to double. This is infinite reconstruction, not a
finite test of positive basin membership. In fact arbitrarily large
powers of two realize each finite unit residue class, so residue
agreement alone is much weaker than global membership.

## 4. Different orders of building the same tree

A finite first-hit prefix is a finite ancestor ideal. Growing at the
frontier embeds it in a parent-first order; precisely the fair infinite
orders exhaust the whole basin. Distance from 1 with any total tie-break
is one such schedule. If priority is by original integer label at the
frontier, it is also fair: once a vertex is eligible only finitely many
smaller labels can precede it. Choosing to exhaust an infinite chain
before any other branch is not a fair enumeration. For a finite tree
of n nodes, the number of parent-first orders is n!/∏ s_v, where
s_v is the size of the subtree in the finite ideal. This follows by
interleaving child subtree orders at each vertex with the multinomial
coefficient. Contracting prescribed connected, contiguous batches gives
a quotient rooted tree, and its hook formula counts batch orders; the
within-batch orders are separate choices.

For an arbitrary finite sample S, not necessarily ancestor-closed,
the induced forest has one component per vertex whose parent lies
outside S (count the root as missing its parent). Such samples
need not respect the temporal order of the dynamics. For example the
numerical order can introduce an odd vertex before its parent.

Arbitrary relabeling has two readings. Transport a bijection and all
arcs and you have the same rooted tree. Alternatively prescribe which
way each tree edge should compare under *new distinct positive
integer labels*. In this second reading the exact condition is that
every vertex have finitely many predecessors in the oriented comparison
order. Necessity is immediate from the finite number of integers
below a label. For sufficiency, use a fixed enumeration to break ties
and repeatedly choose the first vertex whose predecessors have been
emitted. Induction on the finite ancestral history ensures eligibility,
and once eligible a vertex cannot be starved by the finitely many
earlier entries in the tie-break enumeration. Assign consecutive
labels in this order. For a locally finite tree, König's lemma says the
obstruction is exactly an infinite edge chain directed toward a vertex.
Thus every finite pattern of rises and falls extends (orient unmarked
edges away from the root), but an infinite pattern need not. This
does not say that the new labels satisfy the original Collatz formula.

### A tournament where the schedule has an exact price

For an order ρ of the ideal, orient each tree edge child→parent,
against ρ, and every non-tree pair forward in ρ. All true Collatz
arcs in the ideal survive. Write d_v for the number of children and
let t(ρ) be the resulting directed triangle count. Then
```
t(ρ) = Σ_(p→v as parent relation) (ρ(v)−ρ(p)) − Σ_v binom(d_v,2) − d_root.
```
For a triangle listed in increasing order, the cycle pattern is either
a backward outer pair and two forward adjacent pairs, or two backward
adjacent pairs and a forward outer one. For a parent-child outer pair,
all intermediate vertices count except earlier siblings of the child.
The second pattern counts exactly one per vertex of depth at least two.
This proves the formula. In particular swapping two consecutive
incomparable vertices a,b changes t by d_b−d_a. The triangle count
depends on the sum of tree-edge spans. That sum is also the area
under the frontier-size curve: an edge crosses a cut once for
each unit of its span.

**Exact optimum.** A depth-first order visiting smaller child subtrees
before larger ones minimizes this tournament's triangle count over all
parent-first orders. At a vertex, for every two child subtrees A,B,
the cross-area is at least min(|A|,|B|): if A finishes first, every
step in A has at least one still-open edge in B (an unstarted B has
its entering edge). Separate smallest-first traversal attains the
lower bound for every pair simultaneously. Induct inside subtrees.
Equivalently the minimum is
```
n−1−d_root + Σ_v Σ_{pairs of children x,y} (min(s_x,s_y)−1).
```
The optimum is for a fixed finite ideal, not a globally fair
depth-first traversal of the whole infinite tree. Batching preserves
the formula if spans are measured in physical vertices, not just
batch indices. In a connected-batch order a cross edge's span is
the difference of batch start ranks plus the difference of its two
local attachment ranks. So group sizes alone omit relevant data.

The independent verifier checked all 874 increasing parent maps and
59,603 orders through seven nodes, plus four actual finite shortcut
Collatz ideals. At 80 vertices breadth-first gave 882 triangles,
against the exact minimum 168. Same triangle count need not mean
same Hamilton count: the parent vector (−1,0,0,1,1) with orders
(0,1,2,3,4) and (0,1,3,2,4) gives four triangles in both and
Hamilton counts 11 and 13. That adjacent swap tests precisely
why optimizing a local summary is not a Chan–Pak proof.

### Two other completions: trivial spectrum and universal spectrum

1. Compare every pair globally by depth and a total tie-break and
   orient from later to earlier. This contains the basin arcs and is
   transitive: every finite restriction has one Hamilton path.
2. Keep deeper→shallower, but allow arbitrary orientations on
   same-level ties. One can fix a **single** such ambient tournament
   on the basin whose finite ancestor ideals realize the entire
   finite tournament Hamilton-count spectrum. Exponential inverse
   growth gives arbitrarily wide levels arbitrarily deep, in the
   descendants of any live vertex. Enumerate all finite tournaments
   U_j. Choose a live vertex y_j deeper than all previously reserved
   levels, and put a copy of U_j on |U_j| descendants at one common
   level below y_j. Orient only those reserved tie edges exceptionally;
   all other ties follow a fixed total order. In the closure of this
   copy under ancestors, every older exceptional level contains at
   most one vertex, because those levels precede y_j. No future
   exceptional level occurs there, and lower unreserved levels are
   transitive. The closure is therefore the ordinal composition
   of the copy of U_j with a transitive lower part. It has exactly
   H(U_j) Hamilton paths. Every induced ideal is itself a finite
   tournament, giving the reverse spectrum inclusion. These witnessing ideals are not required
   to form one nested growth sequence. Earlier tie gadgets could add
   unintended multiplicative factors on a nested sequence.

Thus the full Collatz basin embeds in completions ranging from H=1 on
every finite ideal to the entire odd-spectrum question. Universality
transfers the problem rather than solving it: the construction inserts
an arbitrary U_j and does not choose it from a prescribed count. The
pattern across the geometry and tournament proofs is not that one
number causes the other; it is that a local picture only determines
the global object when the missing boundary data are supplied.

## 5. The literal local comparison on integer labels

There is a fourth completion, with no basin-membership test needed:
use the shortcut Collatz arcs (cut 1→2) on all positive integers and
compare every other pair from larger to smaller. On [1,N] the only
exceptions to numerical descent are 2j+1→3j+2, with
1≤j≤J=floor((N−2)/3) for N≥2. Each vertex has at most one such
upward edge entering and one leaving. A directed triangle in
increasing numerical order either uses an upward outer edge, with
both interior pairs descending, or uses two consecutive upward edges
and a descending outer edge. An up-arc 2j+1→3j+2 has j possible
middle vertices, all legal since it is the unique up-arc leaving its
source or entering its target. Two consecutive upward steps have
terminal values 8,17,26,... . Therefore
```
directed triangles on [1,N] = J(J+1)/2 + floor((N+1)/9).
```
For J≥1 the unique nonsingleton strongly connected component is
[3,3J+2]. Each up-arc makes its spanning interval strongly connected
by the above triangles, and consecutive intervals overlap; away from
these intervals all arcs point down. Thus H is constant on each
size block N=3J+2,3J+3,3J+4. It increases strictly at the next
block: the newest top vertex y=3j+2 beats all old vertices except
x=2j+1, and the old arc x+1→x is a positive even-path insertion
contribution to the Hamilton count. No comparison of trajectory
depths enters this calculation. It demonstrates a natural batch size
three even though other orderings of the same arcs can erase all
directed triangles.

## Replays

- `python3 check_tree_orders.py`: all small schedules and literal triangle/HP controls.
- `python3 verify_hex.py`: residue and macro-branch checks.
- `python3 check_unordered.py`: sampled finite shape controls for the depth estimate.
- `python3 check_numeric.py`: direct checks of triangles, strongly connected
  components and the first Hamilton-count plateaus for the literal completion.
- `python3 check_lattice.py`: Schläfli shells, macro branching and Catalan frontier checks.
- `python3 check_weights.py`: literal checks of all four length and sterile-ray reconstruction cases.
