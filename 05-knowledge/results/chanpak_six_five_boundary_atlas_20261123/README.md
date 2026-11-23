# One distinguished place and five exits

**PROVED scoped / FINITE-EXACT.** This packet separates exact transports
(tree, polygon, finite field, curve, tournament observation) from analogies
whose maps forget a relevant structure. Chan–Pak remains open.

## The original six, in their natural order

In half-index \(k=(H-1)/2\) the centered digits \(H\mapsto3H-2,3H,3H+2\)
are ordinary ternary digits \(0,1,2\). The two forbidden values are
\(7\leftrightarrow10_3\) and \(21\leftrightarrow101_3\).
Delete those adjacent vertices. Inside the \(10_3\) subtree the exits in
lexicographic order are
\[
 100,\ 1010,\ 1011,\ 1012,\ 102
 \quad\longleftrightarrow\quad 19,\ 61,\ 63,\ 65,\ 23.
\]
The sixth component is the exterior, rooted at 1. The side marked 1
below represents the whole exterior, not the immediate predecessor of
7. The [inherited boundary formula](../six_roots_hidden_bits_20261103/BOUNDARIES.md)
for \(d\) deleted non-root vertices in \(c\) connected groups is
\(1+2d+c\) roots; each deleted vertex has three child stubs, and
\(d-c\) stubs go to another deleted vertex. Two disconnected holes
would make seven, not six. The old uncentered seed set
\(1,23,25,63,65,67\) has the same boundary size but attaches the
second hole at the left child instead of the middle, breaking the
rooted mirror. The old memo's binary/ternary junction \(7\to21\to63\)
is stronger than sharing a cardinality: \(63=2^6-1\) doubles 7's
Mersenne exponent, and the next factor \(65=2^6+1\) is its complement.

The relative exit words are \(0,10,11,12,2\) and their uniform ternary
masses are \(1/3,1/9,1/9,1/9,1/3\). They form a complete prefix code.
If \(I_w\) is the ordinary ternary subinterval for a relative word
\(w\), then \(H_w/3^{|w|}=6+2\,\operatorname{mid}I_w\), by the
centered digit recurrence. For every full half-index digit-length shell
of length \(L\ge4\), the five components together account for exactly
\(3^{L-2}/(2\cdot3^{L-1})=1/6\) of the shell. The old component takes
the other \(5/6\). This is not ordinary natural density: the covered
blocks \([3\cdot3^m,4\cdot3^m)\) give liminf \(1/6\) and limsup
\(3/8\) for the cumulative exposed fraction. Logarithmic density is
\(\log(4/3)/\log3\). These are address statistics, not estimates for
physically realized tournament counts.

## From two quadrilaterals to a hexagon, then a pentagon

Turn each removed node into a quadrilateral, with one incoming side
and three ordered child sides. Glue along the deleted parent-child
edge. The resulting hexagon has cyclic side labels
\[
 1,\ 19,\ 61,\ 63,\ 65,\ 23 .
\]
Reflection fixing the marked side exchanges \(19\leftrightarrow23\)
and \(61\leftrightarrow65\), and fixes 1 and 63. The long seam separates
\(\{1,19,23\}\) from \(\{61,63,65\}\). A connected hole tree with
\(d\) vertices gives a quadrangulated \(2d+2\)-gon with \(2d+1\)
non-incoming sides. This is also the exact edge-count reason a planar
quartic interaction tree with two vertices has six external legs
([quartic Stokes polytopes](https://arxiv.org/abs/1811.05904)).
There is no identification of amplitudes with path counts.

**An exact six-to-five quotient of fillings.** Number the hexagon
vertices \(0,\dots,5\) with marked side \(50\). Contract that side.
The three hexagon quadrangulations, cut by \(03,14,25\), induce
one triangle and one quadrilateral in the pentagon. Their marks
become \(03,14,02\). Triangulate the quadrilateral in either way.
The resulting six **marked fillings** give all five pentagon
triangulations, with exactly one duplicate on forgetting the mark:
\[
 (02,03),\ (02,24),\ (03,13),\ (13,14),\ (14,24).
\]
Only \(02\) and \(03\) are compatible with each other among the
three marks; the two unmarked diagonals \(13,24\) cross.
Each pentagon triangulation is the fan at a unique vertex.
Associate each remaining side to the fan at its opposite vertex.
This gives an explicit, reflection-compatible bijection between the
five exits and the five triangulations. The duplicate fan is at
the contracted vertex, opposite the side labelled **63**.
This is the [Stasheff pentagon](https://arxiv.org/abs/math/0212126)
on the five side labels, but the flips here are diagonal flips.
Indeed the [actual tournament fork](../chanpak_collatz_edge_buckets_20261123/README.md)
has no single arc from 63 to 65. The six boundary sides and the
six marked fillings are *not* the same six-set: the rooted
reflection has two fixed sides but no fixed marked filling.

There is another mismatch of structure. Refining the two ternary
nodes into binary bracketings produces four hexagon triangulations
with the chosen long seam. In each, under fair binary branch choices,
the middle child of the root has mass \(1/4\), never \(1/3\).
Even mixing these refinements leaves error at least \(1/12\) in
total variation, by testing the event of reaching the middle
child. The symmetric mixture achieves the bound: the two outer
masses are \(3/8\), and the middle three sum to \(1/4\) with
each at most \(1/9\). A combinatorial refinement therefore does
not automatically preserve the continuation law.

The same cut is realized in the 4-regular Cayley tree of the free
group on two generators: remove the two vertices of an edge, and
there are six components; root the edge and the other five are
forward directions. A topological dual glues two four-holed spheres
along one cuff to give a six-holed sphere. Its five other peripheral
loops freely generate the fundamental group, and their cyclic product
determines the marked peripheral loop. The common mechanism is gluing
two four-port objects once, not a universal equivalence of groups,
curves, and tournament cores.

## What travels farther

- [Arithmetic branch labels](CURVE.md) are intrinsically
  \(\infty,\pm1,0,\pm1/3\), not ordinary counts mod 5. They give
  all six points of \(\mathbf P^1(\mathbf F_5)\) and a genus-two
  double cover. A signed five-cycle operator is **literally
  Frobenius** on its four nonconstant modes. The marked 1 then
  singles out five of the 15 two-torsion labels; the rest form
  Petersen.
- [Apex transfer in tournaments](APEX-TRANSFER.md) has exact linear
  context ranks **7, 6, 5**, for unrestricted words, first guard O, or both guards O and I
  respectively. The signed part remembers each open end and resets at
  a nonempty intervening block. This is a physical example of what
  an extra port can mean; it is not an identification of six roots
  with six independent generators.
- The old [six-output memo](../six_output_return_20261103/SIX_AND_CATALAN.md)
  supplies other genuinely different sixes: \(3!\) permutations of
  three path slots and a lattice index \(6^t\) from diagonal factors
  \(2^t,3^t\); neither follows from the deleted boundary.
- The [unit-cube hexagon](../hadwiger_six_seven_eight_20261025/metric_spindle/RESULT.md)
  has a true regular six-cycle plus a central point; its distinguished
  object is not a side. Choosing and contracting a rim side does not
  preserve the unit-distance wheel. The [six-colour interface](../six_colour_interface_profile_20261015f.md)
  instead has a genuine palette budget: five exterior colours with
  just one meeting the ports can still allow six-colour extension.
- The [square-sum tripartition](../session147_square_sum_20261103/tripartitions/RESULTS.md)
  has components on 4, 5, 6 vertices with 3, 4, 5 edges, by the tree
  identity. Its size-five branching component costs an extra path
  cover while the size-six component is itself a path. This is a
  useful hostile to counting a component's ports without its internal
  continuation. The [prior square-sum islands](../six_roots_hidden_bits_20261103/BOUNDARIES.md)
  are path orders, not ternary-address cylinders.
- In Lonely Runner the [rank-five](../../../01-canon/theorems/THM-3415-zero-mode-cochain-global-rank-five-support.md)
  and [rank-six](../../../01-canon/theorems/THM-3416-zero-mode-cochain-global-rank-six-support.md)
  phase strata are different modular supports, analogous to the
  apex ranks only in the careful sense that a retained mode is a
  statement about a specified projection. There is no transfer of
  feasibility claims.

Reproduce the finite checks with **python3 replay.py**. The next
tournament problem is to find an operation on the five context reads
that changes a prescribed ternary digit and returns the same marked
ends. The polygon and curve reveal the right distinctions (orientation,
position, weighting), not that physical operation.
