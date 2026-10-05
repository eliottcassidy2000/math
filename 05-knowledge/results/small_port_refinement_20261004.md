# Three ports, six charts, and lossless recursive refinement

2026-10-04. **PROVED:** the pure-family chart-word decoders, compensated port
composition, typed mixed partitions, coarse triangle strata, and guarded
coordinate transport. **REFUTED:** extending single-matrix history recovery
or automatic conforming refinement to the mixed fan/bisection family.
**FINITE-EXACT:** the declared controls. The subdivision itself and the
source-aware monitor are inherited objects; no novelty claim is made for
them. No universal controller payment or Collatz completion follows.

Artifacts: [script](../../04-computation/experiments/small_port_refinement_20261004.py)
and [output](small_port_refinement_20261004.out).

## 1. Recovery: what the small numbers actually count

The owner-specified cumulative map is the invertible shear
\[
 C(x,y,z)=(x,x+y,x+y+z).
\]
The [cumulative-coordinate note](edgewise_cumulative_coordinates_20261004.md)
proves that edgewise adjacency requires a common sign in the cumulative
difference. Componentwise absolute values admit the crossing chord
\((0,2,0)\)--\((1,0,1)\) at degree two. It also gives an exact point decoder
with positive local weights. The
[mixed-refinement note](mixed_triangle_refinement_20261004.md) retains
parent-cell labels and warns that barycentric and edgewise refinement do
not commute. Those results are the closest coordinate mechanisms here.

The [noble flag note](noble_subdivision_flags_20261004.md) supplies the
rank/parent sidecar: an untyped refinement can acquire extra symmetries.
The [four-state tournament codec](tournament_recursive_four_state_20261004.md)
supplies a second inheritance: a finite alphabet can encode an unbounded
ordered word, but forgetting its presentation can destroy a future
operation. No tournament is imposed on tied coordinate values here.

The live board is **ordered coordinates / child charts / port gauges /
boundary faces / original source / exact guards**. The corrected near miss
is treating equal cell counts or equal determinants as equivalent
constructions. The cheapest hostile below is an uncompensated port swap.

For a triangle \(\Delta=\{\lambda_i\ge0:\sum_i\lambda_i=1\}\), one barycentric
subdivision has three original corners, three edge midpoints, and one face
center. Its six triangles are indexed by flags
\[
 \{a\}\subset\{a,b\}\subset\{a,b,c\},
 \qquad (a,b,c)\in S_3 .
\]
Thus three coordinates, six permutations, and seven vertices belong to
one precise object. These numbers count different types. The seven-vertex
graph is a six-cycle with a center joined to every boundary vertex.
Its abstract automorphism group has order twelve; preserving the vertex,
edge, and face ranks leaves order six. Equivalently, the outer cycle must
preserve its alternating corner/midpoint classes. The larger graph group
contains exchanges of a corner and a midpoint, which cannot be affine
symmetries of the original triangle.

This is barycentric subdivision, with counts \((7,12,6)\), not the
degree-two edgewise triangle with counts \((6,9,4)\). The common-sign
edgewise repair remains necessary and is not replaced by these charts.

## 2. A finite alphabet of lossless affine charts

For \(\sigma=(a,b,c)\), let \(B_\sigma\) have the three ordered columns
\[
 e_a,\qquad (e_a+e_b)/2,\qquad (1,1,1)/3 .
\]
It sends local barycentric coordinates \(u\in\Delta\) to the flag triangle
in the chamber \(\lambda_a\ge\lambda_b\ge\lambda_c\). Its inverse is
\[
 B_\sigma^{-1}\lambda
   =\bigl(\lambda_a-\lambda_b,\,
          2(\lambda_b-\lambda_c),\,3\lambda_c\bigr).
 \tag{1}
\]
The inverse coordinates are nonnegative exactly in that chamber and sum
to one. The determinant is \(\operatorname{sgn}(\sigma)/6\). These are
invertible affine embeddings, not six-valued encodings of arbitrary points.
Their rational coefficients retain unbounded precision under iteration.

For a finite chart word \(w=\sigma_1\cdots\sigma_d\), read the outside chart
first and define
\[
 B_w=B_{\sigma_1}\cdots B_{\sigma_d}.
 \tag{2}
\]
Then \(B_{uv}=B_uB_v\). The product maps the parent triangle to the selected
depth-\(d\) barycentric child with its three terminal ports.

**Lossless word theorem.** The exact matrix \(B_w\) determines the entire
finite word \(w\). Moreover, if the terminal ports have been permuted, the
matrix \(B_wP\) determines both \(w\) and the terminal permutation \(P\).

Proof: the absolute determinant is \(6^{-d}\), so it recovers \(d\).
Apply the matrix to the centroid \(c=(1,1,1)/3\). The remaining inner
product sends \(c\) into the strict interior of \(\Delta\). Formula (1)
therefore shows that the three coordinates of \(B_wc\) are strictly ordered,
and their order is exactly \(\sigma_1\). Left-multiply by its inverse and
repeat. After \(d\) steps the remainder is \(I\), or the advertised port
permutation. Since \(Pc=c\), terminal port order does not disturb any
centroid test. This also proves distinct chart words have distinct cells
at each depth: their first differing flag has disjoint relative interior.
No numerical tolerance or general graph-isomorphism procedure is needed.

The word is a discrete refinement history retained in the matrix. A point
alone generally forgets that history: all six charts send their local
third vertex to the same parent center. Even the local coordinate
\((0,0,1)\) is identical in this example, so it does not replace the chart
label. A point can instead be assigned its unique smallest containing
face, as in the inherited cumulative-coordinate decoder; that face
deliberately does not select an arbitrary incident top-dimensional chart.
Neither convention makes ties into oriented pairwise comparisons.

## 3. Interfaces must transport their port frame

If an intermediate triangle is relabelled by a permutation matrix \(P\),
then the interface law is
\[
 (B_uP)(P^{-1}B_v)=B_uB_v.
 \tag{3}
\]
The first factor has the same physical image triangle, but its local port
names changed. Compensating in the second factor preserves the composite.
Dropping \(P^{-1}\) generally changes the finer cell and the represented
point. For example, take \(u=v=(0,1,2)\) and interchange the first two
ports: \(B_uPB_v\ne B_uB_v\).

Thus the unordered image-vertex set in the fixed, labelled parent
coordinate system is sufficient to identify the refinement word.
An abstract triangle alone is not: the ordered original corners supply
the ambient coordinate gauge. Composition using local coordinates still
needs the terminal port correspondence. These are different reconstruction
questions. A directed
acyclic expression can share repeated subwords; flattening to the exact
matrix retains their expanded word, while forgetting the expression
forgets the chosen sharing/factorization. No compression bound for
unrelated words follows.

## 4. Three fan children give a complete adaptive tree carrier

There is also a genuine three-child refinement, not a selection of three
of the six flag triangles. For \(c\in\{0,1,2\}\), put
\(a=c+1,\ b=c+2\) modulo three and define the centroid-fan chart
\[
 A_c=[e_a,\ e_b,\ (1,1,1)/3].
\]
Its image is the triangle where coordinate \(c\) is minimal, and
\[
 A_c^{-1}\lambda
   =(\lambda_a-\lambda_c,\lambda_b-\lambda_c,3\lambda_c),
 \qquad \det A_c=1/3 .
\]
The three images cover the original triangle with disjoint interiors.
Each consists of two barycentric flag cells, distinguished by whether
\(\lambda_a\ge\lambda_b\) or the reverse. Algebraically, if
\[
 E_0=[e_0,(e_0+e_1)/2,e_2],\qquad
 E_1=[e_1,(e_0+e_1)/2,e_2],
\]
then \(A_cE_0=B_{(a,b,c)}\) and \(A_cE_1=B_{(b,a,c)}\).
Thus the three-cell fan, six flag charts, and seven barycentric vertices
have exact refinement maps between them; they are not identical objects.

The earlier decoder works again with base three: the determinant gives
word depth, the unique minimum of the image centroid gives its first slot
\(c\), and multiplication by \(A_c^{-1}\) strips it. A terminal port
permutation is recovered at the end. This gives a lossless address for
every leaf in a finite adaptive fan refinement.

Represent a leaf by the empty tree and an internal node by its ordered
triple of children. Its three children use the charts \(A_0,A_1,A_2\).
The resulting triangles form a conforming subdivision: inserting a
centroid and joining it to its three corners leaves every old boundary
edge whole. Adjacent leaves therefore acquire no hanging boundary nodes.
After \(m\) internal-node insertions,
\[
 (V,E,F)=(m+3,3m+3,2m+1).
\]
Each insertion adds one vertex, three edges, and two faces, proving the
formula. The exact leaf matrices recover their addresses; their
prefix-closed trie recovers the full ordered tree. The fixed original
corner gauge and the terminal port frames remain part of this statement.
The tree does not record which of two independent leaves was refined
first in an implementation.

Now consider the abstract token grammar
\[
 S=\varepsilon\quad\hbox{or}\quad H\,S\,G\,S\,G\,S,
 \qquad H:+2,\quad G:-1.
\]
Its words are exactly those with nonnegative prefix token balance and
final balance zero. A nonempty word starts with \(H\). Its first returns
from balance two to one and then from one to zero supply the two displayed
\(G\)'s; the three intervening/remaining subwords give the unique recursive
parse. This is precisely an ordered full ternary tree, with \(m\) internal
nodes, \(m\) letters \(H\), and \(2m\) letters \(G\). Equivalently it is
the full preorder internal-node/leaf code with its final leaf omitted:
the number of pending tree slots is token balance plus one.

Consequently the grammar parse and a complete adaptive centroid-fan tree
are losslessly interconvertible using exact leaf matrices. This is a
storage statement. If an arithmetic controller uses named operations
\(H,G\), its guards, current/source ports, and terminal proof stay labels
attached to that parse. A grammar separator is not thereby an arithmetic
payment proof; geometric centroid insertion is not an arithmetic update.
The arithmetic operation word is recovered in the specified grammar order,
not by choosing an arbitrary chronological order of geometric insertions.

### Mixed binary and ternary nodes: the exact survivor and its obstruction

Adding a binary node \(P\) gives the balanced grammar
\[
 S=\varepsilon\ \big|\ HSGSGS\ \big|\ PSGS,
 \qquad H:+2,\ P:+1,\ G:-1.
\]
The same pending-slot proof identifies ordered full trees with ternary
\(H\)-nodes and binary \(P\)-nodes. With \(h,p\) such internal nodes there
are \(2h+p+1\) leaves. Geometrically use the three charts \(A_c\) at an
\(H\)-node and the two edge-bisection charts \(E_0,E_1\) at a \(P\)-node.
The typed tree and port correspondences give an exact nested rational
triangle partition. Leaf interiors are disjoint and their union covers
the parent; a path with \(h'\) fan steps and \(p'\) bisections has relative
area \(3^{-h'}2^{-p'}\). This remains a valid hierarchical partition for
every finite tree.

Two additional conclusions would be false.

First, a single product matrix no longer determines its mixed history:
\[
 \boxed{
 A_1A_0A_1E_0E_0
   =E_0E_0A_0A_1A_1
   =\begin{pmatrix}
       4/9&7/12&16/27\\
       1/9&1/12&4/27\\
       4/9&1/3&7/27
     \end{pmatrix}.}
\]
The determinant is \(1/108=3^{-3}2^{-2}\); it recovers the two counts
but not their order. These products have the same marked terminal ports,
so adding only a final permutation does not repair this collision. The
program names \(E_0,E_1\) as B0,B1 in its mixed alphabet. Exhausting every
word in \(\{A_0,A_1,A_2,E_0,E_1\}\) through length five gives no collision
through length four and exactly four two-word collision fibers at length
five: 3,906 words and 3,902 distinct matrices overall. This is the bounded
minimality statement; the displayed identity itself is exact.

The repair is to retain the typed expression or ordered tree and its
local interfaces. The pure-fan and pure-six-chart injectivity theorems
remain valid on their stated alphabets. This witness concerns one
flattened path matrix; it makes no claim that two full geometric leaf
partitions coincide.

Second, arbitrary edge bisections can create hanging boundary vertices.
A fixed-port example is an initial fan split, another fan split inside
its \(A_2\) child, then bisection of the latter's \(A_0\) child. The point
\[
 (1/6,2/3,1/6)
   =\tfrac12\bigl(e_1+(1,1,1)/3\bigr)
\]
becomes a vertex, while the untouched original \(A_0\) neighbor still
has the whole edge from \(e_1\) to the center. The partition is therefore
not face-to-face and is not automatically a simplicial complex.

A finite conforming repair is available: split every shared segment at
all incident vertices of the finite partition, then add an interior
centroid to each leaf and cone its fully subdivided boundary to that
centroid. Adjacent leaves use the same subsegments, so the resulting
triangles form a conforming refinement. Retain the original typed tree
and parent-cell map: these added geometric cells are not extra arithmetic
\(H,P,G\) operations and do not create credits. Native guards and terminal
certificates still belong to the arithmetic interface, not the geometry.

## 5. Seven coarse strata are a separate quotient

Modulo permutations, a continuous nonzero triangle point has seven possible
support/equality types:

| Type, with \(a\ge b\ge c\ge0\) | Defining boundary |
|---|---|
| vertex | \(a>0,\ b=c=0\) |
| equal edge | \(a=b>0,\ c=0\) |
| unequal edge | \(a>b>0,\ c=0\) |
| center | \(a=b=c>0\) |
| two equal large coordinates | \(a=b>c>0\) |
| two equal small coordinates | \(a>b=c>0\) |
| distinct interior | \(a>b>c>0\) |

These are coarse types, not complete orbit addresses. For example,
\((5,1,0)/6\) and \((4,2,0)/6\) have the same type and different unordered
coordinate multisets. Nor are these seven strata canonically the seven
vertices of one barycentric subdivision.

At integer degree six the two-equal-large type is absent: \(2a+c=6\)
and \(a>c>0\) have no solution. All seven first coexist at degree twelve.
Indeed an equal edge requires even degree and the center requires degree
divisible by three, so the first possible degrees are six and twelve.
At twelve, witnesses include
\[
 (12,0,0),(6,6,0),(11,1,0),(4,4,4),
 (5,5,2),(10,1,1),(6,4,2).
 \tag{4}
\]
This distinguishes a continuous list of allowed types from its realization
at a selected finite scale.

## 6. Exact chart transport for a source-aware controller

Use the marked monitor \(z=(x,n,1)^T\), where \(x\) is the current value
and \(n\) is the original source against which payment is compared. Its
normalized triangle point is
\[
 \lambda=\frac{(x,n,1)}{x+n+1}.
 \tag{5}
\]
This is injective: \(x=\lambda_0/\lambda_2\) and
\(n=\lambda_1/\lambda_2\). The original-source inequality is exactly the
coordinate wall
\[
 x<n\quad\Longleftrightarrow\quad\lambda_0<\lambda_1.
 \tag{6}
\]
For \(x,n>1\) with \(x\ne n\), only two of the six strict flag chambers
occur: \((0,1,2)\) when \(x>n\), and \((1,0,2)\) when \(n>x\). Cases
\(x=n\), \(x=1\), or \(n=1\) lie on equality walls. Thus six geometric
charts are not six independent arithmetic branches.

For a positive valuation word \(w\) with affine data
\[
 F_w(x)=(Px+B)/Q,\qquad P=3^r,\quad Q=2^A,
\]
the homogeneous monitor matrix is
\[
 M_w=
 \begin{pmatrix}P&0&B\\0&Q&0\\0&0&Q\end{pmatrix}.
 \tag{7}
\]
After normalization this updates \(x\) and retains \(n\) exactly. These
homogeneous matrices compose in chronological order:
\(M_{uv}=M_vM_u\). They are not the column-stochastic refinement matrices:
the monitor action is projective, with division by the positive sum of
output coordinates.

If the input lies in chart \(\sigma\) and the output in chart \(\tau\),
the local transition is the explicitly marked matrix
\[
 L_{\tau,w,\sigma}=B_\tau^{-1}M_wB_\sigma.
 \tag{8}
\]
Normalize after applying it to local barycentric coordinates. Its entries
need not all be nonnegative, but the actual input guard includes membership
in the input chart and the chosen output chart. On that guarded domain the
local output is nonnegative. Along a sequence of such interfaces,
\[
 L_{\rho,v,\tau}L_{\tau,w,\sigma}
   =B_\rho^{-1}M_vM_wB_\sigma.
 \tag{9}
\]
The intermediate frames cancel exactly, including on shared boundary
faces. A deterministic implementation can break coordinate-order ties
by index, but must retain that chosen chart and its zero local weights.

The arithmetic guard is still required. For a supplied positive odd \(x\),
the word has exactly its listed valuations iff
\[
 Px+B\equiv Q\pmod{2Q}.
 \tag{10}
\]
The inherited carry recursion proves this by induction; it is not merely
formal endpoint integrality. Under a chart, decode
\(x=(B_\sigma u)_0/(B_\sigma u)_2\) before applying this guard. For two
steps the composite domain is
\[
 G_w\cap F_w^{-1}(G_v),
 \tag{11}
\]
together with the chart membership conditions. Changing coordinates does
not weaken or discharge any of these predicates.

Two cheap hostiles separate the losses. Current value five is unpaid
against source three but paid against source seven; both states occur on
actual routes, \(3\to5\) and \(7\to11\to17\to13\to5\).
Forgetting the original source therefore destroys the predicate (6).
Also, the word \(12\) is legal at 27, but \(1212\) is not. An exact matrix
square or a lossless refinement word does not supply its stronger guard.

The monitor is a state representation, not by itself proof that its current
value was reached from its labelled original source. A controller must
retain that actual-prefix evidence. Likewise a smaller dependency still
needs its supplied terminal certificate before it closes the source.
The exact transport above preserves such predicates when supplied; it
neither constructs a universal selector nor proves that every source
eventually crosses the payment wall.

## 7. Reproduction and precise scope

A subsequent [marked-centroid refinement](collatz_carry_interfaces_20261004.md)
shows that a fixed interior centroid's image alone determines a finite word
in either pure chart alphabet. This strengthens the arbitrary-point case
only by retaining the distinguished seed. Final port permutations remain
invisible; a period-three triadic inverse cycle shows why arbitrary rational
points need not have finite addresses. The mixed matrix collision persists.

Run:

    python 04-computation/experiments/small_port_refinement_20261004.py
    python -O 04-computation/experiments/small_port_refinement_20261004.py

The independent standard-library script uses exact Fractions and explicit
exceptions. Its finite universe is:

- all seven subdivision vertices, twelve edges, six cells, and all
  \(7!\) candidate graph automorphisms, with and without rank preservation;
- all 1,555 chart words of depth zero through four; 1,554 terminal-frame
  controls for every word through depth three and all six port permutations;
- 1,849 composition pairs from all words of depth at most two, and all 216
  pairs of single charts with intermediate permutation frames;
- all 364 centroid-fan words through depth five, each with six terminal
  port frames; all 345 full ternary trees with at most five internal nodes,
  giving 3,601 exact leaf addresses, area identities, face counts, and
  independent checks for hanging vertices;
- every binary \(H/G\) word of length at most twelve, independently
  identifying all 72 balanced words and their unique ternary parses;
- all 3,906 mixed chart words through length five, locating the four
  two-word matrix collisions and verifying the displayed exact witness;
  all 79 ordered mixed trees through three internal nodes, checking leaf
  counts and area sums, plus the fixed-port hanging-node counterexample;
- all 454 lattice points of degrees one through twelve, including every
  boundary and the degree-six missing-stratum control;
- all 84 valuation words of lengths one through three with letters one
  through four, two exact source-cylinder representatives each, and four
  independent original-source labels, giving 672 guarded transport checks;
- sixty exact two-step chart telescopings and guard intersections, plus
  ten malformed, incomplete-tree, or undeclared-frame controls.

The transport census deliberately allows arbitrary monitor source labels:
it checks representation and current-word legality, not an asserted path
from every labelled original source. The displayed three/seven hostile has
separate literal prefix evidence.

For either pure chart family, the lossless map sends a chart word with
three terminal ports to its exact rational matrix, and decoding recovers
both word and port frame. The mixed-family collision requires retaining
the expression as well. Point-only, determinant-only, coarse-stratum, and
bare-graph quotients each discard specified data. The controller interface
uses the same compensation rule to preserve a supplied source, guard, and
payment predicate. Its precise consequence is composable representation
fidelity on the stated carrier, not new arithmetic coverage.
