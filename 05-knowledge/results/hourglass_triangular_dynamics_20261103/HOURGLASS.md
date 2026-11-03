# The diagonal-layer hourglass, over every modulus

**PROVED elementary graph statements, with exact controls.** The user
clarified that the intended four-vertex hourglass is the diagonal-layer
relation. It is not the strong four-vertex tournament, nor the seven-cell
magic-square hourglass. We retain those distinct source types.

Inheritance: [HYP-3053, diagonal-layer flow](../../hypotheses/HYP-3053-tournament-diagonal-layer-flow.md)
already gives the GF(2) dimension formulas and H(T)/|Aut(T)| presentation
fibres. Its status is SYNTHESIS, not a proved canon dependency. We prove the
needed statements below, extend to every abelian coefficient group, and give
an explicit integral basis rather than relying on the historical dimension
count. The original [S217 script](../../../04-computation/tournament_diagonal_layer_flow_codex_s217.py)
is provenance. This is standard incidence/cycle-space mathematics in an
explicit layered form; no novelty in the general graph theorem is asserted.

## H1. A rectangle/hourglass basis with no torsion

Take m>=2 nonempty layers V_1,...,V_m, sizes s_i, with every edge between
adjacent layers and no other edges. Orient all edges from layer i to i+1.
Choose anchor a_i in each layer. An abelian-group-valued line word L is of
the form L(u,v)=x_v-x_u if and only if the following relations vanish:

**Local rectangles**, one per u!=a_i, v!=a_(i+1):

    L(u,v)-L(u,a_(i+1))-L(a_i,v)+L(a_i,a_(i+1))=0.

**Hourglasses**, one per internal i=2,...,m-1 and u!=a_i:

    L(a_(i-1),u)+L(u,a_(i+1))
      -L(a_(i-1),a_i)-L(a_i,a_(i+1))=0.              (H1)

The hourglass has FOUR vertices: the left anchor, the right anchor and two
middle-layer vertices. It has four edges, not the six edges of K4.

These displayed cycles form a Z-basis of the integral cycle lattice. There
are respectively

    R=sum_(i=1)^(m-1) (s_i-1)(s_(i+1)-1),
    H=sum_(i=2)^(m-1) (s_i-1)

of them, and R+H=|E|-|V|+1. Thus all relations are needed in this basis;
there are no hidden finite-index or exceptional-modulus constraints.

**Constructive proof.** Use the anchor chain as part of a spanning tree.
Attach each nonanchor in layer 1 to a_2; attach every other nonanchor to
the previous layer's anchor. The nontree edge from a nonanchor in an
internal layer to the next anchor has exactly its hourglass as fundamental
cycle. An edge with both endpoints nonanchors has its rectangle plus the
corresponding hourglass (when the left layer is internal) as fundamental
cycle, with signs set by the orientation. These are triangular changes of
integral basis with diagonal entries +/-1. Consequently the displayed
cycles and fundamental cycles span the identical integral lattice.
Given any abelian group, integrate L along the spanning tree from one
specified root value. The cycle equations then force every remaining edge
to agree. The potentials are unique up to adding the same group element
at every vertex. This proves necessity, sufficiency and the integral claim.

For tournament layers s_i=i, m=n-1:

    |V|=binom(n,2), |E|=2 binom(n,3),
    R=2 binom(n-1,3), H=binom(n-2,2).

Modulo any q>=2 there are exactly q^(|V|-1) realizable line words, and each
has q lifts to potentials. This holds for composite q too. A single anchor
potential is the smallest restoration coordinate: all q choices actually
occur. Adding the next tournament layer adds n potentials and n(n-1)
lines. For n>=3, the rank increment from order n to n+1 is n; the cycle-rank increment is
n(n-2). The latter splits into (n-1)(n-2) rectangles and n-2 hourglasses.

At full order four the layer sizes are 1,2,3: six potential vertices, eight
lines, rank five, two rectangles and ONE hourglass. A line word supported
only on the edge from the first anchor to the nonanchor of the middle
layer satisfies both local rectangles but violates the hourglass. This
is an explicit local-to-global obstruction, not just an extra dimension.

## H2. The fixed-path cube is a different carrier

Fix the Hamilton path 0->1->2->3. Order its three free edges as
(0,3),(0,2),(1,3), giving eight potential words. Its free layers have sizes
1,2, so the line graph is K_(1,2), with two independent XOR lines and no
cycle relation. Its first hourglass appears at tournament order FIVE.
This is why the full order-four hourglass should not be confused with a
relation among the three free bits of the order-four fixed-path cube.

For any tournament class T, Aut(T) acts freely on its ordered Hamilton
paths: an automorphism fixing an ordered path fixes every vertex. Path
orbits are precisely the fixed-path presentations of T. Hence the fibre
size is H(T)/|Aut(T)|. At order four:

| Class | H | Aut order | Fixed-path presentations |
|---|---:|---:|---:|
| transitive | 1 | 1 | 1 |
| source plus directed triangle | 3 | 3 | 1 |
| sink plus directed triangle | 3 | 3 | 1 |
| strongly connected | 5 | 1 | 5 |

Converse merges the two middle classes. Thus the two fibre multisets are
(1,1,1,5) and (1,2,5). They are different from the quadratic folds below.

## H3. An exact quadratic operation on the fixed-path cube

Identify a bounded numerator a in {-7,-5,...,7} with

    j=(a+7)/2 = b0+2b1+4b2,  0<=j<=7.

Use these bits in the free-edge order above. This is a declared coordinate
bijection of presentations, not an intrinsic identification of tournament
classes. Antipodal a->-a is exactly bitwise complement j->7-j. The
quadratic step becomes

    j -> F(j)=(j-3)(j-4)/2,
    [F(0),...,F(7)]=[6,3,1,0,0,1,3,6].              (H2)

Let u=b0 XOR b1, v=b0 XOR b2 be the two K_(1,2) line measurements.
Over GF(2), the three output bits are

    b0'=u,  b1'=1+u+v,  b2'=(1+u)(1+v).             (H3)

Thus both lines determine the next **entire** cube word. In particular,
they form a closed four-state carrier with induced update

    (u,v) -> (1+v, 1+v+uv).                          (H4)

Its cycle is 00->11->01->00, and its fourth state 10 enters 11.
Four is minimal if the observer must report the next exact numerator:
there are four different next values. With only time>=2 observations,
three states suffice and are necessary. With the current signed numerator
as observer, all eight are needed. Minimality is relative to the observer.

The forward-image ladder is literally 8->4->3, followed by a three-cycle.
But its fibre sizes are (2,2,2,2) at one step and (2,2,4) at two steps,
whereas the tournament quotients have (1,1,1,5) and (1,2,5). No bijections
can make those quotient diagrams agree: fibre cardinalities are invariant.
In this chosen coordinate chart the tournament-class quotient does not
even support H2's update: j=0 and j=1 represent the same strong tournament,
but go to j=6 and j=3, respectively a strong class and a triangle with
an extremal vertex. They remain different after merging converse.

More generally, on a connected binary line graph delta(x), an update U
descends to line words iff

    U(x+1)+U(x) is a constant vertex word for every x. (H5)

It predicts the entire next potential word iff that difference is zero.
This follows exactly from ker(delta)={0,1}; it is both necessary and
sufficient. H2 meets the stronger condition. Complement is a kernel
symmetry of F, not a commuting automorphism: F(7-j)=F(j), whereas
7-F(j) generally differs. Future named-port tournament operations need
their own H5 check; the all-size attachment carrier is not thereby reduced
to two bits.

## H4. Literal FRACTRAN, a cellular rule, and the clock

Use state primes p_j=(2,3,5,7,11,13,17,19), clock prime 23, source-tag
prime 29. On promised inputs p_j*29^j the ordered fraction list is

    391/2, 161/3, 69/5, 46/7,
    46/11, 69/13, 161/17, 391/19.

After exactly t FRACTRAN steps its value is

    p_(F^t(j)) * 29^j * 23^t.                        (H6)

Each step has exactly one applicable fraction, since there is one control
prime. This proves the simulation, source preservation, and exact clock
by induction. It is a finite-state encoding with an unbounded clock, not
a simulation of the full unbounded quadratic map or of PRIMEGAME. Outside
the one-control source predicate the semantics differ: input 6 has two
applicable fractions and priority matters.

The sitewise rule (source,current,phase)->(source,F(current),phase+1 mod3)
is a genuine shift-equivariant radius-zero CA on a finite alphabet. It
retains the original source and the three-cycle phase, not unbounded elapsed
time. The repository's [TM/CA/FRACTRAN dictionary](../machines_automata_fractran_tournaments_dictionary_20261101/README.md)
supplies the wider machine context; changing the encoding never by itself
proves an orbit terminates.

For the eight-state transfer matrix A_(F(j),j)=1,

    det(I-zA)=1-z^3,   trace(A^t)=3[3 divides t], t>=1.

The five transient dimensions are visible in initial-state counts but absent
from periodic return traces. This is an exact illustration of why return
spectrum, source multiplicity and clock must remain separate.

## H5. The unique triangular parameter with faithful tournament completions

**PROVED consequence of D1.** Among all integer k with a nonempty bounded
integer state set, k=4 is the ONLY one whose full functional graph can
be retained as a partial tournament orientation. A fixed point forces a
loop and a two-cycle forces antiparallel arcs. Both are forbidden in a
tournament. D1 classifies all cycles and leaves only the three-cycle
parameter k=4. Conversely its eight-state functional graph has no loop
or antiparallel pair, so every unspecified pair may be oriented freely.
This is an intrinsic directed relation with explicitly unspecified pairs,
not a tie-breaking claim that the function itself is a tournament.

The graph has eight prescribed arcs among28 pairs, hence2^20 labelled
tournament completions. The [incoming completion recurrence](../hourglass_quadratic_ports_20261103/COMPLETIONS.md),
already proved there from the unique source SCC, applies without deleting
or reversing a transition. For any partial orientation D it is

    C(U)=2^(binom(|U|,2)-|D[U]|),
    S(U)=C(U)-sum_(empty!=A proper subset U,
                   no D arc U\A -> A) S(A) C(U\A).

Our application gives exactly **982540 strongly connected completions**.
The direct C++ census independently checks all1048576 free orientations
using both forward and backward reachability. The four-state line quotient
has three strong completions out of four, an inherited small control.

The eight-state functional graph has automorphism group C2, generated by
swapping the two deepest leaves, numerator-3 and3. The unique depth-two
branch fixes its parent and hence the cycle phase and every other vertex.
The swap cannot fix a tournament, since it reverses the arc between those
two leaves. Therefore there are491270 isomorphism classes of pairs
(strong tournament completion, retained functional graph), under this
specific decoration-preserving equivalence. This is NOT an ordinary
A000568 count; forgetting the functional graph allows further identifications.
All40320 permutations are checked as an independent finite automorphism
census in `completions.py`.
