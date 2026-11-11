# Every closed tree of the path–cut incidence graph is a star

Date: 2026-11-11. **PROVED:** the all-order star theorem, uniqueness of the
closed tree, and the quantitative domination inequality below. It concerns
the literal incidence graph of [PROOF.md](PROOF.md), including every stub;
removing stubs or retaining only fibre sizes would invalidate the argument.

Let U be a tournament, V=A⊔B, and let the newborn z satisfy z→A and B→z.
Keep the path vertices, dominating-cut vertices, edges and stubs defined
in PROOF.md. A closed component contains no stub. Throughout this note C
is a closed component which is a tree. All degrees below are degrees in C.

**Theorem.** C is a star. More precisely, either its centre is the cut B,
or all its cut vertices have degree one and its centre is an old path.

The two-vertex tree can be viewed in either way when its cut is B.

## 1. Three local facts

### An old path has at most one lower and one upper insertion cut

A lower cut means L⊆B, and an upper cut means L⊇B. A legal lower insertion
slot must occur immediately before the first A vertex of the old path,
after its initial all-B segment. It is therefore unique. Dually, the only
possible upper slot occurs after the last B vertex, before the final
all-A segment. Empty colour classes obey the same statement with the
appropriate endpoint slot.

In particular a path incident with cut B has colour word B*A* and exactly
one legal slot. Thus every such path is a leaf. If B belongs to a closed
component, that whole component is already a star centred at B.

### Branching can occur only below or above B

For a dominating cut L, write R=V\L. Its fibre consists of every pair
of a B-ending Hamilton path p of U[L] and an A-starting Hamilton path q
of U[R], with empty sides assigned their one empty path. Closedness means
every p q is an old Hamilton path; otherwise p z q is a stub at L.

If there are two choices p_1,p_2 and R contains a B vertex, fix q. Since q
starts in A and contains B, it has a second legal insertion slot beyond
its beginning. The two old paths p_1 q,p_2 q share this second cut and L.
If that second cut is not dominating there is a stub; if it is dominating
there is a four-cycle. Both contradict C being a closed tree. Hence

    e_B(U[L])≥2  implies R⊆A.                         (1)

The dual argument gives

    s_A(U[R])≥2  implies L⊆B.                         (2)

The fibre size is the product of these positive endpoint counts.
Consequently every cut of degree at least two is lower or upper.

### Two-sided endpoint rigidity

We use the proved lemma in
[history/ENDPOINT_LEMMA.md](../history/ENDPOINT_LEMMA.md), §1:

> If a nonempty tournament W is partitioned into nonempty A_W,B_W,
> every A_W vertex points to some B_W vertex, and every B_W vertex
> receives from some A_W vertex, then its A_W-starting/B_W-ending
> Hamilton-path count is one exactly when W is transitive.

That proof supplies an explicit second path by three Hamilton-cycle
reroutings in the strong case and by strong-component concatenation in
the remaining cases. It does not assume the present star theorem.

## 2. Canonical neighbours of proper upper and lower cuts

Let L⊃B be a proper upper cut in C. Put R=V\L and define

    D(L)={b∈B : no a∈L∩A points to b},
    X=D(L),   Y=L\D(L).                              (3)

Here D denotes this cut map, not the dominating-set count. We claim:

1. U[X], U[Y], U[R] are transitive;
2. X is a proper lower dominating cut in C;
3. the unique concatenation of their transitive orders is a path vertex
   adjacent to both X and L.

First, U[R] is transitive. Otherwise s_A(U[R])=H(U[R])≥2, and (2) would
force L⊆B, contrary to properness. Empty R has its unique empty order.

Every X vertex beats every A vertex of L. Domination at L implies each
A in Y points to some B in Y: a B out-neighbour cannot lie in X by its
definition. Every B in Y receives from an A in Y, again by definition.
Both parts of Y are nonempty. Apply the two-sided endpoint lemma.

For every Hamilton path p of X and every A-starting/B-ending Hamilton
path q of Y, the concatenation p q is a B-ending path of L; the seam is
X→A. Let r be the unique path of R. Closedness at L makes p q r an old
path. It has legal insertion slots at X and L. Also X is dominating:
the left side is all B, and each right-side B has an incoming A in Y.

Two different choices of (p,q) would give a four-cycle between X and L.
Therefore both choices are unique. The ordinary unique-path theorem
makes X transitive, and endpoint rigidity makes Y transitive. Their
orders concatenate as claimed. Y has an A first and a B last, so X⊊B.

By reversing every arc, swapping A and B, and complementing every cut,
the dual canonical map on proper lower cuts K⊊B is

    C(K)=V\{a∈A : a has no out-neighbour in B\K}.     (4)

It is a proper upper dominating cut in the same component. Its canonical
old path is adjacent to K and C(K), and the three corresponding blocks
are transitive. The maps (3)–(4) are defined from actual arc incidences,
not chosen by endpoint counts.

## 3. Strict progress from an upper branching cut

Suppose the proper upper cut L⊃B has degree at least two. Keep X,Y,R
from §2 and write their transitive orders as

    X=x_1,...,x_p,   Y=y_1,...,y_q,   R=r_1,...,r_s.

Every x_i is B, y_1 is A, and y_q is B. Since the right endpoint factor
is one and L branches, U[L] has at least two B-ending paths and therefore
is not transitive. X and Y are individually transitive, so some backward
arc y_j→x_i exists. Such y_j lies in B because every X beats every A in Y.

Choose j maximal among indices admitting any backward arc to X, and
choose i with y_j→x_i. Then the following is a B-ending path of L:

    x_1,...,x_(i−1),
    y_1,...,y_j,
    x_i,...,x_p,
    y_(j+1),...,y_q.                                 (5)

The first seam, when present, is B→y_1∈A. The middle seam is the selected
arc. The final seam x_p→y_(j+1), when present, follows from maximality of
j: every later Y vertex loses to every X vertex. If j=q, the path ends
at x_p∈B; otherwise it ends at y_q∈B. Thus there is no endpoint exception.

Append the order of R. Closedness at L guarantees the seam from this
B endpoint to r_1, if R is nonempty. We obtain an old path P adjacent
to L. Its initial B segment is exactly

    K={x_1,...,x_(i−1)} ⊊ X=D(L).                    (6)

The insertion slot before y_1 is legal. If K were not dominating it
would be a stub at P. Hence K is a proper lower dominating cut in C,
and P joins K to L.

Moreover its right side contains the directed triangle

    x_i→y_1→y_j→x_i.

The first arc follows because X beats every A in Y, the second follows
from the transitive Y order, and the third is the selected backward arc.
Here j≥2 because y_1∈A while y_j∈B. Thus U[V\K] is nontransitive.
The cut-leaf lemma in PROOF.md now implies that K has degree at least
two. This proves branching directly, before using a canonical neighbour.

Now take the canonical upper neighbour C(K). Every A in L has an
out-neighbour in B\D(L), and this set lies in B\K. Thus L⊆C(K).
For any b∈D(L)\K, domination at K supplies an a∈A with a→b. No A in L
can do this by the definition of D(L). Therefore a∉L, but a∈C(K), and

    C(K)⊃L.                                         (7)

We have proved the exact upper progress rule:

    proper upper branching L
      -> proper lower branching K⊊D(L), with C(K)⊃L. (8)

Every arrow is witnessed by actual old paths in the same closed component.

## 4. The dual progress rule and finite contradiction

The involution

    (U,A,B,L) -> (U^op,B,A,V\L)

preserves domination, path incidences, stubs, closedness, treehood, and
degrees. It interchanges the two canonical maps:

    D_new(V\K)=V\C_old(K),
    C_new(V\J)=V\D_old(J).

Applying (8) after this involution gives the precise lower progress rule:

    proper lower branching K
      -> proper upper branching J⊃C(K), with D(J)⊊K. (9)

Now suppose C has a proper upper branching cut. Choose one L of maximum
cardinality. Rule (8) gives a lower branching K, and (9) then gives an
upper branching J with

    J⊃C(K)⊃L,

contradicting the choice of L. Thus no proper upper branching cut exists.
By (9), no proper lower branching cut can exist either. By (1)–(2), the
only possible branching cut is B.

If B lies in C, §1 makes C a star centred at B. Otherwise every cut has
degree one. A connected bipartite graph with all vertices on one side
of degree one has just one vertex on the other side; here it is an old
path. Thus C is a star centred at that path. QED.

## 5. At most one closed tree, and the quantitative inequality

Here is the uniqueness argument from
[classification/EXCHANGE_IDEAS.md](../classification/EXCHANGE_IDEAS.md),
§4, now combined with the proved antecedent.

There is at most one cut-centred closed star, because its centre must
be the unique cut B. For a path-centred closed star with at least two
leaves, write its old path as P=X Y Z, with X its initial B run, Z its
final A run, and Y starting in A and ending in B. Its earliest and
latest insertion cuts are leaves. The cut-leaf lemma gives

    U[X union Y] and U[Y union Z] transitive in P order.

Thus the orders on A and B are fixed by U, and all backward arcs relative
to P go from Z⊆A to X⊆B. In particular no backward B→A arc is possible.

If two distinct such stars P,Q existed, they would be different shuffles
of the same two shore orders. Choose a∈A,b∈B with a before b in P and b
before a in Q, swapping P,Q if needed. The arc must be a→b, because
b→a would be a backward B→A arc in P. Therefore b∈X_Q and a∈Z_Q.
Let a0 be the first A in Y_Q. Then b→a0, and a0 precedes a in the fixed
A order. In P, a0 precedes a which precedes b, making b→a0 a forbidden
backward B→A arc. This contradiction proves uniqueness of the path star.

It also cannot coexist with a closed cut-B star. Its shores are
transitive, so their unique B-then-A concatenation uses as seam the last
B vertex b0 and the first A vertex a0 of P. Both lie in Y, where a0→b0.
Their newborn path with cut B is therefore deletion-invalid, a stub at
B. A one-edge closed component is already covered by the cut-B case.

Consequently the number tau of closed tree components satisfies

    tau≤1                                             (10)

for every old tournament and every newborn cut. The exact defect identity
and constructive matching in PROOF.md immediately give

    H(U+z)−H(U)−d(G)+1=u+1−tau+E≥0,

and hence

    H(U+z)≥H(U)+d(G)−1.                              (11)

The same proof constructs an injection of all old Hamilton paths and
all but at most one dominating cut into distinct newborn Hamilton paths.
It therefore closes the simultaneous-reservation gap in the earlier
domination injection, rather than only bounding the two fibre totals
separately. The empty old tournament, if included, has H=1 and d=1 and
satisfies (11) directly under the usual empty-path convention.

## Scope and dependency boundary

The proof uses the exact cut fibres and closedness each time a legal
insertion is declared dominating or a deletion seam is supplied. Those
are essential: the maximal-j splice alone gives no theorem after stubs
are discarded. It also keeps both canonical cut maps; a scalar defect
or an arbitrary choice of neighbour would not give strict inclusions.

This proves the quantitative domination inequality for every tournament
and cut. It does not establish that every odd path total other than 7
and 21 is realized; any use for that larger spectrum problem must retain
the actual source tournament and its cut, rather than treating (11) as
a freely selectable count increment.
