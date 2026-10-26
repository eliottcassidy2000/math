# Two composition laws on the marked metagraph carrier

2026-10-26. **PROVED** in the displayed marked carrier; bounded exact
controls are **FINITE-EXACT**. Independent audit is recorded separately.
The contribution is an exact composition interface and its smallest
two-factor obstruction, not a new enumeration of tournament classes.

## Inheritance and the objects that must stay separate

The closest mechanisms are [THM-2183, order-join metric product](../../../../01-canon/theorems/THM-2183-order-join-is-an-exact-tournament-metric-product.md)
and [Session67's relocation carrier](../../session67_metagraph_enumeration_information_20261026/session67-root/COUNT_FREE_ENUMERATION.md).
The former concerns unlabeled **single-arc reversal distance**; the latter
uses **re-marking of a Hamiltonian path**, which changes no tournament class.
[THM-796, three-sorted incidence](../../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md)
already proves why presentations, literal complement lines, and class nodes
must be separate, and gives earlier non-lumpability controls for deletion.
The present operation is concatenation, not deletion.

Let X_n consist of all tournaments on positions 0,...,n-1 with arcs
i->i+1. Put X_0={empty}. Free bits x_(i,j), j>=i+2, are one for i->j.
The involution tau reverses every nonconsecutive arc and preserves the
path. Reflection sigma reverses positions and takes tournament converse;
sigma and tau commute. For n>=3 a literal line is the two-element pair
{x,tau x}; it is blue when sigma x=x. At n<=2 tau is the identity: these
carriers are formal composition inputs, not two-element line pairs.
Ordinary classes identify isomorphism only. A final
converse quotient gives the merged node convention, but is not used to
define the operations below.

The relocation graph R_n has these path presentations as vertices. Retain
each legal rightward move a->after b, b>=a+2, as a named directed move;
the result is the entire relation reindexed in the new order. It is legal
iff that new order is a Hamiltonian path. Include loops and distinct move
instances if they have the same endpoints; leftward moves are inverses.

Anchor: composition across sizes. Niche: a sufficient rooted response
distribution for exact metagraph multiplicities. Wildcard: whether an
unoriented blue loop is enough input for composition. Live board: ordered
joins / relocation moves / complement lines / rooted endpoints / class mass.
The hostile below makes the wildcard decisive.

## M1. Relocation is a Cartesian product on an order-join patch

For x in X_a and y in X_b, let J(x,y)=x->y be their order-join on consecutive
position blocks. All cross arcs point from x to y. For fixed a,b>0,

    J: X_a x X_b -> X_(a+b)

is injective, its image is closed under legal relocation, and its named
relocation multigraph is exactly R_a Cartesian-product R_b. This concerns
relocation, not the complementary-line graph and not the arc-flip graph.

**Proof.** Every Hamiltonian path of a join traverses the first block
entirely before the second: after entering the second it cannot return.
A legal relocation therefore cannot change the two-block order. Within
a block it is legal precisely when the new internal consecutive arcs
form a path. At the cut, any newly exposed cross arc is still forward,
so there is no extra boundary condition. The full reindexed relation is
J of the changed factor and the unchanged factor. The moves, their names,
loops and multiplicities are therefore exactly the Cartesian-product
moves. Restriction to the blocks is the inverse on vertices. QED.

Consequently an ordinary class fibre of A->B is the product of the ordinary
presentation fibres of A and B. This recovers the inherited weight law

    tc(A->B)=tc(A)tc(B),    tc(T)=H(T)/|Aut(T)|.

This proof retains all presentations; a single selected representative
per class would not supply the movement or line multiplicities.

J is associative with empty unit. Every marked tournament has a unique J
factorization into its ordered strong components, which appear as contiguous
intervals of every Hamiltonian path. Thus this graded marked monoid is free
on the strong presentations together with the one-vertex presentation.
This is the classical SCC factorization in this carrier; it is not a
classification of modular-prime tournaments. Its atoms are **strong**, not
necessarily prime in Gallai's modular decomposition.

## M2. Complement transports the product to an associative hinge operation

Define

    x star y = tau_(a+b) J(tau_a x,tau_b y).             (M2)

Inside the two blocks the relations remain x and y. Every cross arc points
from the second block to the first, except the retained path edge

    last(x) -> first(y).                               (M3)

Thus M2 is a reversed order-join with one named edge changed. Equations
M2-M3 are an exact map, not an analogy between a join and a metagraph edge.

**Proof and consequences.** The only consecutive cross pair is the cut
pair, so tau flips exactly the other cross arcs. As tau is involutive,
conjugating associativity of J proves associativity of star, with empty
unit. It also proves

    tau J(x,y) = tau x star tau y,
    sigma J(x,y) = J(sigma y,sigma x),
    sigma(x star y) = sigma y star sigma x.             (M4)

For any ordered list of blocks, the star product has all cross-block arcs
pointing backwards except the named cut edges between consecutive blocks.
There is no dependence on parenthesization. It has a unique star
factorization: take the J-factorization of tau x, and apply tau to each
factor. Hence its atoms are tau-images of strong presentations. These
are not the ordinary strong atoms: the transitive three-tournament is a
star atom because its tau-image is C3, while C3 is the star product of
three singleton atoms.

Transporting the *relocation graph* by tau gives R_n^tau, with edges
x->y iff tau x->tau y is a relocation. In this graph star patches are
Cartesian products, by M1. This does not say that star patches are closed
under ordinary relocation; the move relation was transported too.

If a,b>0 and a+b>=3, x star y is strong: the fixed path, closed by the
cross arc (a+b-1)->0, is a Hamiltonian cycle. The exception at a=b=1 is
the two-vertex transitive tournament. In particular complement sends any
nonstrong presentation of order at least three to a strong presentation.
This is not a test for modular primeness.

## M3. Exact cyclic-triple growth and the missing endpoint data

Writing alpha=last(x), beta=first(y),

    c3(x star y)=c3(x)+c3(y)+indeg_x(alpha)+outdeg_y(beta). (M5)

**Proof.** A cross-block cyclic triple must contain the unique forward
cross edge alpha->beta. Its third vertex is either u in the first block
with u->alpha, or v in the second block with beta->v. These are disjoint
possibilities and both suffice. All other cyclic triples are internal.

For r>=2 nonempty factors x_i of orders n_i, the iterated formula is

    c3(star_i x_i)
      = sum_i c3(x_i)
        + sum_(i<r) [indeg_(x_i)(last)+outdeg_(x_(i+1))(first)]
        + #{i: 1<i<r and n_i=1}.                       (M6)

A triple meeting three blocks is cyclic only when both consecutive
forward hinges occur; their middle endpoints coincide exactly when the
middle block is a singleton. This proves M6 directly, including the
singleton boundary. The n-fold star product of singletons has n-2 cyclic
triples for n>=2. The integers here measure exact edge configurations,
not a transferred number-theoretic invariant.

The degrees in M5 suffice for c3 only. They do not determine the glued
tournament, its Hamiltonian paths, or its metagraph line class.

## M4. A smallest two-factor coarse-composition obstruction

At order four use free-bit order (0,2),(0,3),(1,3), with least significant
bit first. Let x have mask2 and y mask5, so y=tau x. Their adjacency rows
are

    x = 0101, 0010, 1001, 0100;
    y = 0110, 0011, 0001, 1000.

They are isomorphic strong tournaments S4: an isomorphism x->y sends
(0,1,2,3) to (1,3,0,2). Its inverse is (2,0,3,1). Each has H=5 and
trivial automorphism group. Both are reflection-fixed, so {x,y} is **one
literal blue loop**, with both ends at the same ordinary class S4.

Let e be the singleton presentation. The two complements

    z=tau J(e,x),    z'=tau J(e,y)

have respectively H=11 and H=13, so they are not isomorphic, even up to
converse. Both have c3=4, score multiset (1,2,2,2,3), and trivial
automorphism group. Both factor tau-targets have first-vertex outdegree2.
Thus even these data together do not determine composition:

* the source and target ordinary classes of each factor line;
* the literal unoriented four-vertex line itself and its blue colour;
* the size and relevant endpoint degree.

The lost choice is which half-edge is used, equivalently which rooted
vertex of the target class receives the hinge. S4 has trivial automorphisms,
so its distinct degree-two vertices need not be interchangeable. This
failure is earlier and sharper than guessing an all-order coarse recurrence
from the two endpoint class names.

**FINITE-EXACT minimality in the stated operation.** Exhausting every split
a+b=n and every factor presentation for 2<=n<=4 finds no distinct outputs
from the same pair of factor source/target class pairs. The first failure
is n=5, split1+4 (and its reversed counterpart). This is not a universal
minimality claim for all possible metagraph operations.

## M5. A sufficient rooted response and an exact multiplicity convolution

Fix an ordinary tournament class A of order a. For each rooted ordinary
class (C,u), define

    E_A(C,u)=#{x in X_a: [x]=A, [(tau x,last)]=(C,u)},
    B_A(C,u)=#{x in X_a: [x]=A, [(tau x,first)]=(C,u)}.

These are finite, exact integer histograms, with total tc(A). Let
G((C,u),(D,v)) be the unmarked class obtained from two disjoint representatives
by orienting every cross arc D->C except u->v. Rooted isomorphism makes G
well defined. Then the number of join presentations with complement in Z is

    M_(A,B)(Z)
      = sum_(C,u),(D,v) E_A(C,u) B_B(D,v)
          1{G((C,u),(D,v))=Z}.                         (M7)

**Proof.** M1 identifies the entire source fibre with independent choices
of x and y. M2 identifies their complemented join with precisely the
displayed rooted gluing. Group the finite sum by rooted target types.
There is no independence assumption beyond that proven Cartesian product.
This is consistent with the recovered
[odd-cell gluing warning](../../odd_cell_gluing_20261023/README.md): uniform
cross arcs allow independent module choices; partial-bijection transport
requires a separately supplied compatibility relation.
For a+b>=3, the source join is nonstrong and every target is strong, so
each such presentation is the unique join-side half-edge of its literal
line. M7 therefore also counts those literal line instances; no division
by two is appropriate. A final converse quotient can be applied to Z.

M7 is a sufficient interface, not a proof of smallest possible state space.
It does not determine blue/black colour from these rooted types alone;
retain the reflection predicate on full presentations for coloured lines.
It is also not a scalar recurrence for A000568: computing the rooted
histograms still requires structural information about the class fibres.

## Verification and connection contract

`check.py` constructs relations directly. It checks Cartesian relocation
with full reindexing, associativity and reflection, the cyclic-triple
formulas, SCC/star factorizations and the explicit minimal witness. It
computes canonical classes by all permutations only in the small audit
universe, and H independently by subset dynamic programming and direct
Hamiltonian-word enumeration for the witness. Normal/-O receipts agree.

The explicit finite universes are 4,752 factor-presentation pairs with
1<=a,b<=6 and a+b<=8 (10,904 legal named relocation instances); 116
three-factor tuples with each size at most four and total at most seven;
all 33,868 presentations in X_1 through X_7; and 32 ordinary class pairs
(48 presentation pairs) for the rooted convolution, with factor orders at
most four and total at most six. All splits below total order five supply
the positive minimality controls; the one-blue-loop pair supplies the
hostile. Singleton hinge products are checked through order fifteen.
These controls support, and do not replace, the all-order proofs.

Source -> target: marked path presentations -> exact composition and
literal complement-line multiplicities. Maps: J, tau and rooted hinge
gluing. Preserved: the whole binary relation, named endpoint ports and
presentation counts. Lost by class-only projection: which degree-equal
vertex is the endpoint, and which half-edge of a blue loop is selected.
Required sidecar for uncoloured targets: the two rooted histograms in M7.
Cheapest hostile: the single order-four blue loop above.

Board update: join has an exact product law for relocation, while tau
transports both product and movement; the endpoint-rooted target restores
composition multiplicities but ordinary class names and degrees do not.
The open classification obligation is a useful small closed family of
rooted responses for modular-prime substitution and for the coloured
metagraph, not another census of unmarked class counts.

Reproduce from this directory with `python3 -B check.py` and
`python3 -B -O check.py`, comparing `normal.json` and
`optimized.json`. Both are LF files with SHA-256
`7adfd931cbc1c548bc7d78051396ce00352f192517787111b5dd46eccd630644`.
