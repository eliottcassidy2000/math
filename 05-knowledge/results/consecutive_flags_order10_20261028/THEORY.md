# Paired structures and the exact scope of smaller carriers

Local proposition labels below belong to this file, not the global theorem-ID
namespace. **PROVED** means the supplied mathematical argument; finite tables
remain **FINITE-EXACT**. Classical and concurrent dependencies are identified.

## P1. A square matrix is an aligned pair with a nonstandard action

For i<j, read M_ij as the arc i->j of T, M_ji as the arc i->j of W, and
M_ii as the root-outmask s_i extending W to U. This is a labelled bijection

    {0,1}^(n*n) <-> (T_n, W_n, s) <-> (T_n, U_(n+1), distinguished root).

It uses a common identification of the old n vertices. It does not require
W=T. The nested extension problem W=T occupies only 2^(n(n+1)/2) inputs.
Simultaneous relabelling of tournaments induces an affine signed permutation
of the orientation coordinates, not ordinary matrix conjugation. No equivariant
bijection with ordinary matrices can repair this globally: a transposition
fixes the zero matrix but fixes no tournament, since it reverses the arc between
its exchanged vertices. Equal labelled cardinalities do not imply equal orbit
problems. These statements are also proved in the concurrent
[encoding/action packet](../session89_consecutive_tournament_pairs_20261028/encoding/ENCODING_AND_ACTIONS.md).

A class and one node specify a pointed isomorphism type. They do not specify
all labelled realizations or a Hamilton path. For a transitive n-tournament,
there are n! labelled realizations but only n choices of distinguished vertex.
Retain a normalization transport coset whenever later operations address labels.

## P2. Reversing-involution doubles, with the symmetry boundary

Let G_ij=1 iff T and W agree on {i,j}; this is a symmetric graph coordinate.
Use vertices i+,i- and, in the odd construction, infinity. Orient the positive
copy as T and the negative copy as its converse. For i!=j put
i+->j- iff G_ij=1. In the even construction put i+->i- for all i. In the odd
construction put i+->i- iff s_i=1, infinity->i+, and i-->infinity.
Every unordered pair is assigned exactly once. Exchanging + and -, and fixing
infinity when present, reverses every arc. Call this anti-involution pi.

Given a self-converse tournament D of order at least two, Aut(D) has odd order:
an involution would exchange two vertices and contradict their arc orientation.
The group of automorphisms and anti-automorphisms has twice that order. Its
Sylow-2 subgroups give anti-involutions, all conjugate under Aut(D): a negative
conjugator becomes positive after composing with the source anti-involution.
An anti-involution has at most one fixed vertex, hence none in even order and
exactly one in odd order.

For even order choose the winner in each paired arc as i+. For odd order choose
the vertex beaten by the fixed vertex as i+. These normalizations remove the
independent pair-flip choices. Read T, G, and, when present, s; recover W from
T and G. Only simultaneous pair permutation remains. Conjugacy of the marks
proves both surjectivity and injectivity on unmarked isomorphism classes:

    aligned (T_n,W_n) classes <-> self-converse classes of order 2n;
    aligned (T_n,U_(n+1),root) classes <-> self-converse classes of order 2n+1.

The odd structural correspondence is **CITED classical**, after changing graph
coordinates, from Eplett, *Self-Converse Tournaments* (1979), Section 1,
[primary paper](https://doi.org/10.4153/CMB-1979-004-6).
The proof here also supplies the even normalization. Concurrent Session89 uses
a different odd normalization, choosing every paired-arc winner as positive
and putting the diagonal mask on the root arcs; its
[lift and decoder](../session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md)
gives the same kind of class bijection. These are different codecs; do not
compare their raw target codes without applying the gauge change. No external
priority claim is made for either construction.

The source automorphism group is exactly C_Aut(D)(pi). The full target group
can be larger, and the index equals the number of anti-involutions of D.
Thus forgetting pi is a bijection on isomorphism classes but not an equivalence
preserving automorphism groups or inverse-group-order weights. At n=1, s=1
gives C3, whose full group has order 3 although the source group is trivial.
This is the smallest symmetry-loss witness.

`pair_bridge.cpp` exhausts all inputs for n=1..4, both parities. It checks
literal anti-invariance, orbit injectivity, and source-group equality with
the full target centralizer. The even class counts are 1,2,12,176 and the odd
counts 2,8,88,2752. At n=3 an odd source group of order 3 can map to a full
target group of order 21. Counts agree with
[McKay's primary data](https://users.cecs.anu.edu.au/~bdm/data/digraphs.html).

This is **not** a native-complement-closed carrier. Fixed-path rows
(6,4,8,3) at order four are self-converse; their native partner (10,12,9,0)
is not. At order five the analogous witness is (6,4,8,19,7), whose partner is
(26,28,25,16,0). Rows are outgoing-neighbor bitmasks, path 0,1,...,N-1.
The literal-permutation audit verifies both failures.

## P3. The missing parent coordinate is an ordered path forest

Let P be a Hamilton path of U, v its distinguished vertex, and T=U-v. Native
complement tau_P reverses every non-P arc and preserves the P arcs. Deleting
v from P leaves an ordered spanning path forest F, with one component at an
endpoint and two at an interior vertex. Direct restriction gives

    tau_P(U)-v = T^op with precisely the F arcs reversed back.

The root-outmask is complemented except at its one or two P-neighbors.
This proves an exact joint update from (T,F,root mask,insertion position).
There is no need to recompute the two tournaments independently, but F cannot
be replaced by a guessed parent Hamilton path. Already a transitive triple
with its middle vertex deleted has its remaining edge reversed by tau, whereas
native complement of an order-two parent preserves its sole spine edge.
Deleting k vertices similarly leaves at most k+1 path runs.

This specializes the inherited ordered-forest mechanism in
[strong-atom continuation, C1](../strong_atom_continuation_vitali_20261026/README.md).
The checker verifies all 6,504 fixed-path/root states of orders 2..6, including
every root-mask neighbor exception. It checks the actual paired operation, not
just equality of class counts.

## P4. Smallest closed carrier for one declared protocol

Fix N. The literal state is s=(x,j), where x gives the (N-1)(N-2)/2 nonspine
arc bits of a tournament with Hamilton path 0->1->...->N-1 and j is a
distinguished vertex. There are N*2^((N-1)(N-2)/2) states. Define:

- q(s): the complete pointed isomorphism type (T,j).
- o(s): the two ordinary classes ([T],[T-j]).
- tau(s): complement all nonspine bits, leaving j fixed.
- R(s): choose a Hamilton path of T uniformly, retain the actual distinguished
  vertex, and rewrite in its path coordinates.

R is uniform over all normalized states with the same q. Each such state has
the same number |Aut(T,j)| of preimage Hamilton paths, since an automorphism
fixing a fully ordered path fixes every vertex. Thus q determines the reset
distribution, even though o generally does not.

Start with colour o. Repeatedly replace the colour of s by the triple consisting
of its old colour, the old colour of tau(s), and the exact probability vector
of old colours under R(s). Stop when no split occurs. The final partition is
the **coarsest strongly lumpable quotient** preserving o and deterministic tau.
Proof: every admissible quotient refines the initial observation. Inductively,
tau must agree and reset probabilities into every old block must agree, since
old blocks are unions of admissible blocks. Hence every admissible quotient
refines every refinement step. At termination those very tests prove closure.
This proves minimality among these congruences, not arbitrary weak trace or
distributional language equivalence.

| N | Literal states | Pointed types q | Initial o | Smallest closed quotient | (q,q after tau) |
|---:|---:|---:|---:|---:|---:|
| 2 | 2 | 2 | 1 | 1 | 2 |
| 3 | 6 | 4 | 2 | 2 | 6 |
| 4 | 32 | 12 | 7 | 22 | 32 |
| 5 | 320 | 48 | 29 | 293 | 293 |
| 6 | 6,144 | 296 | 215 | 5,442 | 5,442 |
| 7 | 229,376 | 3,040 | 2,591 | 207,070 | 207,070 |

For N=4..7 the stable quotient is reached after three strict refinement steps.
All finite classifiers are saved. Through five an independent implementation
enumerates **literal** permutations and Hamilton reset transitions and obtains
the same partitions. Through seven the larger computation uses the proved
uniform-fibre reduction, exact integer probability ratios, and complete q data.

If the observer is strengthened from o to q, the exact all-size minimum is
(q(s),q(tau(s))). Necessity follows by observing now and after tau; sufficiency
follows because tau swaps the pair and reset depends only on q. This is the
pointed application of the inherited reset/involution theorem
[unnamed dissection, U2](../unnamed_carrier_dissection_20261026/README.md).
For N=5..7 the weaker observer eventually recovers this same partition.
Whether that equality persists for every N>=5 is **OPEN**, not inferred from
the three finite rows.

Equal child and parent classes alone are insufficient: at N=4 free codes 0 and
1 with distinguished vertex 0 have the same o, but one native complement gives
different child classes. The receipt records the rows and exact codes. Conversely
some pointed types genuinely merge at N<=4 under this deliberately weaker
protocol. This is an actual reduction, not a renamed full literal state space.

Adding arbitrary named Hamilton guards changes the answer. The inherited
**THM-4807 +** [full theorem slug](../../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md)
reconstructs the path presentation under that richer language. The 81/27 pair
has the same 67 labelled cyclic triples and 462,915 Hamilton paths, yet different
groups and different answers to a specific named path guard. Neither our weaker
quotient nor the self-converse class bijection overrides that obstruction.

## P5. Signed incidence connects the score map to label equations

Use positive roots e_i-e_j and e_i+e_j for i<j. Put a_ij=1 when i->j in T and
g_ij=G_ij. Define

    b_D = sum_(i<j) ((a_ij-1/2)(e_i-e_j) + (g_ij-1/2)(e_i+e_j));
    b_B = b_D + sum_i (s_i-1/2)e_i.

The even double's positive-vertex **centred outdegree** vector is b_D+1/2,
coordinatewise; the odd double's is b_B-1/2. Negative-vertex scores are their
negatives, and the odd fixed vertex has centred score zero. Count outgoing arcs:
the positive vertex has d_T(i)+d_G(i)+1 wins in the even case and
d_T(i)+d_G(i)+s_i in the odd case. Subtracting (|D|-1)/2 proves the formulas.
All 599 inputs for n=1..3, both parities, pass an independent integer check.

The B_n/D_n root columns are classical signed-incidence coordinates; related
Coxeter tournament embeddings and score-preserving interchange graphs are in
Kolesnik--Mitchell--Przybylowski, *Coxeter Interchange Graphs* (2025), Sections 2
and 6, [primary article](https://doi.org/10.1007/s00026-025-00768-9).
The displayed shifts come from **our** fixed pair/root normalization. They must
not be replaced by that paper's different embedding conventions.

The transpose of this same signed incidence matrix evaluates vertex labels:
(e_i-e_j)^T x=x_i-x_j and (e_i+e_j)^T x=x_i+x_j. Thus graceful edge equations
and square-sum edge equations are two column types of one integral operator;
tournament/graph bits use its forward score map. On a tree, fixed edge equations
x_v=+/-x_u+w determine x_v=epsilon_v*z+t_v. LCA transport recovers every omitted
pair relation; a chord tests a coefficient 0 or +/-2 of z. This precise
conditional solver is inherited from
[relation towers, T4](../relation_tower_frobenius_ports_20261028/THEORY.md).

The source/target map preserves the specified linear equations or score vector.
It destroys the kernel, ordering, and permutation/range constraints. To solve
graceful or square-sum existence one must still choose the edge data and enforce
distinct bounded integer labels. Reducing these signed columns modulo two
identifies + and -, so it loses the sign distinction needed for integral holonomy;
binary Hamming relations alone do not restore them or lost integer precision. This is a useful shared
operator, not a proof of either conjecture or an identification of every
arithmetic dependency with graph cycles.

## Concept-board update

The canonical root supplies a small *enumeration address*. Relative alignment
supplies a *gluing coordinate*. The reversing involution supplies the missing
*symmetry mark*. The ordered forest supplies the *native update coordinate*.
The observation quotient specifies exactly which coordinates may be forgotten.
Output size then determines whether the target should be a scalar, a compressed
bank, or a full transition structure. These are complementary roles; a saving
in one does not automatically constitute a saving in all the others.
