# Paired ports, block repairs, and the two decompositions of ten

**PROVED:** three actual black moves bypass every aligned blue line in
`C3[TT3,V,TT3]`; the main formula is uniform for |V|>=2, with one fixed
certificate for |V|=1. The same formula repairs a nonmodule branch of the
closed `(n-4,3)` endpoint layer. A second pair of block words repairs the
entire adjacent-end-pair A branch at every n>=11, with no chord guards
inside its middle. Two further A-branch repairs retain a whole
prefix block or use the even order of a reflection commutator. Together
with the companion branches, Theorem E2 closes the full endpoint-three
layer at n>=11; see the [packet theorem](README.md). Global black-core
connectivity remains a separate open problem.

## 1. Inheritance: what 3-4-3 and 4-2-4 actually preserve

The closest proved mechanism is the
[two-block endpoint repair](../strong_blue_games_20261020f/detours.md#2-a-two-block-endpoint-repair).
Its sidecar is an endpoint-preserving Hamiltonian path in each reversed
block. The corrected near miss and canonical hostile are the
[oriented port-parity lift and its seven-vertex counterexample](../endpoint_two_repairs_20261020g/lift.md#4-exact-two-port-criterion-for-lifting-a-k-ledger):
macro edge XOR does not determine which physical entry/exit pair is used.
The least-used sidecar below is the order of the product of two positional
reflections. The concept board consists of exceptional ports, reversed-block
paths, prefix parity, reflection commutators, and the different meanings of
the two decompositions of ten.

The requested numerical decompositions already have precise but different
repository meanings.

* The rank levels of the ten proper nontrivial divisors of `p^2 q r` have
  sizes 3,4,3, from `(1+t+t^2)(1+t)^2`. Divisor complement reverses the
  levels. The older split into primes, squarefree composites, and
  square-containing divisors has the same sizes but exchanges `p^2` and
  `pqr`, and does not have this rank-complement property. See
  [the divisor note, section 4.3](../procgen_tcpc_20261001_tournament_clock_prime_collatz.md).
* The explicit address map `(a,b,c,d) -> (a+b,c,d)` folds the Boolean
  four-cube to `{0,1,2} x {0,1}^2`, preserves rank and complement, and has
  fibers of size two at first coordinate one. Removing bottom and top
  gives the rank split 3,4,3. The same dimensions occur in
  `Sym^2(U+V)=Sym^2(U)+(U tensor V)+Sym^2(V)` for two-dimensional U,V.
  These are established elementary maps in
  [the address-lattice note, section 3](../collatz_finite_address_lattice_20261015c.md).
* The actual 4,2,4 decomposition of ten in
  [the boundary-graph note, section 2](../odd_zeta_boundary_graphs_20261019.md)
  allocates marked points to a three-node stable-tree chain. Its vertex
  valences are 5,4,5, giving `Mbar_0,5 x Mbar_0,4 x Mbar_0,5` and thus two
  Petersen factors in that note's boundary model. The same chain with
  mark blocks 3,4,3 instead has valences 4,6,4. This last observation is
  the direct flag-count inference, not an identification of the factors.

Below, 3,4,3 is instead an ordered **tournament substitution**: a cyclic
three-vertex macro tournament, transitive three-vertex end blocks, and a
four-vertex middle tournament. The substitution specifies all cross-arcs.
It preserves the aligned anti-involution when the middle path is blue.
Neither the divisor map nor the stable-tree mark allocation specifies
these arcs. The common integer ten alone supplies no graph isomorphism,
minor, colouring theorem, or Pythagorean invariant.

## 2. A uniform mirrored three-vertex repair

Write `R_Q(T)=T^Q` for reversing the path arcs and `K_Q(T)=(T^Q)^op`.
Three R words Q,U,V correspond to K words Q,reverse(U),V. A word is blue
exactly when its positional reversal is an anti-automorphism.

Let V have r>=2 vertices and a blue directed path M. Form

    A=(a0,a1,a2), Z=(z0,z1,z2),
    A -> V -> Z -> A,
    P=(a0,a1,a2,M,z0,z1,z2),

where A and Z are transitive in their displayed orders. Thus n=r+6 and
the endpoint scores are `(r+2,3)=(n-4,3)`. P is blue. Both endpoints of
its blue edge are strong: P closes in T, and
`(A,Z,M,a0)` is a Hamiltonian cycle in K_P(T).

**Theorem M (PROVED).** Define

    H=(a2,reverse(M),z1,z0,z2),
    Q=(a0,a2,M,z0,z1,z2,a1),
    U=(a0,a1,H),
    W=(a1,reverse(H),a0).                              (1)

These are three black R moves with exact labelled endpoint R_P(T).

**Proof.** In the right block `(a2,M,z0,z1,z2)` after reversing its path,
H is directed and has the required endpoints. The first join goes from
a2 to last(M), an unchanged cross-arc because r>=2. M is traversed
backward; first(M)->z1 is unchanged; z1->z0 is reversed; and z0->z2 is
an unchanged transitive chord. The left block is the singleton a0.
The inherited two-block XOR proof therefore proves all three words are
directed and gives their exact endpoint.

For blackness, let b=first(M). Q's reversal exchanges b,z1 and z2,a2.
But b->z2 holds and a2->z1 does not, contradicting the anti equation.
The middle move has endpoint-score sum n-2 and is black. Just before W,
z2->a0 remains unchanged, whereas a1->a2 has been reversed by U and
was absent from Q. W exchanges z2,a2 and a0,a1, so these two arcs again
contradict the anti equation. This includes r=2 without a score exception.
The score sums, as a consistency check, are `2r+3,r+4,2r+2`. QED.

At r=1 there is one source, with standard-path adjacency rows
`(14,12,8,112,103,71,7)`. The three K words

    (0,2,3,4,5,6,1),
    (3,0,2,4,5,6,1),
    (3,0,1,2,4,5,6)

give the exact labelled target. Their endpoint excesses are -1,0,1;
at the balanced second move the anti equation fails on 3->4 versus
4->1. This fixed certificate was supplied by the companion probe lane
and independently replayed here. Thus the assertion for all r>=1 has
an all-order proof and a single explicit small boundary.

## 3. The same repair without a middle module

Let `P=(a,b,x,W,h,c,z)` be standard-blue, with n>=8 and

    N^-_T(a)={h,c,z},   c->x,
    b->last(W),        x->last(W).                      (2)

Thus h,c are the last two exceptional middle ports, W has n-6>=2
vertices, and blue reflection gives `N^+_T(z)={a,b,x}`.
No uniform cross-arcs or module property of W is assumed.

**Theorem E (PROVED).** There is an exact three-black R detour using

    H=(x,reverse(W),c,h,z),
    Q=(a,x,W,h,c,z,b),
    U=(a,b,H),
    V=(b,reverse(H),a).                                (3)

The same right-block proof applies. The first join of H is supplied by
(2); its join first(W)->c is the blue reflection of b->last(W). The
joins c->h and h->z are respectively a reversed path edge and the
reflection of the unchanged chord a->x. For Q, first(W)->z holds,
while x->c fails; its reversal pairs these two ordered incidences.
The middle move again has sum n-2. For V use unchanged z->a against
the reversed b->x. This proves directedness, blackness, and the literal
labelled target. Theorem M is its module specialization. The two chord
guards are sufficient; their failure is not a claim of disconnection.

There is a structural obstruction to deleting the sidecar from this
particular two-block method. Suppose x beats only first(W) inside the
right block and loses to h. Reversing the right block's displayed path
makes x a sink, so no endpoint-preserving Hamiltonian path starting at x
exists. These conditions are compatible with blue reflection and the
required exceptional pair; the checker records an explicit order-eight
source. Alternatively, if c dominates every vertex of W, its only
inneighbour in the reversed right block is z, which is prescribed as
the final vertex. Again the required H cannot exist. Neither observation
precludes a detour using different blocks or additional moves.

### 3.1 The entire adjacent-end-pair A branch

**Theorem E2 (PROVED).** Suppose n>=11, P=(a,b,x,W,y,c,z) is blue,
`N^-_T(a)={y,c,z}`, and c->x. Then three actual black K moves end
exactly at K_P(T). There is no assumption on b versus c, on the arcs
from b or x to W, or on any module structure.

Write `W=(v3,...,v_(n-4))` in the standard P coordinates. In both cases
use

    Q1=(a,x,W,y,c,z,b).                                (8)

It is directed because a->x and z->b hold. Its reflection pairs
first(W)->z with x->c; the former is true, the latter false. Thus Q1
is black at every order under consideration. Put S=K_Q1(T), B=K_P(T).
The exact neighbourhoods used below are

    N^-_S(z)={a,c,x},    N^-_B(z)={a,b,c,x},
    d_S(a)=d_B(a)=4.                                   (9)

**Case y->x.** Use

    Q2=(W,a,x,y,c,z,b),
    Q3=(W,a,b,x,y,c,z).                                (10)

Their legality and target are the variable-cut identity from section4,
with D=W and U=(y,c). In S, vertex4->a is true, and Q2's reversal sends
vertex4 to z. In B, vertex3->a is true, and Q3's reversal sends vertex3
to z. Both reversals send a to their position5. For n>=12 this position
contains vertex8, outside both neighbourhoods in (9), so both words
are black. At n=11, a itself occupies the fixed middle position5, but
its degree4 differs from the required anti-fixed score5. This handles
the boundary without enumerating its interior.

**Case x->y.** Use

    Q2=(z,W,y,x,b,a,c),
    Q3=(b,z,W,y,x,a,c).                                (11)

In S, the new joins z->first(W), y->x, x->b, b->a, and a->c are
complements of the corresponding original nonpath arcs. The W path
and last(W)->y remain retained. After Q2, the new joins b->z and x->a
are the complements of the previously retained Q1 arcs; all other
Q3 joins remain directed. For the exact endpoint, the added pairs
ax, zb, z-first(W), yx and ac cancel with their repeated occurrences;
Q2 contributes the original pairs ab and bx once. The other boundary
pairs and the internal W edges occur an odd number of times. Thus the
three undirected path-edge sets XOR to E(P), giving a literal labelled
K ledger.

For Q2, its reversal sends a to vertex3 and vertex5 to x. In S the arc
a->vertex5 is false (Q1 complemented it), while x->vertex3 is the
retained original path arc. This contradicts the anti equation. Vertex5
belongs to W once n>=9, so this branch already works from that order.
At the opposite endpoint B, Q3's reversal sends a to z and vertex3 to
x. Here a->vertex3 is false, whereas x->z is true, giving the final
anti contradiction. This proves both words black and finishes the
all-order theorem. QED.

The last-case ledger was found by a finite symbolic search over seven
whole blocks at order12. Its proof above retains only the necessary
boundary arcs and applies at all orders in the stated range. The search
did not enumerate tournaments. The earlier H_R obstruction still
correctly diagnoses why that particular two-block route fails; (11)
changes the route instead of assuming a nonexistent reversed-block path.

## 4. A whole-prefix A repair with one extra port

Number the standard-blue path `P=(0,1,...,n-1)` and put
`a=0,b=1,x=2,c=n-2,z=n-1`. Assume

    N^-_T(a)={h,c,z},   c->x,
    h in {2,...,n-3}.

Reflection gives `N^+_T(z)={a,b,k}`, where k=n-1-h. For a cut t with
4<=t<=n-2 put `D=(3,...,t-1)` and `U=(t,...,c)`. If

    h!=2,   h!=t-1,   t->x,                           (4)

the words

    Q1=(a,x,D,U,z,b),
    Q2=(D,a,x,U,z,b),
    Q3=(D,a,b,x,U,z)                                  (5)

are directed and end exactly at K_P(T). This is the earlier block XOR
identity with its actual additional port guard retained.

For reference, Q1 is black for n>=11 by the cyclic-cut criterion. If
the cut were blue, the square of cyclic successor would be an
automorphism. At odd n score-transitivity would force n=7. At n=2m,
the two parity orbits give `2m-4 <= (m-1)/2+(m-1)`, hence n<=10.
There is also a cheaper A-specific witness: Q1 pairs vertex3 with c
and z with x. Since c->x, it is black whenever h!=n-4, because then
3->z holds. The exceptional equality is a limitation of this witness.

**Theorem A1 (PROVED).** For n>=11, choose t=n-2. All three words in
(5) are black provided

    h not in {2,3,4,n-8,n-3}.                          (6)

Indeed, in S=K_Q1(T),
`N^-_S(z)={a,c,k}`, and 4->a holds when h!=4. Q2 pairs vertex4 with z
and sends a to vertex7. Its being anti would require 7 in {a,c,k},
which is excluded by h!=n-8. At the final endpoint B=K_P(T),
`N^-_B(z)={a,b,c,k}` and 3->a holds when h!=3. Q3 pairs vertex3 with z
and likewise sends a to vertex7, giving the same contradiction.
The remaining two exclusions in (6) are exactly the legality guards.

## 5. Odd prefix blocks: a reflection-order certificate

A tournament automorphism has odd order: an element of even order has
an involution power; a nontrivial involution swaps two vertices and
cannot preserve the orientation of their pair. Consequently, if a
tournament has an anti-involution rho and a proposed blue word has
reflection r, then the permutation r*rho must have odd order.

In B=K_P(T), the standard reflection rho remains anti. For Q3 in (5),
let l=|D|. Its position permutation rotates the first l+3 positions by
three and fixes the rest. When l is odd and `n>=2(l+3)`, the product
`r_Q3*rho` is supported on two disjoint reflected blocks. On each it is
a rotation of an even number l+3 of positions by three. Its cycle length
`(l+3)/gcd(l+3,3)` is even. Thus **Q3 is black independently of all
interior orientations and of h** whenever it is directed.

**Theorem A2 (PROVED).** Let n>=11 and require h!=2,3 and `4->x`.
Choose D=(3), t=4. Then (5) is a three-black detour.
In S, a has outneighbours `{x,h,c,z}` and z has three inneighbours.
Q2 pairs a,z, so the anti condition would require 4=3. Q3 is black
by the disjoint four-cycle argument above. This repairs, for example,
the h=4 and h=n-8 exclusions of A1 whenever this chord guard holds.

More generally take odd l>=3, `n>=max(11,2l+6)`, t=l+3, and require

    t->x,   h not in {2,4,l-1,l+2}.                    (7)

Then (5) is again a three-black detour. The legality exclusions are
2,l+2. Q2 pairs vertex4 with z and sends a to vertex n-l; the latter
is outside `{a,c,k}` unless h=l-1. The remaining exclusion h=4 ensures
4->a in S. The preceding even-order argument handles Q3. In particular
the h=3 branch is repaired whenever the first applicable even-indexed
vertex t>=6 beats x. No claim is made that such a vertex always exists.

## 6. The full-prefix final word and a genuine order-seventeen boundary

The following orbit-capacity improvement was discovered by the root lane
and independently audited here. It concerns the last word of A1, rather
than asserting that all its preceding words are legal or black.

**Theorem A3 (PROVED).** Assume the A-branch hypotheses of section4 and
n>=11. At B=K_P(T), the word

    F=(3,4,...,n-3,a,b,x,c,z)

is directed when h!=n-3 and is black except possibly at n=17 with
h in {3,9}. In particular this final-word test has no exceptions at
any n>=18. The two stated order-seventeen exceptions are real allowances:
the receipt contains a blue F example with h=3.

**Proof.** The only changed joins are (n-3)->a and x->c, supplied by
h!=n-3 and c->x respectively after nonpath complementation. If F were
blue, `g=r_F*r_P` would be an automorphism of B. For n>=11 its complete
permutation is

    g(0)=3, g(1)=4, g(i)=i+6 for 2<=i<=n-9,
    g(n-8)=0, g(n-7)=1, g(n-6)=2,
    g(n-5)=n-2, g(n-4)=n-1,
    g(n-3)=5, g(n-2)=6, g(n-1)=7.

Following these residue-six strings gives the following exact data for
the orbit O containing a. Here the final column concerns only the three
fixed ports; the variable port h may also belong to O.

| n modulo6 | size of O | O intersect {b,c,z} |
| --- | ---: | --- |
| 1,3 | n | {b,c,z} |
| 0,4 | n/2 | {b} |
| 2 | (n+1)/3 | {c} |
| 5 | (n+1)/6 | empty |

For example, at n=6m-1 the orbit is
`{0,3,9,...,6m-9}`; at n=6m+2 its traversal is
`(0,3,9,...,6m-3,6m,6,12,...,6m-6)`.
The displayed permutation similarly traces the other two rows, and in
the first row visits all vertices.

Since `N^+_B(a)={b,h,c,z}`, regularity of the induced orbit tournament
gives respectively n=9, n<=10, n<=14, and n<=17. In the residue-two
case n=14 is the only remaining n>=11; g then has cycle lengths
5,5,2,2, so it cannot be a tournament automorphism. In the residue-five
case n=11 likewise has 2-cycles. At n=17 the a-orbit is {0,3,9}; its
regular outdegree must be one, forcing h=3 or9. This proves the claim.

An exact witness for h=3, supplied by the root lane, has source rows

    (32758,9348,93336,130355,73762,91334,126366,65816,
     75318,118014,67728,78258,92150,18048,101778,81335,8195).

The checker verifies the complete tournament, blue source, endpoint
degrees, exceptional set, directed F, and blue F at B. Thus a blanket
all-order claim for this particular final word would be false. This
does not establish a metagraph bridge. In fact the companion movable-cut
lemma repairs this same source using k=13 and the three K words

    (0,1,...,12,16,13,14,15),
    (1,2,...,12,16,13,14,0,15),
    (1,2,...,14,0,15,16).

All three are independently replayed as black with the exact target.

## 7. Reproduction and limits

`lift_check.py` constructs tournaments directly, imports no producer
module, checks every displayed word in its actual intermediate graph,
checks the complete anti predicate, and compares the literal final arc
set with K_P(T). Its JSON receipt records exhaustive small blue middles,
deterministic larger controls, guarded nonmodule examples, A1/A2 controls,
and the inherited naive block-lift hostile expanded to TT3 end blocks.
It also checks the A3 permutation formula across all six residue classes
and independently replays the fixed order-seventeen blue-word witness.

    python 05-knowledge/results/endpoint_three_repairs_20261020h/lift_check.py
    python -O 05-knowledge/results/endpoint_three_repairs_20261020h/lift_check.py

The source/target/map of the positive connection is a marked cyclic
substitution -> an explicit labelled three-move ledger -> the exact blue
partner. The preserved predicate is directed black reachability. The
discarded macro information is repaired by physical port incidences and
the right-block Hamiltonian path H. The arithmetic and stable-tree
3-4-3/4-2-4 meanings do not supply that sidecar. The exhaustive all-order
case partition for n>=11 is the companion [packet theorem](README.md),
which combines these lemmas with the open-middle and movable-cut repairs.
The local formulas here do not by themselves establish global connectivity.
