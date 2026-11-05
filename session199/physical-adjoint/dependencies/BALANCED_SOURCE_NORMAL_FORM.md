# Balanced sources have a canonical two-vertex inverse

Session 198. **PROVED:** the nontrivial rooted start-one sources with equal
Hamilton-path and reserved-cover counts are exactly the graphs `J(O(W))`.
The source `W` is recovered canonically, at the graph level. Consequently their
scalar image is exactly `D+2U` on actual sources; neither the old endpoint count
nor the old deletion count occurs in this readout.

**FINITE-EXACT:** through 50,001 this balanced scalar image misses exactly the
42 values in Section 6, the largest of which is 3187. Every one of these values
except 7 and 21 already has a physical C1/C2/C3 repair.

**CONDITIONAL / OPEN:** a precise three-coordinate arithmetic coverage lemma
would complete the requested all-odd coverage theorem. The finite computation,
canonical inverse, and positivity do not prove that lemma.

## 1. Inheritance, domain, and the distinction that remains

A source is a nonempty rooted tournament `(V,c)` with exactly one directed
Hamilton path starting at `c`. Write

    (D,E,U,d) = (H(V), e_c(V), R_c(V), H(V-c)).

The reserved-cover count `R_c` counts ordered covers `(P,Q)` in which `Q` is
nonempty and starts at `c`, while `P` may be empty. Set `H(empty)=1`.
The singleton has response `(1,1,1,1)`. A balanced source here means a
**nontrivial** source with `D=U`; the singleton is a separate boundary case.

We use the complete physical grammar and positive coordinates proved in
`../bilinear-descent/PROOF.md`, Sections 1-3, published at commit
`6f3f93a3915ef2a7eba9a611f67a1bb75bd392c9`. Its inherited complete indexed grammar
is `session197/bilinear-c2-inverse/PROOF.md`, Section 1. Every nontrivial rooted
source graph is obtained from the two-vertex transitive source by the moves:

* `O`: append a universal sink, retaining the root;
* `Z`: append a vertex beating only the old root, retaining that root;
* `J`: append the same vertex, marking the newborn as root.

In particular this is a graph grammar, not just an orbit containing some
numerical profiles. Put

    (a,f,e,p) = (D-E-d, U-D+E-d, E, d-E).

These are nonnegative on nontrivial sources, `e+p=d>=1`, and

    D=a+2e+p,  U=a+f+2e+2p,  U-D=f+p.

The moves are

    O(a,f,e,p) = (a+e, a+2f+e+2p, 0, e+p),
    Z(a,f,e,p) = (2a+f+e+p, f, e+p, 0),
    J(a,f,e,p) = (f+p, 0, a+e+p, e).

The closest inherited mechanism is the rank-two identity `Z J = J O = M_2`,
where composition acts from right to left. The new assertion below is its
**converse on every balanced graph**, including a canonical graph inverse.
The hostile boundary is the singleton: it satisfies `D=U` but has no such
inverse. The unresolved sidecar is physical membership of the projected pair
`(D,U)`. We never infer it from inequalities alone.

## 2. Canonical graph normal form

**Lemma 1.** Every nontrivial source with `f=0` has a root of outdegree one.

**Proof.** The two-vertex seed has this property. A last `O` is impossible,
since its new `f=a+2f+e+2p` is strictly positive when `e+p>=1`. A last `J`
creates a root whose sole out-neighbor is the previous root. A last `Z`
preserves `f` and preserves the old root's outdegree. Induction on the graph
grammar proves the claim. This argument concerns the distinguished root of
the given graph, not a different realization of its profile. QED.

A nontrivial source with `E=0` has a universal sink. Indeed, a strong tournament
has a Hamilton cycle and hence a Hamilton path ending at every vertex. A weak
start-one source consists of its initial rooted strong component, or a
singleton, followed by singleton components: a nontrivial strong later
component would contribute at least three choices to the root-start path.
Its last component is therefore a universal sink. This is also the weak-tail
part of the inherited complete source theorem.

**Theorem 2.** The map

    (W,c) -> J(O(W,c))

is a bijection from rooted source isomorphism classes, including the singleton,
to balanced rooted source isomorphism classes. If `W` has response `(D,E,U,d)`,
its image has response

    (h,e,h,e) = (D+2U, D, D+2U, D).                 (1)

Its inverse is canonical: if the distinguished root is `r`, let `c` be its
sole out-neighbor; in `V-r` delete its unique universal sink `z`. The remaining
graph, rooted at `c`, is `W`.

**Proof.** The response laws give

    O(W)=(D,0,2U,d),
    J(O(W))=(D+2U,D,D+2U,D),

so every source produces a balanced source. The construction remains valid on
the singleton and then gives the directed triangle.

Conversely, balance and nonnegativity give `f=p=0`. Lemma 1 gives a sole
out-neighbor `c` of the root `r`. Consider a last graph-grammar step. It cannot
be `O`. If it is `Z`, its newborn is a universal sink after deleting `r`. If it
is `J`, the new value `p=0` says that the predecessor had `E=0`, so that
predecessor has a universal sink by the weak-tail observation. Thus in both
cases `V-r` has a universal sink `z`. It is unique, as a tournament cannot have
two distinct universal sinks.

Every root-start Hamilton path in `V` begins `r,c`. Deleting `r` therefore
identifies these paths with the Hamilton paths of `V-r` starting at `c`, so
there is exactly one such path. The sink `z` is its last vertex and is different
from `c`. Removing `z` gives exactly one `c`-start path in `W=V-{r,z}`. Thus
`W` is an actual source. All its vertices beat `z`; all of `W` except `c`, and
also `z`, beat `r`; and `r` beats only `c`. These are exactly the edges of
`J(O(W))`. The root and the two deleted vertices determine the inverse without
choices, proving the claimed bijection and uniqueness. QED.

In particular every balanced source is strong: a path from `c` through `W`,
then `z`, then `r`, then `c` is a spanning directed cycle. The inverse reduces
`h` to a predecessor path count `D<=h/3`, since every source has `U>=D`.
Equality is possible, for example when `W` is a triangle. This is a descent
for an already supplied balanced graph, not a construction for an arbitrary
integer `h`.

## 3. Exact unsigned scalar projection

**Corollary 3.** The Hamilton-path counts of balanced sources are exactly

    B = { D+2U : (D,E,U,d) is an actual nonempty source }.       (2)

For a prescribed odd `n>=3`, balanced membership is equivalent to existence
of an actual source with

    D odd, 1<=D<=n/3, U=(n-D)/2.                   (3)

The old values `E,d` need not be specified in the target equation, but their
physical existence is still required. The endpoint count of the resulting
balanced source is exactly `D`, not a freely selectable parameter.

One can remove both hidden coordinates from the **arithmetic statement of the
open lemma**, by defining a smaller exact orbit. Let `P` be the least set of
integer triples containing `(1,1,1)` and closed under

    O(D,E,U) = (D,0,2U),
    T_L(D,E,U) = (g(L)D-(L mod 2)E+LU,
                  D-[L=1]E, Q(L)D+B(L)U),  L>=1.             (4)

The scalar functions are explicitly defined by

    (g(1),Q(1),B(1))=(1,1,1),
    (g(2l),Q(2l),B(2l))=(g(l),Q(l),2B(l)),
    (g(2l-1),Q(2l-1),B(2l-1))
       =(g(l)+l-(l mod 2), Q(l)+B(l), 2B(l)),  l>=2.

The complete indexed physical source theorem implies, in both directions,
that `P` is precisely the projection of all physical source responses onto
`(D,E,U)`. Discarding `d` here is legitimate because every transition in (4)
is independent of `d`; every finite arithmetic word has its original graph
recipe. Thus (2) is also exactly `{D+2U:(D,E,U) in P}`. This is not an assertion
that an elementary region in three-space is physical.

## 4. A finite exact inverse and a sufficient arithmetic lemma

Put `G=D+2U`. In positive coordinates it is `(3,2,6,5)` dotted with `(a,f,e,p)`.
After `O,Z,J` the respective rows are

    (5,4,10,9),  (6,5,9,9),  (6,3,11,9).

All dominate `3/2` times the old row; the `O` row even dominates `5/3` times it.
Therefore a forward enumeration pruned at `G<=N` loses no predecessor of any
accepted source. The singleton contributes `G=3`; nontrivial sources start at
the transitive edge with `G=5`. The number of retained letters is at most
`log_(3/2)(N/5)` after the edge seed. These claims do not assert sharpness of
all three individual constants.

There is also a wholly three-coordinate exact enumeration using (4). Its
bound `L<=floor(N/(3U))` is complete, because `B(L)>=L` and
`g(L)D-(L mod2)E>=0` give `G(T_L)>=3LU`. Every physical indexed move increases
`G`, either by compiling it into the one-vertex moves or directly from (4).
The singleton first step is harmless and increases `G` from 3 to at least 5.
Thus both enumerations terminate, return the same scalar image, and retain
physical witnesses. This is a decision procedure, not a positive answer for
every input.

Here is a precise stronger arithmetic statement sufficient for the requested
all-target result:

> **Balanced cofinite lemma (OPEN).** For every odd integer `n>3187`, the orbit
> `P` defined by the explicit recurrences (4) contains a triple with `D+2U=n`.

**Conditional consequence.** If this lemma holds, every odd integer `n>=3`
other than 7 and 21 is realized by a strong tournament with at most three
Hamilton cycles. Above 3187, (1) gives a strong start-one source, and hence a
strong C1 graph (its unique root-start path closes along the universal-sink
construction). At or below 3187, the finite verification in Section 6 and the
40 repair witnesses settle all remaining values. No finite congruence cover,
finite menu, or descent of an existing representation proves the open lemma.

## 5. An all-order physical rebase and the exact residual graph family

Write `F(W)=D-E+U`. The inherited C1 classification, restated in
`../c1-monotone/PROOF.md`, Section 1, says every strong C1 graph is `J(V)` for
an actual nontrivial source `V`: delete an outdegree-one vertex and mark its
sole out-neighbor. This statement concerns actual graphs, not just scalars.

**Theorem 4.** For every source `W`, including the singleton, and every `k>=0`,

    J(Z^k(J(W)))  is isomorphic, after forgetting the root, to
    J(O^(k+1)(W)).                                             (5)

The balanced root in the left graph is the **inner** `J` vertex. Consequently

    F(Z^k(J(W))) = D + 2^(k+1) U.                              (6)

**Proof.** Name the inner `J` vertex `r`, the `k` subsequent `Z` vertices
`z_1,...,z_k`, and the outer `J` vertex `z_(k+1)`. Every vertex of `W` beats
every `z_i`. The `z_i` form a transitive tournament in their displayed order.
Every `z_i` beats `r`; `r` beats only the original root `c` of `W` and loses
to every other vertex of `W`. These are exactly the edges obtained by appending
`k+1` universal sinks to `W` and then `r`. The isomorphism fixes `W`, sends the
ordered sink vertices to `z_1,...,z_(k+1)`, and sends the final root to `r`.
Thus (5) is a literal all-edge statement. Applying (1) after `k` sinks proves
(6). The `k=0` case is the familiar two-`J` rerooting; no nontriviality assumption
on `W` is needed. QED.

**Theorem 5 (exact graph dichotomy).** Every strong C1 graph either admits a
balanced root or has the form

    J(Z^k(O(W))),  k>=1, W a nontrivial source.                 (7)

Every graph in (7) admits **no** balanced root, although its scalar may equal
that of another balanced graph. Its scalar is exactly

    R_k(W) = D + 2(2^(k+1)-1)U - (2^(k+1)-2)d.                (8)

Thus the endpoint count `E` disappears from the remaining rebasing problem.
For a transitive `W` of order `m+1`, `m>=1`, (8) becomes

    (2^(k+1)-1)(2^(m+1)-1)+2.                                (9)

**Proof.** Express the source `V` in a presentation `J(V)` by the complete
`O/Z/J` grammar from the transitive edge. If its last letter is `O`, the graph
is balanced by Theorem 2. If its last non-`Z` letter is `J`, apply Theorem 4.
Otherwise it is `V=Z^k O(W)`. For `k=0` it is again balanced. If no earlier
non-`Z` letter exists, interpret the edge seed as `O` of the singleton. When
`W` is a singleton, the original root `c` is balanced: its sole out-neighbor
is the `O` vertex, and deleting `c` leaves the transitive chain consisting of
that vertex, the `Z` vertices, and the final `J` vertex. Its last vertex is a
universal sink, so Theorem 2 applies, or a direct reconstruction gives
`J(O^(k+1)(singleton))` with a different root. This exhausts the grammar and
leaves precisely (7).

For the converse let `W` be nontrivial, with root `c`; call its appended sink
`z`, and call the final `J` vertex `r`. The root `c` has an out-neighbor in `W`
and also beats `z`, so its outdegree in the final graph is at least two.
Every other vertex of `W` beats `z` and `r`. The vertex `z` beats the `k` new
`Z` vertices and `r`. Each `Z` vertex beats `c` and `r`. Thus `r` is the unique
outdegree-one vertex. Moreover its deletion has no universal sink: every
vertex of `W` has an outgoing edge to `z`, except possibly `c`, which also
beats `z`; `z` beats the first `Z` vertex; every `Z` vertex beats `c`.
By Theorem 2 no root can be balanced.

Iterating `Z` gives, for any input `(D_0,E_0,U_0,d_0)` and `k>=1`,

    D_k=D_0-E_0+(2^k-1)U_0-(2^k-2)d_0,
    E_k=d_0, U_k=2^k U_0-(2^k-1)d_0, d_k=d_0.

Substitute `O(W)=(D,0,2U,d)` and take `F=D_k-E_k+U_k`; this is (8).
Substitution `(D,U,d)=(1,2^m,1)` proves (9). QED.

The graph dichotomy does not say that the scalar images of its two sides are
disjoint. Indeed the graph in (7) with transitive `W` of order two and `k=4`
has scalar 95, and 95 also occurs in the balanced image. A target-dependent
replacement of (7) by a different source remains a concrete positive route.

### Simple balanced return maps do not yet suffice

For a balanced response `(h,e,h,e)`, two simple returns are

    Z:       (h,e) -> (2h-e,e),
    M_2 J:   (h,e) -> (6h-e,2h-e).                (10)

The latter means apply `J` first. Another return is
`M_2:(h,e)->(3h,h)`. These maps are physical and follow directly from (1) and
the source laws. They are useful possible steps of a future arithmetic descent,
but are not a complete balanced grammar. For instance the orbit of `(3,1)`
under these three maps misses the balanced scalar 19. A direct source
`O Z` from the transitive edge has `(D,E,U,d)=(5,1,7,1)`, so (1) gives a
balanced source of scalar 19. The verifier checks both this witness and the
restricted three-map omission.

The intended role of (10) is therefore a return mechanism, not a claim that
these low-dimensional maps already prove coverage.

## 6. Finite evidence and reproducibility

The verifier independently enumerates the four-coordinate `O/Z/J` orbit and
the three-coordinate indexed orbit, each pruned by `G<=50,001`, and compares
their entire projected profile sets and scalar sets. It also checks the graph
bijection on every labelled rooted source through order five, reconstructing
every edge under the canonical deletion map, and uses literal path/cover
enumeration as a separate small-order control.

The balanced gaps among odd integers from 3 through 50,001 are exactly

    7,11,13,21,23,41,43,47,55,67,77,83,149,163,173,231,233,247,251,
    277,299,317,391,421,433,443,449,503,519,533,569,583,617,637,821,
    911,1067,1339,1861,1903,2647,3187.

The 29 previously certified C1 gaps are contained in this list. The additional
13 balanced gaps are

    11,23,47,149,317,391,433,443,519,617,911,1903,3187.

`REPAIRS.json` gives actual C1 or C2/C3 source recipes for the 40 listed values
other than 7 and 21. Its C2/C3 records are inherited from the audited packet
`session198/c1-monotone`, commit `70a6809847b098f41d4da9296a807cbeb47e5c25`, and re-evaluated here. Its 13 C1
records are generated and checked here. This paper does not assert a new
nonexistence theorem for 7 or 21.

Run `python3 verify.py` and `python3 -O verify.py` in this directory. Correctness
checks use explicit exceptions, not removable Python assertions. The two
receipts must agree. The independent audit is distributed separately.
