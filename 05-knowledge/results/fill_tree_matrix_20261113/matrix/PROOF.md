# Marked matrices, spanning trees, and what a perfect fill can certify

**PROVED elementary identities / CITED inherited OCF / FINITE-EXACT controls.**
The square-zero matrix carrier and odd-cycle formula are existing mathematics,
recovered here with an endpoint-preserving tree bridge and a concrete hostile
to confusing subset coverage with coverage of integer path counts. This note
does not prove Chan–Pak or Collatz.

## Inheritance and scope

Closest proved mechanism: [THM-002, OCF](../../../../01-canon/theorems/THM-002-ocf.md).
Closest tree mechanism: [THM-1460](../../../../01-canon/theorems/THM-1460-arborescences-are-the-determinantal-relaxation-of-H.md)
and [THM-1745](../../../../01-canon/theorems/THM-1745-leaf-graded-arborescence-filtration-721-shadow.md).
Corrected near miss: [THM-1870](../../../../01-canon/theorems/THM-1870-the-cycle-count-spectral-boundary-is-the-hamiltonian-length.md),
whose false sub-Hamiltonian spectral extrapolation is overridden by its warning.
Least-used relevant sidecar: the entire matrix of starting and ending counts,
rather than its scalar sum. The [October square-zero packet](../../squarefree_cycle_inventory_20261019/SQUAREFREE_CYCLE_BRIDGES.md)
already inventories simple cycles with vertex marks; no rediscovery claim is made.

The broader comparison is already developed in
[THM-4986](../../../../01-canon/theorems/THM-4986-the-7-21-omission-proof-has-an-exact-spanning-tree-twin-block-floor-crevasse-lemma-and-the-first-creator-contrast.md):
connected simple-graph spanning-tree counts omit only2, because cycles realize
every integer at least3. The analogous tournament completeness statement stays
open. [HYP-9640](../../../hypotheses/HYP-9640-chan-pak-as-a-matrix-semigroup-orbit-and-one-defect-digraph-spectra.md)
already treats matrix-semigroup coverage and one-missing-pair spectra; its
finite bounds are not all-order coverage theorems.

## 1. The adjacency-power repair

For a loopless digraph with adjacency matrix A, `(A^k)_{ij}` counts directed
walks of k edges, allowing revisits. An n-vertex Hamiltonian path has n-1
edges; a Hamiltonian cycle has n edges. Neither can generally be recognized
from an unmarked power alone.

The transitive three-vertex tournament has rows `(6,4,0)`, where row bits
encode out-neighbors. Its Hamiltonian count is1 but `A^3=0`. A directed
triangle together with three isolated vertices has `tr(A^6)=3` but no
spanning path or cycle. These are exact hostile controls. For oriented
graphs, `tr(A^k)/k` does count simple cycles when k=3,4,5; repetition first
becomes possible at length6, through two traversals of a triangle.

Let `R=Z[x_1,...,x_n]/(x_1^2,...,x_n^2)` and `X=diag(x_1,...,x_n)`.
The variables commute; this is not ordinary exterior degree-one algebra.
Then the full endpoint matrix is

    K_rs = [x_1...x_n] (X(AX)^(n-1))_rs.                 (1)

Every walk contributes its vertex-visit monomial. A repeated vertex kills
the term; a surviving degree-n term visits every vertex exactly once.
Thus `H(D)=sum_rs K_rs`. Likewise, for n>=2,

    HC(D) = (1/n)[x_1...x_n] tr((XA)^n).                (2)

Division by n removes only the choice of cyclic starting point. These are
the standard marked-matrix identities, also discussed in
[Irving–Omar, §2](https://arxiv.org/html/2412.10572v3#S2), not new algorithms
with polynomial-size scalar state.

An ordinary-matrix version retains states `(S,v)` with `v in S` and arcs
`(S,v)->(S union {w},w)` when `w notin S` and `v->w`. This is an acyclic
digraph graded by `|S|`. Its paths from singleton states to full-set states
count (1). It has `n*2^(n-1)` states: the missing information in A is the
visited set. Hamiltonian paths become ordinary paths after that lift.

## 2. An endpoint-preserving weighted tree identity

Assign weight `x_i` to every arc `i->j`. Let

    W_ij=x_i A_ij,
    L=diag(column sums of W)-W,
    tau_r(x)=det L with row/column r deleted.

The directed matrix-tree theorem makes `tau_r` the sum over spanning
out-arborescences rooted at r, with monomial `product_i x_i^(outdegree(i))`.
Consequently the exact endpoint refinement of the tree relaxation is

    K_rs = [product_(i != s) x_i] tau_r(x).             (3)

**Proof.** Every spanning tree has n-1 arcs. The squarefree monomial in (3)
forces outdegree1 at every vertex other than s and outdegree0 at s.
A finite rooted out-tree with these degrees is exactly a directed path
starting at r and ending at s. Conversely every such path contributes1.
This also handles n=1 with the empty determinant equal to1. QED.

Equivalently reduce each cofactor in R: all genuinely branching trees
vanish, and the coefficient of the unique missing vertex records the path
endpoint. Setting all variables to1 in the ordinary polynomial instead
counts every spanning out-tree. That operation is not a homomorphism on R;
it must not be performed after forgetting the square-zero semantics.

This makes the relationship to a spanning-tree model precise. Existence of
an out-arborescence only requires reachability from its root; Hamiltonian
existence additionally forbids branching. The out-star on four vertices
has one spanning out-tree, three leaves, and zero Hamiltonian paths.
For the leaf polynomial `sum_A z^(number of leaves)`, inherited THM-1745
already identifies its coefficient of z with H. Equation (3) retains
both endpoints, which the scalar leaf polynomial forgets.

The same construction recovers all unmarked path-cover bin counts, not just H.
Expand `det(zI+L(x))` in the square-zero ring. For each coefficient of z,
sum the coefficients of its unique multilinear representative, and call
the result `F_D(z)`. Then

    F_D(z)=sum_k c_k(D) z^k,
    P_D(b)=sum_k b*(b-1)*...*(b-k+1) c_k(D).             (3a)

Here c_k counts unordered covers by k nonempty directed paths, and P_D(b)
counts covers in b labelled bins with empty bins allowed, using
`H(empty)=1`. This includes the unmarked two- and three-bin counts P2 and
P3 used by growing-defect constructions. Root-specific starting or ending
constraints such as S2, E2 and S3 require additional marks; they are not
recovered by the unmarked polynomial F_D(z) alone.

**Proof.** The directed matrix-forest theorem expands the determinant
over spanning rooted out-forests, with z per component and the product
of outgoing vertex weights. Square-zero reduction deletes exactly the
forests with branching, leaving disjoint directed paths. Each path has
one forced root, its initial vertex, so it is counted once. Assign its k
components injectively to b labelled bins in the falling-factorial number
of ways; this proves the second formula. QED. This is an elementary
marked refinement of the inherited matrix-tree/forest mechanism, with
no priority claim. It carries exponentially many subset coefficients;
an n-by-n display is not a polynomial-size scalar shortcut.

## 3. A genuine subset-cell fill, and its growth law

For a tournament T on vertex set V define its entire induced-count array

    P_T(x)=sum_(S subset V) H(T[S]) x_S,
    x_S=product_(v in S) x_v,       H(empty)=1.

In the commutative square-zero ring,

    P_T(x)=product_(v in V)(1+x_v)
           product_(C an odd directed cycle of T)(1+2x_V(C)).   (4)

Different cyclic orders on the same support are different factors.
**Proof.** Products of overlapping cycles vanish. A disjoint cycle
collection contributes weight2 per chosen cycle, and unused vertices
can be supplied by singleton factors. The coefficient at S is exactly
the inherited OCF for T[S]. QED.

This gives a precise spatial interpretation: cells of the Boolean cube
are labelled by subsets S; each carries the positive odd integer H(T[S]).
Modulo2 the array is exactly `product_v(1+x_v)`, so every subset-cell is
filled once modulo2. Over the integers additional occupation comes in
pairs, with the complete cycle inventory tracking the overlaps.

For one new vertex v, let U be the enlarged tournament. All new cycles
contain v, so no two can multiply nontrivially. Equation (4) gives

    P_U=P_T*(1+x_v+2 sum_(C in U, v in C, |C| odd) x_V(C)).     (5)

The old face is unchanged, the new face is filled by positive odd values,
and the exact increase at any cell is governed by induced-count weights
on complements of cycles. At the full cell this is the usual deletion
identity `H(U)=H(T)+2 sum_(C through v) H(U-V(C))`.

Thus there really is a fill maintained at every growth step. Its cells
are vertex subsets of one tournament. It does **not** imply that the
values on those cells, or across all tournaments, contain every allowed
odd integer. The following five-vertex hostile isolates the difference.

## 4. The three-page fake fill and its forced missing cycle

Take five labels0,1,2,3,4 and declare just three odd-cycle supports

    {0,1,2}, {0,1,3}, {0,1,4}.

Insert those three factors in (4). All32 subset coefficients are positive
and odd; multiplication respects the same nonnegative growth rule. Since
all three cycles intersect, the full coefficient is `1+2*3=7`.

But that array is not realizable by a tournament. Reverse every arc if
needed to make `0->1`. The three prescribed triangles then force
`1->j->0` for every page j=2,3,4. The page tournament has a directed
Hamiltonian path `a->b->c`, hence

    0 -> 1 -> a -> b -> c -> 0

is an unavoidable directed five-cycle. With transitive pages there is
exactly one such cycle and no other odd cycles, so the actual count is9.
The transitive-page completion has rows `(2,28,25,17,1)`.

More generally, with m transitive pages between the same two hubs, every
odd page subset supplies one odd cycle, all intersect at the hubs, and

    H = 1 + 2 * sum_(j odd) binom(m,j) = 2^m+1, m>=1.       (6)

The controls reproduce this for m=1..8. This is a familiar cyclic blow-up
family, not a newly discovered count family. What matters here is the
specific missing constraint: local cycle declarations force further
cycles when orientations are completed. Positive odd occupancy and
cycle-packing positivity alone cannot certify tournament realization.

The same concrete family supplies5 at m=2 and17 at m=4; its first values
are3,5,9,17,33. Thus5 and17 already occur in a consistent cycle-completion
geometry, while the attempted three-page value7 is corrected to9 by the
forced five-cycle. This is a structural comparison, not an implication
from Fermat-prime arithmetic to tournament realizability.

This argument uses a *shared edge*, not the invalid implication that
pairwise intersecting cycles share one common vertex. See the repaired
[THM-029](../../../../01-canon/theorems/THM-029-alpha1-gap.md) for that
historical correction and the all-order H!=7 theorem.

## 5. Missing arcs expose the tournament-specific constraint

The directed cycle C_m has exactly m Hamiltonian paths for every m>=3.
Thus both7 and21 are attained by strongly connected oriented graphs.
Even the slightest departure from tournaments suffices in small orders:

| Missing-pair digraph rows | Deleted arc from a tournament | H(D) | Two tournament completion counts |
|---|---|---:|---|
| `(16,1,2,7,14)` | `2->0` |7|9 and9|
| `(48,1,3,7,14,28)` | `5->1` |21|23 and41|

These are literal independently recounted witnesses, not an all-integer
one-missing-pair theorem. That latter coverage question is HYP-9640(D),
still open in its current source.

For an arbitrary missing ordered arc u->v, the exact repair cost is

    H(D+uv)-H(D)
      =sum_(S: u in S, v notin S) E_D(S,u) S_D(V\S,v),       (7)

where E counts paths ending at u and S counts paths starting at v.
Cut a new Hamiltonian path at its unique added arc to prove the formula;
concatenation proves the converse. Equation (7) is nonnegative and keeps
the endpoint sidecar needed for exact repairs. When both completions are
tournaments and H(D) is odd, each repair cost is even. In the7 witness
both costs are2; interpolating between its two tournament orientations
does not pass through7 at all.

## 6. What to transfer to Chan–Pak

The useful map is `actual tournament -> marked path/tree polynomials ->
subset coefficients -> scalar H`. The first two descriptions retain
vertex labels and endpoints; the scalar projection loses both. Perfect
subset fill, arborescence existence, and matrix-cone membership therefore
do not provide an inverse realizing a prescribed scalar.

A productive return construction must supply the missing inverse: actual
orientations or legal repairs, plus a proof that their reachable scalar
counts cover the required integers. Equations (3),(5),(7) are exact
tools for certifying such a construction; the three-page hostile is a
cheap test against an insufficient proposed closure law.

## Exact validation

Run from the repository root:

    python -B 05-knowledge/results/fill_tree_matrix_20261113/matrix/check.py
    python -O -B 05-knowledge/results/fill_tree_matrix_20261113/matrix/check.py

The program compares literal Hamiltonian permutations, marked walk
products and weighted tree cofactors on all4,165 loopless labelled
digraphs through order4. It also compares their weighted forest minors
with literal path permutations and all cut positions, then checks (3a)
against direct vertex-to-bin assignments for b=0,1,2,3. On all1,099
labelled tournaments through order5
it checks all33,866 subset coefficients against independent path DP and
all1,098 available birth factorizations. It also checks the fake/real
three-page inventory, the book family, both missing-pair witnesses and
the walk/tree hostile examples. No removable assertions are used.
The source hashes and all resulting finite data are saved in matching
normal/optimized receipts. These controls support the proofs above;
they are not a global spectrum census.
