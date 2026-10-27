# Parent-only root-switch operator for ordinary tournament enumeration

Status: **PROVED** elementary statements below; independently audited.
**FINITE-EXACT** executable checks are restricted to orders 1–6.
No novelty or practical asymptotic speedup claim is made. This is an operator
packaging and implementable refinement of the exact root-orbit formula already
in the earlier [weighted-core transport packet, §3](https://github.com/eliottcassidy2000/math/blob/b4234f728c0516b5fc521c001d69e26f3d2508d5/05-knowledge/results/weighted_core_transport_20261020/README.md), not a replacement discovery of that formula.

## 1. Inherited bridge and the computation it still needs

Let `U_n` be ordinary (not converse-merged) isomorphism classes of tournaments
on `n` vertices, `a_n=|U_n|`, and `a_0=1`. For a representative parent `T`, let
an attachment mask `S⊆V(T)` mean that the new root `x` beats exactly `S`.
Put

    R_N = disjoint union over T in U_(N-1) of 2^V(T)/Aut(T).

These are exactly pointed tournament classes `(U,x)` of order `N`. The proof
is deletion of the root: changing an identification of its deletion parent
is precisely the parent automorphism action. Consequently

    p_N := |R_N| = sum_T (1/|Aut(T)|) sum_(g in Aut(T)) 2^c(g),

where `c(g)` is the number of cycles on parent vertices, including fixed
points. This is a pointed count, not `N a_N`.

The inherited exact unpointing formula is

    a_N = (1/N) sum_[(U,x) in R_N] |Orb_Aut(U)(x)|.             (1)

The present contribution is an explicit local procedure for its orbit-size
summands using only parent-size rooted-deletion transport, together with a
particularly clean binary-rank certificate.

## 2. Construct Q without child canonicalization

For states `i,j in R_N`, take a representative `i=(U,x)` and define

    Q[i,j] = #{v in V(U) : [(U,v)] = j}.                       (2)

This is independent of representatives. Here is a noncircular way to evaluate
it from the complete `(N-1)`-parent bank:

1. Construct `U=T+x_S` explicitly, with root label `x`.
2. For each vertex `v`, form the `(N-1)`-vertex deletion matrix `D=U-v`
   and its attachment mask `S_v={w : v beats w}`.
3. Canonicalize the **pair** `(D,S_v)` by minimizing `(code(D),mask(S_v))`
   lexicographically over relabelings of its `N-1` vertices. Equivalently,
   canonicalize `D`, then minimize the transported mask under `Aut(D)`.
4. Find that key in `R_N`, and increment its column once.

No isomorphism or automorphism computation on `N` vertices occurs in these
steps. The elementary implementation first orders deletion vertices by
outdegree, then searches permutations within equal-degree groups. This is a
complete invariant because every isomorphism preserves degree, and every
remaining allowed order is searched. Lexicographic minimization of parent
code followed by mask must be joint; minimizing the two independently would
lose their alignment.

The complete state bank may itself be obtained by applying the same key to
all `2^(N-1)` masks on every parent. Thus no parent automorphism oracle is
required, although one could save work with one.

## 3. Exact block theorem

Fix a child class `U`. Let its vertex orbits be `O_1,...,O_r`, with sizes
`m_1,...,m_r`, all positive and summing to `N`. Its pointed states are exactly
these orbits. Regardless of the initially chosen root, rerooting at every
vertex visits state `j` exactly `m_j` times. Therefore its block in `Q` is

    Q_U = 1_r (m_1,...,m_r).

Different child classes have disjoint blocks. It follows, over the integers,

    Q 1 = N 1,       Q² = N Q,       tr(Q) = N a_N.            (3)

Over characteristic zero, each nonzero rank-one block contributes one to
rank, so

    a_N = tr(Q)/N = rank_Q(Q) = number of distinct rows of Q.  (4)

Its support blocks give the ordinary child equivalence relation. Choosing
one pointed state per distinct row supplies one representative child per
class for the next iteration. The representative is retained as its raw
adjacency matrix, without canonicalizing that child. Recursion therefore
needs only canonicalization on the current parents, never on the larger
children to recognize their isomorphism classes. This is not a complexity
miracle: parent canonicalization is real work, and the state bank is large.

The diagonal `Q[i,i]=m_i` is the missing weight in (1). To count only `a_N`,
one can compute the diagonal directly by comparing rerooted keys against
state `i`, without storing the whole matrix. The full matrix adds an exact
certificate and generates the next representative bank.

## 4. Binary field: no lost blocks, including even N

A tournament automorphism group has odd order. Indeed, an involution on
vertices would exchange some pair and reverse its arc, impossible; Cauchy's
theorem excludes even group order. Every vertex orbit consequently has odd
size. Reducing `Q_U` modulo two therefore gives the all-ones matrix `J_r`,
which has rank one for every positive `r`. Hence the all-order identity

    a_N = rank_F2(Q mod 2)                                   (5)

holds for **every** `N`, not just odd orders. This is not a graph Laplacian
formula. The known binary graph-Laplacian failure in the Hamilton relocation
carrier does not apply to this different, weighted rerooting operator.

Because `m_j` is odd, `r≡N mod2`. If `N` is odd, `Q mod2` is an idempotent
projection. If `N` is even, it is a square-zero differential. In the latter
case `im Q⊆ker Q`, and

    dim_F2(ker Q / im Q) = p_N - 2 a_N >= 0.                 (6)

In particular every even-order tournament has at least two vertex orbits.
This homology dimension is an exact symmetry-deficit statistic, not a closed
scalar recurrence for `a_N`.

For an arbitrary prime `ell`, the block rank is zero exactly when every
vertex-orbit size is divisible by `ell`. Thus

    rank_Fell Q = a_N - #{U : ell divides every |O|, O a vertex orbit}. (7)

If `ell` does not divide `N`, no block can vanish. If it divides `N`, loss
is possible but not mandatory. At `N=3`, the cyclic tournament contributes
`Q_U=[3]` and disappears mod3, whereas the transitive block `J_3` survives.
For `ell=2`, no block ever vanishes. This stronger oddness argument is needed;
the polynomial `Q²=NQ` alone does not establish binary rank when `N` is even.

## 5. Hostile controls and the state-compression boundary

* At `N=3`, the unique two-vertex parent has four attachment states. Three
  belong to the rigid transitive child and one to the directed triangle.
  Thus `p_3=4`, `a_3=2`, and dividing `p_3` by three gives the false `4/3`.
* Even **the exact same parent plus attachment Hamming weight** does not
  determine the needed diagonal weight. The two masks of weight one on the
  two-vertex parent produce roots of outdegree one in a transitive triangle
  and a cyclic triangle. Their `Q[i,i]` are respectively one and three.
  Lost information: which parent vertex the new root beats, not just how many.
* Selecting a minimum unrooted deletion card is not a root-orbit selector:
  at order two both roots delete to the singleton, but source and sink are
  distinct pointed classes. Keep the attachment mask aligned with the card.
* The already-known infinite family obstruction in session67 shows that
  current parent core, SCC sizes, and automorphism cycle index do not determine
  future symmetry-corrected family transitions. No contrary claim is made here.

Within this construction the exact row-support equivalence is the necessary
and sufficient gluing datum for identifying ordinary children. The number of
components/rank alone can count a level but cannot regenerate its child bank.
No assertion of a universally smallest encoding or impossibility of smaller
analytic recurrences is justified.

## 6. Adjacent-level attachment metagraph (secondary observation)

Let `B_n` be the bipartite multigraph with one shore `U_n`, the other `U_(n+1)`,
and one edge for each pointed child, joining its deletion parent to its child.
Parallel edges are retained: different root orbits can have isomorphic deletion
parents. The graph is connected. Starting from any parent, add a new source,
then delete one old vertex; this is a two-edge walk through a child. Repeating
until all original vertices have been replaced yields the transitive parent.
Every child has at least one deletion-parent neighbor. The empty/singleton
case is connected as well. Therefore its cycle-space dimension satisfies

    beta_1(B_n) = p_(n+1) - a_n - a_(n+1) + 1.

Equivalently `a_(n+1)=p_(n+1)-a_n+1-beta_1(B_n)`. This usefully identifies a
global gluing correction but is not advertised as a computationally closed
recurrence: constructing the right shore already identifies child classes.
The root-switch algorithm above avoids that circularity by working solely
on pointed states and parent-size transport.

## 7. Reproducibility and bounded evidence

Run from this directory:

    python3 check_root_switch.py
    python3 -O check_root_switch.py

The standalone script uses nonremovable validation checks and constructs every
row using only deletion parents. Full-child canonicalization is exclusively
a separately marked audit oracle, not a dependency of row construction or the
next-level representative bank. Independent audit used a separately implemented
full-permutation engine, checked every state and matrix entry through six, and
compared the entire labelled tournament universe through five. It also disabled
full-child canonicalization and successfully ran the production generation
through six after removing its explicitly separated audit block. Normal and
optimized replays agreed. No larger exhaustive census is claimed.

| N | a_N | p_N | tr Q | rank mod2 | rank mod3 | rank mod5 |
|---:|---:|---:|---:|---:|---:|---:|
|1|1|1|1|1|1|1|
|2|1|2|2|1|1|1|
|3|2|4|6|2|1|2|
|4|4|12|16|4|4|4|
|5|12|48|60|12|12|11|
|6|56|296|336|56|53|56|

The corresponding adjacent-level multigraph cycle dimensions are
`0,1,2,7,33,229`. At even orders 2,4,6 the binary homology dimensions are
`0,4,184`. All block identities are checked on every pointed state, not only
by comparing final totals. Ordinary counts agree with A000568 in this finite
range; the sequence itself is not an input to generation.
