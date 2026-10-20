# A finite defect grammar and eventual polynomial tournament counts

**PROVED in the scopes below. FINITE-EXACT:** the complete template compiler at
budgets zero, one and two, independently compared with all tournament classes
through order seven. No claim of external novelty. An efficient all-budget
recurrence for A000568 remains **OPEN**.

Let `f(T)` be the minimum number of backward arcs over all total orders of the
vertices of a tournament. Let `b(T)` be the order of its canonical quotient by
maximal transitive clone chains. A pair is homogeneous when every third vertex
sees both members identically. Such pairs form directed chains; their maximal
chains give the unique quotient and positive weights, up to quotient
isomorphism. This is inherited from
[THM-4463, homogeneous tournament contacts](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md)
and [the previous weighted-core packet](../weighted_core_transport_20261020/README.md).
“Core” below always means this quotient, not a modularly indecomposable tournament.

## 1. Sparse defects limit every automorphism

If `f(T)<=k`, every automorphism of `T` moves at most `3k` vertices.
More precisely, if its nontrivial vertex cycles have lengths `l_i`, then

\[
k\ \geq\ \sum_i\frac{l_i^2-1}{8}\ \geq\ \frac{\#\text{moved vertices}}3.
\tag{1}
\]

An even cycle in a tournament automorphism is impossible: its half-turn
reverses the pair of opposite vertices. On an odd cycle `l=2s+1`, the induced
tournament is regular, of outdegree `s`. Restrict any proposed total order to
these vertices. Its last `s` vertices require at least `s,s-1,...,1` backward
outgoing arcs, respectively. Thus the cycle contributes at least `s(s+1)/2`
distinct backward arcs. Internal arcs of different cycles are disjoint.
Finally `(l^2-1)/8 >= l/3` for odd `l>=3`, with equality only at three.

Ordered joins of `k` directed triangles attain support `3k` at budget `k`.
Consequently a tournament in this class has **no fixed-point-free automorphism
when `n>3k`**. For a fixed budget, the rootless cycle-index boundary only uses
nontrivial cycle mass at most `3k`, with
`sum m_l(l^2-1)/8 <= k`. This is a necessary bound, not a sufficient test for a
particular prescribed tournament or weighted core.

## 2. The core-size bound is also sharp

\[
\boxed{\quad f(T)\leq k\quad\Longrightarrow\quad b(T)\leq4k+1.\quad}
\tag{2}
\]

Choose an order with at most `k` backward arcs. Mark their at most `2k`
endpoints. Every interval of unmarked vertices is a homogeneous transitive
chain: all of its incident arcs agree with the order. There are at most `2k+1`
such intervals. Compressing these intervals, and retaining each marked vertex,
leaves at most `4k+1` vertices. The canonical quotient can only compress further.

The ordered join
`K1 -> C3 -> K1 -> C3 -> ... -> C3 -> K1`, with `k` triangles, attains equality.
It is pair free: distinct strong components separated by a triangle cannot
merge, and each triangle distinguishes its own vertices. Each triangle needs
one backward arc, while the displayed order achieves exactly `k`.

This does not give a constant-size core for the full tribonacci family.
That family has an economical grammar with an **unbounded** defect budget.
“Bounded grammar” and “bounded number of defects” are different parameters.

## 3. Weighted feedback straightens exactly

Write `Q[w]` for the transitive inflation of a pair-free core `Q` with positive
integer chain lengths `w_i`. Then

\[
f(Q[w])=\min_{\pi\text{ an order of }V(Q)}
       \sum_{\substack{i\to j\text{ in }Q\\j\text{ before }i\text{ in }\pi}}
                w_iw_j.
\tag{3}
\]

The upper bound comes from ordering the blocks contiguously. For the reverse
bound, take **any** order of all inflated vertices, including interleavings.
Choose one uniformly random representative from each block and order blocks
by these representatives. For a pair of blocks, its expected weighted
backward cost equals the number of backward cross-arcs in the original order.
Summing gives expected block cost equal to original cross cost, at most its
total cost. Some sampled block order therefore costs no more. This proves (3)
without assuming a priori that optimal orders keep modules together.

For fixed `k`, put `K=k+1` and cap each weight:

\[
c_i=\min(w_i,K).
\]

Then

\[
f(Q[w])\leq k\quad\Longleftrightarrow\quad f(Q[c])\leq k.
\tag{4}
\]

One direction follows by reducing weights in (3). In a capped order of cost at
most `k`, no backward arc can meet a coordinate of weight `K`. Increasing any
such coordinate has zero cost in that same order, proving the other direction.

Call coordinates of capped weight `K` **heavy**. A weight-preserving
automorphism of an admissible capped core fixes every heavy vertex. Otherwise
its nontrivial odd cycle permutes at least three whole chains of size `K`;
its lift to `Q[c]` moves at least `3K>3k` vertices, contradicting (1).
This pointwise fixing is why the following denominator has only `1-x`, rather
than factors `1-x^d` arising from permutations of unbounded weights.

## 4. An exact finite compiler, without an input core census

For a fixed budget `k`, generate the following **finite** objects.

1. An ordered endpoint set of size `t<=2k`.
2. A simple defect graph `D` on it, with no isolated vertices and at most `k`
   edges. These are the arcs to reverse in the reference total order.
3. Endpoint weights `1,...,k`, with
   `sum_{ij in D}w_i w_j <= k`.
4. For each of the `t+1` gaps, a weight in `0,...,K`. Zero omits the gap;
   `K` stands for an arbitrary length at least `K`.

Construct the ordered tournament, normalize its transitive clone chains, cap
the resulting weights at `K`, and deduplicate by weighted-core isomorphism.
Include the transitive cases separately. Call the finite set of resulting
nonempty capped weighted-core orbits `S_k`.

**Completeness.** Given any admissible capped core, choose an optimal order
from (3). Every endpoint of a backward arc has weight at most `k`; there are at
most `2k` endpoints, and at most `k` backward arcs. Each intervening unmarked
interval is homogeneous. In the pair-free core it contains at most one
vertex. Thus the core is represented by one of the templates above.

**Soundness.** Each template has a displayed order of cost at most `k`.
Clone normalization preserves the tournament; reducing any weight by capping
preserves admissibility by (3). Formula (4) then permits all expansions of its
heavy coordinates. Multiple ordered descriptions are removed by exact
weighted-core isomorphism, not by dividing by a generic symmetry factor.

For `s in S_k`, let `h(s)` be its number of heavy coordinates and `v(s)` the
sum of its capped weights. Since the cap stabilizer fixes each heavy
coordinate separately, there are no further identifications among its
independent expansions. Consequently, including the empty tournament,

\[
\boxed{\quad A_k(x)=\sum_{n\geq0}a_{n,\leq k}x^n
       =1+\sum_{s\in S_k}\frac{x^{v(s)}}{(1-x)^{h(s)}}.\quad}
\tag{5}
\]

This is an exact structural counting algorithm. It canonicalizes finitely many
bounded templates once for each `k`, and never generates all cores at the
target order `n`. It is not a claim that the finite compilation is efficient
as `k` grows.

## 5. Polynomial degree, leading coefficient, and a recurrence

Every heavy coordinate avoids the backward arcs in a good core order.
The at most `2k` marked endpoints leave at most `2k+1` gaps; pair freeness
allows at most one heavy vertex in each gap. Thus `h(s)<=2k+1`.
Equation (5) has denominator dividing `(1-x)^(2k+1)`, and its coefficients
are eventually polynomial in `n`, of degree at most `2k`.

The top-degree terms admit a complete classification. To have `2k+1` heavy
vertices, all gaps must occur, there must be exactly `2k` light endpoints,
and exactly `k` backward arcs with disjoint endpoints. Every endpoint weight
must be one: their `k` positive products have total cost at most `k`.
The defect graph is therefore a perfect matching on `2k` ordered points.
Conversely each such matching, with a heavy vertex in every gap, is pair free
and admissible. Two light vertices are distinguished by an intervening heavy;
two heavy vertices by an intervening light; an adjacent light/heavy pair by
the light vertex's matching partner. Other pairs have an intervening vertex
that distinguishes them.

These capped cores are pairwise nonisomorphic. An isomorphism preserves the
heavy vertices and their transitive order, hence fixes each heavy vertex.
Each light vertex is then identified by its unique gap among the heavies,
so the matching itself is recovered. There are `(2k-1)!!` matchings.
This also covers `k=0` with the usual empty-matching convention.

It follows that

\[
\boxed{\quad
a_{n,\leq k}=\frac{n^{2k}}{2^k k!}+O_k(n^{2k-1})
\quad(k\geq1),\qquad a_{n,\leq0}=1.
\quad}
\tag{6}
\]

In particular, for all sufficiently large `n`,

\[
\sum_{j=0}^{2k+1}(-1)^j\binom{2k+1}{j}a_{n-j,\leq k}=0.
\tag{7}
\]

An explicit, deliberately loose onset follows by bringing (5) to the common
denominator. A term on a core of order `q` contributes degree at most
`v+2k+1-h <= k(q-h)+(k+1)h+2k+1-h = kq+2k+1 <= 4k^2+3k+1`.
Thus (7) holds for `n>4k^2+3k+1`. The actual onset can be much smaller.

## 6. Two exact all-order outputs of the compiler

The checks in [sparse_count.py](sparse_count.py) produce

\[
A_0(x)=\frac1{1-x},\qquad
A_1(x)=\frac{1-2x+x^2+x^3}{(1-x)^3},
\]

and

\[
\boxed{\quad
A_2(x)=
\frac{1-4x+6x^2-3x^3-x^4+5x^5+2x^6-4x^7+x^8}
     {(1-x)^5}.\quad}
\tag{8}
\]

These formulas are **FINITE-EXACT finite compilations plus the all-order
completeness proof**, rather than fits to initial coefficients. Budget one
uses 27 raw templates and 20 capped core orbits; budget two uses 4,032 raw
templates and 2,691 capped core orbits. The latter count includes all finite
small-weight cases, not just distinct underlying cores.

For `n>=1`, `a_(n,<=1)=1+binom(n-1,2)`. Budget two begins, at `n=0`,

`1,1,1,2,4,11,33,82,173,324,556,893,1362,1993,2819,...`.

It satisfies the fifth-order recurrence
`a_n=5a_(n-1)-10a_(n-2)+10a_(n-3)-5a_(n-4)+a_(n-5)`
for `n>=9`, and more explicitly

\[
a_{n,\leq2}=\frac{n^4-6n^3-n^2+62n-72}{8}\qquad(n\geq4).
\]

This follows by expanding the binomial coefficients of (8); it is not an
interpolation argument. An independent complete class augmentation through order seven,
with subset-DP feedback minima, agrees at every budget `0,1,2`.
An independent inflated-vertex DP also checks (3) in 1,098 cases.

## 7. What this does and does not give for A000568

The filtration is exhaustive: every `n`-vertex tournament has feedback at most
`floor(n(n-1)/4)`, since the two opposite orders have complementary costs.
Thus the ordinary A000568 entry is exactly the coefficient of `x^n` in
`A_floor(n(n-1)/4)(x)`.

This supplies a recurrence in each fixed-budget layer and an exact diagonal
recovery of the full sequence. It does **not** supply a bounded scalar
recurrence for the full sequence, or prove that growing the budget is cheap.
The sharp remaining target is a recurrence on the finite endpoint templates
under insertion of one new defect, with duplicate normalization handled
structurally. Perfect matchings control its top-degree boundary; shared
endpoints, omitted gaps and clone mergers are its lower-degree corrections.

Connection contract: the source is an ordered tournament with bounded backward
cost; the target is a finite capped weighted-core orbit; the preserved predicate
is `f<=k` and the exact count of all positive-weight expansions. The chosen
order and redundant gap descriptions are destroyed. Canonical clone chains,
the cap vector and the pointwise-fixed heavy coordinates are the necessary
sidecars. The independent all-class feedback census is the cheapest decisive
control. No Collatz implication is asserted by this encoding.
