# The first missing-gap layer of zero-heavy tournament cores

2026-10-20. **PROVED:** an exact structural classification and class count
at mass deficit $k$, for every budget $k\ge1$, including the strong case.
**FINITE-EXACT:** complete record and weighted-isomorphism controls through
$k=4$, with independent small-budget comparisons. No fitted polynomial or
unrestricted tournament enumeration claim is used.

Artifacts: [checker](verify.py), [receipt](receipt.json).

## 1. Inheritance, correction, and connection contract

The closest mechanism is the endpoint/gap mass identity in
[the previous zero-heavy mass theorem](../charge_palettes_boundary_20261020g/zero_heavy.md).
It proves the sharp maxima and counts mass deficits $0\le\delta\le k-1$.
The canonical hostile is the next value at budget two: extrapolation gives
45 records, but the complete compiler gives 42 classes. The
[full-position owner](../kernels_games_detours_20261020f/general_fixed_deficiency.md)
retains complete weights and external defect partners. The
[zero-heavy port-word theorem](../zero_heavy_ports_20261020g/zero_heavy.md)
supplies the complementary warning that arbitrary long unit interfaces
cannot be replaced by a fixed small separator.

The board is **weighted mass / distinct backward edges / endpoint overlap /
occupied gaps / unit cells / exact external partners**. The operation here
is a classification of the next mass layer, followed by the actual
weighted-isomorphism quotient. It is not clone contraction, deletion of
an unobserved port, or a transfer between feedback-vertex and feedback-arc
parameters.

Source: a pair-free weighted tournament core with weights in $\{1,\ldots,k\}$.
Target: an ordered perfect matching and a full numerical gap word.
Map: choose a good order, retain every reversed endpoint pair and every
gap weight, and use the intrinsic large-vertex order to recover the
record. Preserved predicates: the complete tournament, all weights,
pair freeness, weighted feedback at most $k$, and strongness when required.
Vertex names are quotiented only by weight-preserving isomorphism. The
needed sidecars are the complete gap word and every external matching
partner, including a partner outside a unit cell. The cheapest test is
the three rejected missing-gap records among the 45 budget-two records.

**Correction of a boundary description.** The earlier note's statement
that shared defect endpoints can enter already at $\delta=k$ was too
early. Missing gaps do enter there, but shared endpoints require
$\delta\ge k+1$. The sharper identity below proves this. A pair-free
budget-two example at deficit three is retained in the receipt, so the
distinction has an actual next-layer witness.

## 2. Exact statement

Use the inherited unique transitive-chain quotient, with a homogeneous
pair meaning two vertices seen identically by every third vertex. The
feedback number is the exact weighted minimum

$$
f(Q[w])=\min_\pi\sum_{ij\text{ backward in }\pi} w_iw_j,
$$

as proved in [the sparse weighted-core theorem](../defect_boundary_20261020/sparse_theory.md).
The weights here are all at most $k$, so there are no coordinates heavy
at the inherited cap $k+1$. They are nevertheless retained as integers;
weight one and weight greater than one play different roles below.

Set $M_k=2k^2+3k$ and $S_k=2k^2+k$. Let $Z_k$ count the nonempty,
pair-free weighted core classes with $f\le k$ and mass $M_k-k$.
Let $Z_k^{\rm s}$ count their nontrivial strong counterparts at mass
$S_k-k$. Write $I_k$ for the number of perfect matchings on $2k$ ordered
points crossing every proper initial cut:

$$
I_1=1,\qquad
I_k=(2k-1)!!-\sum_{j=1}^{k-1}I_j(2k-2j-1)!!.
$$

**The next-layer theorem.** For every $k\ge1$,

$$
\boxed{Z_k=(2k-1)!!\left[\binom{3k}{k}-1\right].}       \tag{1}
$$

For every $k\ge2$,

$$
\boxed{Z_k^{\rm s}=I_k\binom{3k-2}{k}-(2k-3)I_{k-1},} \tag{2}
$$

and $Z_1^{\rm s}=0$. Every admitted core in either statement has weighted
feedback number exactly $k$.

| $k$ | all classes $Z_k$ | strong classes $Z_k^{\rm s}$ |
|---:|---:|---:|
|1|2|0|
|2|42|11|
|3|1245|344|
|4|51870|15490|

## 3. Completeness: only two record types survive

For any good order write $r$ for the number of distinct backward edges,
$t$ for their endpoint count, $W$ for total endpoint weight, and $g$ for
the number of nonendpoint vertices. Pair freeness permits at most one
nonendpoint vertex per endpoint gap. Let
$U=\sum_{\rm occupied\ gaps}(k-w_i)$ and $q=k-r$.
The inherited exact identity is

$$
M_k-v=(2k+1)q+(k+r-W)+k(2r-t)+k(t+1-g)+U.              \tag{3}
$$

All terms are nonnegative. For strong cores replace $M_k$ by $S_k$ and
$t+1-g$ by $t-1-g$; the first and last ordered vertices must be endpoints.

At deficit $k$, (3) first forces $r=k$, because one unit of $q$ costs
$2k+1$. There are $k$ positive integer backward products with sum at most
$k$, so every product is one. Every endpoint has weight one, and $W=t$.
With $d=2k-t$ the identity sharpens to

$$
k=(k+1)d+k(t+1-g)+U,                                  \tag{4}
$$

or $k=(k+1)d+k(t-1-g)+U$ in the strong case. The gap terms
remain nonnegative before setting $d=0$, so (4) forces $d=0$:
the backward graph is a perfect matching on $2k$ unit endpoints. There
are precisely two possibilities:

1. Every allowed gap is occupied, with deficits summing to $k$.
2. Exactly one allowed gap is missing, and every other gap has weight $k$.

For the strong construction the allowed gaps are the $2k-1$ interior
ones; there are no exterior gap vertices. For the unrestricted
construction there are $2k+1$ allowed gaps.

These conclusions hold for **every** good order of the class. In
particular a good order of cost less than $k$ would have $r<k$ and
contradict (3). This proves exact feedback $k$ and shows that every
weight-greater-than-one vertex has all incident arcs forward in every
good order.

## 4. Pair freeness and the intrinsic owner

With every gap occupied, the matching construction is pair free, as in
the inherited theorem. If one gap is missing, it fails to be pair free
**if and only if** that gap is interior and its two adjacent endpoints
are matched to one another.

For the forward implication of the bad case, the two adjacent endpoints
have identical incidences to everything outside them. Conversely,
two gap vertices are distinguished by an endpoint between them; two
endpoints are distinguished by an occupied intervening gap unless they
are the adjacent pair at the omitted gap. If these adjacent endpoints
have different external partners, a partner distinguishes them. A gap
vertex and an adjacent endpoint are distinguished by the endpoint's
partner. For a nonadjacent gap/endpoint pair, an intervening vertex
distinguishes them unless the only such vertex is that partner; this is
exactly the same bad omitted-gap configuration. This covers all pairs.

It remains to count classes rather than numerical presentations. In any
good order all vertices of weight greater than one are untouched by
backward arcs. Their induced tournament is transitive with an intrinsic
order. Incidences to these ordered vertices locate each remaining unit
vertex in a fixed unit cell. No outside matching partner is dropped.

For $k\ge3$, a full-gap record has at most one unit gap, because two such
gaps would cost $2(k-1)>k$. Its nontrivial unit cell has at most three
vertices, exactly as at the inherited $k-1$ boundary. Outside partners
identify all roles unless the two local endpoints are paired together.
In that case the cell is a directed triangle; its three rotations change
the chosen gap vertex but preserve the **same numerical matching and
weight record**. Hence no division by three is warranted.

A valid missing-gap record has all occupied gaps of weight $k\ge2$.
It has at most one two-unit cell, at an interior omission. The arc in
that cell fixes its order; if its vertices were paired together the
record would already have been rejected. Exterior omissions introduce
only a singleton extreme cell. Thus these records also have unique
numerical owners. Full and missing records have different vertex counts
and cannot become isomorphic to one another.

### The two-unit-gap case at $k=2$

In a full record the two unit gaps may be separated, giving only the
preceding singleton/three-vertex cells. If two such cells have cross-paired
endpoints, each local triple is transitive and its intrinsic order fixes
the roles; otherwise the preceding internal-triangle argument applies.
If adjacent, there are just two
new local shapes.

At an exterior end the four-unit cell consists of two endpoint roles and
two gap roles. The endpoints are either distinguished by fixed external
partners, or paired internally. In the latter case the cell is
$1\to C_3$ or $C_3\to1$; its unique extreme vertex fixes the gap location,
and triangle rotation preserves the numerical record.

An interior cell, also possible at a strong end, has the ordered shape
$a,g,b,h,c$. Exactly one of the three endpoints has the fourth endpoint
as its fixed external partner. If that pinned endpoint is $a$ or $c$,
the other two endpoints and the intervening gap form a triangle module;
the remaining ordered singleton distinguishes the placement. If the
pinned endpoint is $b$, its in/out split separates the two pairs
$\{a,g\}$ and $\{h,c\}$; the unique reversed cross edge $c\to a$
identifies $a,c$ and hence the gap roles. These exhaust the local cases.
Again distinct records are not identified. The checker independently
tests every weight-preserving isomorphism class in this small case.

At $k=1$ all weights are one. The three possible missing gaps give
$1\to C_3$, a transitive four-vertex clone record, and $C_3\to1$.
The first and last are distinct pair-free classes. The sole strong
record is the rejected middle omission, so there are no strong classes.

## 5. Counting the rejected records, including strongness

There are $g$ gap deficits in $\{0,\ldots,k\}$ with sum $k$.
Stars and bars counts $\binom{g+k-1}{k}$, because the cap cannot be
exceeded at this total. A deficit equal to $k$ means the unique missing
gap. All other records have every gap occupied.

In the unrestricted case, a rejected record is a matching with one marked
adjacent matched endpoint pair and its intervening gap omitted. There
are $2k-1$ locations and $(2k-3)!!$ remaining matchings. Thus exactly
$(2k-1)!!$ records are rejected. Taking $g=2k+1$ gives (1).

In a valid strong record every consecutive arc of the displayed full
vertex order is forward. A tournament with this Hamilton path is strong
exactly when every proper cut has a backward crossing: its strong
components must occur contiguously along the path. Here this is precisely
indecomposability of the endpoint matching. Hence the raw strong count
is $I_k\binom{3k-2}{k}$.

For $k\ge2$, deleting a marked adjacent pair from an indecomposable
matching leaves an indecomposable $(k-1)$-matching. The pair could not
have been exterior; conversely it can be inserted into any of the
$2k-3$ internal gaps of an indecomposable $(k-1)$-matching. This is an
inverse pair of maps, so the rejected strong count is
$(2k-3)I_{k-1}$, proving (2). It retains the marked location, rather than
dividing by a guessed automorphism factor.

At $k=2$, the raw all/strong counts are $45/12$ and the clone rejections
are $3/1$, giving $42/11$. The discrepancy is therefore a specific
missing-gap clone obstruction, not an unexplained failure of a fit.

## 6. Exact tests, scopes, and the next boundary

The checker generates every structural record at $k=1,2,3,4$, including
all gap compositions, matchings and rejected clone records. It compares
the missing-gap criterion against every literal homogeneous pair. The
weighted quotient fixes the intrinsic transitive large-vertex order and
exhausts all permutations inside each remaining unit cell. An independent
quotient at $k\le2$ instead uses full outdegree/weight cells and exhausts
all permutations there. These give the same class counts. Reversed and
cyclic relabellings test the previously corrected mutable-sort-key failure.
Small exact subset dynamic programs independently test weighted feedback.

The complete older budget compiler's mass coefficients agree at $k=1,2$.
They are independent controls, not assumptions in the new completeness
proof. The receipt also records exact-type rejection controls, the raw
$45$ versus $42$ hostile, and a first shared-endpoint example at the next
layer: on seven ordered vertices of weights $(2,1,2,1,2,1,2)$, reverse
edges $(1,3),(3,5)$ using zero-based full positions. It is pair free,
has exact feedback two, and mass eleven, hence deficit three.

Thus the new classification closes exactly $\delta=k$. At $\delta=k+1$
shared endpoints really become possible; this note makes no formula
claim there. No result here counts all zero-heavy cores, all tournaments,
or all feedback-vertex classes. The full port and clone information is
precisely what makes this next layer tractable.

```
python -B 05-knowledge/results/zero_heavy_next_20261020h/verify.py
python -B -O 05-knowledge/results/zero_heavy_next_20261020h/verify.py --verify
```

Normal and optimized verification agree on **224,169 active checks** and
the saved receipt, LF SHA256
`ab4172650a0c387a39fbf4f6c450105e20a264896079c66c4b130aca16ec5391`.

Prepublication proof correction: equation (4) initially replaced $t$ by
$2k$ in its gap term before proving $d=2k-t=0$. Keeping $t+1-g$
(respectively $t-1-g$) gives the exact identity and the stated
nonnegativity argument. The checker already used the correct lower bound;
its formulas, counts, and receipt are unchanged.

Independent structural review of the owner proof, including all two-unit-gap
cases at $k=2$, passed. The root agent independently identified the same
equation-(4) transcription issue before promotion; its correction is
recorded above.
