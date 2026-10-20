# Exact local ownership at deficiency three, and at every fixed deficiency

**PROVED below:** the deficiency-three menu; an exact symmetric-difference
rule for all clone mergers; an intrinsic local ownership algorithm with at
most `7!` candidate orders at deficiency three; and the corresponding
`(2D+1)!` theorem for arbitrary fixed deficiency. The resulting enumeration
does not require a core bank, a global set of previously seen classes, or a
general graph-isomorphism canonicalizer. It still enumerates ordered defect
patterns and retains their complete external incidences.

**FINITE-EXACT:** the accompanying independent controls give third-stratum
counts `5, 756, 36473` at budgets `1, 2, 3`. These values are checks of the
ownership theorem, not evidence for a fitted formula. An aggregate formula
for this stratum and a small scalar recurrence for unrestricted A000568
remain **OPEN**. No external novelty claim is made.

## Inheritance and the missing coordinate

The closest proved mechanisms are the canonical chain core, weighted
feedback straightening and capping in
[`sparse_theory.md`](../defect_boundary_20261020/sparse_theory.md), and the
deficiency identity and bounded exceptional motifs in
[`strata.md`](../defect_transitions_20261020/strata.md).
The earlier exact two-deficiency ownership proof is
[`strata_audit.md`](../defect_transitions_20261020/strata_audit.md).
Its canonical hostile example is the unit directed triangle with one external
defect port: one matching presentation and two shared-endpoint presentations
describe the same weighted core. The port must remain visible.

The correction to avoid is to regard a raw presentation's deficiency as
invariant under clone normalization. Mergers can both destroy old heavy
coordinates and create a new heavy coordinate by summing light weights.
The least-used relevant sidecar is the common symmetry-support theorem in
[`continuation_width_symmetry_20261020`](../continuation_width_symmetry_20261020/README.md).
Here a sharper location statement is available for isomorphism search:
the intrinsically ordered heavy vertices fix every light gap.

The live concepts are weighted chain cores, intrinsic heavy gaps, backward
edge neighborhoods, external ports, and canonical owners. The map from an
ordered defect template to its weighted core preserves the actual weighted
tournament and feedback budget, but destroys the chosen order and endpoint
roles. Retaining all gap vertices and all external defect incidences, then
minimizing over their admissible local orders, repairs exactly that loss.

Our parameter is **feedback deficiency** `D=2k+1-h`. It is distinct from
the other packet's **symmetry slack** `delta=3f-s`, and also from the incoming
center-deficit series `D(z)=z^4`. The symbol `delta` below means only the
graph defect `2r-t`.

## 1. Objects and the exact deficiency-three menu

Fix a feedback budget `k>=0` and cap `K=k+1`. A state is a nonempty,
pair-free tournament `Q` with positive vertex weights `w<=K`, satisfying
weighted feedback `f(Q[w])<=k`, modulo weight-preserving isomorphism.
A vertex of weight `K` is heavy; let their number be `h`.

Choose a good order whose backward weighted cost is at most `k`.
Let `E` be its undirected backward-edge graph, `r=|E|`, and `t` its number
of incident vertices. These are the defect endpoints. Every endpoint gap
contains at most one other core vertex, since two consecutive vertices
unincident to all defects would be a homogeneous pair.
Let `ell` be the number of light occupied endpoint gaps and `u` the number
of endpoint gaps that are not heavy, including empty gaps. Thus

    h=t+1-u,   0<=ell<=u,
    q=k-r>=0,   delta=2r-t>=0,
    D=2k+1-h=2q+delta+u.

Moreover `sum_{ab in E}(w_a w_b-1)<=q`. In particular every endpoint
weight is at most `q+1`.

At `D=3` there are exactly six rows:

| `(q,delta,u)` | Defect graph, besides isolated matching edges | Endpoint weights |
|---|---|---|
| `(0,0,3)` | matching with `k` edges | all one |
| `(0,1,2)` | one `P3` | all one |
| `(0,2,1)` | one `P4`, one `K1,3`, or two `P3` | all one |
| `(0,3,0)` | a tree on five vertices; a triangle; `(P4 or K1,3)+P3`; or three `P3` | all one |
| `(1,0,1)` | matching with `k-1` edges | all one, or exactly one endpoint has weight two |
| `(1,1,0)` | one `P3`, with `k-1` total edges | all one, or one degree-one endpoint has weight two |

There are three unlabelled trees on five vertices: `P5`, `K1,4`, and the
tree with degree sequence `(3,2,1,1,1)`. Every exceptional endpoint gap
has weight `0,1,...,k`; all other gaps have weight `K`. Infeasible rows
at small `k` are simply absent. For `k=1`, the fifth row includes the
single light vertex of weight one, with zero defect endpoints. The empty
template is excluded because the empty tournament is counted separately.

**Proof of completeness.** Solve `2q+delta+u=3` in nonnegative integers.
When `q=0`, all endpoint weights are one. For a connected nonmatching
defect component on `v` vertices and `e` edges, `2e-v>=v-2`.
At `delta=1` only `P3` occurs; at `delta=2` the connected possibilities
are the two trees on four vertices. At `delta=3` a connected component
has odd `v<=5`: `v=3,e=3` gives a triangle, and `v=5,e=4` gives a tree.
Partitions of the remaining defect into components give the displayed
disjoint unions. With `q=1`, weighted excess is at most one, so only one
degree-one endpoint may be doubled. This proves both necessity and the
claimed exhaustive generation rule. Soundness follows from the displayed
order's weighted cost. A generated template may need normalization; that
issue is handled next. QED.

## 2. One rule covers every clone merger

Take any tournament represented by an order `0,...,b-1` and its set `E`
of backward pairs. Write `N_E(i)` for the undirected defect neighborhood
of `i`, regardless of vertex weights.

**Exact merger theorem.** For `i<j`, the pair `{i,j}` is homogeneous iff

    (N_E(i) symmetric_difference N_E(j)) minus {i,j}
          = {i+1,...,j-1}.                                      (1)

The homogeneous pair is directed `i -> j` iff `ij` is not in `E`.
The resulting directed homogeneous-pair paths are precisely the canonical
transitive chains. Contract each whole path, sum its weights, and cap at
`K`; the result is the canonical weighted state.

**Proof.** Outside the interval `[i,j]`, the base order gives `i` and `j`
the same comparison to any third vertex, so equality after flips requires
equal defect incidences. Inside the interval their base comparisons are
opposite, so equality requires exactly one defect incidence. These are
exactly (1). The orientation of the pair is immediate. The inherited
canonical homogeneous-pair-chain theorem now supplies the contraction;
weights affect only the sum on each contracted chain. QED.

In particular a homogeneous pair need not be adjacent in the displayed
order. On three ordered vertices with defects `{02,12}`, the canonical
chain is `(2,0,1)`. The rule also permits an endpoint to merge with an
unmarked neighbor across intervening endpoints. Conversely, an unmarked
vertex strictly inside `(i,j)` obstructs the merger, because it belongs
to neither defect neighborhood. This is an exact external-port test,
not a claim based on the induced local tournament alone.

The checker compares (1) with the separate adjacency-based chain algorithm
on all 1099 labelled tournaments of orders one through five, and again on
every generated template. Both paths agree.

The smallest decreasing-deficiency witness has `k=1`, cap `K=2`, endpoint
weights `(1,1)`, the one defect `{01}`, and gaps `(0,0,0)`. Its raw heavy
count is zero, so its source deficiency is three. Its two vertices form
one chain of weight two, giving one heavy core vertex and final deficiency
two. Any claim that normalization only increases deficiency already fails
on this two-vertex example.

## 3. Intrinsic heavy gaps bound all ownership ambiguities

**Local ownership theorem.** For any canonical capped state of budget `k`
and deficiency `D`, let `m_0,...,m_h` be the sizes of the light gaps in the
intrinsic heavy order. Then

    e := sum_j max(m_j-1,0) <= 2D,
    product_j m_j! <= (e+1)! <= (2D+1)!.                       (2)

Every good order, and every weighted isomorphism, preserves the heavy
order and the light-gap partition. Thus at most `(2D+1)!` candidate
orders suffice to select a canonical owner. At `D=3` this is at most
5040; it does not grow with `k`.

**Proof.** A backward edge meeting a heavy vertex would cost at least
`K>k`, so it cannot occur in a good order. Hence the heavy vertices induce
a transitive tournament, whose order is unique, and every light vertex
has one intrinsic cut in that order. Its comparisons to the heavies
determine its gap. Weight-preserving isomorphisms preserve those cuts.

Use any good endpoint presentation from section 1. There are `L=t+ell`
light vertices. If `h>=1`, every one of the `h-1` internal heavy gaps is
nonempty: two consecutive heavy vertices with no intervening light vertex
would have identical comparisons to all others and be homogeneous. If
`c` is the number of nonempty light gaps, then

    e=L-c <= t+ell-(h-1) = u+ell <= 2u <= 2D.

When `h=0`, there is one nonempty light gap and `u=t+1`, so directly
`e=t+ell-1<=2u-2<=2D`. This handles the small-budget boundary with no
heavy vertices. Repeatedly applying
`(a+1)!(b+1)! <= (a+b+1)!` for positive `a,b` proves the factorial bound
from `sum(m_j-1)_+=e`. QED.

The bound `e` is intrinsic and computable without choosing a good
presentation. The stronger presentation-dependent bound `e<=2u` is
sometimes useful, but no choice of minimum `u` is required. The total
number of vertices in non-singleton light gaps is at most `2e<=4D`.
Singleton light gaps and the arbitrarily long matching background retain
their exact positions and cross-gap incidences.

For a faster implementation, the backward cost between different light
gaps is independent of their internal orders. Compute that crossing cost
once. Then list the orders of each small gap by its internal weighted
cost, and combine only those whose cost sum fits the residual budget.
This is an exact pruning rule, not a probabilistic filter.

## 4. Canonical ownership without a global duplicate table

For each candidate order inside the intrinsic light gaps:

1. Interleave it with the fixed heavy order and compute its weighted
   backward cost. Discard it if the cost exceeds `k`.
2. Recover its defect endpoint set, ordered endpoint weights, endpoint-gap
   weights, and backward pairs. Canonical pair-freeness guarantees that
   every endpoint gap contains at most one unmarked vertex.
3. Form the tuple

       (q, delta, u, endpoint_weights, gap_weights, backward_pairs),

   with vertex positions numbered in the recovered endpoint order.

The **owner** is the lexicographically least such tuple. The precise tuple
ordering is a convention; its invariance and reconstructibility are what
matter. A weighted automorphism may give several vertex orders with the
same tuple. Those are the same presentation, not extra owners.

**Exact acceptance rule.** Generate each ordered menu record once. Reject
it if (1) finds any homogeneous pair. For every remaining record, compute
its owner. Accept iff the record is its own owner.

**Proof of exactly one acceptance per class.** An intrinsic weighted
isomorphism preserves the complete collection of candidate records, so
it preserves their minimum. Conversely, a record reconstructs all
vertices, weights and arcs; equal owners therefore imply weighted
isomorphism. A pair-free state of the target deficiency has a direct good
presentation in the exact menu, and its least good presentation is also
in that menu. The raw menu enumerates that numerical record exactly once.
Thus one and only one record is accepted for each class. QED.

This argument explains why rejecting all non-pair-free templates loses
nothing. Their normalized outputs have their own direct pair-free menu
presentations. One need not distribute a merged presentation among other
rows, guess its new deficiency, or retain a global `seen` set. The checker
does retain sets and uses the old general canonicalizer **as independent
controls**; the direct acceptance counts are separately checked against
them.

In particular the D2 ported-triangle collision and every D3 shared-endpoint
or changed-gap collision are handled by the same rule. This is a bounded
symbolic ownership grammar, not a manually enumerated all-k collision
count. It deliberately retains all external ports; replacing a gap by
only its induced unmarked isomorphism class would not justify the rule.

Here are two explicit D3 hostile controls at `k=2`.

* **Changing endpoint roles and defect count.** The record with endpoint
  weights `(1,1,1)`, gaps `(0,0,3,3)`, and defects `{01,02}` has signature
  `(q,delta,u)=(0,1,2)`. Its displayed weighted order has five vertices.
  Interchanging its first two vertices yields endpoint weights `(1,1)`,
  gaps `(1,3,3)`, and the single defect `{01}`, of signature `(1,0,1)`.
  Both are pair-free presentations of the same state. The first order
  spends two backward edges and the second spends one; the budget is an
  upper bound, so both are legitimate. An endpoint becomes an occupied gap.
* **Identical local gap types, different external ports.** Use four unit
  endpoints and endpoint gaps `(1,3,1,3,1)`. The matchings `{01,23}` and
  `{02,13}` give two pair-free states. Both have light-gap tournaments
  `TT2,TT3,TT2`, with identical weights, but the cross-gap arcs differ.
  Those transitive gap orders are rigid, so the states are not isomorphic.
  The complete owner keys distinguish them. Retaining only local gap
  isomorphism types would incorrectly identify them.

## 5. All-deficiency enumeration and the unrestricted sequence

The proof works for every `D`, not merely three. Its complete ordered menu
is obtained as follows:

- Choose nonnegative `q,delta,u` with `2q+delta+u=D`, set `r=k-q` and
  `t=2r-delta`, and generate every simple graph on the ordered endpoint
  set `{0,...,t-1}` with exactly `r` edges and no isolated vertex.
- Assign positive endpoint weights with
  `sum_{ab in E}(w_a w_b-1)<=q`. Each is bounded by `q+1`.
- Choose exactly `u` of the `t+1` endpoint gaps to have weight in
  `0,...,k`; every other gap has weight `K`. Omit empty reconstructed
  tournaments and infeasible parameter choices.
- Apply the pair-free filter and local-owner acceptance rule.

All endpoint vertices and all weighted gap positions are ordered, so a
specified numerical record is generated once. Define `O_(k,D)` to be
the accepted records, and let `mass(R)` be the sum of all capped weights
reconstructed from `R`. Then the exact ordinary generating function is

    A_k(x) = 1 + sum_{D=0}^{2k+1}
                   sum_{R in O_(k,D)}
                   x^mass(R) / (1-x)^(2k+1-D).                (3)

Each heavy coordinate can be inflated independently to any actual weight
at least `K`. Heavy vertices are individually fixed by every weighted
automorphism, because their intrinsic transitive order is fixed. Hence
their generating factors are independent and have no additional orbit
division. Light weights are exact. This proves (3).

For `n>=0`, every tournament on `n` vertices has feedback at most
`floor(n(n-1)/4)`: compare any order with its reversal. Consequently

    A000568(n) = [x^n] A_floor(n(n-1)/4)(x).                  (4)

This is an exact local-ownership enumeration of the unrestricted sequence,
not an efficient scalar recurrence. The ordered endpoint graphs and their
long matching pairings still grow with `k`; the factorial locality bound
depends on `D`, which itself can grow with `k`. The achieved reduction is
specific: quotient ownership needs only bounded local order search at each
fixed deficiency, rather than a global table of previously generated
isomorphism classes or an independently enumerated pair-free core bank.
The accepted templates still represent every class being counted; this is
not a claim that the total enumeration has become small.

## 6. Reproduction and finite evidence

Run from the repository root:

```text
python 05-knowledge/results/deficiency_three_20261020/ownership_check.py
python -O 05-knowledge/results/deficiency_three_20261020/ownership_check.py
```

The explicit universe is every menu template of source deficiency at most
three at budgets one through three. Templates satisfy the inherited
weighted cost filter and are normalized before the independent class-key
comparison; final deficiency is recomputed. Separately, the direct
pair-free owner test counts classes without deduplication.
Normal and optimized (`-O`) runs pass and produce byte-identical 3496-byte
output in [`ownership_check.out`](ownership_check.out). The measured
optimized replay took about 23 seconds on this workspace.

| Budget | `D=0` | `D=1` | `D=2` | `D=3` |
|---|---:|---:|---:|---:|
| 1 | 1 | 5 | 9 | 5 |
| 2 | 3 | 45 | 264 | 756 |
| 3 | 15 | 435 | 5365 | 36473 |

The k=1 and k=2 sets agree with the independent complete capped-template
compiler. Every local owner through k=3 is crosschecked bijectively against
the prior weighted degree-cell canonicalizer. The D3 states include light
gaps of size seven, so the `7!` bound is attained as a candidate search
space. The independent normalization and all-class comparisons test the
actual quotient, not just the number of raw presentations.

For k=3 the owner signatures have counts
`30113, 5823, 272, 1, 255, 9` in the table's six-row order.
These are finite evidence only; this note does not extrapolate their
dependence on k.
