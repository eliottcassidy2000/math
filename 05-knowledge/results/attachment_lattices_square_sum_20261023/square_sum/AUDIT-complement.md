# Independent audit of the complement/observer packet

23 October 2026, square-sum lane. **PASS for the two all-size theorems, their
quantifiers, and the stated bounded census.** One minor merged-state indicator
clarification was requested and incorporated by the author; it does not change
either theorem or any computed value. The subsequent addressed score-vector
transition was also audited and passes.

Audited sources:
[complement README](../complement/README.md),
[complement checker](../complement/check.py), and its frozen receipt.
No complement file was edited by this auditor.

## Proof audit

1. **Adjacent-swap delta.** Write `e=G_uv`, `a=G_uw`, `b=G_vw`. If w is before
   the adjacent ordered pair `(u,v)`, the triple is cyclic exactly when
   `a=e` and `a!=b`. If w is after the pair, it is cyclic exactly when
   `a!=e` and `a!=b`. Therefore placing every unequal-incidence w on the
   stipulated side makes the triple initially transitive and cyclic after
   the swap. Equal-incidence w never contributes; triples not containing both
   u,v do not change. The delta is exactly `d_G(u,v)`, with no cancellation
   or assumption on the rest of the order.

2. **All-order single-class boundary.** Order-independent isomorphism class,
   even after merging converses, implies constant directed-triangle count.
   The delta theorem then makes every pair twins. On each triple this forces
   all three graph bits equal; overlapping triples propagate that equality
   to all pairs. For n=2 every graph is already complete or empty. Conversely
   complete/empty graphs yield transitive tournaments for every observer.
   The four-way equivalence is sound for every n>=2.

3. **Twin compression.** A twin-class permutation is an actual graph
   automorphism, including permutations of twins occupying nonadjacent slots
   in the order. Equivariance supplies the tournament isomorphism. The word
   count `n!/prod m_i!` is consequently a safe order-state upper bound, not a
   census of distinct tournament types. Fixing marked ports is necessary.

4. **Constant-transfer theorem.** For any G and edge coordinate uv, choose an
   order with u,v adjacent and fix `T=T(G,P)`. Swapping only those two order
   positions changes the backward graph of this same T to `G xor {uv}`.
   Thus the assumed unqualified transfer equates f at every adjacent pair
   of the graph Boolean cube. Cube connectivity proves constancy, including
   after the complement quotient. No assumption that a canonical observer
   exists is hidden in the proof.

5. **Incidence balance.** Pairs `(graph,total order)` and `(tournament,total
   order)` are equivariantly bijective. Automorphisms act freely on total
   orders, so division by the respective automorphism orders gives integer
   orbit counts and exactly the displayed two-sided incidence identity.
   This is not a bijection of unmarked graph/tournament classes.

6. **Triangle expectation.** A mixed graph triple has exactly two vertices
   incident with one edge and one nonedge, and is cyclic for exactly two of
   six orders. Thus the expectation is
   `sum_v d_v(n-1-d_v)/6`. Empty/complete induced triples never contribute.

## Minor merged-type precision

The raw definition of `mu_C` uses converse-merged tournament types. An
indicator on that space is directly defined for a converse-invariant
property. For an arbitrary tournament property, first lift the orbit mass
equally to its two types (or average their indicators). This lift is canonical:
reversing P pairs every occurrence of T with an occurrence of its converse.
Thus the averaged-statistic claim is recoverable in full, but the word
"indicator" should not silently identify two different property values.
This clarification was sent to the author before publication and is now
incorporated explicitly in the README. There is no outstanding correction.

## Subsequent vertex-insertion formulas

The added formula `O=(L intersect N) union (R minus N)` follows directly from
reversing graph edges relative to the observer order. A new directed triangle
is exactly `w->u->v->w`, with `u in O`, `v outside O`. Consequently

```
Delta c_3 = e_T(O,V minus O)
          = sum_(u in O) d_T^+(u) - binom(|O|,2).
```

Every internal pair of O contributes exactly one outarc to the sum, which
justifies the subtraction. Old vertices outside O gain one outarc to w;
vertices inside O do not. The new vertex has outdegree `|O|`. These are
complete exact updates for the addressed score observable, not reconstruction
of all tournament arcs. The isolated/universal cases specialize to crossing
edges/nonedges of G, respectively, and the Laplacian cut identity is correct.
The stated two-owner hostile correctly produces increments1 and0 with the
same unindexed score multiset and set cardinality.

The added `complement/cut_response.py` was rerun and passed389,180 checks,
including every graph, observer, cut, and neighborhood through four old
vertices. This audit found no mathematical error in the extension.

## Subsequent incidence-flow synthesis

[INCIDENCE_FLOW.md](../complement/INCIDENCE_FLOW.md) also passes independent
algebraic review. The two-step transition is
`D_r^(-1) B D_c^(-1) B^T`, with stationary weights `r/sum r`. Multiplication
by `D_r^(1/2)` on the left and its inverse on the right gives the stated
positive semidefinite Gram matrix. The conductance/conditional-variance
identity follows by expanding the pairwise square; the factor1/2 and total
mass normalization are correct.

Connectedness of the unmerged incidence carrier also follows directly from
the audited constant-transfer construction: every graph and its single-edge
toggle share one intermediate tournament type under two adjacent observers.
Thus the two-step graph chain connects all graph types. Positive self-loops
give aperiodicity. The conclusions of a unique stationary law and convergence
are per finite n and do not assert a uniform rate.

The exact mixture relation
`mu_star=(1/4)mu_empty+(3/4)mu_matching` is correct coordinatewise. Since rows
are probability vectors, distinct rows cannot have a two-row linear
dependence. The explicit three-row relation is therefore support-minimal in
the stated order-four census. More generally a left-null vector sums tozero,
so its normalized positive and negative parts give genuinely different input
mixtures with the same observation. Individual-row separation does not
contradict this mixture failure.

`complement/check_flow.py` was rerun in normal and optimized modes; both
passed363 checks and reproduced the same exact rational receipt. No error or
outstanding correction was found in this synthesis.

## Computational audit

The original checker was rerun unmodified through order7 and passed
**864,044 checks** including the frozen-receipt comparison. Its universe is
accurately described: every labelled mask is assigned its full symmetric-group
orbit; all5040 orders of each of1044 graph isotypes at7 are inspected. The
converse quotient has272 tournament types, and the522 graph-complement pairs
have522 distinct merged distributions and522 distinct triangle histograms.
These are finite statements; neither reconstruction target was promoted to
an all-size theorem.

For an independent implementation,
[audit_complement.py](audit_complement.py) uses ordinary adjacency matrices,
direct orientation comparisons, and direct directed-triangle tests, rather
than the original XOR observer implementation. It canonizes every graph and
tournament through n=5 under all vertex permutations and recomputes:

- the adjacent-swap formula on every pair of each graph isotype;
- the empty/complete independence boundary;
- both directions of the incidence count and its integral divisors;
- the triangle expectation and histogram separation of all18 complement
  pairs at5.

Normal and optimized Python each passed **997 checks**, with identical output.
The independent implementation's horizon is5, not7. Order7 was reproduced
using the original exact algorithm and independently audited structurally.

Reproduction:

```
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/audit_complement.py
python -B -O 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/audit_complement.py
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/check.py
```
