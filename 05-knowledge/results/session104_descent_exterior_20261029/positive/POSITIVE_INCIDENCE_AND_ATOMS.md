# Positive evidence, incidence, and well-founded clocks

**PROVED elementary synthesis; independent audit passed.** The graph alternative below is classical. The contribution here is an explicit interface among the repository's target atoms, K4 edge coordinates, and the auxiliary DAG in Session 103. Numerical equality of dimensions is not an isomorphism claim.

## 1. Signed compression versus positive detection

Let L be a real d-by-N matrix, N finite and positive. The following are equivalent:

1. Lz≠0 for every nonzero z∈R^N with z≥0.
2. Zero is not in the convex hull of the N columns of L.
3. There is a row vector y with every entry of yL strictly positive.

The first two are equivalent by normalizing a nonzero nonnegative kernel vector to have coordinate sum one. For 2⇒3 choose the point p of least Euclidean norm in the compact convex hull. Since p≠0, minimality along the segment from p to any column c gives p·c≥||p||²>0. Take y=p^t. The converse follows by multiplying a putative nonnegative kernel vector by yL.

For rational L, strictness allows y to be chosen rational by density, then scaled so yL has positive integer entries. Thus a finite linear detector that never erases positive evidence can be replaced, for the **zero/nonzero question alone**, by a single positive scalar detector. This does not reconstruct z, preserve signs, preserve individual sources, or give a closed representation for later operations. A sum of all coordinates already illustrates the distinction.

Conversely, signed injectivity is much stronger: it requires ker L=0. Calling a positive zero-detector faithful without specifying its cone and predicate would conceal this difference.

## 2. Exact graph alternative

Orient each edge of a finite simple graph and let B have column e_v−e_u for an arc u→v. The following are equivalent:

- ker B intersects the nonnegative orthant only in zero;
- the oriented graph has no directed cycle;
- an integer vertex potential f satisfies f(v)−f(u)≥1 on every arc u→v.

A directed cycle gives a positive vector of unit edge flows with boundary zero. If a nonzero nonnegative circulation exists, choose an edge with positive flow and follow positive outgoing flow: conservation prevents stopping, and finiteness forces a repeated vertex and a directed cycle. A DAG has a topological ordering, furnishing f. Finally f^tB has entries at least one, so f^tBz≥Σ_e z_e>0 for every nonzero z≥0.

For K4, B has rank three and signed kernel dimension three for **every orientation**. Nevertheless it detects all positive edge evidence exactly in the transitive orientations. A directed triangle is the smallest hostile: its three positive charges disappear under B. Thus three coordinates can detect positivity without retaining the six edge coordinates, but cyclic cancellation destroys that property.

The gradient B^t f is the list of potential differences, represented in exterior notation by 1∧f after fixing edge signs. This is the genuine interface to wedge products; no cross-product convention or Hodge metric is needed for the graph alternative.

## 3. The Session 103 observer is an instance of this mechanism

For a tournament with marked Hamilton cycle successor s, form D_C with arc s(u)→v for every noncycle arc u→v. If D_C is a DAG, choose f as above. A directed path cover has excess

    Σ_used arcs [f(v)−f(s(u))]
       = Σ_ends f(s(t_i)) − Σ_starts f(a_i).

Cycle arcs have zero excess; all noncycle arcs have strictly positive excess. Thus the sum cannot silently cancel an off-cycle contribution. Equality detects precisely the canonical cycle-cut cover. Ordering rows and columns by their endpoint potential makes every ordered-k observer unit triangular.

This is inherited from [Session 103](../../session103_words_necklaces_20261029/triangularity/AVOIDING_WORD_TRIANGULARITY.md), not a new proof of its classification. The connection identifies *why* the positive grading works: an acyclic incidence carrier admits a separating positive functional. If D_C is cyclic, a positive circulation defeats this grading. Singularity of M2 then still needs the separate bridge–sink argument in the actual grammar; cyclicity alone does not prove singularity for arbitrary digraphs.

## 4. Finite DAGs versus infinite arithmetic dynamics

Reverse the potential direction to obtain a decreasing rank. For a deterministic transition F on a set X, with a terminal set A where execution stops, existence of r:X→N satisfying r(F(x))<r(x) for all x∉A proves every orbit reaches A. More precisely its hitting time is at most r(x). Conversely, if every orbit reaches A, its hitting time itself is such a rank. This converse is a characterization, not a method of proving Collatz: calculating an unproved hitting-time function assumes the conclusion.

Acyclicity alone is insufficient on an infinite graph. The transitions n→n+1 on positive integers have no directed cycle but never terminate. They have a strictly increasing integer potential, and a decreasing real potential 1/n, but no decreasing natural-number rank. Finite acyclic residue charts likewise do not suffice if transitions escape charts or reset the clock. A global well-founded rank needs those carries and resets included.

The graph distinction also clarifies finite algorithms. Testing N candidate LRC endpoints one by one has the decreasing search rank N,N−1,...,0, regardless of whether a safe endpoint exists. This proves the **test terminates**, not that its answer is positive. A proof of the mathematical target needs an additional invariant forbidding exhaustion without a witness.

## 5. Target atoms and sequential endpoint insertion

The inherited endpoint measure for a finite nonempty set V of distinct positive integer speeds and rational 0<δ<1/2 has tags (v,k,ε), 0≤k<v, ε=±1, located at t=(k+εδ)/v modulo one. Repeated locations retain their owner tags. Let Safe(V) mean ||vt||≥δ for every v∈V. Its accepted atomic count is positive exactly when Safe(V) is nonempty, including isolated boundary points.

Adding a new speed w has an exact birth/filter decomposition of accepted tags:

    old accepted tags surviving the w test
       disjoint union
    new w-owned endpoint tags passing every old speed test.

The union is disjoint **by tags**, not necessarily by locations. Hence an incremental implementation must retain the old accepted locations (or an equivalent future-test carrier), not only their count. All previously rejected old tags remain rejected because constraints only increase. The birth term is essential: candidate locations change when a runner is added. There is no monotone decreasing accepted-count rank in general.

For V={1,...,13}, δ=1/14, the inherited exact safe set is {j/14:gcd(j,14)=1}; its six locations have twelve accepted owner tags among 182 candidates, despite having zero Lebesgue measure. This endpoint certificate is finite and source-visible. It is not a uniform LRC(14) positivity proof.

The corresponding Collatz source measure is countable rather than finite: positive weight on every integer source converts zero limiting survival mass to pointwise termination, but supplies no decay estimate. In contrast to the finite matrix theorem, its atom weights have no positive uniform lower bound. The current descent lane examines quantitative, source-preserving conditions that would yield an actual integer rank.

## Scope and provenance

Classical ingredients: strict separation of a finite convex hull; directed-cycle/circulation/topological-order alternative; well-founded ranking of terminating deterministic transitions. Repository inputs: operation_closed_synthesis_20261029/THEORY.md §III (owner-tag endpoint and Elias source measures), and Session103 all-k potential theorem. New packaging here: the common positive-cone test, its signed/positive distinction, explicit K4 hostile census, and exact endpoint birth/filter implementation.

No claim is made that the six-state legality machine is Λ² of its four-state acceptance machine. The exterior lane supplies a rank obstruction to that stronger identification. No Collatz, LRC(14), or abc conjecture is proved by these interfaces.

## Independent audit

See [descent/positive-interface audit](../audit_descent/REPORT.md). Its independent checker includes all K4 orientations, endpoint cell and birth/filter controls, native integer-rank edges, exact source cylinders, and paid-reset comparisons.
