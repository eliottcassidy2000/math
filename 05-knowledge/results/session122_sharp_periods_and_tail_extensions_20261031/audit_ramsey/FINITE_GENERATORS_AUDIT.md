# Independent audit of finite capacity generators

2026-10-31. **PASS**, scoped to the nonempty singleton/join/union/full-C5
class D. This is a separate automated mathematical review, not a formal
proof-assistant build or external human review.

Audited source: `../ramsey/FINITE_GENERATORS.md`.

## Conclusions

For every fixed positive clique budget A, an explicit recursively computable
M_A makes the exact capacity row the unbounded-knapsack closure of its first
M_A values. A finite residue graph then computes all eventual penalties,
a certified tail threshold, and the finite exceptional prefix. This replaces
the previous abstract Presburger-effectivity step with a concrete finite
algorithm. It does not establish efficient complexity or unrestricted Ramsey
bounds. The general arbitrary-A compiler has not been implemented here.

## Proof audit

**A global deficit bound is available before computing the row.** The sharp
L-attainer and r copies of K_A yield F(A,tL+r)>=w(A)tL+Ar. Hence every optimum
has deficit at most D=(w(A)-A)(L-1). This is a uniform bound with no empirical
stabilization premise.

**The optimum has exact independence budget.** Any purported optimum with
independence below B could gain a vertex by disjoint union with a singleton.
This stays in D and within the budget. Consequently every capacity attainer
has actual independence equal to its specified positive argument.

**Connected components remain in D.** Flattening top-level unions leaves
singletons, joins, or full nonempty C5 roots. Each is connected. The
nonnegative component deficits sum to the whole optimum's deficit, so each
is at most D. No assertion about arbitrary induced subgraphs is needed.

**Enlarging clique budgets does not invent an optimum.** At a C5 component,
actual clique budgets can be enlarged to a maximal feasible integer pattern.
Replacing each child by an F(a_i,b_i)-attainer preserves its exact
independence number and the parent's clique upper bound. If this increased
order, it would improve the whole maximizing graph, a contradiction.
Thus the old child orders equal the values of the enlarged-budget capacities.
This step is what licenses use of the previously computed child rows.

**Bounded connected-root exceptions.** Joins and nonsharp C5 patterns have
strictly smaller slopes than w(A), by the audited strict-superadditivity and
full polytope-vertex arguments. Their actual independence budgets are at most
floor(D/gap). The floor, rather than the strict-violation ceiling-minus-one
in the rows5..8 theorem, is correct: a genuine optimum may have deficit D.

**Sharp duals are positive at every A.** The half vertex is the unique
maximizer, and its five active constraints are independent. Varying just one
edge slack slightly keeps all ports positive. A zero or negative dual
coefficient would give another maximizer or an improvement. Thus all five
dual coefficients are positive. Equivalently each coefficient is the sharp
half-objective minus a two-coordinate indicator objective, already known to
be strictly positive.

**The child tail threshold is retained.** If S=floor(D/min dual) and
N=max N_(a_i), port slacks total at most S. The odd-cycle inverse gives
2b_i=b-t_i with |t_i|<=S. Reducing b modulo2L into
[2N+S,2N+S+2L-1] leaves every b'_i at least N and preserves integrality.
Each child shift is a multiple of L, which is divisible by its own
L(a_i). The certified child tail is therefore valid in both directions of
this reduction. One cannot omit N or replace it by an observed plateau.

**Actual graph realization survives reduction.** Capacity attainers at the
reduced child arguments have exact independence b'_i. Since at least one
original port slack was zero, the reduced C5 root has independence exactly
b'. Its order decreases by exactly w(A)(b-b'). Disjoint union with
(b-b')/L sharp L-attainers restores both order and independence, with clique
number at most A. The reduced root need not have the original exact clique
number; the upper budget is the required invariant.

**Finite generation follows in both directions.** Every component of an
optimal graph is either bounded already or replaced by pieces of independence
at most M_A. Replacing each such piece by the corresponding F(A,b)-attainer
cannot improve the optimum. Conversely every sum of prefix attainers is an
actual legal disjoint union. Thus the equality with finite knapsack closure
holds for every B, including the exceptional prefix.

**The recursion is not circular.** Computing M_A uses only child rows with
clique budget below A. The finite F(A,b) prefix is computed in increasing b:
union children decrease b, while join and C5 children decrease the clique
budget. A=1 is the explicit initial row. Hence no uncomputed current-A tail
threshold appears in the input.

**Residue paths supply an explicit threshold.** Edge costs are integer
scaled deficits and nonnegative; secondary weights are positive lengths.
A lexicographic shortest path has no repeated vertex, since deleting a cycle
either lowers cost or preserves cost and lowers length. It has at most L-1
edges and weight at most (L-1)M_A. Its recorded item list is an actual union
construction. The zero-cost L-item extends that realization to every larger
weight in the same residue. Conversely finite generation shows no graph can
have smaller deficit. This proves the stated threshold and allows a finite
scan for the least threshold. The all-A least period remains exactly L by
the separately audited sharpness theorem.

## Independent exact controls

`check_generators.py` imports only this auditor's `check.py`. It recomputes
all bounds for certified rows2..8 using a digit-sum weight formula, exact
polytope-vertex reconstruction, and independently solved sharp duals.
It uses L-1 rounds of lexicographic Bellman-Ford relaxation, unlike the
producer's heap-based Dijkstra implementation, and then verifies every edge
satisfies the optimality inequalities.

The recomputed M-values are

    A=2,3,4,5,6,7,8: 6,15,14,55,74,49,35.

All match the producer. Sufficient shortest-path thresholds are
1,3,3,7,7,7,7, and the finite check recovers the already-proved least
threshold1 in each row. All path weights and profits are checked as actual
item lists. The independent coin-DP replays a total of1715 positive prefix
values, through at least the coarse threshold plus one period for rows5..8.
These controls instantiate the theorem; they do not replace its all-A proof.

For late-generator hostiles M=3,11,101, items(1,1),(2,4),(M,2M) have slope2,
nominal modulus2, and first zero-deficit odd realization exactly at M.
The exact knapsack row is 2B-1 for odd B<M and2B otherwise. Its least tail
threshold is M-1. This rejects the inference that an early plateau locates
the final tail or that shortest residue cost alone gives its first arrival.
These are abstract knapsack hostiles, not claimed members of D.

Reproduce:

    python3 -B check_generators.py > generators.result.json
    python3 -B -O check_generators.py > generators.result.optimized.json
    cmp generators.result.json generators.result.optimized.json

Both runs passed with byte-identical outputs. The completed rows1..8 remain
inherited/certified inputs to these finite controls; no general arbitrary-A
capacity compiler or numerical row9 result is claimed.
