# A quadratic-time exact cut test for labelled blue feasibility

Take a proposed integer score vector $s$ at the labelled Hamiltonian backbone $0\to1\to\cdots\to n-1$. Require the mirror symmetry $s_{n-1-i}=n-1-s_i$. Then it has a blue presentation **if and only if**
$$
\sum_{i\in S}s_i\ \ge\ \binom{|S|}{2}
 +|\{i<n-1:i\in S,\ i+1\notin S\}|
\qquad(S\subseteq\{0,\ldots,n-1\}).
$$
The total is already $\binom n2$ by symmetry. The extra term counts prescribed backbone wins leaving the subset. These inequalities give ordinary feasibility with the backbone fixed: assign each remaining unordered pair to its winner with vertex capacities equal to the required residual wins. The bipartite matching capacity tests reduce exactly to the displayed cuts (count the free edges entirely in $S$); equality of total demand and total free edges makes all vertex capacities exact.

## Why imposing the mirror does not add an obstruction

Average an ordinary feasible orientation and its mirrored converse. This gives a fractional mirror-compatible orientation with the same scores. Use the [blue half-score variables](OBSTRUCTIONS.md) for its free mirror orbits, all in $[0,1]$. Their columns are signed incidences (one $\pm1$ at each of two distinct coordinates) or unary columns (one $-1$). Choose an extreme point of the feasible box slice. Its strictly fractional columns must be linearly independent, or a kernel direction allows small perturbations in both senses.

Draw the fractional two-coordinate columns as edges of a multigraph; retain unary columns at their vertices. In a connected fractional component of $v$ vertices there are at most $v$ fractional columns by independence and at least $v-1$ pair edges by connectedness. If it is a tree, it has at most one unary. A leaf without a unary then has only one fractional term in an integral score equation, impossible. An isolated fractional unary is impossible for the same reason. The remaining possibility is a unicyclic component with no fractional unary. There can be no tree hanging off its cycle, again by the leaf argument. Parallel edges with different signs may give a two-edge cycle.

At every cycle vertex, integrality says the two fractional edges must agree or sum to one, depending on the incidence signs. On a dependent cycle the relation returns to the original value; independence forces a return to its complement. Thus all its edge variables are $1/2$. Every vertex except the innermost vertex at even order has a unary orbit in the original box, so the cycle (of length at least two) has a vertex with a unary pinned at 0 or 1. Round the cycle along either of its two alternating choices, cancelling at every other vertex. At the chosen vertex the change is $+1$ or $-1$, with either sign possible by reversing all rounding choices. Choose the sign cancelled by toggling its pinned unary. The cycles are disjoint, so the roundings do not interfere. This constructs an integral blue presentation.

## Test without listing cuts

Scan the backbone. Store, for each subset size and last membership bit, the minimum cost among prefixes of that type. Selecting a vertex adds its score; a selected-to-unselected transition subtracts one. At the end the minimum for each cardinality must be at least its binomial bound. There are $O(n^2)$ states and constant work per transition in integer arithmetic. The [exact replay](replay_cut_feasibility.py) matches this test with Laurent support through order nine.

This is a support theorem for **labelled** blue scores, including prescribed first score. Its partner's scores follow by reversing the path and shifting the two endpoints. It does not solve which ordering exists for an unlabelled one-box score pair, and does not identify the unquotiented edge endpoints; those still use the [alignment recurrence](ALIGNMENT_RECURRENCE.md).
