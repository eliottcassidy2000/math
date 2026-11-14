# Literal 2ᵃ3ᵇ endpoint fibres, and what prevents arbitrary count assembly

**Proved:** positive smooth fibre expansion; a marked trijection fixing each fibre; an all-integer marked endpoint decoder; and an explicit unmarked single-edge gap constructor. None identifies an arbitrary odd target with one unmarked Hamilton count.

## Integer sums, and the missing graph step

Every positive N is a sum of 2^a3^b with a strictly increasing and b strictly decreasing. Remove powers of two until the remaining number is odd, take its largest pure 3-power, and subtract it. The remainder is even and less than twice that power, so the residual after removing at least one factor two is smaller than the just-chosen 3-power. Induction gives the two strict orders, with no division by a graph count. Examples are 34=18+16, 45=27+18, 69=27+18+24. Consecutive ternary exponents, as in the inherited Collatz carry words of THM-4640, are a stronger and false demand for arbitrary N. For 17, a maximal exponent b=2 with all three exponents present costs at least 9+6+4=19; b=1 would require 17−3 to be a power of two, and b≥3 exceeds 17. The 17=9+8 representation skips b=1.

## The positive smooth fibre theorem

Delete k marked vertices to a transitive chain. In a full Hamiltonian path, fix the permutation of these marks and the endpoints of each of the k+1 intervening increasing chain runs, allowing each run to be empty or a singleton. Different nonempty runs must have disjoint endpoint sets. Check the finitely many edges between the runs and the marks; when an internal run is empty, check the corresponding direct mark-to-mark arc. These contacts no longer depend on the interior assignment.

Each unanchored chain vertex may then be assigned independently to any run whose open interval contains it. Within a run its internal order is forced. Thus the fibre size is the product of the numbers of open runs available at each unanchored vertex; if any has zero, the fibre is empty. Distinct endpoint data and mark orders give disjoint fibres and count every Hamiltonian path once. With k=2 this is exactly 2^a3^b, where a counts vertices in exactly two open runs and b in all three. With k=3 the same holds with a counting exactly-two vertices plus twice the exactly-four vertices, and b exactly-three. Five-way membership is possible for k=4, so the two-prime alphabet does not persist just by allowing more apices: choose five distinct starts and five distinct ends in front of and behind one common central vertex, and set the required apex contact edges. That fibre has cardinality five.

For any a,b≥0 every monomial 2^a3^b is physically available as one literal fibre: choose six distinct run endpoints in the order

    a₀, a₁, [a double-choice vertices], a₂, [b triple-choice vertices], b₀, b₁, b₂,

with apex order x,y and runs (aᵢ,bᵢ). Choose the four necessary arcs b₀→x→a₁ and b₁→y→a₂. Choosing also x→a₀ and y→a₁ both gives the closers below and makes it easy to keep the tournament strong: the first apex reaches the first chain vertex and the final chain vertex can beat both. Each double/triple membership is literal.

This can be counted without scanning endpoints. For each permutation of the k marks, read the transitive chain in order. The state of each run is empty, last vertex admissible at its exit, or last vertex inadmissible. When appending a vertex to an empty run, check its entering edge from the preceding mark if that mark exists. When appending to a nonempty run, its order is legal automatically. Update the last-vertex flag and sum over choices. At the end forbid inadmissible last flags and check the direct mark arc for each empty internal run. This finite positive transfer counts H in O(k!(k+1)3^(k+1)n) additions. THM-4224 gave two- and three-run terms for a particular one-defect family; this argument uses arbitrary apex incidence, arbitrary endpoint fibre and standard Hamilton paths.

## A controlled disconnected carrier at exactly n−1 edges

For a three-run fibre R₀→x→R₁→y→R₂ with all runs nonempty, suppose the first two have distinct endpoints and extra arcs x→a₀ and y→a₁. Replace the departures x→a₁ and y→a₂ by these two closers. The carrier is now (R₀+x cycle) ⊔ (R₁+y cycle) ⊔ R₂ path, with n vertices, n−1 edges and three components. For just the first closer one gets one cycle and a remaining path. The marked apex cycles and marked remaining path recover the original word uniquely by cutting the closers and restoring the departures. Thus this is an actual weight-preserving bijection on the specified endpoint fibre. It does not extend to all graphs with n−1 edges; the extra arcs are genuine gates. In any graph with n−1 edges, the total cycle rank is its component count minus one, so disconnected carriers are not automatically forests. This is the restricted constructive bridge to the repository's line/cycle-plus-solo/star trijection, not a conflation of their quotient spaces.

## Why the fibres do not freely add as ordinary H

Even the obvious common-apex attempt injects extra cycles. Glue k one-apex graphs at the shared apex v, with the outside blocks ordered and all interblock arcs forward. If block i has p_i out-neighbors of v and block j has q_j in-neighbors, every i<j gives p_i q_j distinct mixed directed triangles through v. Every cycle uses v, so odd-cycle packing gives

    H_glue ≥ 1+Σ_i(H_i−1)+2Σ_(i<j) p_i q_j.

Longer mixed odd cycles may add further mass. An exact positive-sum compiler must account for these genuinely extra fibres, not just write an antichain expansion of the target.

## Arbitrary exact marked endings, then exact unmarked single-edge gaps

Given 0≤N<2^r, take a transitive chain y,c₀,c₁,…,c_r,t and an apex x. Make x beat y,c₀, make t beat x, and let c_i beat x precisely at a 1 bit of N (weight 2^(i−1)). The graph is strong since t→x→y→…→t. A full path *from c₀ to t* must use x,y consecutively because y only has x as an in-neighbor. If c_i immediately precedes x, earlier vertices c₁,…,c_(i−1) each choose before or after the fixed x,y pair, and later vertices must follow it. This gives exactly 2^(i−1) paths at each selected bit, no others. Thus the prescribed marked count is N in O(log(N+1)) vertices. It is not the ordinary H.

For a scalar gap, insert a spacer z between c₀ and c₁, set x→z, and compare the cases x→c₀ and c₀→x (the latter is paired, not substituted into the preceding marked count proof). Both graphs remain strong. All odd cycles use x. In a transitive carrier whose endpoints i<j give x→i and j→x, the odd-cycle packing contribution to H is 2^d when there are d≥1 interior vertices and 2 when d=0: choose an even subset of the d vertices. Reversing only x→c₀ removes contributions 2^i from each selected c_i and 2^(r+1) from t, and adds just the adjacent pair y,c₀ of weight 2. Therefore the exact positive gap is

    H(x→c₀)−H(c₀→x) = 2N+2^(r+1)−2.

As N varies through [0,2^r), these adjacent step-two intervals tile all nonnegative even numbers as r varies (r≥1 covers all positive evens). This is an explicit single-arc constructor with r+5 vertices, not a claim of priority over other all-even-gap constructions. It still cannot anchor an endpoint at will: from the compulsory outgoing vertices y,z to the selected in-neighbors and t, the lower endpoint already satisfies H_lower≥(5/2) gap +8. The gap mechanism is absolute, but choosing the absolute endpoint remains a separate condition.

The companion replay checks the endpoint transfer and literal graphs, including actual ternary fibres and the three-component marking. Numerical examples of 69 and 45 are typed in the separate rail note.
