# Compact arrow rules: five, ten and fifty

**Explicit and replayed.** These are smaller *ordinary tournaments from which the original graph can be read*, not smaller versions of the exact same observer Gram. Their full unlabelled classes are characterized by the arrow rules and the intrinsic decoder below. Blue means there is a Hamiltonian spine whose reversal is an anti-automorphism.

## Pentagonal field rules

For the pentagon itself use Z/5 with i→i+1,i+2. The skew row Gram equals 1 at cyclic distance one and −3 at distance two. Thus ordinary tournament structure recovers the pentagon on five vertices. Reflection makes it self-converse; the spine 0,1,2,3,4 is blue.

For Petersen and Hoffman–Singleton use q=1 and q=5 respectively. There are two halves U and W, each with q cohorts of five points, with coordinates U_(i,x), W_(j,y). Here i and j range from 0 to q−1, and x,y∈F5. The complete arrow rule is:

- Within a cohort use the cyclic five-tournament x→x+1,x+2.
- Between U cohorts orient in increasing cohort index, and between W cohorts in decreasing index.
- Cross arcs normally go U→W. Reverse just the affine contacts y=x+ij mod 5.

The graph decoder uses score groups, not external vertex names. At q=1, U has score 6 and W score 3. At q=5, U_i scores are 42−5i and W_j scores are 7+5j. Only scores 22 and 27 are ambiguous: an ambiguous W vertex beats all five members of the score-7 group, while an ambiguous U vertex beats four of them. This identifies the halves and then their cohorts. Read graph edges in each U cohort at Gram value 1 and each W cohort at value −3. Read cross graph edges from the exceptional W→U arcs. This is the familiar two pentagons and matching model of Petersen and the [Robertson five-pentagon/five-pentagram model](https://doi.org/10.1016/0095-8956(71)90015-3) of Hoffman–Singleton. Its parameters are checked directly in compact.py.

The involution U_(i,x)↔W_(i,−x) reverses all arrows. It reverses the cohort orders and their cyclic rules; the affine exception transforms to −x=ij−y, equivalent to the original. Read all U cohorts in sequence, taking consecutive cyclic increments +1 within a cohort. Choose the last residue x so that 2x≠−i² at the final cohort; append the reflected order in reverse. This is a Hamilton spine, with the central arc forward by that condition. These hosts are also strong: any U can go through W to any other U since their affine contact rows differ (at q=5 they agree on at most one slope, at q=1 on none). Every W is accessible from some U and reaches one. Thus both compressed plain tournaments are strong and blue on exactly the number of labels in the graph they recover.

In the displayed spine the backward tiles consist of three in each of the 2q cyclic cohorts, plus 5q² affine exceptions, for 6q+5q² total. Of these exactly q are axis-fixed (one solution of 2x=−i² per cohort); the rest are in pairs. So changing one reflection orbit at a time in the **edit cube** reaches these masks in exactly 6 moves for Petersen and 80 for Hoffman–Singleton. This is a blue edit-cube path, **not** a path in the native whole-complement metagraph.

## General unmarked compression by binary addresses

Let G have m Moore labels, and choose any total order v_0,...,v_(m−1). Orient the tournament among those m sites forward in that order except on graph edges, which go backward. To make the order recoverable from the *plain tournament*, add r=ceil(log2(m+1)) beacon sites (a transitive chain), and two markers a→b. All m sites beat both markers; both markers beat all beacons. Site i beats beacon j iff bit j of i is 1, with low bit first in the beacon chain.

For m≥5 the two markers are the unique twin pair (identical directions to all vertices except one another). Two original sites differ on some beacon. Two beacon columns differ on some original site (the powers of two separate the nonconstant columns, and at most the highest is constantly zero). Marker twins cannot spill into either layer: each original address has a zero, and every beacon has at least one zero in its column. The marker pair distinguishes original sites from beacons, and the beacon chain gives bit order. The unmarked tournament therefore recovers every address and then G by reading the backward arcs.

At degree 57 this is conditional and uses **3,264 vertices**, not 6,500. Only 12 beacons and two markers are overhead. For a marked-spine presentation without those overhead vertices choose a Hamilton path in the complement as the original order and read all backward tiles; this uses 3,250 vertices with a retained spine. Such a complement path exists for d≥3 by the minimum degree bound for a Hamilton cycle: the complement has degree m−1−d≥m/2. Equivalently the longest-path proof finds a cycle through its vertices by the intersection of the endpoint-neighbor index sets, then extends if there are vertices outside. Appending markers and beacons to that spine works as well.

The unmarked addressed construction is not blue: an anti-automorphism would preserve its unique twin pair while interchanging the layers of different sizes. The separate [partial observer method](moore_observer_compression.md) retains blue and unmarked decoding at at most **6,097** vertices for a degree-57 input. These are different contracts, not claims that the 3,264 construction satisfies the original blue Gram contract.

## If the exact binary Gram is required

The [overlap theorem](moore-gram-overlap.md) says at most seven physical sites can be both sources and probes, even allowing a different binary factorization of that Gram. Its induced-cycle construction works for both degree 7 and degree 57 (a greedy induced path of length five only excludes at most five neighbors at any step). All four graphs' Gram has full rank m, so the optimal physical totals, with logical self bits allowed, are **5, 13, 93, 6,493** for degrees 2, 3, 7, 57. The first two cases use the [small overlap certificates](gram-small.py). This is an optimum for the much stricter exact binary coefficient/Gram model, not a lower bound for the field or beacon decoder.
