# Session 304: peelable triangles, mixed Collatz, and a missing ternary digit

[Arithmetic proof and conditional digit tree](ARITHMETIC.md). Mixed signed odd-part Collatz has a sharp three-step decreasing macro avoiding 7 and 21. This is an exact arithmetic tree, not a Chan–Pak construction or a result for either single-sign Collatz map. Triangle peeling cannot simulate the tree: it never changes the 3-free part. A different conditional ternary compiler needs exact 3H−2, 3H, 3H+2 operations with compatible renewable marks.

## Physical negative digit, in two infinite families

A word describes a transitive spine with a single apex that beats precisely its 1 positions. In both lines below insert two letters without moving any old vertex; all old arcs are preserved and the same apex is retained. For k≥1:

- 100(1^k)00 → 1000(1^k)010 gives H=15·2^k+3 → 45·2^k+7.
- 101(0^k)10 → 1001(0^k)110 gives H=15·2^k+1 → 45·2^k+1.

Each output is 3H−2, and all these tournaments are strong because the apex beats the first spine vertex and loses to the last. These are not yet repeatable ready transitions.

Proof: all cycles meet the apex. The odd-cycle packing identity therefore gives H−1 as weighted inversions of the word, with weight W(1)=2 and W(d)=2^(d−1) for d≥2. Indeed, choose the first and last spine vertices of a cycle at distance d; the interior subset must be even, and each cycle contributes 2. In streaming form the current binary prefix integer is X and its last bit is b (initially both 0). Appending 0 adds X+b to H; appending 1 adds zero. Updating X←2X+e and b←e proves the displayed formulas, including the adjacency correction. For a run of k≥1 ones, X becomes (X+1)2^k−1; for a run of zeros X becomes 2^k X and H increases by (2^k−1)X+b.

## Why the two insertions are essential, and why they cannot simply repeat

Fix an active spine (first bit 1 and last bit 0). A single insertion in or immediately adjacent to this active spine cannot realize 3H−2. By reverse-complement, take the new bit to be 1. Every old inversion gains at most a factor two. Map each new inversion to the old inversion formed by the first old 1 and the same zero; it has no larger weight if the new 1 is after the first old 1. Equality in the resulting bound 3(H−1) would require the mapped inversions to exhaust all old inversions, hence only one old 1 (every old 1 has an inversion with the final zero). It would also require all old inversions to double, so the cut is between the first 1 and its immediately following 0. That adjacent inversion has W(1)=W(2) and does not double. Inserting immediately before the active spine leaves all old weights fixed and each new weight at most twice its mapped old weight, with strict inequality for any old inversion not hit, or for a mapped old adjacent inversion. If neither occurs, the word has one 1 and starts 10, which does have an adjacent inversion. Insertion immediately after the spine of a 1 has no effect. Inactive padding farther away is not covered by this local no-go.

For any active span L≥3, its endpoints alone contribute 2^(L−2) to H−1, while summing weights over all pairs gives H−1<2^L. Hence any infinite one-apex trajectory with exact negative digits must satisfy L_j=j log₂3+O(1), by H_j−1=3^j(H_0−1). A strategy that adds two active positions every time is impossible. More generally an eventually periodic active-width increment sequence is impossible since its average is rational, whereas log₂3 is irrational. If widths grow by one or two, two-increment stages must have bounded discrepancy from frequency log₂3−1, the same quantitative criterion that motivates rotation and Christoffel-word comparisons. This does not give any physical itinerary. In particular an infinite rule confined to one or two insertions inside the current active spine cannot succeed: the one-insertion option is ruled out, while two insertions force width to grow too quickly.

A different impossibility concerns cross-port edits. Keeping the internal arcs of induced blocks A,B fixed, any orientation between them has H≥H(A)H(B): all internal odd-cycle packings still contribute, and mixed packings add nonnegative weight. Thus starting from the product baseline T→C3, no number of cross-port edits alone yields 3H−2. The physical charts above instead change the internal spine/apex geometry. Equivalently a genuine negative digit triples the nonempty packing weight H−1 without tripling the empty contribution.

## Reproduction and status

The companion check independently enumerates literal directed paths on small physical charts, checks the exact arithmetic macro over a bounded interval and verifies its sharp example by direct search. None of these finite checks supplies the all-order proof; the arguments above do. Chan–Pak, the one-sign Collatz problems, and a common renewable ternary-ready state remain open.
