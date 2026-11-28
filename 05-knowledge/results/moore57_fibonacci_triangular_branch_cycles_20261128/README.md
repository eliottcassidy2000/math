# Branch-ring transfer, Fibonacci–triangular schedules, and a 63-cycle lower bound

2026-11-28. This route is separate from the [21-to-63 construction by detouring edges of a physical 21-cycle](../moore_57_heptagon_63_detours_20261127): here the fixed seed is a cyclic order of *branches*, and there need not be any 21-cycle over that order. Results about a full degree-57 Moore graph are conditional. The replay also works on the [actual shared-edge chart](../moore57_shared_edge_boards_20261128), which is not a Moore graph.

## All 22 coefficients, not just the two endpoints

Root at a site r. Its d branches consist of its neighbors (the hubs) and their q=d−1 leaves. In a diameter-two Moore graph every pair of branches is joined leaf-to-leaf by a perfect matching. The next argument requires only the complete matchings along the specified ring of m distinct branches, so applies equally to a partial chart with those matchings installed.

Write the permutation matrices for crossing successive branch pairs as P_i and their product around the ring as P. Write f=tr(P). In each branch choose an outgoing sheet. Its matching edge determines the incoming sheet of the next branch. Within a branch, either the incoming and outgoing sheets coincide (one leaf, no internal path), or connect the distinct leaves by their two-edge path via the hub. The branches are distinct, so every choice gives a simple cycle. Conversely the cycle identifies the outgoing sheets in the fixed orientation; counting with a single representative of the cyclic order already counts physical unoriented cycles once. Extra chords do not matter.

Using a variable t to mark the two-edge hub detours, put B(t)=(1−t)I+tJ. All permutation matrices commute with B(t), so
$$
G_m(t)=\operatorname{tr}(B(t)^mP)
=(1+(q-1)t)^m+(f-1)(1-t)^m .
$$
The decomposition is into the all-ones line and the centered subspace; the trace of P on the latter is f−1. For *each specified subset* of k detour positions, even with separate contact variables, the number of physical cycles is
$$
w_k=(q-1)^k+(-1)^k(f-1).
$$
Indeed the multivariate generating function is the product of (1+(q−1)t_i) plus (f−1) times the product of (1−t_i). Therefore the coefficient of t^k in G_m is the count of this family at length m+2k; also G_m(1)=q^m. For odd m, the two extremes add to (q−1)^m+1. For even m the difference N_{3m}−N_m equals (q−1)^m−1 instead. In particular at d=57 and m=21,
$$
N_{21}=f,\qquad N_{63}=55^{21}+1-f.
$$
One detour already recovers the same hidden count: N_{23}=21(56−f). These identities give no existence claim for a degree-57 Moore graph.

## Fibonacci–triangular numbers arise from a restricted edit grammar

This uses the repo's [triangular–Fibonacci/boat experiment](../../../04-computation/triangular_fibonacci_heron_boat_kps.py), which juxtaposes 21=3·7=T_6 and 55=5·11=T_10, and the prior [Fibonacci state-clock work](../../../04-computation/sturmian_fibonacci_triangular_kakeya_sixty_clock_thm4035.py). The bridge here is a local rule on the actual branch-ring edit positions.

Require that no two chosen hub detours be cyclically adjacent. Let I_m(z) be the independent-set polynomial of the cycle on the positions. The same physical transfer then gives
$$
G_m^{\rm sep}(t)=I_m((q-1)t)+(f-1)I_m(-t),
\qquad
I_m(z)=\operatorname{tr}\left(\begin{pmatrix}1&z\\1&0\end{pmatrix}^{m}\right).
$$
The two states record whether the preceding position was edited; the zero entry forbids two consecutive edits. With one marked position, the schedules omitting it are an independent-set problem on a path of m−1 sites (F_{m+1} schedules at z=1); those using it leave m−3 sites (F_{m−1}). Thus the projection to schedules has the golden-ratio transfer, while actual physical counts carry an extra ground-state weight and the signed monodromy term. This is not a claim that the cycle count is Fibonacci.

There are m(m−3)/2 nonadjacent two-position schedules. Allow exactly those and the empty schedule: the count is
$$
1+\frac{m(m-3)}2=\binom{m-1}2=T_{m-2}.
$$
In this *formal schedule count*, [Luo Ming's classification](https://www.fq.math.ca/Scanned/27-2/ming.pdf) says that the Fibonacci hits occur at exactly m=3,4,8,12, with counts 1,3,21,55. The empty schedule can have zero physical realizations if f=0. The corresponding physical tally is f+(T_{m-2}−1)((q−1)^2+f−1), not T_{m-2}. At m=21 there are 190 schedules of these two kinds. That is also the number of branch-pair strips installed in the previous session; we have no map identifying the two sets.

## One open path gives a simultaneous insertion constraint

Fix an ordered path through m−1 distinct branches from c to a, using only cross-matchings. For each starting sheet x in c let y be its endpoint in a. The pair (y,x) has exactly one return of length at most two in a Moore graph. It is either the direct matching edge a→c or a two-edge matching route through one of the other d−2 branches. The root and hubs do not supply such routes between different leaves' branches. Summing over starting sheets gives a q-slot ledger. In particular,
$$
f(\text{shortcut ring})+
\sum_{b\text{ fresh}} f(\text{ring inserting }b) \le q .
$$
The sum is an exact equality if one also includes the two-edge returns through all internal branches of the path, counted as closed walks. A repeated branch can either use a second leaf, giving a simple cycle, or revisit the identical leaf, giving a nonsimple walk. On a partial chart the inequality still holds for every known return, because two distinct returns give a triangle or square. It is a simultaneous constraint, not an independent cap on each ring.

For m=21, each fixed 20-branch path has 37 fresh-branch candidates in a Moore graph. In particular a candidate with f=54 leaves at most two sheets for all the other candidates and the shortcut together (and f=55 is impossible for a permutation of 56 sheets), and at most two candidates can have f≥21. Summing the inequality over ordered paths gives an average f≤56/37 over ordered 21-branch rings (write the inserted branch last). Every unoriented cyclic order is counted exactly 42 times in that average, so the bound transfers without any sampling hypothesis.

## A larger conditional 63-cycle bound without symmetry assumptions

Fix a root. Every resulting all-detour 63-cycle reveals which 21 vertices are hubs: they are exactly the vertices adjacent to that root. Reading these hubs in order gives the unique cyclic branch order up to reversal. Thus there is no repeated count at a given root. Summing over the (57)_21/42 unoriented branch orders and using the average above counts at least
$$
\frac{(57)_{21}}{42}\left(55^{21}+1-\frac{56}{37}\right)
$$
distinct physical 63-cycles at that root. Here (57)_21 is the falling factorial. Globally one cycle can be counted from at most three roots: the hubs must be one of the three every-third-vertex sets on the cycle, and any two hubs determine at most one common neighbor. Multiplying the bound by 3250 and dividing by 3 gives **more than 9.91×10^72 simple 63-cycles if a degree-57 Moore graph exists**. This improves the earlier heptagon-seeded lower bound, but does not prove that the graph exists or that any of the counted cycles is induced.

## Replayed transport on the current partial graph

The shared-edge certificate has 355 complete branch pairs. A seeded sample of 10,000 distinct 21-branch orders found f ranging from 0 to 6, including Fibonacci–triangular values 1 and 3. This is not an exclusion of 21 or 55. Those sampled orders alone certify more than 3.52×10^40 distinct simple 63-cycles in the *partial* graph, by the exact transfer and the unique root projection, without enumerating individual cycles.

One sampled 20-branch path had complete transport from both endpoint branches to all the others. Its 56-sheet ledger split into 46 fresh-branch returns, 7 simple returns repeating a branch but not a vertex, 1 repeated-vertex return, and 2 direct shortcuts. Thus the chart already supplies an observed gap of eight between the fresh-only insertion inequality and the actual full ledger.

Run python3 replay.py for brute-force transfer controls. To reproduce the chart sample, replay the preceding [shared-edge witness](../moore57_shared_edge_boards_20261128/README.md) to generate /tmp/moore-shared, then run python3 replay.py /tmp/moore-shared. The transcript is replay.out.
