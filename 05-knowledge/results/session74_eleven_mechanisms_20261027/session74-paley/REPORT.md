# Paley11 is exactly six native path complements from a forest frame

**Status:** PROVED (cycle-rank lower bound and the exact six-move distance, combining the proof with the independently checked explicit certificate); FINITE-EXACT (the six named transitions, their feedback optima, and bounded hostile checks). No literature novelty claim. No LRC, Collatz source-payment, or arithmetic identification is inferred from the shared integer 11.

## 1. Inheritance and what changes

The input is THM-4810, `01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md`, its S63 source `05-knowledge/results/h_spectrum_atlas_complement_budgets_paley_clocks_20261027.md`, and the independent audit beside that source, pinned at commit `ca5a9f7bee48a08bad9d1f471ece85f8b085c345`. The mathematical baseline worktree was `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`. Exact file hashes are in `DEPENDENCIES.json`.

Inherited, not new here:

- For sorted outdegrees `s_(0)<=...<=s_(n-1)`, `L(T)=sum_i |s_(i)-i|/2` is an integer, `L<=f`, and a legal native path complement changes L by at most one.
- A native move `tau_P` reverses every arc OFF the named directed Hamilton path P, preserving all arcs of P. Repeating the same named move is the identity.
- Paley11 has L=15, exact feedback number f=20 and intrinsic backward-graph cycle rank beta=10. S63 gave distance at least five to a forest frame and at least fifteen to a transitive tournament.
- S63 already established exact minimum f=14,13,12 at depths one, two and three. Our witness does not improve those feedback minima.

New deduction here: a forest-framed tournament has **L<=n-2**, not just the inherited L<=n-1. More generally every tournament satisfies

    beta(T) >= max(0, L(T)-n+2),                  n>=2.

This sharpens the Paley11 forest lower bound from five to six. An explicit six-move legal word attains it. Thus the open endpoint is settled:

    d_native(Paley11, forest-framed class) = 6.

The existing score/forest results were checked for overlap: `04-computation/experiments/lattice_bridges_20261021/score_vectors_forests.py` and `05-knowledge/results/session56_rmq_lca_connections_20261025/score-forest-counts/{REPORT,COROLLARIES}.md` count labelled score fibres via graphical-zonotope forests and describe fixed-path complement pairings. Their forests are counting objects, not the actual backward graph whose cycle rank is minimized here. They do not state the present bound. This is a targeted repository comparison, not a claim of global priority.

Concept board: **anchor**—native Paley11 forest distance; **niche**—actual shared-cycle rank rather than feedback-edge count; **wildcard**—regular staircase family exposing sharpness and information loss. The useful connection is an inequality on these same labelled tournament objects, not a numerical analogy to a different occurrence of 11.

## 2. Definitions and exact scope

Fix n named vertices, n>=2. A rank order defines a transitive baseline. Its **backward graph** D is the undirected graph on ALL n vertices of the pairs on which T disagrees with that baseline. Write m=|E(D)| and c=c(D), counting isolated vertices as components. Its cycle rank is

    b(D)=m-n+c.

Define `f(T)=min_order m` and `beta(T)=min_order b(D)` separately. Minimizing f and minimizing beta need not select the same order. A tournament is forest-framed iff beta=0.

For the proof, index the baseline so vertex i has baseline outdegree i (higher indices beat lower indices). For certificate orders, we use the familiar opposite listing convention: **earlier vertices beat later vertices** in the baseline. Reversing the whole list converts between conventions; L, m, c and beta do not depend on which convention describes the same frame.

Only actual native complements are charged. Legal path remarking, relabelling, and full converse may be interspersed freely: each preserves L and beta. The order stays n, and no substitution, copying, vertex insertion, port erasure, or free reconstruction is supplied. The named paths and all vertex labels are retained in the certificate. No source weights or probability masses are invented or changed.

## 3. Connected-frame strictness

**Lemma (PROVED).** If the backward graph D is connected and n>=2, then `L(T)<=m-1`.

For each i, let p_i be the number of backward edges from i to a larger index and q_i the number to a smaller index. The actual outdegree is

    s_i=i+p_i-q_i,
    sum_i p_i=sum_i q_i=m.

The assignment taking vertex i to baseline score i gives

    L(T) <= (1/2)sum_i |p_i-q_i| <= m.

The quantity on the middle is an integer. If it is less than m, it is at most m-1, proving the claim. Otherwise equality forces `min(p_i,q_i)=0` for every vertex. Connectivity and n>=2 imply each vertex has an incident edge, so each displacement `p_i-q_i` is strictly positive or strictly negative. At i=0 it is positive; at i=n-1 it is negative. There is consequently some adjacent pair i,i+1 whose displacements have signs +,-. Thus

    s_i >= i+1,    s_(i+1) <= i.

In the score-to-target assignment, exchange just the assigned targets i and i+1. Its total absolute cost decreases by exactly two. Since sorted matching minimizes absolute assignment cost, `2L<=2m-2`. This proves the lemma. Connectivity is used only to rule out isolated zero-displacement vertices in the equality case; the graph need not be acyclic.

## 4. Cycle-rank bound and native clock

**Theorem (PROVED).** For every n>=2 tournament and every rank order,

    L(T) <= n-2+b(D).

If c>=2, the inherited elementary inequality L<=m gives

    L <= m=n-c+b(D) <= n-2+b(D).

If c=1, the lemma gives

    L <= m-1=n-2+b(D).

Minimizing ONLY b(D) over orders gives the claimed intrinsic bound

    beta(T) >= max(0,L(T)-n+2).

In particular, beta=0 implies L<=n-2. Notice the statement is about L; the stronger-looking assertion `f<=n-2` for every forest frame is FALSE. Our final eleven-vertex forest witness has f=10=n-1, while L=9=n-2.

Let a word contain k actual native complements, with only the free operations listed in Section 2 between them. The inherited score Lipschitz bound gives `L(T_k)>=L(T_0)-k`. Consequently

    beta(T_k) >= max(0,L(T_0)-k-n+2),
    d_native(T,{beta<=b}) >= max(0,L(T)-n+2-b).

This is a lower bound for the cycle rank of EVERY final frame. It does **not** assert that beta itself is 1-Lipschitz. Indeed the first move below sends Paley11 from inherited beta=10 to exactly beta=5 in one move.

## 5. Six explicit moves, all labels retained

Start with vertices 0,...,10 and `x->y` iff `(y-x) mod 11` is in `{1,3,4,5,9}`. The following are directed Hamilton paths in their respective CURRENT tournaments, not six paths asserted to be simultaneously legal in the initial tournament:

| Move | Named path P_k |
|---:|---|
| 1 | 0,1,2,3,4,5,6,7,8,9,10 |
| 2 | 1,3,9,0,10,5,2,8,4,6,7 |
| 3 | 0,3,2,1,6,5,9,8,4,7,10 |
| 4 | 1,6,2,4,7,8,3,5,9,10,0 |
| 5 | 0,4,3,5,8,6,9,10,2,7,1 |
| 6 | 1,8,7,3,10,6,5,0,4,2,9 |

After each move, this certificate records an actual rank frame; m=f is verified by an independent exact subset dynamic program:

| k | L | exact f | frame components c | frame b | intrinsic beta |
|---:|---:|---:|---:|---:|---:|
| 0 | 15 | 20 | inherited | inherited | 10, inherited |
| 1 | 14 | 14 | 2 | 5 | 5 |
| 2 | 13 | 14 | 1 | 4 | 4 |
| 3 | 12 | 12 | 2 | 3 | 3 |
| 4 | 11 | 12 | 1 | 2 | 2 |
| 5 | 10 | 11 | 1 | 1 | 1 |
| 6 | 9 | 10 | 1 | 0 | 0 |

The intrinsic beta column follows by matching the explicit frame upper bound with the theorem's lower bound. No enumeration over all ranks is required for this equality. In particular the minimum beta over ALL words of exactly k native complements from P11 equals `6-k` for k=1,...,6. The case k=0 is different: beta=10, not 6. The depth-two f value in this particular word is 14; it is not the already-known depth-two global minimum 13.

The final forward rank order is

    1,10,4,6,8,7,5,2,3,9,0.

Its actual ten backward pairs and connected component are listed in `normal.json` and independently regenerated by `check.py`. The final graph is a connected tree. All intermediate row bitmasks, forward ranks, paths, and search counts are in `certificate.json`.

**Exact-distance corollaries.** Six moves suffice by this witness, and at least six are needed since L(P11)=15 while all forest frames have L<=9. Every prefix has length k and drops L by exactly k, so each prefix is a geodesic from the initial labelled P11. Each suffix is also shortest from its own stage to the forest-framed class. More generally, for integer b=0,...,5,

    d_native(P11,{beta<=b})=6-b.

The same conclusions hold if free converse and relabelling are allowed, because the lower-bound invariant survives them and the given word does not need them. This does not settle the distance from P11 to the transitive class: the inherited lower bound fifteen remains, and no matching transitive witness is asserted.

## 6. Why eleven is useful, and where the bound is sharp

Every regular tournament on odd n=2r+1 vertices has `L=r(r+1)/2`, so

    beta(T) >= max(0,(r-1)(r-2)/2),
    d_native(T,forest frames) >= max(0,(r-1)(r-2)/2).

For Paley clocks q=4t+3 this is `t(2t-1)=(q-3)(q-5)/8`: at q=3,7,11,19,23,31 the lower budgets are 0,1,6,28,45,91. The mechanism is regularity plus the native score update, not a unique numerical power of Paley11. Eleven is the first prime Paley order after seven; at eleven the budget jumps from the known one-move P7 forest bridge to an exactly six-move P11 bridge.

The cycle-rank lower bound for regular tournaments is sharp for every odd n>=3, although not for P11 itself. Set n=2r+1, take the baseline j->i for i<j, and reverse exactly the staircase pairs

    0<=i<=r-1,    i+r+1<=j<=2r.

Every vertex then has score r. There are `r(r+1)/2` reversed edges. Vertex r is isolated in the backward graph; every other vertex belongs to one connected component (through the extreme left and right vertices). Thus c=2 and

    b = r(r+1)/2-(2r+1)+2 = (r-1)(r-2)/2.

This attains the lower bound for intrinsic beta. It does not prove that every regular tournament attains the distance lower bound. At n=11 the staircase has beta=6, while inherited beta(P11)=10: the degree profile loses four units of Paley-specific cycle debt.

The forest score ceiling is sharp at n=11 even without using the six-move search: reverse the edges of the undirected path

    4--7--0--8--1--9--2--10--3--5--6

in the ascending-score baseline. Its scores are `[2,3,4,5,5,5,5,5,6,7,8]`, so L=9=n-2. This explicitly prevents any universally stronger order-eleven ceiling based solely on L.

**Hostile controls.** The staircase at n=7 has m=6=n-1 but contains a cycle, and L=6>n-2=5. Therefore 'at most n-1 backward edges' cannot replace the forest predicate. The first P11 move drops beta by five, so replacing the score Lipschitz theorem by an alleged cycle-rank Lipschitz theorem is also invalid. Isolated vertices are counted throughout; omitting them breaks the staircase equality calculation.

## 7. Construction, verification, and reproduction

`greedy.cpp` is a witness finder, not a breadth-first search. At stage k it enumerates directed Hamilton paths of the one current literal tournament in lexicographic order, retains score-decreasing candidates, computes an exact optimal feedback rank, and stops at the first frame of cycle rank 6-k. The numbers of paths inspected at the six stages are

    1,12369,48,8357,44,10123   (total 30942),

and 2662 candidates receive subset-DP evaluation. The search is heuristic in its selection of a single branch, but its successful output is a complete exact witness; its counts are not claimed as class censuses or distance-search completeness.

`check.py` is independent of the C++ implementation. It reconstructs P11 directly from quadratic residues; validates every named path at its own stage; implements complement by an independent unordered-pair loop; checks every literal output row, same-path involution and endpoint score update; reconstructs backward graphs and their connected components; and verifies f by a first-vertex DP, whereas the finder uses a last-vertex DP. It also checks:

- EVERY ordered backward graph, equivalently every labelled tournament relative to the fixed baseline, for n=2,...,6 (2,8,64,1024,32768 graphs), with no inherited filters;
- EVERY legal Hamilton-path pair on EVERY labelled tournament at n=2,...,5;
- the connected-frame strict bound and general cycle bound on the entire graph universe;
- sorted absolute score transport against its independent cumulative-count formula, plus brute-force assignment minima through n=4;
- regular staircase equality for r=1,...,15, the independent eleven-vertex sharp forest, and the edge-count/forest hostile.

The all-size theorems follow from Sections 3-4, not from these finite checks. All checker guards use explicit exceptions, not removable assertions. Both normal and optimized Python runs produce byte-identical PASS receipts.

    g++ -std=c++17 -O2 greedy.cpp -o /tmp/session74-paley-greedy
    /tmp/session74-paley-greedy > /tmp/session74-paley-greedy.out
    cmp greedy.out /tmp/session74-paley-greedy.out
    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json
    sha256sum -c SHA256SUMS

An independent agent separately checked the proof, all six transitions, feedback DPs and final tree without importing producer code, and enumerated all ordered backward forests through n=9. Its report is separately supplied to the integrating session. This packet changes neither the inherited source files nor the historical claims: the old lower bound five remains a correct weaker bound, now sharpened to the exact value six.
