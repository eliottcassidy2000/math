# A one-flip family and its actual blue metagraph neighbor

This provides an infinite family where the geometry can be described exactly, then asks which members are genuine first appearances. The key warning from the toy examples is that three notions differ: small feedback arc set, close to the transitive class in the fixed-order flip graph, and early birth for a numerical count.

## 1. One arc from transitive: class, count, and deletion

Use the ordinary **one-arc graph** on isomorphism classes at fixed order: an edge means some labelled representatives differ in exactly one orientation. Its distance to the transitive class is exactly the minimum number of backward arcs over all linear orders of the vertices (the feedback arc number). One direction is to flip every backward arc in an optimal order. The other follows by lifting any class walk to labelled representatives successively, then comparing with its transitive endpoint. An order minimizing backward arcs cannot have a backward adjacent arc, since swapping that adjacent pair would remove one and not change the others. Thus the minimum can be taken over Hamiltonian frames.

At distance one, every class is

    TT_a → E_c → TT_b,       a+b+c=n, c≥3,

where E_c is transitive in the order 0,…,c−1 except for c−1→0. Indeed reversing one arc of a transitive tournament is either still transitive (adjacent reversal), or gives precisely this interval of cyclic activity. The ends outside the interval are singleton SCCs. Thus the classes at distance one are indexed uniquely by (a,c,b), amounting to (n−1)(n−2)/2. Converse swaps a and b.

The unique Hamiltonian path not using the back arc is the increasing path. One that uses the back arc chooses any subset of the interior to appear before c−1; the remainder comes after 0, each segment increasing. This is the familiar H(E_c)=2^(c−2)+1 and it is unaffected by padding. Remove an endpoint of the exceptional arc and the graph becomes transitive. Remove a core interior vertex and get E_(c−1) with the same pads (transitive for c=3). Outside deletions preserve the count. This gives an exact, tiny recursive family of order-to-order deletion profiles, not a proof of minimal order for its H-values.

For example the literal c=6 graph first realizes 17 at order six. But E_7 reads 33 and E_8 reads 65; the atlas has already realized them at orders six and seven. In particular an induced chain of simple strong graphs with strictly growing H need not be a chain of first births. These are claims about numerical birth, not about which class is easiest to build.

## 2. What the same family does in the actual merged tiling metagraph

Here switch graphs: a [tiling metagraph line](../../01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md) complements **all nonconsecutive chord bits** in one marked Hamiltonian frame, not one arc. Blue means the chord word is fixed by frame reversal.

E_c is self-converse by the index reflection j↦c−1−j. It is rigid for c≥4: all scores are distinct except the two extreme pairs; within either tied pair a swap is impossible because it reverses their common arc. It has exactly one spanning cycle (every cycle uses the unique back arc, and a spanning one must take all remaining vertices in increasing order). Thus its class has H=2^(c−2)+1 tilings, with c closed and H−c open.

The increasing frame is blue. Any other frame uses c−1→0 and, to be blue, must use it in the center, where the frame equals its reverse reflected in the index order. That is possible only for c even. In that case each mirror pair of interior vertices chooses which half it occupies independently. Hence the blue incidence count for c≥4 is

    B(E_c) = 1 + [c even] 2^((c−2)/2).

Exactly 1+[c even] of these are closed, because closed frames are rotations of the unique spanning cycle. For the closed blue frames the increasing frame always works, and the rotation with the exceptional arc in the middle works just when c is even. For c=3 divide the frames into the single automorphism orbit: there is one tiling, which is blue and closed.

Padding on the left and right replaces closure by openness whenever a+b>0. The node is self-converse precisely when a=b. If a=b>0, all core blue frame patterns survive with the same number of blue tilings; if a≠b the two different converse classes merge and the resulting node has twice the tilings of one class, all black. In either padded case the node has no self-line, by the established closed/open law. This distinguishes symmetry of a core, symmetry of its padded placement, and actual blue incidence.

## 3. The exact path-count sequence on the other end of a blue line

In the increasing frame of E_n its only backward chord is the extreme one. Complementing all chords gives Q_n: every nonconsecutive chord is backward **except** the extreme chord, which goes forward. The marked pair E_n–Q_n is a literal blue line, closed at its E_n end and open at its Q_n end. This is not the unique blue line from the transitive node. The latter goes to the full backward-chord tournament F_n, whose counts h_n are the [inherited tribonacci law](../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md), h_0=h_1=h_2=1, h_3=3.

**Neighbor formula.** For n≥3,

    H(Q_n) = h_n − (n−1) + [z^(n−2)] (1−z)/[(1+z)(1−2z−z²−2z³)].

Proof: F_n and Q_n differ by one exceptional arc between the two extreme vertices. Count the affected paths by removing this adjacent pair, leaving the interior F_m, m=n−2, split into two ordered optional paths P,Q. An induced subset of F_m splits by its runs of *consecutive indices*: between runs the later one uniformly beats the earlier one, so the path counts are products of h for their run lengths.

For the old direction (top→bottom), P must either be empty or end at the largest interior index, and Q must either be empty or start at the smallest. The consecutive-run law forces P to be a suffix and Q a prefix. Both prescribed endpoint paths are unique (begin at smallest forces the whole increasing path; apply converse for ending at largest). Thus exactly m+1 paths use the old arc, including the two empty-side splits.

For the new direction (bottom→top), exclude P ending at the smallest or Q starting at the largest. Let A(z)=Σ_(r≥1)h_r z^r and E(z)=Σ_(r≥1)e_r z^r where e_r is the path count ending at the smallest of F_r. The inherited endpoint recurrence is e_r=h_(r−1)−e_(r−1) for r≥2 and e_1=1, hence E=z(1+A)/(1+z). Coloring the interior with colors P,Q gives alternating runs; total split GF is (1+A)/(1−A). The first excluded event replaces the first run of a word starting with P by E, giving E/(1−A); the second gives the same count on the last Q run. Their intersection starts P and ends Q, giving E²/(1−A²). Therefore the new-arc GF is

    (1+A)/(1−A) − 2E/(1−A) + E²/(1−A²)
      = (1+A)/((1−A)(1+z)²)
      = (1−z)/[(1+z)(1−2z−z²−2z³)].

This counts the actual partner and derives its recurrence from the decomposition rather than fitting numerical values.

| n | H(E_n) | H(F_n), transitive's partner | H(Q_n), E_n's marked partner |
|---:|---:|---:|---:|
| 3 | 3 | 3 | 1 |
| 4 | 5 | 5 | 5 |
| 5 | 9 | 9 | 11 |
| 6 | 17 | 17 | 29 |
| 7 | 33 | 31 | 69 |
| 8 | 65 | 57 | 169 |

At n=4 the E_n–Q_n line is a blue self-line: E_4 rows (6,12,8,1) map to Q_4 rows (10,4,9,2) by label images (2,0,3,1). At n=5 and 6 the marked blue line joins two different count fibers born at that same size. At n=7 it goes from an old value (33) to a new value (69), and at n=8 from an old value (65) to a new value (169). Thus even one fixed infinite line construction changes its birth role. A marked blue line is not in general a chronological link between first creators.

This is an explicit case where the best way to understand a family is to retain all three layers: the one-arc edit radius, the order-changing deletion skeleton, and the complement of a fixed Hamiltonian frame. `python verify.py` checks each count, mirror frame and closure with literal graph or path DP through order 16; the recurrences above prove the claims at every order.
