# Middle-out simultaneous cut inversion for a marked physical seam

**Proved.** An exact best-first decoder finds the common apex tags that give a prescribed signed gap and land nearest a chosen center on either side. It searches the paired cut, not two independent lists of scalar counts. Its intermediate certificates are the smallest anchored affine integer lattice of the two reads and intersected directional ranges. This does not give a lift at every exponent.

## Physical input and the two reads

Let \(T_-,T_+\) on the same labelled old vertices differ at one marked arc, with counts \(h,h+2\). Give a new vertex the **same** set \(S\) of beaten old vertices on both sides. The preceding [cut-gauge packet](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md) gives
\[
G_T(S)=H(T+x_S)-H(T)
      =\sum_{u\to v,\ u\in S,\ v\notin S} W_T(u,v).
\]
Here \(W_T(u,v)\) is the sum of the two complementary ordered-path-cover counts with reversed seam endpoints. It is symmetric and even. Precompute the covers for every endpoint pair with a subset table \(E(A,u)\) of paths ending at \(u\) and its converse table \(B(\bar A,v)\) of paths starting at \(v\); sum \(E(A,u)B(\bar A,v)\) over \(A\), adding the swapped-endpoint version. This computes all weights in \(O(n^2 2^n)\) arithmetic operations, not one Hamilton DP per tag set.

Put
\[
a(S)=G_{T_-}(S)/2,\qquad d(S)=(G_{T_+}(S)-G_{T_-}(S))/2.
\]
The two new counts in old order are \(h+2a\) and \(h+2+2a+2d\). Thus \(d=0\) preserves the sign, and \(d=-2\) also yields a literal unordered unit seam, reversing the sign. The center at fixed \(d\) is \(h+1+2a+d\). Mixing the two signs before adjusting this shift gives an off-by-two target error.

For Boolean tags \(x_i=[i\in S]\), the pair \(F(x)=(a,d)\) is quadratic: an arc contributes a vector multiple of \(x_i(1-x_j)\). Fix any tags (the **anchors**). The residual variables \(U\) have a polynomial
\[
F(x)=c+\sum_{i\in U}v_i x_i+\sum_{i<j\in U}v_{ij}x_ix_j.
\]
Fixing tag \(i\) to zero removes it; fixing it to one adds \(v_i\) to \(c\) and \(v_{ij}\) to each corresponding remaining unary vector. The pair coefficients on the still-unfixed vertices do not change.

## Exact affine lattice at the anchor

Let \(\Gamma\) be the integer span of the current unary and quadratic vectors. Then \(c+\Gamma\) is the **smallest affine integer lattice containing all completions' two-read values**. Containment follows from the polynomial. Conversely, evaluate at the empty, singleton and pair supports of \(U\): their differences and second differences give the generators. This is stronger than separate divisibility in the two reads. For example, generators \((1,1),(1,-1)\) impose no marginal divisibility restriction, but their lattice requires the coordinates to agree modulo 2.

For prescribed final defect \(d_\ast\), set \(\delta=d_\ast-c_d\). Write each generator as \((b_j,e_j)\) and let \(g=\gcd_j e_j\ge0\). If \(g>0\), extended Euclid gives an actual lattice vector \((b_\ast,g)\). The slice is empty unless \(g\mid\delta\), and otherwise its first coordinates are exactly the progression
\[
c_a+(\delta/g)b_\ast\pmod{k},\qquad
k=\gcd_j\bigl(b_j-(e_j/g)b_\ast\bigr).
\]
Modulus 0 denotes a singleton. Indeed the vectors \((b_j,e_j)-(e_j/g)(b_\ast,g)\) span the zero-defect kernel, and subtracting \((\delta/g)(b_\ast,g)\) from any vector at defect \(\delta\) reduces to that kernel. If all \(e_j=0\), require \(\delta=0\) and use modulus \(\gcd_j b_j\). This is linear in the number of generators, aside from extended-gcd arithmetic. An integer in the progression need not come from a Boolean completion.

## Joint interval and best-first stopping

For any direction \((\lambda,\mu)\), take the anchored scalar polynomial \(\lambda a+\mu d\). Select disjoint pairs of unassigned vertices, for example by descending absolute interaction coefficient. For each pair, take the exact min and max of its two unary terms and internal interaction over its four assignments. Bound leftover unaries individually and each unused term \(w x_ix_j\) between \(0\) and \(w\). Adding gives a certified interval. Larger disjoint clusters work with the same proof.

Use \((0,1)\) to reject impossible defect. Every direction with \(\lambda=1\) bounds \(a\) after subtracting \(\mu d_\ast\). The implementation uses \(\mu=0,\pm1,\pm R,\pm2R\), where \(R\) is the ratio of total absolute sizes of the current first-read and defect generators, rounded down and at least 1. Both signs matter because there can be cancellations between the two unsigned growth terms.

Intersect the joint interval with the anchored progression and, optionally, with \(a\le K\) or \(a\ge K\). Prune empty intersections. Put each surviving node in a queue keyed by distance from \(K\) to its nearest admissible point. Branch on one further tag, preferably one with high leverage in both reads. The first complete node popped is an exact nearest witness: every unresolved queue key is a lower bound for every completion in that node, and each not-yet-popped feasible assignment belongs to some remaining node. Exhausting the queue proves infeasibility on the selected side. The same decoder can target two final counts from an intermediate pair not presently at gap two. No subexponential worst-case bound is claimed.

## One-apex physical evidence

Write \(E_p=2^{p-1}(2^p-1)\). Rows are outgoing bit masks; tags index the same labels. The tagged reads passed full mask comparison and fresh path DP of the actual tournaments. One old marked arc is preserved throughout each query.

| Step | Source search | Next center | Nearest at or below | Nearest at or above |
|---|---|---:|---:|---:|
| 3 to 4 | Complete six-vertex marked birth bank, 6 marks | 120 | 116 (sign changes) | none |
| 5 to 6 | Complete eight-vertex marked birth bank, 18 marks | 2016 | 1952 (sign retained) | 2034 (sign changes) |
| 7 to 8 | Named ten-vertex seam from prior packet | 32640 | 28984 (sign retained) | 36152 (sign retained) |

Direct witnesses:
- 116: rows (32,9,3,5,15,30), mark (1,5), tag 21. The old lower side becomes 117; the old upper side becomes 115.
- 1952: rows (132,65,74,163,15,23,57,118), mark (1,2), tag 134.
- 2034: rows (12,145,82,22,225,15,43,109), mark (2,5), tag 53. The old lower side becomes 2035; the other becomes 2033.
- 28984 and 36152: same rows (988,433,74,34,236,645,938,780,60,286), same mark (8,9), and tags 942 and 952 respectively.

In the named ten-vertex source, only 5 of 1024 tags preserve the sign, including the two uniform tags; the two exact nearest sign-preserving queries popped 61 and 74 nodes. Across all eighteen order-eight marks each sign-preserving query direction has 4608 raw assignments, versus 562 and 503 live queue pops. Bound evaluations at children and cut preprocessing are additional work. These small cases show pruning, not a large-instance runtime comparison.

## Two and three common apex vertices

For a small complete bank it is efficient to enumerate common tags on all but the final apex, then evaluate all final apex tags from the two cut polynomials. A subset recurrence adds the least set bit \(i\) to a tag \(S\): add its unary coefficient and the pair coefficients against all previous bits. Thus a single table of endpoint cover weights supplies every final Hamilton count in \(O(n2^n)\) arithmetic operations after the weight preprocessing. Intermediate gaps need not remain unit. Compare the two resulting physical counts in both orders and independently recount every retained extremal witness.

- All 18 order-eight marks, two extra vertices: 2,359,296 labelled tag sequences. Around 2016 the nearest attainable centers are **1960** and **2034**. The lower uses rows (132,65,74,163,15,23,57,118), mark (2,4), tags (1,1), retaining sign. The upper is the one-apex witness with one padding vertex, hence still exchanges signs.
- Both marked six-vertex orbits, three extra vertices: 4,194,304 labelled tag sequences. Around 120 the nearest centers are **118** and **128**, both retaining sign. Both use old rows (32,9,3,5,15,30) and mark (1,5), with tag sequences (0,109,239) and (0,65,231) respectively. The two marked orbits suffice because each eligible six-vertex lower class has one orbit of its three eligible marked arcs.

These searches exclude an exact next-center hit only in the stated scopes. One-sided sign checks alone would not exclude unordered lifts.

## No induced \(p=3\) to \(p=4\) lift at any order

This first step is globally excluded **for the operation of retaining the old seam by induced extension**, even when its sign changes. It does not exclude independently built \(119\leftrightarrow121\) pairs.

Use the [primitive unit-seam grammar](../chanpak_support_birth_grammar_20261122/PRIMITIVE-SEAMS.md): after common ordinal cuts, every one-arc gap-two pair of lower count at least 5 is a primitive both-strong pair with transitive singleton padding around it. Since the old marked arc is not separated by any common ordinal cut of the final pair, the final primitive block contains the entire old strong primitive block: a cut dividing that old block would restrict to an ordinal cut of a strong tournament. Thus any induced step gives a final strong primitive pair containing an old strong primitive pair with the same marked arc.

Busch's [strong minimum](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) is 125 at order ten and increasing thereafter, so the final block at \(119,121\) has at most nine vertices. The old \(27,29\) block has only six or seven vertices: the order-five maximum is 15, while the minimum at order eight is 45. Consequently the six-vertex source requires at most three common extensions, and the seven-vertex source at most two. A smaller final pair can be padded with common ordinal sinks to reach these scan sizes without changing either count. The order-six tests above check both signed gaps. The complete 14-mark order-seven bank was checked over 458,752 two-apex tag sequences, again in both signed gaps. Its nearest final centers about 120 are also 118 and 128. Neither source order permits center 120, proving the obstruction.

## Replay

The Python replay checks the one-apex decoder against all tags and recounts returned rows directly; the small birth-bank file supplies all marks from the order-six and order-eight class lists. The C++ final-apex cut replay checks larger extensions and recounts extremal witnesses by ordinary subset DP.

For the published bank convention, compile the C++ replay and feed it the row records:
    c++ -std=c++17 -O3 check_nearest.cpp -o /tmp/check_nearest
    /tmp/check_nearest 5 2 < midout-audit-deck.rows
    /tmp/check_nearest 3 3 < midout-audit-deck.rows
    /tmp/check_nearest 3 2 7 < midout-p3-seven.rows

The six-vertex deck has six marked rows, three relabellings for each of its two marked isomorphism types. The replay drops duplicates for its three-apex order-six scan. The standalone global replay in [the exact-seam note](EXACT-SEAMS.md) independently derives canonical source classes and does not need the deck files.
