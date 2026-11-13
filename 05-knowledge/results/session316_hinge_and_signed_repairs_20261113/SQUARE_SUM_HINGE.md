# Square-sum hinges: literal edges, composite reflections, and feasible ports

## The requested arrows have two different types

Write \(Q_N\) for the loopless graph on \(1,\ldots,N\) with \(x\sim y\) if \(x+y\) is square, and \(R_c(x)=c-x\). For a square \(c\), its graph on admissible pairs is a matching. This is the inherited [reflection-orbit calculus](../procgen_sumgraph_20260926_reflection_orbits.md).

The arrows \(4\to5\) and \(2\to7\) are single moves of \(R_9\). The arrow \(8\to6\) is not a square-sum edge: its sum is 14. It is the **composite** reflection
\[
R_9R_4R_9=R_{14},
\qquad 8\xrightarrow{9}1\xrightarrow{4}3\xrightarrow{9}6.
\]
At order eight these three motifs are not just selected examples: \(Q_8\) is exactly the disjoint union of this \(P_4\) and the two edges \(2-7,4-5\). Indeed distinct sums are at most 15, so only 4 and 9 can occur. The 9-reflection gives the full matching, and \(1-3\) joins two of its pairs. The distinction between an edge and a composite port is already visible in the first eight labels.

More generally \(R_bR_cR_a=R_{a+b-c}\). For any fixed allowed central edge \(u-v\) of sum \(c\), outer squares \(a,b\) attach the ports \(a-u\) and \(b-v\), provided all four vertices are admissible and distinct. One cannot insert this composite as a free edge in a spanning path: it consumes the two actual internal vertices \(u,v\). Different composite moves using the same hinge therefore cannot be used independently. Compare [the marked ear at the 15/16 seam](../seam_threshold_20260925.md), where suppressing a degree-two vertex similarly gives an edge outside the original arithmetic category.

## A complete arithmetic classifier for four-cycles on the hinge

For integers \(s,t\ge3\), the path
\[
(s^2-1)-1-3-(t^2-3)
\]
is induced except possibly for its endpoint edge. The cross chords have sums \(s^2+2\) and \(t^2-2\), neither square. The outer vertices are distinct since two squares cannot differ by 2 mod 4. Thus it is a \(C_4\) exactly when
\[
s^2+t^2-4=w^2.
\]
Conversely every four-cycle containing \(1-3\) is of this form, by reading its two outer edges from the hinge. For fixed \(s\), the allowed \(t\) are completely determined by positive factor pairs
\[
de=s^2-4,\quad e>d,\quad e\equiv d\pmod2,\qquad
t=(e-d)/2\ge3,\quad w=(e+d)/2.
\]
The graph-order cutoff is \(\max(s^2-1,t^2-3)\). This gives an exact enumerator for the cycles through the hinge, with no Hamilton search and no double-counting by reversing the cycle (the neighbors of 1 and 3 specify the roles of \(s,t\)). The factors of 5 and 12 show that outer squares 9 and 16 can never close such a four-cycle. For orders below 46 both \(s,t\le6\). For \(s=5\) only \(t=10\) qualifies and for \(s=6\) only \(t=7\) qualifies. Hence the first hinge four-cycle is
\[
35-1-3-46-35
\]
at order 46 with targets \(36,4,49,81\), agreeing with the [older global four-cycle threshold](../crossroads_family_20260926_squares.md) but obtained here by factoring the hinge relation.

For equal outer squares \(s=t\), a closure is possible exactly for \(s=2k\) with \(j^2-2k^2=-1\). Indeed \(w^2=2s^2-4\) forces both \(s,w\) even, then divide by four. The Pell solutions are \(j+k\sqrt2=(1+\sqrt2)^{2r+1}\), \(r\ge0\): a larger positive solution can be reduced by multiplying by \(3-2\sqrt2\), which stays positive and reduces \(k\) for \(k\ge3\), and \(k=2\) is not a solution. The first \(k=1\) degenerates to the hinge itself; nondegenerate closures begin at \(s=10,58,338,\ldots\). In particular the two 9 sides of the requested motif give a genuine induced \(P_4\), not a cyclic completion.

## Local support is not Hamilton support

The [September first-window calculation](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md) already gives the unique \(Q_{15}\) path and its extensions through \(Q_{17}\). In \(Q_{15}\), leaves 8 and 9 force the endpoints; the leaf 8 and the certified-interior degree-two vertex 15 saturate vertex 1 with edges \(1-8,1-15\). Thus \(1-3\), though the hinge at order eight, is the one unused edge of \(Q_{15}\). The direct 9-reflection edges \(4-5\) and \(2-7\) are both used. At order 23 the hinge reappears in all three Hamilton paths, but with *mixed* outer targets 9 and 25, not the 9-and-9 detour. It appears in seven of ten paths at 25 and none of twelve at 26; hence even the single-edge support is nonmonotone with order. Counts are up to reversal, as independently forced by the leaf 18 in this range.

For the three equal-side hinges, the first spanning path containing the prescribed four-vertex segment appears at these orders:

| Outer square | Prescribed segment up to reversal | First \(N\) |
|---:|---|---:|
| 9 | \(8,1,3,6\) | 32 |
| 16 | \(15,1,3,13\) | 27 |
| 25 | \(24,1,3,22\) | 25 |

The delay is a real marked-interface effect, not delayed edge birth. In the last row \(Q_{24}\) already contains the entire segment but no Hamilton path. The first use of the numerically smallest ports is the latest of the three.

An exact endpoint-aware certificate for the order-32 threshold works without Hamilton path branching. Keep the three requested edges, and delete every other edge incident with the saturated vertices 1 and 3. For \(8\le N\le20\) the result is disconnected. For \(N=21,\ldots,26\) the following sets \(S\) leave more than \(|S|+1\) components, impossible for a spanning path:

| \(N\) | \(S\) | Components after deletion |
|---:|---|---:|
|21|4|3|
|22|4|3|
|23|4,7,14|5|
|24|2,5,12|5|
|25|5,7,12,14|6|
|26|7,10,14,23|6|

At orders 27 through 31, branch only on the two endpoint labels, including forced leaf 18 where it exists. Assign degree demand one to those two and two elsewhere. Repeatedly reject if a demand is exceeded by forced edges or cannot be met by available edges; delete unforced edges at any vertex whose demand is filled, and force all available edges at any vertex with exactly its demand available. Also reject a forced cycle. At 27 this rejects 18 endpoint pairs by degree and the remaining 8 by cycle. At 28,29,30,31 respectively it rejects all 27,28,29,465 pairs by degree. This is a finite exact propagation certificate, not a claim that degree propagation decides the general problem. It deliberately avoids the [retracted unconditional degree-two rule](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md), which would wrongly treat a potential endpoint as interior.

At 32 a Hamilton witness, with the requested segment reversed, is
\[
2,7,18,31,5,4,32,17,19,30,6,3,1,8,28,21,15,10,26,23,13,12,24,25,11,14,22,27,9,16,20,29.
\]
The replay also finds the two other first orders by exact prescribed-segment traversal. In that traversal a segment is not an unmarked contracted vertex: when first entering it, the only permissible entry is one end and the state must leave at the other, marking all segment vertices used. This is a bijection with paths that contain that segment, and retains precisely the orientation and incidence information lost by an ordinary contraction. It is a template for the metagraph: paired roots and legal path order, not a free added edge, must survive any reduction.

## A caution from the adjacent arithmetic models

The [summand closure lane](../collatz_mod6_20260922_w6_summand_closure_square_filter.md) already separates square reflections from its additive closure \(\{2,3\}\) with three holes \(\{1,4,6\}\). The [triangular-ten seam](../seam_threshold_20260925.md) likewise identifies a true two-edge ear operation shared between two graph-growth problems, while distinguishing their first thresholds and edge labels. The present hinge is another literal example: a composition law and a rooted transfer can be exact even though neither connectivity nor Hamilton feasibility is determined by an unmarked scalar. No square-sum, Collatz, or Chan–Pak implication follows from matching the labels 4, 5, 8, 9 alone.



## A new marked 25-fold lift: the small labels themselves survive

The [earlier Gerbicz construction](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md#6-gerbiczs-25-fold-blow-up-s6-proved-by-a-finite-junction-check) takes any odd Hamilton chain \(a_1=1,\ldots,a_n=3\) and glues 25 shifted copies of it with singleton labels 1 through 12. Its usual glue does not prescribe the neighbors of 1 and 3. Replacing the glue by the following exact port word preserves **all three requested motifs on their original small labels**, at arbitrarily large orders:
\[
8-1-3-6,\qquad 2-7,\qquad 4-5 .
\]
The small edge \(1-3\) will be the closing edge of the lifted cycle, while \(1-8\) and \(6-3\) will be its next edges. For reproducibility, the full port word is:

~~~
1,8,28,78,66,16,84,34,2,7,9,72,22,27,77,4,5,31,81,19,69,
75,25,11,14,64,80,30,70,20,29,79,21,71,73,23,13,63,18,68,
76,26,74,24,12,37,87,82,32,17,67,33,83,86,36,85,35,65,15,
10,6,3
~~~

Read singleton symbols \(1,\ldots,12\) literally. Each consecutive port pair \((25+c,75+c)\) for \(-12\le c\le12\) stands for the block \(T_c=(25a_i+(-1)^{i+1}c)_i\); the reversed port pair stands for the reversed block. Because \(n\) is odd and \(a_1=1,a_n=3\), these are the exact endpoint labels. The word uses every port pair exactly once as a pair and every singleton once. The 36 junctions between tokens all have square sum, checked directly from the displayed word. Within a block the consecutive sums are 25 times old square sums, and, as \(c\) varies, the labels around a fixed \(a_i\) fill \([25a_i-12,25a_i+12]\). Thus the output visits every label in \([1,25n+12]\) once. Its endpoints are again 1 and 3, so the closing edge gives the asserted cycle. In particular \(8,2,7,4,5,6\) remain literal singleton labels in the required positions; this does not confuse a scaled copy of the motif with preservation of the small labels.

Starting at the inherited order-35 chain, iteration gives Hamilton cycles with all three motifs at every \(N=(71\cdot25^m-1)/2\), \(m\ge1\). Deleting an edge outside the motifs gives a Hamilton path containing them as well. The port word was found by a finite degree-constrained search, but its validity is certified by the 36 literal junctions and the all-order block argument, not by assuming the search is complete. This is an honest constructive analogue for a metagraph transfer: preserve marked endpoints and boundary edges by changing the gluing word, rather than treating an unmarked count as a sufficient input.

