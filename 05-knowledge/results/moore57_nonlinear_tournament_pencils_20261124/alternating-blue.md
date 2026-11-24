# Two alternating blue motifs detect exact Moore saturation

Let \(G\) be a \(d\)-regular simple graph on \(m\) vertices, with adjacency \(A\). Write \(N=J-I-A\) for its nonedge adjacency and \(k=m-1-d\). Form the usual self-converse observer tournament: its two copies of label \(i\) are a source \(u_i\) and observer \(w_i\), the two internal orders are opposite, and \(u_i\to w_j\) exactly for a graph edge or \(i=j\). The anti-automorphism \(\alpha\) swaps the copies. A **partial blue word** here is a palindromic Hamilton path on a selected set of mate-pairs, with respect to this specified \(\alpha\), even if a smaller induced tournament has extra symmetries. Count actual signed words, not quotient classes.

## A scalar tournament objective equivalent to the degree-57 problem

Let \(F_4\) count partial blue words on four mate-pairs with alternating first-half layers \(W,U,W,U\). These are closed and correspond precisely to injective walks with edge types \(N,A,N\) on the four labels. Let \(F_5\) use five pairs with first-half layers \(U,W,U,W,U\). These are open and correspond to injective walks of types \(A,N,A,N\). Neither statistic depends on the order of labels used inside a layer, because no transition stays in one layer.

For \(F_4\), the only possible repeated labels are the first and fourth: adjacent equality is excluded and gap-two repeats would demand one edge be both an edge and a nonedge. For \(F_5\) the possible collisions are \(1=4\), \(2=5\), or \(1=5\). These events are disjoint: simultaneous \(1=4,2=5\) would demand opposite types on the same unordered pair, and any other overlap produces adjacent equality. It follows that
\[
 F_4=mdk^2-\operatorname{tr}(NAN),\qquad
 F_5=md^2k^2-k\operatorname{tr}(ANA)-d\operatorname{tr}(NAN)
       -\operatorname{tr}((AN)^2).
\]
The first term in each counts unrestricted walks by regularity. For the gap-three collisions in \(F_5\), close a triangle word and extend one step, with \(k\) or \(d\) choices, respectively. The gap-four collisions give exactly the last trace.

With \(t,c\) the numbers of undirected triangles and (not necessarily induced) four-cycles in \(G\), respectively, regularity gives
\[
\begin{aligned}
 \operatorname{tr}(ANA)&=md(d-1)-6t,\\
 \operatorname{tr}(NAN)&=md(m-2d)+6t,\\
 \operatorname{tr}((AN)^2)&=m^2d^2-2dm(d+d^2)+2md^2+12t+8c.
\end{aligned}
\]
Indeed \(AJ=dJ\), \(AN=dJ-A-A^2\), and the standard closed walk decompositions give \(\operatorname{tr}(A^2)=md,\operatorname{tr}(A^3)=6t,\operatorname{tr}(A^4)=md(2d-1)+8c\). The \(2d-1\) part counts immediate out-and-back repetitions, leaving eight directed rooted traversals per four-cycle.

Consequently **lexicographically maximizing the two counts to their theoretical targets** means first requiring \(t=0\) by \(F_4\), then \(c=0\) by \(F_5\). One nonnegative scalar works as well: for \(d\ge3,m=d^2+1\), put \(\gamma=m-2d-2\). Then
\[
 F_5+\gamma F_4 = (F_5+\gamma F_4)_{t=c=0}-6t-8c .
\]
For \(d=2,m=5\), \(F_5\) alone has both coefficients strictly negative. At order \(m=d^2+1\), the condition \(t=c=0\) is equivalent to Moore diameter two: the radius-two list \(1+d+d(d-1)\) then has no duplicates and exhausts all vertices. Thus the explicit scalar target at \(d=57\) is **6,021,051,175,776,000**, where \(\gamma=3134\), \(F_4=1,886,906,112,000\) and \(F_5=107,487,420,768,000\). This is an exact tournament-motif rephrasing and an objective for search, not proof that the target is attained.

The local count changes can be scored without summing five-label words. Adding or deleting edge \(ab\) changes triangle count by the signed number of common neighbors. Adding a nonedge changes four-cycle count by \((A^3)_{ab}\), all length-three walks then being simple. For deleting an existing edge subtract \((A^3)_{ab}-\deg(a)-\deg(b)+1\), which removes the repeated-vertex walks. A degree-preserving graph 2-switch acts by eight literal observer-arc reversals (two mirrored arcs per undirected edge toggle), leaving the canonical spine intact. Reversing each mirror pair consecutively passes through a black position between blue endpoints. These are auxiliary literal free-arc walks, not the antipodal involution itself. This is a natural move for attempting to reach the theoretical value; no monotonic ascent to it is asserted.

The [independent replay](moore-alternating-blue-check.py) gets these known-case targets by direct injective-walk enumeration:
 
| Graph | \(F_4\) | \(F_5\) |
|---|---:|---:|
| Pentagon | 30 | 20 |
| Petersen | 960 | 2160 |
| Hoffman–Singleton | 604800 | 4057200 |

It additionally checks regular examples with triangles or squares, so the collision corrections are exercised away from a Moore equality case.

## Mate-arc degree separates spanning graph paths from blue paths

Weight each of the \(m\) arcs between mates by \(z\), and form the weighted directed Hamilton path and cycle polynomials in the observer tournament. In the original \(U\to W\) diagonal gauge their top coefficient \([z^m]\) counts oriented Hamilton paths and cycles of the **complement** of \(G\). In fact every mate-pair must be an adjacent block \(u_i,w_i\), and consecutive blocks connect precisely at nonedges. If all diagonal arcs alone are reversed, the top coefficient instead counts Hamilton paths and cycles of \(G\) itself, by blocks \(w_i,u_i\). The reversed-diagonal gauge still admits the same anti-automorphism \(\alpha\). Cycle counts here are modulo rotation, not reversal.

Every \(\alpha\)-palindromic spanning path uses exactly one mate arc: the middle arc has that form, and a nonmiddle mate arc would repeat at its mirror position. Thus the top graph-reading mode is entirely outside the \(\alpha\)-blue stratum for \(m>1\), even though the short alternating blue motif detects Moore equality. When this anti-automorphism is unique, the top mode consists solely of black rooted tilings.

The two pentagon gauges give 10 graph Hamilton paths and two oriented graph cycles. Petersen gives 240 paths and zero cycles in the graph, versus 95,520 paths and 6,432 cycles in its complement. Hoffman–Singleton has an explicit graph Hamilton cycle in Robertson coordinates: alternate pentagon \(P_i\) and pentagram \(Q_i\) for \(i=0,\ldots,4\), starting the pentagon at label \(j_i\), traverse the whole pentagon by steps \(+1\), and the whole pentagram by steps \(+2\). The matching \(P_{i,j}\sim Q_{k,ik+j}\) forces \(j_{i+1}=j_i-3-i\bmod5\); it closes after five rounds, since the total shift is \(-25\). The replay verifies this literal cycle. Thus any conjecture about the top graph-Hamilton-cycle coefficient at degree57 should distinguish the degree-three obstruction from the degree-seven positive case.

## Calibrated switch landscape: strict descent fails already at degree seven

Write the short-motif defect as \(E=6t+8c\), so decreasing \(E\) means improving the tournament objective. At order ten, [the exact cubic replay](moore-objective-cubic.py) finds all 21 ordinary graph classes using exact isomorphism checks and legal 2-switches. Completeness uses the [fixed-degree switch connectivity theorem](https://arxiv.org/abs/1503.06073). Every non-Petersen class has a strict improving switch.

That low-degree behavior does **not** generalize to degree seven. The [static 50-vertex fixture and exhaustive local check](moore-objective-d7-fixture.py) give a 7-regular graph with \(t=29,c=43,E=518\). Of its 22,223 legal 2-switches, none decreases \(E\), and exactly two are neutral. A neutral switch followed by an improvement is possible:
\[
\begin{array}{rcl}
 \{36\,46,\ 6\,31\}&\mapsto&\{6\,36,\ 31\,46\},\\
 \{28\,46,\ 24\,39\}&\mapsto&\{28\,39,\ 24\,46\}.
\end{array}
\]
The labels are zero-based; each two-number item is an unordered edge. The exact \((E,t,c)\) trace is \((518,29,43)\to(518,33,40)\to(516,34,39)\). Thus an actual degree57 strategy based on the alternating-blue objective should distinguish single-switch strict descent from neutral or compound moves. A natural sharply testable next question is whether bounded-length plateau blocks suffice in degree seven, before hypothesizing any such behavior at 57. Neither that question nor the existence of the degree57 equality state is settled here.
