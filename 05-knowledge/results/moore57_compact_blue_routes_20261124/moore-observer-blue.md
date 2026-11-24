# Hamilton roots and blue lines of the Moore observer tournaments

For a degree-\(d\) diameter-two Moore graph on \(m=d^2+1\) labels, the observer tournament has sources \(u_0,\ldots,u_{m-1}\) in increasing order and observers \(w_{m-1},\ldots,w_0\) in decreasing order. Across the layers, \(u_i\to w_j\) iff \(i=j\) or \(i,j\) are graph-neighbors. Its involutive anti-automorphism \(\alpha\) exchanges \(u_i,w_i\). We work with the rigid, canonically decoded charts, so \(\alpha\) is its unique anti-automorphism. For degree 57 every statement remains conditional.

## Exact signed fold and cycle split

Every path-rooted representation is one rooted tiling since the tournament is rigid. Its tiling is blue iff the path is fixed by simultaneous reversal and application of \(\alpha\), namely iff it has form
\[
 p_1,\ldots,p_m,\alpha(p_m),\ldots,\alpha(p_1).
\]
Its first half chooses one sign per label. Consecutive sources must increase in label, consecutive observers decrease; a source-to-observer transition requires a Moore edge and an observer-to-source transition requires a Moore nonedge (labels are distinct). The final sign must be source for the central mate arc to agree. The initial sign distinguishes cycle closure: the path closes precisely when \(p_1\) is an observer. These conditions are necessary and sufficient, since the remainder of the spine is forced by \(\alpha\). Thus both blue counts are obtained by subset DP on the \(m\) labels with a signed last label, not the \(2m\)-vertex set. This gives blue open and closed incidences individually, not just parity.

## A local obstruction to a blue line between observer charts

Let \(P\) be a Hamilton path from \(x\) to \(y\). Reversing all nonspine arcs complements vertex scores \(s\mapsto2m-1-s\), but adds 1 at \(x\) and subtracts 1 at \(y\). All observer charts of given \(m,d\) have the same reversal-symmetric score multiset. Hence its physical mate has any chance of belonging to this observer family only if \(s(x)=s(y)+1\). For a blue path \(y=\alpha(x)\), this is \(s(x)=m\), restricting the first label to \(u_d\) (open) or \(w_{d+1}\) (closed).

There is a stronger all-order test. Work in the converse gauge \(T^P\), reversing only the spine arcs; it is the converse of the actual mate and belongs to an observer class precisely when the mate does. For a blue path its anti-automorphism is still \(\alpha\). If the target is a rigid observer chart, this is its unique mate-pairing as well. The orientation on each mate-pair identifies its source, and exactly one mate-pair orientation changed, the central one. If \(u_k\) is the end of the first half, the proposed target source layer is \(U-\{u_k\}+\{w_k\}\).

The observer mate of the new source \(w_k\) is \(u_k\). Its cross-layer graph degree is the number of retained sources beating \(u_k\): this is \(k-1\) if the preceding first-half vertex was a source (one arc from a lower index was reversed), or \(k\) if it was an observer. Therefore \(k\in\{d,d+1\}\). On \(U-\{u_k\}\), the reversed arcs come from adjacent same-sign runs in the first half, taking the mate image for observer runs. Each run increases in the original source order after reading it in the appropriate direction. Transitivity of the retained set forces these reversed arcs to be a matching of consecutive pairs in the retained order. Indeed any two consecutive reversals sharing an endpoint yield a directed triangle with their unflipped outer pair. If an isolated reversed arc spans a third retained index, neither subordinate arc to that index is reversed, again creating a triangle. Conversely disjoint consecutive reversals are just disjoint adjacent swaps.

For the admissible \(d\ge7\) charts whose late anchor set \(F=\{2d+2,\ldots,m-1\}\) dominates the Moore graph, this rules out **every** internal blue line to a rigid observer chart. The central label \(k\) has a neighbor \(z\in F\). This vertex \(u_z\) still beats \(w_k\) in \(T^P\): the only possible other flipped spine arc touching \(w_k\) and a retained source comes from an original observer-to-source nonedge transition, not this edge. More generally \(w_k\) has at most \(d+1\) incoming arcs from retained sources (its \(d\) graph-neighbors and possibly the one flipped nonedge). In a transitive source layer it must therefore be inserted after at most \(d+1\) retained vertices, but the late neighbor \(u_z\), even after removing \(u_k\) and the disjoint adjacent swaps, has zero-based rank at least \(2d\). That is impossible. The same proof applies to Petersen charts with both potential central labels \(3,4\) dominated by the two late anchors. This can be achieved by taking nonadjacent anchors and placing the three undominated labels in the three other dual index-pairs. It does not assert an obstruction to black lines in all orders or to nonrigid observer targets with several possible anti-automorphisms.

## Exact small class incidence census

The [Python replay](moore-observer-blue.py) enumerates the eight admissible pentagon observer classes, using degree-refined exhaustive canonicalization. It computes:
 
| Pentagon order | \(H\) | Blue open | Blue closed | Distinct blue mate classes |
|---|---:|---:|---:|---:|
| 02341 | 8517 | 57 | 74 | 129 |
| 02314 | 8513 | 57 | 74 | 126 |
| 02431 | 8097 | 57 | 70 | 126 |
| 02134 | 8173 | 59 | 70 | 127 |
| 01432 | 6461 | 57 | 56 | 113 |
| 01342 | 6821 | 59 | 60 | 119 |
| 01243 | 6821 | 57 | 60 | 117 |
| 01423 | 6473 | 57 | 56 | 113 |

There are 980 physical blue incidences, all external, meeting 857 distinct ordinary self-converse classes. In fact the eight source classes have **no** internal physical line, of either color: of all Hamilton path roots, 10,851 pass the necessary endpoint-score filter, and exhaustive canonicalization of those mates gives no target in the eight classes. These are counts of physical views and ordinary isomorphism classes, not a direct one-arc-edit graph.

For a Petersen example, order the two-subset labels as \(02,03,12,13,04,14,24,34,01,23\). The all-source-first spine gives one literal external blue line. Both of its endpoint tournaments are rigid; the following two rows are checked by the folded recurrence and, for total path and cycle counts, the [compiled subset replay](moore-observer-height.cpp). Directed cycles are identified by rotation.

| State | \(H\) | Hamilton cycles | Blue open | Blue closed |
|---|---:|---:|---:|---:|
| Observer | 2,065,541,352,033 | 48,704,047,433 | 960,507 | 1,278,014 |
| Canonical spine mate | 2,197,946,981,193 | 51,264,500,761 | 1,195,661 | 1,270,678 |

The closed ordinary fibre masses are twenty times the listed cycle counts. Subtracting the blue columns from the corresponding open and closed masses gives the black incidences. The original all-source-first spine is open; its specified mate is closed. The upper and lower states have substantially different path counts, so literal complement here is emphatically not a height-preserving symmetry.
