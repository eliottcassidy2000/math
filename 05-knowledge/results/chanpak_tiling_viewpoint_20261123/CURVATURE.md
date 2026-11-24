# Which balanced unit squares survive pair contraction?

**PROVED.** Suppose the editable arcs are disjoint pairs \(X=\{a,b\}\) and \(Y=\{c,d\}\). Assume **each** vertex outside the four sees the two members of \(X\) in the same direction, and likewise for \(Y\). Its direction to \(X\) may differ from its direction to \(Y\). Then all four heights form a balanced unit square (zero mixed difference, both edge increments of magnitude two) **if and only if** the four vertices form the square of heights \(1,3,3,5\), with an arbitrary transitive prefix and suffix. In particular every balanced unit square with minimum greater than one has an outside vertex which distinguishes the endpoints of at least one editable pair. The hypothesis is essential. The five-vertex square of heights \(11,13,13,15\) in [LINES.md](LINES.md) is an example.

## The contraction identity

Contract \(X,Y\) to singletons \(x,y\). Leave \(xy\) unspecified and let \(U_+,U_-\) denote its two orientations, with \(x\to y\) in \(U_+\). Write \(M_{ij}=[i\to j]\) for \(i\in X,j\in Y\). Let \(H_{st}\) count directed Hamilton paths when \(s=1\) means \(a\to b\) and \(t=1\) means \(c\to d\). Then
\[
 H_{11}-H_{10}-H_{01}+H_{00}
 =(M_{bc}-M_{bd}-M_{ac}+M_{ad})\,[H(U_+)-H(U_-)].
\]

For the mixed difference a word contributes zero unless it uses both editable pairs consecutively within each pair. Collapse the pairs inside such a word and sum with alternating signs over their internal orientations. If at least one other vertex separates the two pairs, flipping the internal order of either pair preserves every outside incidence, so its signed contributions cancel. When the two pairs touch in order \(X,Y\), summing their four possible orders reads the cross entries with signs \(+,-,-,+\), giving the displayed factor \(K\). The outside prefix and suffix admit exactly \(W_{xy}\) choices, the number of paths in \(U_+\) using \(xy\). With the pairs in order \(Y,X\) the factor is \(-K\), since the reverse cross indicators are \(1-M_{ij}\). All paths avoiding the contracted pair are common to the two orientations, hence \(W_{xy}-W_{yx}=H(U_+)-H(U_-)\).

For binary \(M\), \(K=0\) implies the two rows agree or the two columns agree. Transposing the corresponding marked vertices gives an isomorphism on both parallel flips, whose side length is then zero. If \(|K|=2\), the cross matrix is a checkerboard. Simultaneously transposing the two pairs exchanges opposite corners isomorphically, so flatness forces all four heights equal. The sole remaining type has one exceptional cross entry among three agreeing entries. By converse and relabeling use \(X\to Y\) except for \(d\to b\). Flatness is now equivalent to \(W_{xy}=W_{yx}\).

## The positive witnesses in each side

Fix \(d\to c\). Reverse \(a\to b\) to \(b\to a\). Outside \(d\), the endpoints \(a,b\) have identical incidences to the rest. All paths avoiding this editable adjacency cancel between the counts. Among the remaining words, interchanging the internal order of the pair is a bijection unless \(d\) touches it. Only the \(b\to a\) side gains words in this case, with contiguous triple \(dba\) or \(bad\). Two disjoint sets of gained words give the bound
\[
 H_{00}-H_{10}\ \ge\ W_{xy}+H(U-x).
\]
Expand an \(xy\) adjacency to \(badc\) for the first set. Replace \(y\) in any path on \(U-x\) by \(dbac\) for the second set. The end types of the inserted words are respectively \((X,Y)\) and \((Y,Y)\). This proves validity against every outside neighbor, and contraction recovers the inputs.

In a flat square, the common value \(W_{xy}=W_{yx}\) is positive. Indeed, take a Hamilton path of the vertices other than \(x,y\). Insert the pair before the first vertex which fails to dominate both, orienting the pair so its last endpoint dominates that vertex. Its predecessor, when present, dominates both endpoints. If no such vertex exists, append the pair. The empty remainder also allows either pair. Thus at least one orientation has positive marked adjacency mass. A side length of two now forces \(H(U-x)=1\). Taking converse and exchanging the roles of the two pairs gives \(H(U-y)=1\) from the other unit side.

A tournament has height one precisely when it is transitive. Indeed a cyclic triple has three paths and each extends injectively when adjoining more vertices, by inserting each new vertex at its first available position. Hence the unmarked remainder is a transitive order, and both \(x,y\) have an insertion cut in this order because deleting either leaves a transitive tournament. If their cuts differ, one orientation of \(xy\) is transitive, and the other creates a directed triangle through an intermediate vertex. Their heights differ, contradicting flatness. Thus both cuts agree, placing all four original vertices in one interval in a transitive order. The four local orientations have heights \(3,1,5,3\) in order \((0,0),(1,0),(0,1),(1,1)\). Ordinal transitive prefix and suffix preserve each height. This also proves the converse.

## A single flip with exactly one distinguishing vertex

**PROVED.** Suppose that outside \(a,b\) precisely one vertex \(d\) sees the two endpoints differently, with \(a\to d\to b\). All remaining incidences may be arbitrary. Contract \(a,b\) to one vertex \(x\), leaving the pair \(xd\) unspecified. Then
\[
 H(b\to a)-H(a\to b)=W_{dx}+W_{xd}.
\]
Indeed, in the positive-witness argument above, a gained word contains exactly one of the consecutive triples \(dba,bad\). Contract \(ba\) to \(x\). This is a bijection with words traversing respectively \(dx\) or \(xd\) in the appropriate contracted orientation. In particular this flip always strictly raises the count. It is a unit increment precisely when the two smaller-board adjacency masses sum to two. In an arbitrary environment this criterion needs only the two orientations on one fewer vertex.

## Incident editable pairs

**PROVED.** Every flat square with two nonzero increments uses disjoint editable pairs. In fact for two incident pairs \(vi,vj\), take their Boolean variables to indicate arcs directed *out* of \(v\). The mixed difference is strictly negative, always an even integer at most \(-2\). Only words traversing both editable adjacencies survive the mixed difference. They traverse consecutive triple \(ivj\) or \(jvi\), in one of the two mixed orientations, contributing with a minus sign. Contract the middle vertex out of the triple. The resulting multiplicities are precisely the two marked adjacency masses for the pair \(ij\) in the two orientations on the board without \(v\). Their sum is positive by the pair insertion argument above, and even because their difference is the difference of two odd tournament heights. This proves that incident edits never propagate an increment unchanged around a physical square, with no outside symmetry assumption.

Run `python3 -B curvature.py` for a direct finite check of the disjoint formula and classification over all marked exterior-homogeneous boards through order six and the incident sign law over all labeled tournaments through order five. Higher orders in both statements follow from the proofs.
