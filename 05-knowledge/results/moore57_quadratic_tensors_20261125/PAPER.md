# Testing Moore graphs by the tensor of routes through a matched pair

**PROVED, conditional on existence for the numerical restrictions.** A degree-57 Moore graph exists if and only if there is a 112-regular simple graph on 1,625 vertices with a binary label on each edge incidence satisfying the local table below. Its underlying 1,625-vertex graph must have eigenvalue 12 at least 104 times and cannot be strongly regular. The test has no distinguished vertex, coordinate field, or group-action assumption.

## What comes from the attached paper

The supplied 111-page preprint *The matrix multiplication exponent is nine fourths* (OpenAI, November 24, 2026) claims in Theorem 1.1 that the formal arithmetic exponent is \(9/4\) over every field. Section 1.1 gives the matrix multiplication tensor, §1.3 compares shared pair labels with a shared triple label, and §2 uses independent list-distinctness tests in its lower-bound argument. The all-field **upper** bound has a separate exposition by Amanov, van den Berg and Zuiddam. We attribute the matching lower bound to the new preprint. The graph theorem below uses just the support indices \((i,k),(k,j),(j,i)\) of ordinary matrix multiplication, and requires no complexity estimate or lower-bound claim. [Supplied preprint and source](https://github.com/openai/math/blob/main/preprints/The-matrix-multiplication-exponent-is-nine-fourths-November-24-2026/README.md) · [Independent upper-bound exposition](https://arxiv.org/abs/2610.30054)

The earlier rooted model retained the intermediate fiber when composing two permutation rows. The older Petersen contraction gave one smaller graph, \(K_5\), after forgetting the ends of the matched edges. The missing information in that contraction is one bit **per incidence**, since each inter-pair edge touches a particular end of each matched edge. Requiring two local labels at a common pair to agree recovers an actual two-step route. The four endpoint sign pairs then have to be tested for distinctness. [Rooted reciprocal Latin boards](../moore57_thirteen_20261125/LATIN.md) · [Petersen contraction](../planar_shadows_petersen_switching_20261116/README.md)

## Root-free exact construction

Let \(A\) be the adjacency of a hypothetical graph. Its standard Moore identity is \(A^2+A=56I+J\), with nonprincipal eigenvalues 7 and \(-8\). It has a perfect matching. Indeed, any cut with smaller side of size \(3\le s\le1625\) has at least \((57-7)s(3250-s)/3250\ge25s\) edges by the spectral bound on its centered indicator. A two-vertex side has at least 112 and a one-vertex side has 57. For any nonempty deletion set \(U\), each odd component of the remainder therefore uses at least 57 edges to \(U\); summing gives at most \(|U|\) such components. Connectivity and even order handle the empty deletion set. Apply Tutte’s perfect-matching theorem. The construction only needs one matching. [Tutte’s theorem and matching extendability in distance-regular graphs](https://arxiv.org/abs/1507.06254)

Contract any perfect matching. Mark the two endpoints of pair \(i\) by signs \(+,-\). There is at most one edge between two pairs: two with a common end give a triangle, and two with distinct ends give a square. Write \(R_{ij}=1\) for a cross edge and \(s_{ij}\in\{\pm1\}\) for the endpoint it uses at \(i\). At each pair 56 signs are positive and 56 negative, so \(R\) is simple and 112-regular. At a common neighbor \(k\) of \(i,j\), the wedge \(i,k,j\) goes through one physical vertex precisely when \(s_{ki}=s_{kj}\). It then delivers endpoint label \((s_{ik},s_{jk})\).

| Pairs \(i,j\) | Required multiset of delivered labels at coherent common neighbors |
|---|---|
| \(R_{ij}=0\) | \((+,+),(+,-),(-,+),(-,-)\) each exactly once |
| \(R_{ij}=1\) | only \((-s_{ij},-s_{ji})\), exactly once |

For nonadjacent pairs, all four pairs of physical endpoints need a unique midpoint. For adjacent pairs the cross edge takes care of its own endpoints. The two mates individually supply two further length-two routes. The remaining endpoint pair needs one external midpoint. In both cases any second delivery to the same endpoints creates a forbidden triangle or square. Conversely, given the table and balanced signs, split all 1,625 vertices into joined pairs and add the specified cross edges. Every pair of distinct vertices in the split graph has either its direct edge and no midpoint, or exactly one midpoint and no direct edge. Within a matched pair no outside vertex can touch both ends. Thus the result is 57-regular on 3,250 vertices and has diameter two and girth five.

This is an integer contraction of two incidence supports at \(k\), with the agreement factor \((1+s_{ki}s_{kj})/2\). Counting common neighbors without this factor retains routes whose two halves stop at opposite ends of the matched edge. Even the correct count four at a nonedge loses the joint condition on the two outer labels.

A search can reject a quotient immediately if a nonedge has fewer than four common neighbors or an edge has none. A nonedge with exactly four common neighbors forces agreement of the two incidence bits at each of the four middle pairs. At its two outer vertices the four incidences each split two and two, and the two resulting bipartitions cross in all four ways. For any fixed middle pair collect these compulsory equalities among its 112 incidences. A component larger than 56, or an inability to group the components into two sets of size 56, rejects the quotient before orienting any other incidence. Flipping all signs at one vertex just relabels a matched edge. Every automorphism preserving the matching induces a faithful action on the quotient, since fixing each pair must also fix its unique contact to every neighboring pair.

## Two small matrices retain the spectral obstruction

Define \(D_{ij}=R_{ij}s_{ij}\), and let \(S_{ij}=R_{ij}s_{ij}s_{ji}\), using zeros off the support. The signed adjacency \(S\) is symmetric while \(D\) need not be. In the orthonormal constant/difference basis within each matched pair, \(2A+I\) becomes
\[
 \begin{pmatrix} R+3I&D^t\\ D&S-I\end{pmatrix}.
\]
On a cross edge the factor 2 cancels the two basis denominators. The matched edge and the identity give the two diagonal shifts. Squaring the Moore identity yields exactly
\[
 \begin{split}
 (R+3I)^2+D^tD&=225I+8J,\\
 DD^t+(S-I)^2&=225I,\\
 D(R+3I)+(S-I)D&=0.
 \end{split}
\]
Conversely these three identities, with the same simple support and the incidence sign rule, give \((2A+I)^2=225I+4J\) after splitting the pairs and changing basis back. The integer local table and the matrix tests are equivalent for this chart.

Write \(f_+,f_-\) for the multiplicities of 12 and \(-18\) in \(R\); write \(g_+,g_-\) for the multiplicities of 16 and \(-14\) in \(S\). Then
\[
 f_+-g_-=104,\qquad g_+-f_-=105,\qquad \operatorname{rank}D\le1520.
\]
Indeed \(D\mathbf1=0\). On the 1,624-dimensional nonconstant pair-sum space put \(C=R+3I\), and on the full 1,625-dimensional pair-difference space put \(V=S-I\). The squared equations give spectra in \([-15,15]\) and \(DC=-VD\). Every \(C\)-eigenvector with eigenvalue \(c\) strictly in this interval is paired via \(D/\sqrt{225-c^2}\) with a \(V\)-eigenvector of eigenvalue \(-c\), and conversely via the transpose. The two-dimensional blocks contribute one eigenvalue of each sign \(\pm15\) to the full operator. The unpaired spaces are exactly its extremal eigenspaces. The full multiplicities of \(+15,-15\), after the constant eigenvalue 115 is removed, are 1,729 and 1,520 by the trace and dimension of \(A\). Thus \(f_++g_+-f_--g_-=209\), while \(g_++g_- - f_+-f_-=1\) by the dimension difference; these give the two displayed excesses. The number of interior pairs equals \(\operatorname{rank}D\) and cannot exceed 1,520.

Consequently the 1,625-vertex quotient must have at least four distinct eigenvalues. It is connected, has principal eigenvalue 112 and has eigenvalue 12. If there were just one other nonprincipal eigenvalue \(t\), the graph would be strongly regular. Its common neighbor number on nonedges would satisfy \(\mu=112+12t\), and evaluation at the constant vector gives \((112-12)(112-t)=1625\mu\). These force \(t=-61/7\), impossible for an algebraic integer. A connected graph with just two eigenvalues would be complete, also incompatible with the degree. Prescribing a strongly regular quotient before solving the bit constraints thus excludes every completion.

## The five-cycle census remembers the matching

A triangle in the quotient can have zero or one coherent corner. Indeed following the three cross edges and using an internal matching edge at each incoherent corner gives a simple cycle of length three plus the number of incoherent corners. Two or three coherent corners would give a square or triangle. Each quotient edge has exactly one triangle coherent at the opposite corner by the adjacent row of the table. Therefore there are exactly \(1625\cdot112/2=91,000\) quotient triangles with one coherent corner. These correspond bijectively to the pentagons in the original graph containing two matched edges. A zero-coherent triangle lifts to an alternating hexagon. Its number can depend on the matching.

For each matched pair \(i\), its 112 incident quotient edges give exactly 112 pentagons using its matching edge and one other: it is an incoherent corner exactly when it lies on the edge opposite the coherent corner. Each edge of the original graph belongs to \(56^2\) pentagons. Choose one other neighbor of each endpoint; the two chosen vertices are nonadjacent by girth and have a unique further common neighbor. Hence every perfect matching partitions the 58,094,400 pentagons by their number of matched edges as follows.

| Number of matching edges | 0 | 1 | 2 |
|---|---:|---:|---:|
| Pentagons | 53,089,400 | 4,914,000 | 91,000 |

The count with one matched edge follows from \(1625(56^2-112)\). This provides a global completion test without choosing any root.

## Small-graph controls and repair after an erased pair

In the 50-vertex Hoffman–Singleton graph use the five point rows and five line rows over the field of five elements. For every permutation of the five row indices, match every point to its unique neighbor in its assigned line row. All 120 resulting perfect matchings pass the joint endpoint table. Each has 150 one-coherent triangles in its 25-site quotient. Twenty have no zero-coherent triangles and 100 have exactly forty. They give two distinct characteristic polynomials for the ordinary quotient, with the same signed boundary multiplicity differences 3 and 4 in degree seven. All six Petersen matchings give the inherited quotient `K5` and its ten one-coherent triangles. Run `python3 -B matching_tensor.py` in this directory. [Exact replay on both known graphs](matching_tensor.py)

A missing incidence partition has an intrinsic recovery rule. Delete the two endpoints of one matched pair while remembering its 112 external contact vertices. These sites do not meet by an edge: two on the same side would make a triangle, and two on opposite sides a square. Two contacts on the same side already had the deleted endpoint as their unique midpoint. Two on opposite sides have a unique midpoint outside the deleted pair. Joining sites with an external common neighbor thus produces exactly the complete bipartite graph with parts of size 56, recovering the two sides up to exchange. This distinguishes a true pair completion from an arbitrary binary labeling of its contact sites.

### Recovering 28 erased incidence rows

**PROVED.** Keep the ordinary quotient \(R\) and erase the entire incidence partitions on any \(t\le28\) quotient vertices. All sign rows are uniquely recoverable up to swapping the two endpoints of each erased pair. For \(t\le55\), every affected incidence between an erased and an unerased pair is still recoverable. Both assertions assume a valid completion and give necessary consistency tests on partial data.

Fix an erased pair \(i=\{a,b\}\), let \(l=t-1\), and let \(h\) other erased pairs be neighbors of \(i\) in \(R\). The external endpoints on unerased pairs contacting \(i\) are known from their rows. Their true sides have sizes \(A,B\ge56-h\) and \(A+B=112-h\). On them, join two if they share a physical midpoint in an unerased pair. This graph is obtained from \(K_{A,B}\) by some deletions. An erased neighbor of \(i\) causes no deleted midpoint edges among these remaining contacts. Its contact endpoint cannot touch another contact, and its mate cannot touch a contact on that same side without creating a square. Of the other \(q=l-h\) erased pairs, each physical vertex can touch at most one contact on either side by girth, so at most \(2q\) edges disappear. The losses at a single contact have degree at most \(q\), since a contact can touch at most one endpoint of any erased pair.

For \(l\le54\), we have \(q<\min(A,B)\), excluding isolated vertices, and \(2q<A+B-2\). Disconnecting a complete bipartite graph of part sizes at least two without isolating a vertex requires at least \(A+B-2\) deletions: a component containing \(u,v\) sites on the two sides loses \(u(B-v)+v(A-u)\) edges, minimized on the corners of \([1,A-1]\times[1,B-1]\). Consequently the surviving graph is connected; its unique bipartition recovers both visible sides for each erased pair.

It remains to recover an edge between two erased pairs when \(t\le28\). Fix a recovered visible side \(C\) of pair \(i\), with hidden endpoint \(a\), and an endpoint \(y\) in another erased pair. Test which contacts \(c\in C\) have either a direct edge to \(y\) or a two-step route via an unerased pair. All these tests now use known incidences. If \(y\) contacts \(a\), every such contact already has its unique route through \(a\), so none is seen. Otherwise, every unseen contact needs a midpoint in an erased pair distinct from \(a\). Such a midpoint can touch at most one contact of \(C\). Vertex \(y\) has at most one neighboring vertex in each erased pair, including its own mate. At most \(t\) contacts can therefore be unseen. Here \(|C|\ge56-(t-1)>t\); a total absence of visible routes characterizes adjacency to \(a\). Repeating this for both sides recovers every remaining cross edge.
