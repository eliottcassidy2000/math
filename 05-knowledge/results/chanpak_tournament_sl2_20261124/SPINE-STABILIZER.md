# The native height-preserving symmetry of the fixed-spine path cube

A vertex of the chart is a labelled tournament with fixed Hamilton spine \(0\to1\to\cdots\to n-1\). The free coordinates are its nonconsecutive pairs \((i,j)\), \(i<j-1\). Edges in this cube reverse one free arc. Write \(H\) for the number of Hamilton paths, and \(\tau\) for reflecting the spine labels while reversing every arc. The latter induces reflection on the free coordinates, so it is a cube isometry preserving \(H\).

**Theorem.** For \(n\ge4\), the only isometries of the fixed-spine cube preserving \(H\) on every vertex are the identity and \(\tau\). Moreover \(H\) at distance at most two from the transitive tiling already certifies this: any isometry fixing the transitive tiling and preserving \(H\) on that ball lies in the same two-element group. For \(n=3\) the group is trivial.

### Recovering the spine from its short probes

An isometry preserving \(H\) fixes the transitive tiling: if any tournament is nontransitive it has an induced directed triangle, whose three paths extend injectively by inserting successive outside vertices. A transitive tournament with the prescribed spine has one orientation, with every arc forward. Any cube isometry fixing this anchor is a coordinate permutation (its action on the neighbors fixes its action on all distance profiles).

Reversing the sole chord \((a,b)\) of span \(d=b-a\ge2\) gives \(H=1+2^{d-1}\). Without the reversed arc \(b\to a\) the only path is the increasing spine. Paths using it consist of two increasing runs separated by \(b\to a\); the \(d-1\) vertices strictly between \(a,b\) may go in either run, whereas the vertices below \(a\) and above \(b\) are forced to the initial and final runs respectively. Thus any proposed coordinate symmetry preserves span.

Write \(f_i=(i,i+2)\) for the shortest chords. The simultaneous reversal of \(f_i,f_j\) gives \(H=5\) when \(|i-j|=1\), and \(H=9\) when \(|i-j|\ge2\). For neighboring starts the four involved vertices have the five paths
\[
 (0,1,2,3),\\ (0,3,1,2),\\ (1,2,0,3),\\ (2,0,3,1),\\ (3,1,2,0).
\]
 For disjoint supports use ordinal product. If the intervals share just an endpoint, count separately no descent (one path), just one descent (two per chord), and both descents (four paths, from the two freely assignable midpoints). The coordinate permutation on the \(f_i\)'s therefore acts on a path graph, and is identity or reversal.

### The longer chords have no hidden permutations

Fix a chord \((a,b)\) of span \(d\ge3\) and probe with each \(f_i\). The minimizing probes are **exactly** \(i=a,b-2\). The details below classify all relative placements, up to reflecting the spine. Every path is counted by its increasing runs and its chosen descents, of which only the two reversed arcs are available.

- **End-aligned:** \(i=a\) or \(i=b-2\). By reflection assume \(i=a\). Both descents cannot occur, since they enter the same vertex. Among \(2^{d-1}\) long-descent paths, the short increasing arc \(a\to a+2\) occurs and is now forbidden exactly when \(a+1\) is in the first run and \(a+2\) is in the second. Remove \(2^{d-3}\). Of the two short-descent paths, both are permitted if \(d\ge4\); at \(d=3\), exactly one is forbidden by its use of \(a\to b\). Including the increasing path, the total is \(3+3\cdot2^{d-3}\) for \(d\ge4\), and 5 for \(d=3\).
- **Strictly inside:** \(a<i<i+2<b\). The short increasing chord occurs in long-only paths if its endpoints fall in the same run and its midpoint in the other; exactly a quarter of the assignments do so. The two short-only paths are still allowed. Paths using both descents contribute \(3\cdot2^{d-4}\) for each of the two orders. In the long-then-short order, the midpoint \(i+1\) of the short interval may lie in any of the three runs. All other vertices strictly between \(a,b\), apart from the two short endpoints and that midpoint, have two options (first or second below \(i\), first or third above \(i+2\)). The opposite order follows by reflection. Thus the total is \(3+3\cdot2^{d-2}\), again strictly larger.
- **Proper crossing:** \(i=a-1\) or \(i=b-1\). Assume \(i=a-1\). Among long-only paths, exactly \(2^{d-2}\) avoid the forbidden forward short arc. Both short-only paths are allowed. Paths using both descents are possible in either order, with \(2^{d-2}\) choices in each: assign the free vertices \(a+2,\ldots,b-1\) to the two runs that can contain them. In the long-then-short order those are the first and third runs, and in the short-then-long order the second and third. This gives \(3+3\cdot2^{d-2}\), strictly larger.
- **One endpoint shared, or disjoint:** the value is \(3(1+2^{d-1})\). Disjoint intervals give the ordinal product. If \([i,i+2]\) ends at \(a\), paths with both descents have \(b\to a\to i\) consecutively, with the \(d\) remaining intermediate vertices freely allocated to the first or last runs. The long-only and short-only contributions are unchanged, and the total is \(1+2^{d-1}+2+2^d\). The other touch case is its reflection.

Once the shortest chords are fixed individually, the unordered pair of minimizing indices determines \(a\) as its smaller index and \(b\) as its larger index plus two. Hence every coordinate is fixed. If the initial action reverses the shortest chords, first compose with \(\tau\).

### Scope for finite and infinite group encodings

This obstruction concerns ambient **native arc geometry plus the exact count on every rooted state**. It does not preclude a rich action on a specially selected subfamily, on a multigraph quotient without an isometric lift, or on port-transfer matrices. For example, the Golay code can support \(M_{24}\) as an edit-metric action inside a restricted framed family, and can be used with a constant-height auxiliary factor, but it cannot extend to an isometry preserving all ambient heights. Similarly, no faithful nonabelian subgroup of \(\mathrm{SL}_2(\mathbb R)\) acts in this ambient form. Under any collection of compatible level-preserving ambient actions the only possible level observations remain the independent two-element mirrors; requiring outer-shell restriction equivariance forces their choices to agree once there are nontrivial mirror pairs. This does not assert that all abstract metagraph automorphisms lift to this action.

## A positive \(\mathrm{SL}_2\) action after retaining switching gauges

There is a different, exact representation in which the full group does act. Take \(F=\mathbb R\) with \(\chi=\operatorname{sgn}\), or \(F=\mathbb F_q\), \(q\equiv3\bmod4\), with \(\chi\) the quadratic character. Choose a nonzero representative \(v_x\) for each point \(x\) of \(\mathbb P^1(F)\). For distinct points define the orientation sign \(S_{xy}=\chi(\det(v_x,v_y))\). The sign is antisymmetric because \(\chi(-1)=-1\). Changing representatives multiplies \(S_{xy}\) by a product \(s_xs_y\), precisely the operation of reversing all arcs across the cut between negative and positive gauges. For \(A\in\mathrm{SL}_2(F)\) write \(Av_x=\lambda_A(x)v_{Ax}\); then
\[
 S_{Ax,Ay}=\chi(\lambda_A(x))\chi(\lambda_A(y))\,S_{xy},\qquad
 \lambda_{AB}(x)=\lambda_B(x)\lambda_A(Bx).
\]
Thus projective transformations act exactly on the switching class, preserving the alternating triple sign \(S_{xy}S_{yz}S_{zx}\). The triple signs determine the switching class: to compare two sign systems with the same triple signs, fix a point \(x_0\), set the other gauges from their sign ratios on edges incident to \(x_0\), and use the triangle equations for every remaining edge. Over \(\mathbb R\) these are the cyclic orientation signs on the projective line. The projective action is faithful modulo the center, and the original \(\mathrm{SL}_2(F)\) acts faithfully on signed projective representatives \(v\) modulo multiplication by scalars with \(\chi=1\), with \(-I\) interchanging the two sheets. Its determinant-sign encoding is a **switching** action; it does not preserve Hamilton height across gauges or act by native pair relabellings on the fixed-spine cube.

Over the reals, a finite chart on affine points starts with a transitive order. Switch one vertex \(x\) with \(a\) predecessors and \(b\) successors, both positive. The resulting tournament is the directed triangle of blocks \(x\to A\to B\to x\), with \(A,B\) transitive of widths \(a,b\). A spanning path with \(x\) at either endpoint is forced for each endpoint. Otherwise choose a nonempty subset of \(B\) before \(x\) and a nonempty subset of \(A\) after it, freely; inside each segment the old transitive order is forced. Consequently
\[
 H=(2^a-1)(2^b-1)+2.
\]
There is exactly one directed Hamilton cycle: start at \(x\), use the whole \(A\) block and then the whole \(B\) block. The examples \((a,b)=(1,1),(1,2),(1,3),(2,2),(2,3)\) give heights \(3,5,9,11,23\) exactly. In particular the \(H=11\) realization has five cycle-closing and six open roots. The \(H=23\) realization is a projective-switch representative with that height, not the two-chain rail template: its 15 possible marked pairs give second-width heights among \(95,107,393,509\), never 79.

A more general real chart uses any binary gauge word on the transitive order. Within a gauge color, arcs advance in that order; between colors they go backward. A Hamilton path is therefore an alternating sequence of nonempty monochromatic runs that partition the vertices, with each run sorted increasingly and with each run's last position larger than the next run's first. This is an exact positive-word grammar for this switching class. It does not by itself cover the Chan–Pak spectrum.

Replay: run spine-stabilizer-check.py with Python 3.
