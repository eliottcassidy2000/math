# A projective switching action on tournament orientations

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
There is exactly one directed Hamilton cycle: start at \(x\), use the whole \(A\) block and then the whole \(B\) block. The examples \((a,b)=(1,1),(1,2),(1,3),(2,2),(2,3)\) give heights \(3,5,9,11,23\) exactly. In particular the \(H=11\) realization has five cycle-closing and six open roots. The \(H=23\) realization is a projective-switch representative with that height, not the prime two-chain rail template: the transitive blocks are proper nontrivial modules.

A more general real chart uses any binary gauge word on the transitive order. Within a gauge color, arcs advance in that order; between colors they go backward. A Hamilton path is therefore an alternating sequence of nonempty monochromatic runs that partition the vertices, with each run sorted increasingly and with each run's last position larger than the next run's first. This is an exact positive-word grammar for this switching class. It does not by itself cover the Chan–Pak spectrum.

For a fixed ordered spanning word \(p=(x_1,\ldots,x_n)\), define its **signed** path weight \(W(p)=\prod_{i<n}S_{x_i x_{i+1}}\), whether or not it is directed. Switching multiplies this by \(s_{x_1}s_{x_n}\), since all internal gauges cancel. Hence the signed endpoint sum matrix
\[
 B_{xy}=\sum_{p:\,x_1=x,\ x_n=y}W(p)
\]
transforms by \(B\mapsto DBD\), with \(D=\operatorname{diag}(s)\), in addition to the projective relabeling. Closing the same word takes \(W(p)S_{x_nx_1}\), which is genuinely invariant under switches. This makes the distinction between paths and cycles exact at the signed level. Neither the signed weight of a path nor that of a cycle decides whether every edge is positive; ordinary directed counts do not inherit this invariance.


On a fixed finite set of \(n\ge3\) real projective points, the stabilizing permutations of the entire switching class are exactly the rotations of the cyclic order, \(C_n\). The triple invariant recovers that cyclic order: fix a point, read off the linear order of all remaining points from their ordered triples with that point, and move the anchor. An order-preserving cyclic permutation preserves every triple sign, hence gives a switching equivalence. This says nothing about whether each rotation is induced by an actual real projectivity of that finite configuration; four marked points already carry cross-ratio information.

Compare [weighted and framed SL₂ transfers](SL2-TRANSFERS.md) and the [ordinary spine height stabilizer](SPINE-STABILIZER.md). Marking one gauge, retaining only the switching class, and tracking exact ordinary counts give three different symmetry questions.
