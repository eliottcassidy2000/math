# Conditional near-Hamiltonicity of Moore57

Assume a Moore graph on \(n=3250\) vertices with degree 57 exists. It contains a path through at least **3138 vertices** and a cycle through at least **2388 vertices**. These are finite guarantees, not a claim of Hamiltonicity.

Let \(A\) be its adjacency matrix and \(W=A+\tfrac12 I\). The exact Moore identity gives
\[
W^2=\tfrac{225}{4}I+J,\qquad W\mathbf1=\tfrac{115}{2}\mathbf1.
\]
Thus its operator norm on the centered space is \(15/2\). For any vertex set \(R\) of size \(r\), the nonnegative vector \(W\mathbf1_R\) is supported precisely on the closed neighborhood \(N[R]\). Cauchy–Schwarz gives
\[
|N[R]|\ge\frac{(\tfrac{115}{2}r)^2}{\mathbf1_R^\top W^2\mathbf1_R}
=\frac{13225\,r}{4r+225}.
\]

Choose a longest path, fix one of its ends, and let \(R\) be the set of free endpoints reachable by Pósa rotations fixing that end. Pósa's rotation lemma bounds the external neighborhood by \(2r-1\), so \(|N[R]|\le3r-1\). The two inequalities force \(r\ge1047\): the difference \((3r-1)(4r+225)-13225r\) is \(12r^2-12554r-225\), negative for positive integral \(r\le1046\).

For disjoint anticomplete sets \(R,T\), centered mixing with \(W\) gives
\[
529\,rt \le 9(3250-r)(3250-t).
\]
Indeed \(\mathbf1_R^\top W\mathbf1_T=0\); subtracting the constant terms gives a centered bilinear form of magnitude \((115/2)rt/n\), bounded by \((15/2)\sqrt{r(1-r/n)t(1-t/n)}\). Every rotated path uses the same vertex set as the longest path, so \(R\) is anticomplete to its omitted vertices \(T\). The right-hand bound on \(t\) decreases with \(r\); at \(r=1047\) it is \(165225/1471<113\). Thus the path omits at most 112 vertices.

The same mixing inequality shows that disjoint sets each of size 376 cannot be anticomplete: equal anticomplete sets have size at most \((3/(23+3))3250=375\). Use an edge between the first and last 376 vertices of the path to close a cycle, losing at most \(375+375\) vertices.

## Literature guardrail

[Brouwer and Haemers, *Hamiltonian strongly regular graphs*](https://aeb.win.tue.nl/preprints/hamilton2g.pdf) conjecture that Petersen is the only connected non-Hamiltonian strongly regular graph. Their applicable sufficient criterion (nonintegral least eigenvalue, or minus the least eigenvalue at most \(\mu+1\)) does not cover Moore57: the least eigenvalue is \(-8\) and \(\mu=1\). The [constant-expander Hamiltonicity theorem](https://arxiv.org/abs/2402.06603) assumes a sufficiently large expansion constant, so does not provide a certified numerical conclusion here. In particular these results are not substitutes for a Moore57 Hamilton-path or Hamilton-cycle proof.

The [large strongly regular graph theorem of Pyber](https://arxiv.org/abs/1409.3041) also gives only a sufficiently-large-order statement here, not an explicit bound covering 3250. On the other hand, a full graph Hamilton spine is explicit for d=2,3,7 in the accompanying replay; the Petersen example is a path and not a cycle. A full spine in the *complement* is an easier, different assertion.
