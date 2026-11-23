# Unrestricted flip squares as exact planar charts

**PROVED** for the identities, inverse and incidence bound. **FINITE-EXACT** for the marked examples. Start with the [exact but sparse one-apex inverse](ONE-APEX-INVERSE.md), the [joint apex cut](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md), and the prior [planar chart](../chanpak_planar_prime_ports_20261123/README.md). The new plane is not an arbitrary interpolation of vertex names: its axes are two **literal arc reversals**, allowed anywhere. The two planar charts below use different observables. Neither chart identifies modular, real or algebraic invertibility with Boolean target admission.

## 1. Two measured gains form a planar gradient

On one labelled tournament mark arcs \(e,f\), and let \(x,y\in\{0,1\}\) record their reversals. Interpolate the height by the unique multilinear polynomial
\[
 h(x,y)=h_{00}+ax+by+cxy,\qquad c=h_{11}+h_{00}-h_{10}-h_{01}.
\]
This is also its expectation when the two bits are independent, with probabilities \(x,y\) in the real square. All four heights are odd, so \(g=(h-h_{00})/2\) is integral. Put \(d=c/2\). The two incremental readings are
\[
 \nabla g=(a/2+dy,\ b/2+dx),\qquad \det D\nabla g=-d^2 .
\]
These are the measured gains (divided by two) of toggling each bit **from zero to one**, with the other bit held fixed. For \(d\ne0\) the chart is an affine automorphism over \(\mathbb Q\); over a field of prime characteristic it is invertible iff \(d\) is nonzero there. Over \(\mathbb Z\), its image is the translate of \(d\mathbb Z^2\) and it is unimodular exactly when \(d=\pm1\). The inverse is the literal swap and division by \(d\); a physical witness must still decode to bits 0 and 1. Over real weights this is genuine Jacobian invertibility, not an assertion about every lattice point being an endpoint count.

**Integral sevenfold square.** The row masks \((202,76,9,64,239,207,132,14)\) with marks \((3,4),(0,6)\) give \((h_{00},h_{10},h_{01},h_{11})=(31,209,37,217)\). The normalized gain map is \((89+y,3+x)\), with determinant \(-1\); both orders are increasing. The counts of Hamilton paths using **both** marked pairs are \((0,42,0,44)\): the unit is a cancellation of two large terms. On the other hand, at \(p=103\), the [one-apex contact chart](ONE-APEX-INVERSE.md#literal-checks) has no monotone route at any order. Keep a one-apex transitive lower tournament but let the edits hit **spine arcs**: masks \((254,508,504,240,224,448,128,256,89)\), flips \((1,7),(3,6)\), give \((103,411,221,721)\). All four tournaments are strong; \(c=192\) and the gain map has full rank even mod 7. A second strong disjoint square at 149 appears in the [exact checker](jacobian-squares-check.js).

## 2. Incident marks have an unconditional gap in the Jacobian

A directed Hamilton order contributes to \(c\) only if it uses both marked pairs; its contribution is the checkerboard sign of their two orientations. If the marks are incident, write them as consecutive edges of an undirected block \(u,v,w\). Gauge the bits so that zero on both means \(u\to v\to w\); then both admitted block orders \(u,v,w\) and \(w,v,u\) have positive sign. There is no cancellation, and
\[
 c\ge2H(T-\{u,v,w\}),
\]
with \(H(\varnothing)=1\). Reversing one gauge changes the sign of \(c\), not its magnitude.

In fact, take a Hamilton order \(x_1,\dots,x_k\) of the outside tournament. Write \(a_i=[x_i\to u]\), \(b_i=[x_i\to w]\), with \(a_0=b_0=1\) and \(a_{k+1}=b_{k+1}=0\). Inserting one of the two directed blocks at gap \(i\) gives
\[
 a_i(1-b_{i+1})+b_i(1-a_{i+1})\ \ge\ a_i+b_i-a_{i+1}-b_{i+1}.
\]
Summing telescopes to 2, and all these insertions give distinct full paths. Full paths whose deleted block leaves a backward splice only add to the count. The constant is sharp by making the outside tournament ordinally precede the whole block. Consequently the planar gradient of **every** incident face is an automorphism over \(\mathbb Q\), with \(|d|\ge H(T-\{u,v,w\})\). This survives independent probabilistic interpolation of other arcs because the gauge sign is fixed for every orientation of those arcs.

The same proof applies to a marked simple path of any **positive even** length: a mixed derivative counts the two block orders with the same sign, and its magnitude is at least twice the outside height. For odd length the block orders have opposite signs; empty outside already gives cancellation. For a set of marks containing a cycle or a vertex of degree at least three, the full mixed derivative vanishes, because no Hamilton path can contain all the marks. These are structural restrictions on the high-degree orientation polynomials.

## 3. The joint two-height chart: rank from a literal apex experiment

Let the original lower and upper hosts \(L,U\) differ by one or two arcs. Attach a common apex with contact set \(S\) to each. Use the joint map \(F=(H(L+x_S),H(U+x_S))\) and the [forward weighted-cut formula](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md). Each value is its base plus a Boolean quadratic; its linear coefficients are positive-weighted outdegrees \(A_i^L,A_i^U\) (zero only for a sink), and its \(i,j\) term is the negative symmetric slot weight \(W_{ij}^L,W_{ij}^U\). The weights are all positive.

**Physical local rank.** There is a two-bit face through \(S=\varnothing\) with nonzero Jacobian at a Boolean corner of weight at most one. To prove it, put \(A_i=(A_i^L,A_i^U)\) and \(M_{ij}=-(W_{ij}^L,W_{ij}^U)\). If two \(A_i\)'s are independent, use the sink corner. Otherwise all \(A_i\)'s lie on a line of positive slope \(\lambda\) (their sum has two positive coordinates). For each pair \(i,j\), at least one \(A_i,A_j\) is nonzero since a tournament has at most one sink. If \(M_{ij}\) leaves this line, the Jacobian at one of the two singleton corners is nonzero: the column containing the nonzero \(A\) can be kept fixed and the other is translated by \(M_{ij}\). If all \(M_{ij}\)'s lie on the line, then \(W^U=\lambda W^L\); equality of the scaled outdegrees would force the changed arcs in their \(L\) orientations to carry zero divergence with their positive weights. A nonempty one- or two-edge digraph cannot carry such a circulation. If both hosts are strong, some nonempty corner on that face has rank two and both extended tournaments stay strong: the determinant is affine and cannot vanish at all three nonempty corners but not the fourth.

**Sevenfold persistence obstruction.** If \(H(U)=7H(L)\), the common-contact difference \(H(U+x_S)-7H(L+x_S)\) is a nonzero quadratic. If it vanished identically, mixed coefficients would require \(W^U=7W^L\); singleton coefficients then force the same impossible circulation. On a Boolean cube any nonzero polynomial of degree at most two is nonzero on at least one quarter of the vertices: pin the complement of a highest-degree monomial, whose coefficient survives every pinning. Therefore at least a quarter of common attachments break the exact ratio. This does not prohibit carefully chosen attachments, as the [scalable exact families](TWO-APEX-SCALABLE.md) demonstrate in a different construction.

## 4. What the Jacobian analogy yields and where it stops

On a fixed two-contact face the joint map has form \(F=A+Bx+Cy+Dxy\). Its determinant is affine; the exact quadratic secant identity is
\[
 F(u)-F(v)=DF((u+v)/2)(u-v).
\]
So it is injective on the **open real square** iff its determinant has no interior zero. At an interior zero, choose two nearby points symmetric along a kernel direction. For the **closed** square, a zero on an edge need only cause a collision if its kernel has a direction tangent to that edge; a corner alone does not suffice. This is a precise boundary condition, not merely sampling. The secant identity belongs to the degree-two setting of [Wang's quadratic Jacobian theorem](https://doi.org/10.1016/0021-8693(80)90233-1); here it yields a physical-face criterion.

Over characteristic zero, a dominant chart with \(D=0\) is affine. If \(D\ne0\), project off \(D\) to get a nonzero linear output \(\alpha x+\beta y\) after deleting constants. If \(\alpha\beta\ne0\), take that output as \(u\), substitute for one variable, and complete the square in the other; affine input and an output shear give \((u,v^2)\). If exactly one of \(\alpha,\beta\) is nonzero, the second output is linear in the other input with a nonconstant affine coefficient in \(u\); translations and an output shear reduce it to \((u,uv)\). Thus the only nonconstant-Jacobian failures in this chart are a fold or a contracted line. The real quadratic context is treated more generally by [Nien, Peckham and McGehee](https://arxiv.org/abs/1507.02732). Neither of these normal forms transports the Boolean square to itself in general, so admission must use the original bits.

The same planar rigidity extends beyond quadratic total degree if **one input enters only affinely**. For a polynomial map \(F(x,y)=xA(y)+B(y)\) over characteristic zero, nonzero constant determinant implies \(\det(A,A')=0\). Ratios of nonzero components of \(A\) are then constant, so \(A=a(y)v\) for a fixed vector \(v\); \(a(y)\det(v,B')\) a nonzero constant makes \(a\) constant. In linear coordinates adapted to \(v\), the second output is affine in \(y\) with nonzero slope and the first is a constant multiple of \(x\) plus a polynomial in \(y\). Hence the map is a triangular automorphism. This works for Hamilton-count charts with one free arc and a second parameter tying any number of other arc probabilities, provided both readouts are linear path counts and the free arc occurs in neither of the tied assignments. It is a scoped Keller rigidity theorem, not a resolution of the unrestricted planar Jacobian conjecture.
