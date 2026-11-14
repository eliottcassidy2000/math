# Hidden composition: what three probes recover

Write \(Q_T(z)\) for the ordered path-cover polynomial with \(z\) optional slots, \(H_T=Q_T(1)\), \(q_T=Q_T(2)\), and \(K_T=q_T-H_T\). This note uses two operations: ordinal stacking, which multiplies \(Q\), and the strong triangle wrapper with modules \(U\to W\to x\to U\) and \(x\) singleton.

## 1. An exact 4+7 to 5+6 alias and the depth that breaks it

There is an 11-vertex ordinal pair \(A,B\) with strong-component orders \((4,7)\) and \((5,6)\), respectively. It has the **same** \((H,q)=(615,49920)\), with the exact hidden difference
\[
Q_A(z)-Q_B(z)=2z^3(z^2-1)(z^2-4)(z^2-2).
\]
The four component profiles, in order, are \((5,32),(123,1560),(15,120),(41,416)\). The replay regenerates these from the literal outgoing row masks below, checks the factorization and verifies from the full small-order catalogue that this is the only common \((H,q)\) across the two strong-component-size splits. The ordered SCC sizes also prove the aggregates nonisomorphic, even after converse. It is **not** the same as the unpadded cyclic cuts with the two sizes on opposite sides of a singleton: their \(H\)-value lists are actually disjoint.

| order | outgoing row masks |
|---:|---|
|4| \(2,4,9,3\) |
|7| \(6,36,104,19,7,89,27\) |
|5| \(6,12,24,17,3\) |
|6| \(6,12,24,49,35,7\) |

The cyclic-wrapper polynomial, denoted \(\mathcal W_W(Q_U)\), is
\[
z\big[Q_U(z)Q_W(z)+
 (Q_U(z+1)-Q_U(z))(Q_W(z+1)-Q_W(z))+
 (Q_U(z-1)-Q_U(z))(Q_W(z-1)-Q_W(z))\big].
\]
Proof: in the inherited [odd-cycle packing formula](../../../01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md), a packing with no cross-cycle contributes \(zQ_UQ_W\). Any cross-cycle uses the unique \(x\) and nonempty directed paths in \(U,W\), hence only one such cycle occurs. Its oddness requires the two paths' total vertex count even. The unweighted subset sum over a nonempty distinguished path is \(Q(z+1)-Q(z)\) by convolution in the number of slots; its signed variant is \(Q(z-1)-Q(z)\) because \(Q_S(-1)=(-1)^{|S|}Q_S(1)\). Averaging the parity signs cancels the cycle's multiplicity two. This proves the displayed identity, extending the [two-slot specialization from Session 329](../session329_activation_and_regeneration_20261114/TWO_SIDED_ACTIVATION.md).

**Observation-depth theorem.** Let \(D=Q_A-Q_B\) for nonempty modules and let \(r\ge2\) be its *first different positive integer evaluation*. Wrap each once with the same arbitrary nonempty \(W\). For \(z<r-1\) the difference stays zero; at \(r-1\) it is exactly
\[
(r-1)\big[Q_W(r)-Q_W(r-1)\big]D(r)\ne0.
\]
Here \(D(0)=0\), and the slot difference in brackets is strictly positive (there are new covers using the last slot). Thus the distinguishing slot moves down exactly one level at each wrapper. For the 11-vertex pair, the first wrap keeps \(H\) equal for **every** \(W\), but its \(Q(2)\) gap is exactly \(30240[Q_W(3)-Q_W(2)]\). After a second common wrap the scalar Hamiltonian count detects it. In particular a single fixed wrapper, iterated and observed only through \(H\), recovers all positive integer evaluations of \(Q\) successively. The proof is triangular, not an assumption that one count identifies a graph. If two module polynomials agree identically, this observer cannot distinguish them.

At the singleton companion the literal first two strong depths (orders 13 and 15) have \((H_A,H_B)=(50535,50535)\) then \((2632275,2602035)\). These values are also tested by the full-path replay.

## 2. Three measurements identify the rank-one boundary, not the graph

For \(U,W\) of arbitrary positive order, write \(h=H_U,q=q_U,k=H_W,r=q_W\). A neutral singleton appended within \(U\) doubles \(q\) and leaves \(h\) fixed; similarly in \(W\). With independent scale factors \(x,y\in\{1,2,4,\ldots\}\), the count is
\[
F(x,y)=qrxy-qkx-hry+3hk.
\]
This is the exact response to two hidden ports, not a guess from the total order. Its matrix of unsigned products is rank one. Measure just
\[
N=F(1,1),\quad a=F(2,1)-N,\quad b=F(1,2)-N.
\]
**Theorem.** These three totals determine all future \(F(x,y)\) and the two ordered ratios \(h/q,k/r\). They need not identify the modules.

Indeed put \(p=qr\). The other products are \(qk=p-a\), \(hr=p-b\), \(hk=(N+p-a-b)/3\), so the rank-one relation forces
\[
2p^2-[N+2(a+b)]p+3ab=0.
\]
The quadratic appears to leave a second branch. If the true ratios are \(u=h/q,\ v=k/r\), then by the product of roots the rival would require ratios \((T(v),T(u))\) with
\[
T(t)=1-\frac{2}{3(1-t)}.
\]
The physical set of ratios lies in \((0,1/6]\cup\{1/4,1/2\}\): weak cuts give \(q\ge(n+1)H\) for \(n\ge5\); the four-vertex classes have maximum \(5/32<1/6\); smaller cases are explicit. The value \(1/9\) is impossible since \(H\) is odd and \(q\) even. Now \(T(1/2)<0\), \(T(1/4)=1/9\); and for \(0<t\le1/6\), \(T(t)\in[1/5,1/3)\), where the only possible physical ratio is \(1/4\), attained only if \(t=1/9\). Hence the rival cannot be physical. Recover the four products by the formulas above. The discarded information is reciprocal rescaling of the two profiles and the interior tournament classes, not an extra positive branch. The argument also proves the quadratic's two roots cannot coincide in the physical case, since the involution exchanges ratios between two disjoint domains.

One same-size example of why to probe: with a single added neutral vertex, the **cyclic** 4+7 and 5+6 cuts can both give \(H=50739\), though there is no unpadded overlap. Their extra-left and extra-right probes are respectively \((105203,104075)\) and \((104523,105411)\); the equal scalar hides different future responses.

This is an exact transfer of the *method* in [the affine Jacobian clutch calculation](../../../01-canon/theorems/THM-3225-affine-jacobian-clutch-resultant-and-two-boundary-no-escape.md): eliminate the hidden coordinate, but remove branches using genuine boundary data rather than declare every resultant root physical. Different proof ingredients apply in the two problems. The tournament reconstruction neither supplies a Keller mate nor proves the planar Jacobian conjecture.
