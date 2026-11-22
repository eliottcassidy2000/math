# The reversed arc certifies a planar Jacobian

**Proved for every single-arc tournament seam.** A two-point integer bracket certifies joint rank for the two normalized apex cuts. It is a polynomial contact result, not a proof that every target cut has a Boolean realization.

Let \(T_-\) have \(a\to b\) and \(T_+\) be the reversal of that arc only. Put
\[
 F(S)=(Q_{T_-}(S)/2,\ Q_{T_+}(S)/2).
\]
In the forward-arcs gauge of [CUT-GAUGE](../chanpak_mirimanoff_unit_lifts_20261122/CUT-GAUGE.md), both coordinates are integer multilinear quadratic cut polynomials. Write \(w^\pm_{ij}\) for their normalized two-slot pair weights. Crucially \(w^-_{ab}=w^+_{ab}=w>0\), since the arc \(ab\) is unused in a contributing two-slot cover with \(a,b\) on separate slots.

## 1. Marked two-point identity

Fix all contacts other than \(x=x_a,y=x_b\) at a tape \(z\). Write
\[
 F(x,y;z)=D+A x+B y+Mxy,\qquad M=(-w,-w)^{\mathsf T}.
\]
Let \(J(x,y;z)=\det(F_x,F_y)\) and \(d_a=A_+-A_-\), \(d_b=B_+-B_-\). Direct expansion gives
\[
 J(x,y;z)=\det(A,B)+w d_a x-w d_b y.                 \tag{1}
\]
In particular the bracket is affine, not quadratic, in the marked variables.

There is a canonical outside tape \(z_*\): give contact one to every outneighbor of \(a\) other than \(b\) (the same in both tournaments), and zero to its inneighbors. In either cut polynomial, an unchanged arc \(a\to v\) is then suppressed regardless of \(x\); an unchanged arc \(v\to a\) is suppressed in the same way. Hence \(A_- =w\), \(A_+=0\) and \(d_a=-w\) on this tape. Equivalently \(F_x=(w(1-y),-wy)\) there. It follows that
\[
 \boxed{J(1,1;z_*)-J(0,1;z_*)=-w^2.}                \tag{2}
\]
The same difference holds at \(y=0\). For strong sources, the two contact sets in (2) are nonempty and proper (they contain \(b\) and omit an inneighbor of \(a\)), so both give strong apex extensions of both sources. To see this elementary preservation, an old strong graph has a directed path from every old vertex to some predecessor of the apex and from some successor of the apex to every old vertex.

Consequences:
- The minor in the two *marked* directions is nonzero somewhere, not merely some unspecified minor of the full gradient. The exact same two Boolean points certify this over every field where \(w\ne0\).
- The gcd of the two integer bracket values divides \(w^2\). At every prime \(\ell\), at least one of the points has bracket valuation at most \(2v_\ell(w)\). Primes dividing \(w\) are thus bounded-depth singular primes, not necessarily singular everywhere to unbounded precision.
- At every characteristic not dividing \(w\), both cut polynomials are algebraically independent: the two-variable specialization with nonzero Jacobian has independent formal differentials and is dominant. One can see dominance directly by taking the formal expansion at a nonsingular point; its linear term is invertible, so a nonzero polynomial relation cannot have vanishing lowest homogeneous component after substitution.

**Uniform corner density.** On the outside Boolean cube \(d_a(z)\) is an affine polynomial which is nonzero at \(z_*\), also after reduction at any prime not dividing \(w\). Hence it is nonzero on at least half those tapes: if nonconstant, pin all variables except one with nonzero coefficient, giving at most one zero in each pair; if constant, all tapes survive. On each surviving tape and each fixed \(y\), (1) gives at least one nonvanishing value among \(x=0,1\). Therefore at least one quarter of all Boolean contacts have nonsingular marked minor over the field. A second finite certificate uses sparse pins: if \(d_a(0)=d_a(e_k)=0\) for every outside vertex \(k\), all coefficients of \(d_a\) vanish, contradicting \(d_a(z_*)=-w\).

## 2. Exact target inverse on the canonical marked face

Let \(g\) be the desired difference of the two cuts and \(L\) the desired lower cut. On \(z_*\) the gap is
\[
 F_+-F_-=g_0-wx+\eta y,\qquad \eta=d_b.
\]
Hence solve \(x=(g_0+\eta y-g)/w\) and substitute into \(F_-=L\). The lower read has exact form \(D_-+B_-y+w(1-y)x\); the fibre equation is
\[
 (D_--L)+B_-y+(1-y)(g_0+\eta y-g)=0.                \tag{3}
\]
This is degree two if \(\eta\ne0\), with leading coefficient \(-\eta\). If \(\eta=0\) it is linear except precisely at the contracted line \(x=B_-/w\), which maps to one target. These statements hold over fields with \(w\ne0\); over the integers one must additionally enforce division by \(w\) exactly. In the actual tournament chamber one must enforce \(x,y\in\{0,1\}\); noninteger or other integral roots do not supply contacts. In particular (3) is an exact inverse of this marked face, not of the entire contact space.

A related control makes the connection to the prior planar Jacobian packets precise: a bilinear map \(D+Ax+By+Mxy\) with \(M\ne0\) cannot have constant **nonzero** Jacobian. Constancy makes both \(\det(A,M)\) and \(\det(M,B)\) vanish, leaving \(A,B\) collinear with \(M\), so the constant is zero. Thus the marked-face argument uses a moving evaluation point, not a claim that a cut map is a Keller map. Compare the repository's Session267 planar proof (bracket and top-term cancellation) and Session269 moving-centre proof (avoid committing to one exceptional evaluation point); neither requires nor implies the present integer identity.

## 3. Formal prime-adic local inverse, and its limitation

For a prime \(\ell\), let \(P\) be the better of the two points in (2) and \(s=v_\ell(J(P))\le 2v_\ell(w)\). For every formal target \(\mathbf t\equiv F(P)\bmod\ell^{2s+1}\), there is exactly one formal marked pair in \(P+\ell^{s+1}\mathbb Z_\ell^2\) mapping to \(\mathbf t\). Indeed, write \(x=P+\ell^{s+1}u\), \(D_*=DF(P)\) and let \(H(u)\in\mathbb Z_\ell^2\) be the degree-two remainder. The equation for \(u\) becomes
\[
 u=\ell^{-s-1}D_*^{-1}(\mathbf t-F(P))
        -\ell^{s+1}D_*^{-1}H(u).
\]
Since \(\ell^sD_*^{-1}\) is integral by the adjugate formula, the first term is integral and the nonlinear term is an \(\ell\)-adic contraction on the full integral ball. Iteration gives unique existence. In particular if \(\ell\nmid w\), every target over the selected residue ball lifts in this formal sense, with unique lift in the input ball. Neither (2) nor this inverse says that an arbitrary such \(\ell\)-adic parameter is a physical Boolean contact or that a width gadget implements it. That physical parameter realization is a separate obstruction.

## 4. A quantitative statement for actual Boolean outcomes

On \(z_*\), the four physical cut vectors contain a horizontal segment of length \(w\) at \(y=0\) and a vertical segment of length \(w\) at \(y=1\): the endpoints are the four corner outputs, as follows from \(F_x=(w(1-y),-wy)\). Their convex hull has area at least \(w^2/2\). Use the horizontal segment as a base. If the vertical endpoints lie on opposite sides of its line, the two triangles on that base have disjoint interiors and their heights add to \(w\). If they lie on the same side, the farther endpoint is at height at least \(w\), and its one triangle suffices. The bound concerns the convex hull of four genuine outcomes, not a claim of lattice saturation or literal realization of intermediate counts. It says that correlated random choice among these four physical contacts alone already gives a two-dimensional set of expected joint cuts of quantitatively positive area.
