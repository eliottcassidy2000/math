# Fixed rational cylinders and finite period values

Date: 2026-10-20. Status: **PROVED; independently reviewed.**

This generalizes the frozen nodal-cylinder argument. It does not require
that every intrinsic automorphism of the curve extend to the ambient plane.

## Theorem

Let \(A=\mathbb C[f(y),g(y)]\subset\mathbb C[y]\) have fraction field
\(\mathbb C(y)\), with \(f,g\) nonconstant as a pair. Assume:

1. its affine plane curve has a unique embedding up to polynomial ambient
   automorphism;
2. there exist distinct \(r,s\in\mathbb C\) with
   \(f(r)=f(s)\) and \(g(r)=g(s)\).

Put \(R=\mathbb C[x]\) and \(E=R\otimes_{\mathbb C}A\subset R[y]\).
For every \(P,Q\) with \(R[P,Q]=E\), the polynomial
\[
                       \int_r^s P(x,y)Q_y(x,y)\,dy
 \tag{1}
\]
is constant in \(x\). Consequently, if the \((x,y)\)-Jacobian belongs to
\(R\), it is zero. In particular \(E\) admits no generating Keller pair.

The same proof works if \(A\) has finitely many plane embeddings up to
ambient equivalence, provided that finiteness is independently supplied.
No general finiteness theorem for embeddings is asserted here.

## Proof

### A. The intrinsic automorphism group is finite

The normalization of \(A\) is \(\mathbb C[y]\): this ring is integral over
\(\mathbb C[f]\) for any nonconstant \(f\), finite over \(A\), normal, and
has the same fraction field. Its conductor is a nonzero ideal
\(q(y)\mathbb C[y]\). Both \(r,s\) belong to its zero set.
Indeed, if \(h\) is in the conductor, then \(h,yh\in A\), so
\[
 h(r)=h(s),\qquad r h(r)=s h(s),
\]
which force \(h(r)=h(s)=0\).

Every automorphism of \(A\) extends to its normalization and has form
\(y\mapsto\alpha y+\beta\). It preserves the finite conductor support.
An affine map is determined by the images of two distinct points in that
support, so only finitely many such maps are possible. Denote the resulting
finite group by \(\Gamma\).

### B. A global conormal determinant

Let \(p_0=(f(r),g(r))\) be the chosen singular point of the curve. It is
singular because two normalization points lie over it; a smooth local
one-dimensional ring is normal and has only one normalization preimage.
For the maximal ideal
\[
                    \mathfrak m=(f-f(r),g-g(r))\subset A,
\]
the displayed classes form a basis of the two-dimensional conormal
\(\mathfrak m/\mathfrak m^2\): a singular plane equation has no linear term.

In \(E\), let \(I=R\mathfrak m\), and let \(p(x),q(x)\) be the values of
\(P,Q\) at \(r\), equivalently at \(s\). Generation gives
\[
                           I=(P-p,Q-q).
 \tag{2}
\]
Thus the matrix \(L_r(x)\) of these two classes in the fixed conormal basis
is invertible over \(R\), and
\[
                          \Delta_r=\det L_r\in\mathbb C^*.
 \tag{3}
\]

### C. Fiberwise equivalence gives only finitely many possible areas

Write \(\gamma(y)=(f(y),g(y))\). For each \(\xi\in\mathbb C\), uniqueness
of plane embedding gives
\[
                (P(\xi,y),Q(\xi,y))
                     =\Phi_\xi\bigl(\gamma(\sigma_\xi(y))\bigr),
 \tag{4}
\]
where \(\Phi_\xi\) is an ambient polynomial automorphism and
\(\sigma_\xi\in\Gamma\). Unlike in the node proof, no extension of
\(\sigma_\xi\) to an ambient automorphism is needed.

For each fixed \(\sigma\in\Gamma\), its pullback is an isomorphism between
the conormal at \(\gamma(\sigma(r))\) and the conormal at \(p_0\).
Using the bases induced by \(f,g\), let its matrix be \(T_\sigma\).
It is an invertible constant complex matrix. Equation (4), reduced to the
conormal at \(p_0\), gives
\[
                 \det D\Phi_\xi
                    =\Delta_r/\det T_{\sigma_\xi}.            \tag{5}
\]
Here \(\det D\Phi_\xi\) is the constant Jacobian of the ambient automorphism,
so evaluation at its indicated singular point loses no information.

For any ambient automorphism with determinant \(\delta\), the form
\(\Phi^*(U\,dV)-\delta U\,dV\) is a closed polynomial form and is exact.
Integrating along a normalization path with coincident image endpoints
annihilates the exact term. Thus (4)–(5) give
\[
 \int_r^s P(\xi,y)Q_y(\xi,y)\,dy
     =\frac{\Delta_r}{\det T_{\sigma_\xi}}
         \int_{\sigma_\xi(r)}^{\sigma_\xi(s)} f(y)g'(y)\,dy.
 \tag{6}
\]
The right side belongs to a finite set of complex numbers, since
\(\Gamma\) is finite. The left side is a polynomial in \(\xi\).
A polynomial over \(\mathbb C\) with finite image is constant.
This proves (1). No nonvanishing of any individual period is required.

Finally, fixed-endpoint integration by parts gives
\[
 \frac d{dx}\int_r^s P Q_y\,dy
       =\int_r^s(P_xQ_y-P_yQ_x)\,dy+[P Q_x]_r^s.
 \tag{7}
\]
The boundary term is zero because \(P,Q\in E\) have equal values at \(r,s\).
If the Jacobian equals \(j(x)\), the right side is \((s-r)j(x)\).
Since \(r\ne s\), \(j=0\). \(\square\)

## Examples supplied by the cited prime-exponent theorem

The Abhyankar–Sathaye theorem, with its source and exact scope recorded in
[the moving-algebra note](MOVING_HYPERELLIPTIC_GENERATOR_OBSTRUCTION.md),
applies to
\[
 A=\mathbb C[y^2,\;y\,g(y^2)]
\]
whenever \(g(0)\ne0\), \(g\) is nonconstant, and \(2\deg g+1\) is prime.
This includes:

- the ordinary node, \(\deg g=1\);
- two nodes or a tacnode, \(\deg g=2\);
- the corresponding fixed degree-six conductor models, \(\deg g=3\).

All these constant cylinders admit no generating Keller pair. This statement
does not classify moving families or arbitrary actual-source algebras.

## Failure boundary

For the standard node, \(P=y^2-1,\ Q=x\,y(y^2-1)\) generates only after
inverting \(x\), not over \(R\). Its period is \(8x/15\), so generic
generation cannot replace the entire algebra equality.

For a moving curve, even the reference periods can vary, and the finite
set in (6) need not exist. The family in Section 7 of the moving-algebra
note explicitly exhibits this. Its obstruction requires the additional
joint-node scalar descent and coefficient argument, not this lemma alone.
