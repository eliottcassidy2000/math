# A nodal-cylinder period obstruction to a Keller generating pair

Date: 2026-10-20. Status: **PROVED; independently reviewed.**

This is a separate continuation of
[the actual-source conductor packet](../source_conductor_collision_multiplicity_20261020/ACTUAL_SOURCE_CONDUCTOR_AND_COLLISION_MULTIPLICITY.md),
especially its minimum-defect nodal-cylinder reduction. That earlier packet is
frozen and is not changed here. This note excludes its remaining generic
conductor-degree-two case. It does not prove the complex Jacobian conjecture.

## 1. Exact algebraic statement

Put
\[
 R=\mathbb C[x],\qquad B=R[y],\qquad
 w=y^2-1,\quad z=y(y^2-1),\quad E=R+wB.
\]
Then
\[
 E=R[w,z]\cong R[W,Z]/(h),\qquad h=Z^2-W^3-W^2.
\]
For example, the equality follows by writing any polynomial in \(y\) as an
even polynomial plus \(y\) times an even polynomial.

**Theorem 1 (period rigidity for generating pairs).** Suppose \(P,Q\in E\) and
\[
                         R[P,Q]=E.                         \tag{1}
\]
Let \(I=(w,z)\subset E\), set \(p(x)=P(x,1)=P(x,-1)\) and
\(q(x)=Q(x,1)=Q(x,-1)\), and define the matrix \(L(x)\) by
\[
 \binom{P-p}{Q-q}
       \equiv L(x)\binom wz\pmod {I^2}.
\]
Then \(L\in\operatorname{GL}_2(R)\), its determinant
\(\lambda\) belongs to \(\mathbb C^*\), and the polynomial integral satisfies
\[
 \boxed{\quad
 A(x):=\int_{-1}^{1}P(x,y)Q_y(x,y)\,dy=\frac{8}{15}\lambda.
 \quad}                                                     \tag{2}
\]
Here and below integration means polynomial antidifferentiation followed by
evaluation at the indicated endpoints. No analytic path choice is involved.
In particular,
\[
 \boxed{\quad
 \int_{-1}^{1}(P_xQ_y-P_yQ_x)\,dy=0.
 \quad}                                                     \tag{3}
\]
Consequently no pair satisfying (1) has
\(P_xQ_y-P_yQ_x\in\mathbb C^*\). More strongly, if this Jacobian belongs to
\(\mathbb C[x]\), it is zero.

The distinction between membership \(P,Q\in E\) and generation (1) is essential.
The theorem does not assert that every pair in \(E\) is subject to (2).

## 2. Precisely the classical embedding input

We use S. S. Abhyankar and A. Sathaye, *Uniqueness of plane embeddings of
special curves*, Proc. Amer. Math. Soc. **124** (1996), 1061–1069,
[DOI](https://doi.org/10.1090/S0002-9939-96-03254-6),
[primary paper PDF mirror](https://scispace.com/pdf/uniqueness-of-plane-embeddings-of-special-curves-omujaliiec.pdf).

Their Theorem (1.3), p. 1062, applies over an algebraically closed field of
characteristic zero to
\[
 Y^n-X^m+\sum_{in+jm<mn}a_{ij}X^iY^j
\]
when \(m,n\) are distinct primes. Its affine curve has a unique plane embedding:
every alternative embedded image with isomorphic coordinate ring is ambiently
equivalent by a polynomial automorphism. We use equivalence of curves/defining
ideals, not exact equality of arbitrarily rescaled defining equations. This
quantifies over all alternate plane embeddings, not only bounded-degree ones,
and has no smoothness hypothesis. For \(h\), take \(n=2,m=3\);
the term \(W^2\) has weighted degree \(4<6\).

Only this single-curve theorem over \(\mathbb C\) is invoked. We do **not**
assume that fiberwise ambient automorphisms assemble into an automorphism over
\(\mathbb C[x]\), or that uniqueness of embeddings holds for arbitrary curves.

## 3. Proof of Theorem 1

### 3.1 The conormal determinant is a global constant

The evaluation homomorphism \(E\longrightarrow R\) at \(y=1\) has kernel
\(I=(w,z)\). Equation (1) implies
\[
                         I=(P-p,Q-q).
 \tag{4}
\]
Indeed, every element of \(E\) is a polynomial in \(P,Q\) with coefficients in
\(R\), and a polynomial minus its value at \((p,q)\) belongs to
\((P-p,Q-q)\).

Since \(h\in(W,Z)^2\), the classes of \(w,z\) form a free \(R\)-basis of
\(I/I^2\). By (4), the classes of \(P-p,Q-q\) generate that module.
Thus \(L\) is a surjective square matrix over \(R\), has a right inverse, and
\[
                      \det L=\lambda\in R^*=\mathbb C^*.
 \tag{5}
\]
This is the only relative rigidity needed in the proof.

### 3.2 The prescribed parametrization also lifts, one fiber at a time

Fix \(a\in\mathbb C\). Reducing (1) modulo \(x-a\) gives
\[
 \mathbb C[P(a,y),Q(a,y)]
          =E_0:=\mathbb C[y^2-1,y(y^2-1)].
 \tag{6}
\]
In particular the pair in (6) defines a closed embedding of the abstract
nodal curve \(C_0=\operatorname{Spec}E_0\) in the affine plane. The classical
theorem supplies an ambient polynomial automorphism carrying the standard
embedded curve \(C_0\) to this image. Its restriction may initially differ
from the prescribed isomorphism in (6) by an automorphism of \(C_0\).

Every automorphism of \(E_0\) extends uniquely to its integral closure
\(\mathbb C[y]\): the fraction field is \(\mathbb C(y)\), because \(y=z/w\);
\(y\) is integral, because \(y^2=w+1\); and \(\mathbb C[y]\) is normal.
The conductor is \((y^2-1)\mathbb C[y]\). For completeness, \(E_0\) consists
exactly of polynomials \(f\) with \(f(1)=f(-1)\); requiring both \(f\) and
\(yf\) to have this equality forces \(f(1)=f(-1)=0\), which proves the
conductor assertion.

An automorphism of \(\mathbb C[y]\) has form \(y\mapsto by+d\), \(b\ne0\).
Preservation of the conductor forces preservation of the unordered pair
\(\{-1,1\}\). Taking the sum and difference of its two images gives
\(d=0,b=\pm1\). Both possibilities lift to ambient automorphisms
\[
                         (W,Z)\mapsto(W,\pm Z).
 \tag{7}
\]
Correcting the initially chosen ambient equivalence by (7), we obtain an
ambient automorphism
\[
 \Phi_a=(F_a(W,Z),G_a(W,Z))
\]
whose restriction is exactly
\[
 F_a(w(y),z(y))=P(a,y),\qquad
 G_a(w(y),z(y))=Q(a,y).
 \tag{8}
\]
No consistency of the choices \(\Phi_a\) as \(a\) varies is required.

The Jacobian determinant of any polynomial automorphism is a nonzero
constant. Reducing the expansions of (8) at \((W,Z)=(0,0)\) modulo
\((W,Z)^2\), the absence of a linear term in \(h\) identifies its derivative
matrix with \(L(a)\). Therefore
\[
                  \det D\Phi_a=\det L(a)=\lambda.
 \tag{9}
\]
The translation \(\Phi_a(0,0)=(p(a),q(a))\) does not change this calculation.
It is important to correct the parametrization in (7): otherwise an
untracked sign could reverse the node loop.

### 3.3 Polynomial area and the contradiction

The polynomial one-form
\[
           F_a\,dG_a-\lambda W\,dZ
 \tag{10}
\]
is closed: its exterior derivative is
\((\det D\Phi_a-\lambda)dW\wedge dZ=0\).
Every closed polynomial one-form in two variables over \(\mathbb C\) is
exact. One elementary proof is to integrate its \(dW\)-coefficient with
respect to \(W\), then integrate the remaining coefficient, which by
closedness depends only on \(Z\).

Pulling (10) back along \((w(y),z(y))\) and integrating from \(-1\) to \(1\)
annihilates its exact part, because both endpoints map to \((0,0)\).
Using (8) yields
\[
 A(a)=\lambda\int_{-1}^{1}(y^2-1)(3y^2-1)\,dy
      =\lambda\left(\frac65-\frac83+2\right)
      =\frac{8}{15}\lambda.
 \tag{11}
\]
This holds for every \(a\), so it is the polynomial identity (2).

Finally, differentiating \(A\) and integrating one term by parts gives
\[
 A'(x)
 =\int_{-1}^{1}(P_xQ_y-P_yQ_x)\,dy
       +[P Q_x]_{y=-1}^{y=1}.
 \tag{12}
\]
Membership in \(E\) gives equal endpoint values for both \(P\) and \(Q\),
and differentiating the latter equality gives equal endpoint values for
\(Q_x\). Thus the boundary term vanishes. Equations (2) and (12) prove (3).
If the Jacobian is a polynomial \(j(x)\) independent of \(y\), equation (3)
reads \(2j(x)=0\). This proves the theorem. \(\square\)

## 4. Consequence for the actual-source conductor

Use exactly the hypotheses of the frozen conductor packet:
\[
 F=(P,Q):\mathbb C^2\to\mathbb C^2,\quad
 J\in\mathbb C^*,\quad P\text{ monic of degree }n\ge1\text{ in }y,\quad
 \mathbb C(x,P,Q)=\mathbb C(x,y).
 \tag{13}
\]
For the monic-in-\(v\) irreducible equation \(R(u,v,Z)\) of its primitive
source graph, put \(c=R_Z(P,Q,x)/J\). The earlier packet proves that
\(cB\) is the actual conductor and that
\[
 \deg_y c(a,y)=2\dim_{\mathbb C}
       \bigl(\mathbb C[y]/\mathbb C[P(a,y),Q(a,y)]\bigr).
 \tag{14}
\]
If \(F\) is not invertible, every source line has a collision, by
Gwozdziewicz's injectivity-on-one-line theorem as used there. If the generic
\(y\)-degree of \(c\) were two, that packet's Section 6, Theorem D proves:

1. every slice has degree-two conductor and a unique ordinary node;
2. its leading coefficient and discriminant are nonzero constants;
3. a polynomial shear and nonzero constant scaling of the source \(y\)
   coordinate makes the conductor ideal \((y^2-1)\);
4. in those coordinates, the **entire actual source algebra** is
   \(\mathbb C[x,P,Q]=\mathbb C[x]+(y^2-1)\mathbb C[x,y]\).

The coordinate change has constant nonzero Jacobian, so the transformed
pair remains a Keller pair. Theorem 1 now contradicts item 4.

**Corollary 2.** Under (13), a noninvertible Keller map must have generic
source-conductor degree at least four:
\[
                  \deg_y c\in\{4,6,8,\ldots\}.
 \tag{15}
\]
Equivalently, its generic source-line normalization defect \(\delta\) is
at least two. If the conductor has generic degree at most two, the Keller
map is invertible.

This is not a statement that every exceptional slice has defect at least
two: the argument does not exclude a degree drop to two inside a family
whose generic conductor degree is four or more. Nor is conductor degree
the generic sheet number of \(F\), or its total polynomial degree.

## 5. Hostile controls and what was not proved

1. **The surviving non-Keller node.** \(P=w,\ Q=x+z\) generates \(E\),
   has \(\lambda=1\), and has area \(8/15\). Its Jacobian is \(-2y\), whose
   integral over the node interval is zero. Thus this familiar hostile
   example agrees exactly with the obstruction rather than disproving it.

2. **Generic generation is insufficient.** \(P=w,\ Q=xz\) generates
   \(E\) after adjoining \(1/x\), but not over \(\mathbb C[x]\): at \(x=0\)
   the second generator vanishes. Its conormal determinant is \(x\),
   and its area is \(8x/15\), not constant.

3. **Even an invertible conormal matrix is insufficient without generation.**
   \(P=w,\ Q=z+xw^2z\) belongs to \(E\) and has conormal matrix the identity,
   but
   \[
       A(x)=\frac8{15}+\frac{64}{315}x.
   \]
   It is not a generating pair. For instance at \(x=-1\), both parameters
   \(y=\pm\sqrt2\) give \((P,Q)=(1,0)\), although they are distinct points
   of the standard nodal curve. Therefore this pair does not give a closed
   embedding of \(C_0\).

4. **There is no asserted relative ambient lifting theorem.**
   We have not proved that every pair in (1) is obtained from \(w,z\) by an
   automorphism of \(\mathbb C[x,W,Z]\) fixing \(x\). Pointwise lifting plus
   the global conormal determinant is enough for (2). The stronger relative
   question is unnecessary here and remains unclaimed.

5. **No broad embedding generalization.** The cited theorem has its precise
   prime-exponent hypothesis. The proof uses the standard nodal cubic and
   its two-element automorphism group, not an assertion that all rational
   one-place curves have unique plane embeddings.

## 6. Reproduction and next target

`verify_nodal_period.py` checks exact polynomial identities for a declared
finite collection of ambient polynomial automorphisms, including nonlinear
\(x\)-dependent shears, reversals, scalings, and translations. It separately
checks the hostile examples above and the integration-by-parts identity.
These checks validate signs, normalizations, and boundary behavior, not the
classical embedding theorem or its all-embedding quantifier.

The next actual-source target is generic conductor degree four (\(\delta=2\)),
where immersed slices may contain two ordinary nodes or one two-branch
tangency of intersection multiplicity two. The minimum-defect proof does not
automatically extend: collision pairs may move and the closed-path periods
can couple. A useful next question is whether the rank-two collision/period
module and its conormal determinants impose a nonconstant-Jacobian
obstruction under the exact actual-source algebra hypotheses.
