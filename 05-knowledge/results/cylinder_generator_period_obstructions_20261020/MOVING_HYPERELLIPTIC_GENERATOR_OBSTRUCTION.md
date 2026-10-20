# Moving hyperelliptic source algebras: a polynomial-unit obstruction

Date: 2026-10-20. Status: **PROVED; independently reviewed.**

This is separate from the frozen nodal-cylinder packet. It treats an explicit
class of entire source algebras, not all Keller maps or all algebras of a given
conductor degree.

## 1. Setup and statements

Write \(R=\mathbb C[x]\), \(B=R[y]\), and
\[
 g(x,w)=c+\sum_{i=1}^{d}a_i(x)w^i,\qquad
 c\in\mathbb C^*,\quad a_d\ne0,\quad d\ge1.
 \tag{1}
\]
Set
\[
 w=y^2,\quad z=y g(x,y^2),\qquad
 E_g=R[w,z]=R[w]\oplus y g(x,w)R[w]\subset B.
 \tag{2}
\]
The direct sum in (2) is the even/odd decomposition. The plane presentation is
\[
 E_g\cong R[W,Z]/\bigl(Z^2-Wg(x,W)^2\bigr).
 \tag{3}
\]

We call a fiber **nonempty in collisions** if \(g(\xi,w)\) is nonconstant.
Equivalently, not all of \(a_1(\xi),\ldots,a_d(\xi)\) vanish.

**Theorem A (pure period obstruction).** Suppose:

1. every fiber is nonempty in collisions;
2. for every \(\xi\) in some nonempty Zariski-open subset of \(\mathbb C\),
   the affine curve \(Z^2-Wg(\xi,W)^2=0\) has a unique plane embedding, in
   the sense of ambient equivalence of every isomorphic embedded image.

Then for every \(P,Q\in E_g\) satisfying the **entire algebra equality**
\[
                         R[P,Q]=E_g,                         \tag{4}
\]
the Jacobian \(P_xQ_y-P_yQ_x\), if it belongs to \(R\), is zero.

A sufficient condition for hypothesis 2 is that \(2d+1\) is prime: every
fiber of full degree has leading exponents \((2,2d+1)\). No constancy of
the leading coefficient is required. In particular all \(1\le d\le3\)
qualify, corresponding to exponents \(3,5,7\).

**Theorem B (prime-odd-degree-\(g\) Keller sector).** Let \(g\) be as in (1),
and suppose \(2d+1\) is prime. No pair generating \(E_g\) over \(R\) has
nonzero constant \((x,y)\)-Jacobian. Common zeros of the nonconstant
coefficients \(a_i\) are allowed. In particular this excludes every
nontrivial cubic-or-lower-\(g\) model.

If \(g=c\) is constant, \(E_g=B\); this excluded trivial case
does admit generating Keller pairs, for example \(P=x,Q=y\).

Theorem A uses only the stated single-curve embedding input and elementary
algebra/calculus. Theorem B's additional common-zero case uses the classical
injectivity-on-one-line theorem, separately identified below.

## 2. Classical inputs and their exact scope

The embedding input is Abhyankar–Sathaye, *Uniqueness of plane embeddings of
special curves*, Proc. Amer. Math. Soc. 124 (1996), 1061–1069,
Theorem (1.3), p. 1062:
[DOI](https://doi.org/10.1090/S0002-9939-96-03254-6).
For distinct primes \(m,n\), an equation \(V^n-U^m\) plus strictly lower
\((n,m)\)-weighted terms has a unique plane embedding over an algebraically
closed characteristic-zero field. This includes singular curves.
If \(\deg_w g(\xi,w)=e\) and its leading coefficient is \(l\ne0\), replace
\(Z\) by \(lZ\) in (3); the resulting equation has leading terms
\(Z^2-W^{2e+1}\) and all other terms have lower weight.
Thus hypothesis 2 holds on the full-degree open set when \(2d+1\) is prime.
Special fibers may drop to degrees whose leading odd exponent is composite;
no embedding theorem is invoked for those fibers.

The source check, including indexed original pp. 1061–1063 and the
curve/defining-ideal convention, is recorded in the independently reviewed
[nodal packet](../nodal_cylinder_period_obstruction_20261020/SOURCE_VERIFICATION.md).
No relative ambient lifting over \(R\) is assumed.

The separate one-line input is Gwozdziewicz, *Injectivity on one line*,
Theorem 1.1:
[primary paper](https://arxiv.org/pdf/alg-geom/9305008).
A complex polynomial Keller map injective on one affine source line is an
automorphism. It is used only in Section 8.

## 3. Fibers, normalization, conductor, and parametrized lifts

Fix a fiber and abbreviate its nonconstant polynomial by \(g(w)\), with
\(g(0)=c\ne0\). Its curve ring is
\[
 A=\mathbb C[y^2,y g(y^2)]
   =\mathbb C[w]\oplus y g(w)\mathbb C[w].
 \tag{5}
\]
Its fraction field is \(\mathbb C(y)\), because \(y=z/g(w)\). Its normalization
is \(\mathbb C[y]\): \(y^2=w\) is integral and \(\mathbb C[y]\) is normal.
Its conductor is
\[
                         g(y^2)\mathbb C[y].                 \tag{6}
\]
Indeed, for \(f=e(w)+y o(w)\), the two conditions \(f,yf\in A\) are exactly
\(g\mid o\) and \(g\mid e\). These two conditions also imply
\(f\mathbb C[y]\subset A\), by the even/odd decomposition.

All normalization branches are immersed: \(w_y=2y\) is nonzero away from
zero, and \(z_y(0)=c\ne0\). Every root \(r\) of \(g\) is nonzero and gives
the collision \(y=\pm\sqrt r\) over \((W,Z)=(r,0)\). These are precisely
the singular points of (3), since its two partial derivatives are
\[
 -g^2-2Wgg_W,\qquad 2Z.
 \tag{7}
\]
At \(W=0,Z=0\), the first derivative is \(-c^2\ne0\).

Every automorphism of \(A\) extends to \(y\mapsto\alpha y+\beta\) on its
normalization and preserves the conductor divisor. The support of this
divisor is a nonempty symmetric finite set with at least two elements,
whose sum is zero. Preservation of its average forces \(\beta=0\).
Preservation of multiplicities, followed by evaluation at \(y=0\), gives
\[
                         g(\alpha^2 W)=g(W).                \tag{8}
\]
Consequently every abstract automorphism lifts ambiently by
\[
                         (W,Z)\mapsto(\alpha^2W,\alpha Z).
 \tag{9}
\]

By (4), every specialized pair gives a closed embedding of the same abstract
curve. On the open set supplied by hypothesis 2, ambient equivalence of its
image is available. Composing by (9) corrects any discrepancy in the chosen
abstract parametrization. Therefore, for each \(\xi\) in that open set, there
is an ambient polynomial automorphism
\(\Phi_\xi\) whose restriction is exactly
\[
       (P(\xi,y),Q(\xi,y))
          =\Phi_\xi\bigl(y^2,y g(\xi,y^2)\bigr).
 \tag{10}
\]
No dependence of \(\Phi_\xi\) on \(\xi\) is required.

## 4. A repeated generic root already contradicts a nonzero base Jacobian

Assume that \(g\) has a repeated root \(r\) over an algebraic closure of
\(K=\mathbb C(x)\). At that root,
\[
                         g=g_W=g_x=0.                       \tag{11}
\]
For the last identity, differentiate \(g(x,r)=0\), using \(g_W(x,r)=0\);
in characteristic zero the derivation on \(K\) extends to algebraic
extensions. Equivalently a repeated irreducible factor of \(g\) divides
both partial derivatives. Also \(r\ne0\), since \(c\ne0\).

Choose polynomial representatives \(F,G\in R[W,Z]\) of \(P,Q\). At
\(y=\pm\sqrt r\) one has \(Z=Z_x=Z_y=0\), so
\[
 P_xQ_y-P_yQ_x
   =2y(F_xG_W-F_WG_x)(x,r,0).                                \tag{12}
\]
The two values are opposites. If the Jacobian is \(j(x)\in K\), they are
also equal, so \(j=0\).

This argument needs only membership \(P,Q\in E_g\), not generation or
embedding uniqueness. In proving Theorem A for \(j\ne0\), we may therefore
assume that \(g\) is squarefree over \(K\).

## 5. The ambient Jacobian scalar descends to a constant

Continue with squarefree generic \(g\). Using the global representatives
above, put
\[
 D(x,W)=\bigl(F_WG_Z-F_ZG_W\bigr)(x,W,0)\in R[W].
 \tag{13}
\]
For fibers in the open set of hypothesis 2, the derivatives of these
representatives at each singular point agree with those of \(\Phi_\xi\).
Their difference as functions on the ambient plane is a multiple of the
defining equation, whose value and gradient vanish there. Thus \(D(\xi,r)\)
has the same nonzero value for **all** roots \(r\) of \(g(\xi,W)\): namely
the constant Jacobian of \(\Phi_\xi\).

Divide \(D\) by \(g\) over \(K\). At infinitely many fibers of full degree
with distinct roots, its remainder has degree less than \(d\) and is
constant on all \(d\) roots. Hence that remainder is a scalar
\(\lambda\in K\):
\[
                            D\equiv\lambda\pmod g.           \tag{14}
\]
For \(d=1\) this conclusion holds directly by division.

Write \(\lambda=U/V\) in coprime terms in \(R\). The polynomial \(g\) is
primitive over the UFD \(R\), since its constant coefficient \(c\) is a
unit. Gauss's lemma applied to (14) yields
\[
                         g\mid VD-U\quad\hbox{in }R[W].
 \tag{15}
\]
If \(V(\xi)=0\), coprimality gives \(U(\xi)\ne0\). But specializing (15)
would make the nonconstant polynomial \(g(\xi,W)\) divide the nonzero
constant \(-U(\xi)\), impossible. Thus \(V\) is constant and \(\lambda\in R\).
Applying Gauss again gives \(g\mid D-\lambda\) over \(R[W]\).
For **every** fiber, including those outside the uniqueness open set,
generation makes the cotangent matrix at any singular point invertible:
the classes of \(W-r,Z\) form its two-dimensional conormal basis, and the
classes of \(P-P(r),Q-Q(r)\) generate the same maximal ideal modulo its
square. This matrix is precisely the derivative matrix of \(F,G\) there.
Consequently \(D(\xi,r)\ne0\), without any ambient lifting assertion for
that exceptional fiber.
Every fiber has at least one root, so \(\lambda(\xi)\ne0\) for every \(\xi).
Therefore
\[
                              \lambda\in\mathbb C^*.          \tag{16}
\]

This is the crucial descent step. Individual functions on the finite root
cover can be nonconstant units; it is their **common value across every
node of a fiber** that descends and is forced to be constant.

## 6. Period transport and an elementary coefficient contradiction

For every \(\xi\) in the uniqueness open set, the polynomial one-form
\[
                \Phi_\xi^*(U\,dV)-\lambda W\,dZ
 \tag{17}
\]
is closed, hence exact. Pullback along the normalization and integration
between the two preimages of any singular point kills the exact part.

Work on a small analytic open set inside the uniqueness open set, where
all roots of \(g\) are distinct and of full degree. Choose one analytic
root \(r(x)\) and a square root
\(t(x)\), so \(t^2=r\ne0\). Equation (17) implies
\[
 \int_{-t(x)}^{t(x)}P Q_y\,dy
       =\lambda\int_{-t(x)}^{t(x)}y^2
                       \frac{\partial}{\partial y}\bigl(y g(x,y^2)\bigr)\,dy.
 \tag{18}
\]
For any pair with equal target values at the two moving endpoints,
Leibniz's rule and integration by parts give
\[
 \frac d{dx}\int_{u(x)}^{v(x)}P Q_y\,dy
       =\int_{u(x)}^{v(x)}(P_xQ_y-P_yQ_x)\,dy.
 \tag{19}
\]
The boundary term is the common value of \(P\) times the difference of
the total \(x\)-derivatives of the equal endpoint values of \(Q\); it is zero.

Suppose that the Jacobian of \(P,Q\) is \(j(x)\in R\). The Jacobian of the
standard pair \(w=y^2,z=y g(x,y^2)\) is
\[
                      -2\sum_{i=1}^d a_i'(x)y^{2i+2}.
 \tag{20}
\]
Differentiate (18), use (19) on both sides, and divide by \(2t\). This gives,
for every generic root \(r\),
\[
 j=-2\lambda H(r),\qquad
 H(W)=\sum_{i=1}^d\frac{a_i'}{2i+3}W^{i+1}.                  \tag{21}
\]
Assume \(j\ne0\) and put \(\kappa=-j/(2\lambda)\in R\setminus\{0\}\).
Since \(g\) is squarefree,
\[
                       H-\kappa=(A W+B)g                   \tag{22}
\]
over \(K[W]\), with \(A,B\in K\); the degree bound is \(d+1\).
Comparing the constant and linear coefficients gives
\[
                         B=-\kappa/c,\qquad
                         A=\kappa a_1/c^2.                 \tag{23}
\]
Comparing leading coefficients then gives
\[
                    a_d'=(2d+3)\kappa a_1a_d/c^2.            \tag{24}
\]
If \(a_1\ne0\), the right side is a nonzero polynomial of degree at least
\(\deg a_d\), while a nonzero derivative \(a_d'\) has smaller degree.
If \(a_d\) is constant, the left side is zero and the right side is not.
Either way (24) is impossible.

If \(a_1=0\), (23) gives \(A=0\), \(B\ne0\). The coefficient of \(W^2\)
in (22) forces \(a_2=0\). Inductively the coefficient of \(W^i\) is
\(a_{i-1}'/(2i+1)=Ba_i\), so it forces every remaining \(a_i=0\), contrary
to \(a_d\ne0\). For \(d=1\), \(a_1=0\) already contradicts the definition
of \(d\).

This proves Theorem A. \(\square\)

## 7. The quadratic case and a hostile moving family

For \(g=aW^2+bW+c\), equation (22) implies, if \(a,b\ne0\),
\[
 5a'(b^2-ac)=7ab b',\qquad
 j=-\frac{2\lambda c^2}{7}\frac{a'}{ab}.
 \tag{25}
\]
The leading-coefficient argument above proves the contradiction without
dividing by \(a b\), and therefore makes the degeneracies easier to track.

The moving family
\[
        w=y^2,\qquad z=y(y^4-s(x)y^2+c),\quad c\ne0
 \tag{26}
\]
is immersed on every slice. Its conductor is
\((y^4-s(x)y^2+c)B\), and its singular root cover
\(W^2-s(x)W+c=0\) has a nonconstant unit \(W\) whenever the cover genuinely
varies. Thus a unit on that cover is not automatically constant.

If \(a,b\) are its roots and a local square root of \(a\) is chosen, its
standard loop period is
\[
 \frac8{105}a^{5/2}(3a-7b),\qquad
 \frac d{dx}(\text{period})=\frac45s'(x)a^{5/2}.
 \tag{27}
\]
Its Jacobian is \(2s'(x)y^4\), not a nonzero base polynomial. This is a
control showing that moving periods really exist, and explaining why
the common-scalar descent in Section 5 is necessary.

## 8. Common-zero fibers: the separate one-line shortcut

In the setting of Theorem B, suppose that all nonconstant coefficients
vanish at some \(\xi\). Then
\[
                (E_g)_\xi=\mathbb C[y^2,c y]=\mathbb C[y].
 \tag{28}
\]
If (4) holds, \(y\) is a polynomial in \(P(\xi,y),Q(\xi,y)\), so the source
line \(x=\xi\) is mapped injectively.

If the Jacobian is a nonzero constant, Gwozdziewicz's cited theorem makes
\((P,Q)\) a polynomial automorphism. It would follow that
\[
                      B=\mathbb C[P,Q]\subset R[P,Q]=E_g.
 \tag{29}
\]
But (2) shows \(y\notin E_g\) whenever \(g\) has positive \(W\)-degree:
its odd coefficient would require \(1\in gR[W]\). This contradicts (29).

If there is no common-zero fiber, Abhyankar–Sathaye supplies hypothesis 2
on the full-degree open set because \(2d+1\) is prime, and Theorem A applies.
This proves Theorem B. \(\square\)

## 9. An exact recognition criterion for actual source algebras

**Proposition C (quadratic-coordinate recognition).** Let
\(S=R[P,Q]\subset B\) contain \(y^2\). Suppose its conductor in \(B\) is
a nonzero principal ideal \(C(x,y)B\). Then
\[
             S=R[y^2,\;y g(x,y^2)],\qquad C(x,y)=g(x,y^2)
 \tag{30}
\]
after multiplication of \(C\) by a nonzero constant. If \((P,Q)\) is
Keller, \(g(x,0)\in\mathbb C^*\).

To prove this, put \(A_0=R[y^2]\). Every element of \(B\) has a unique
even/odd decomposition \(e(y^2)+y o(y^2)\). Since \(A_0\subset S\),
\[
             S=A_0\oplus yI,\qquad
 I=\{o\in A_0:y o\in S\},
 \tag{31}
\]
where \(I\) is an ideal of \(A_0\). Testing membership of \(f\) and \(yf\)
in \(S\) shows that the conductor is \(IB=I\oplus yI\).

The involution \(y\mapsto-y\) preserves \(IB=CB\). The units of \(B\) are
\(\mathbb C^*\), so \(C(x,-y)=\epsilon C(x,y)\), with \(\epsilon=\pm1\).
If \(C\) were odd, write \(C=y h(x,y^2)\). Its membership in
\(IB=I\oplus yI\) would imply \(h\in I\), hence
\(h\in IB=yhB\), impossible in the polynomial domain. Thus \(C\) is even,
say \(C=g(x,y^2)\). Intersecting \(IB=gB\) with \(A_0\) gives \(I=gA_0\),
which proves (30).

Finally, every element of \(S\) has \(y\)-derivative at \(y=0\) divisible
by \(g(x,0)\). If this polynomial vanished at \(x=\xi\), both \(P_y\) and
\(Q_y\) would vanish at \((\xi,0)\), contradicting a nonzero constant
Jacobian. Therefore \(g(x,0)\) has no zero and is a nonzero constant.
\(\square\)

The condition can equivalently be imposed by asking \(S\) to contain
a polynomial
\[
                    y^2+b(x)y+d(x).
 \tag{32}
\]
The source shear \(Y=y+b(x)/2\), followed by subtraction of an element of
\(R\), puts \(Y^2\) in \(S\). A quadratic element with any nonzero constant
leading coefficient can first be made monic. A nonunit leading coefficient
is not covered by this reduction.

**Corollary D (conditional actual-source conductor bound).** Under the
actual-source hypotheses of the frozen conductor packet, assume in addition
that \(S=\mathbb C[x,P,Q]\) contains a monic quadratic polynomial in \(y\).
If the Keller map is noninvertible, its generic source-conductor degree is
at least eight.

Indeed, the frozen packet supplies a nonzero principal conductor. After (32),
Proposition C identifies the entire source algebra with \(E_g\), not just
its conductor divisor. If the conductor has generic degree \(2d\le6\),
then \(d\le3\), and Theorem B excludes the nontrivial case.
If \(d=0\), the conductor and \(g(x,0)\) are units, so every source-line
algebra is \(\mathbb C[y]\); the cited one-line theorem makes the map
invertible. The even-degree conductor rule in the frozen packet gives the
stated lower bound.

More generally, under the same extra quadratic-element hypothesis,
noninvertibility forces the positive even generic conductor degree \(m\)
to have \(m+1\) composite, by Theorem B and Proposition C.

The extra quadratic-element hypothesis is **not** asserted for arbitrary
Keller maps. In particular this is not an unconditional degree-eight
conductor bound.

## 10. Scope, exact controls, and the next missing classification

- **PROVED from the stated CITED inputs:** the generating-pair obstruction
  for the explicit rings (2) whenever \(2d+1\) is prime, including all
  nontrivial cubic-or-lower-\(g\) models.
- **PROVED recognition:** a principal-conductor source algebra containing a
  monic quadratic normalization coordinate has presentation (2).
- **Not proved:** that an arbitrary actual source algebra with generic
  conductor degree four or six contains such a quadratic element or has
  presentation (2) after permissible coordinate changes.
- **Not proved:** the arbitrary-pair problem \(P,Q\in E_g\) without (4),
  except for the repeated-generic-root obstruction of Section 4.
- **Not inferred:** that arbitrary conductor-cover units are constant;
  (26) shows the opposite. The scalar matching at every singular point,
  forced by global fiber embedding, is the necessary additional information.
- **Not inferred:** any resolution of the complex Jacobian conjecture.

Finite exact checks are supplied separately for coefficient identities,
period derivatives, root collisions, repeated-factor controls, and
representative changes by the defining equation. They do not prove the
all-embedding classical input or classify all actual source algebras.

The next target is to detect, prove, or obstruct a monic quadratic element
inside actual source algebras of generic conductor degree four or six.
Any further criterion must preserve the *entire algebra*, not just the
conductor divisor, and must track specialization. A conductor polynomial
with the same degree or root multiplicities by itself does not supply
the required presentation.
