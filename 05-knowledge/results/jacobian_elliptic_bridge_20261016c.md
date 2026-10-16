# A marked elliptic bridge to the Lyness five-cycle

2026-10-16. **PROVED:** the explicit curve maps, finite moduli matching,
degree-twelve golden comparison, and the specified differential obstruction.
**FINITE-EXACT:** the source-bound controls in the companion checker.
**CITED / CLASSICAL:** Lyness integrability and its elliptic translation
interpretation. **OPEN:** planar JC, source-normal chart entry, and any
all-row polynomial lift. No new Keller pair is supplied.

## 1. Inheritance and the precise target

The closest actual source is
[THM-4385, the row-ten elliptic sign quotient](../../01-canon/theorems/THM-4385-source-normal-row-ten-elliptic-sign-quotient.md),
corrected in scope by
[THM-4389, the non-isotrivial weight-thirteen pencil](../../01-canon/theorems/THM-4389-source-normal-weight-thirteen-row-ten-nonisotrivial-elliptic-pencil.md).
The latter retains a parameter which the fixed-curve description loses.
Both are finite-row carriers, not fibres of a supplied polynomial Keller map.

The current inverse-clock result is
[the source-preserving Keller cokernel operator](keller_inverse_clock_20261015f.md).
It does not turn factorially small representatives into zero classes.
The closest differential exclusion is
[THM-3981, centered-cusp nonexact quadrature](../../01-canon/theorems/THM-3981-centered-cusp-quadrature-genus-two-transcendence.md):
finite algebraic covers do not repair that nonexact time form. The positive
source-polynomial alternative remains
[the finite dual-closure reconstruction criterion](jc_bilinear_reconstruction_20261013a.md).

The live concepts here are a marked source point, its elliptic modulus, the
chosen torsion translation, the differential, and the omitted boundary.
The main positive operation is an exact finite moduli compiler. The hostile
is that a golden fixed point can lie on a singular fibre, while a different
golden presentation is smooth. The least-used sidecar is the *scale and pole
divisor of the differential*, rather than another value of the invariant.

For the classical input, see Andrew Hone,
[Efficient ECM factorization in parallel with the Lyness map](https://arxiv.org/html/2002.03811),
§§1–2, especially equations (2)–(5), (15), and the five-torsion exception in
§3. His general map is $(x,y)\mapsto(y,(ay+b)/x)$; we specialize to
$a=b=1$. His conserved quantity $K$ is our $h-3$. The five-cycle, logarithmic
symplectic form, and elliptic-translation interpretation are classical.
The formulas below are independently derived and checked for this marked
specialization; no claim of historical novelty is made for them.

## 2. Exact dynamics, including the missing boundary

Set

\[
 L(x,y)=\left(y,\frac{1+y}{x}\right),\qquad
 F=(x+1)(y+1)(x+y+1),\qquad H=\frac{F}{xy}.
\]

The successive scalar coordinates are

\[
 x,\ y,\ \frac{1+y}{x},\ \frac{1+x+y}{xy},\ \frac{1+x}{y},\ x,\ y.
 \tag{1}
\]

Thus $L^5$ is the identity as a rational map. A literal five-step evaluation
in this affine chart is defined **if and only if**

\[
                  xy(x+1)(y+1)(x+y+1)\ne0.                \tag{2}
\]

This guard cannot be removed by cancelling the final expression. At
$(-1,2)$ the first five states are
$(-1,2),(2,-3),(-3,-1),(-1,0),(0,-1)$; the fifth application requires $0/0$.
An extension to a smooth projective fibre is a different, explicitly
completed interface.

Direct substitution gives

\[
 H\circ L=H,\qquad \det DL=\frac{1+y}{x^2},\qquad
 L^*\Omega=\Omega,\qquad \Omega=\frac{dx\wedge dy}{xy}.
 \tag{3}
\]

The product of the five actual determinants is one on (2), but the individual
determinant is not constant; at $(1,2)$ it is three. These are rational,
log-symplectic dynamics, not a polynomial Keller map.

## 3. A literal elliptic curve map and its differential scale

The projective completion of $C_h:F-hxy=0$ is smooth exactly when

\[
                    h\ne0,\qquad h^2-11h-1\ne0.          \tag{4}
\]

Here is an exact derivation, retaining all coordinate choices. Put

\[
\begin{split}
 B_h(x)&=x^2+(3-h)x+2,\\
 v&=2(x+1)y+B_h(x),\\
 q_h(x)&=x^4+(2-2h)x^3+(h^2-6h+1)x^2-4hx.
\end{split}
\]

Then $v^2-q_h(x)=4(x+1)C_h$. This gives the inverse
$y=(v-B_h(x))/(2(x+1))$ when $x\ne-1$. For $h\ne0$, the only affine point
lost on $x=-1$ is $(-1,0)$, where $\partial_y C_h=h\ne0$.
All three points at infinity $[1:0:0],[0:1:0],[1:-1:0]$ are smooth.
Consequently the quartic discriminant detects all singular parameters.

For a binary quartic $(a,b,c,d,e)$ in descending coefficient order, use

\[
 I=12ae-3bd+c^2,\quad
 J=72ace+9bcd-27ad^2-27b^2e-2c^3.
\]

In this family these are

\[
\begin{split}
 I(h)&=h^4-12h^3+14h^2+12h+1,\\
 J(h)&=-2(h^2+1)(h^4-18h^3+74h^2+18h+1),\\
 4I^3-J^2&=6912h^5(h^2-11h-1).
\end{split}                                                     \tag{5}
\]

For example, $h=0$ gives three lines. The other two excluded values have
the singular fixed points $(\phi,\phi)$ and $(\phi',\phi')$, where
$\phi^2=\phi+1$ and $\phi'=1-\phi$:

\[
             H(\phi,\phi)=3+5\phi=\phi^5
                        =\frac{11+5\sqrt5}{2}.            \tag{6}
\]

So the positive golden fixed point is **not** a point on a smooth member
of this elliptic pencil.

There is also a completely literal short-Weierstrass map. Write
$c=h^2-6h+1$ and $d=-4h$. On $x\ne0$, put

\[
             W=\frac d x+\frac c3,\qquad Y=\frac{dv}{x^2}.
\]

Then

\[
       Y^2=W^3-\frac I3 W-\frac J{27},\qquad
       \frac{dW}{2Y}=-\frac12\frac{dx}{v}.                 \tag{7}
\]

The inverse is $x=d/(W-c/3)$, $v=Yx^2/d$, followed by the preceding
$y$-formula, with its stated denominators retained. Equation (7) identifies
the point and differential, not merely the modulus.

On a smooth compact fibre,

\[
                 \tau_h=\frac{dx}{\partial_y C_h}
                         =\frac{dx}{v}                   \tag{8}
\]

is a nonzero holomorphic differential. In the ambient plane, substituting
$h=H(x,y)$ gives $\tau\wedge dH=\Omega$. Thus (3) implies fibrewise
preservation of $\tau_h$. The order-five automorphism is nonidentity, since
its affine fixed-point equations have only the two points in (6), on
singular fibres. An automorphism of an elliptic curve preserving its
nonzero differential is a translation: after subtracting its value at the
origin, the remaining group automorphism has derivative one and is the
identity in characteristic zero. Hence $L$ acts on each smooth fibre as
translation by a nonzero point of order five, as in the classical description.

## 4. A finite moduli compiler and its exact golden cost

Equations (5) and (7) give

\[
 j(h)=\frac{I(h)^3}{h^5(h^2-11h-1)}.
\]

Put $b=h-h^{-1}-11$. Then

\[
       j(h)=\frac{(b^2+10b+5)^3}{b},\qquad
       h^2-(b+11)h-1=0.                                  \tag{9}
\]

**Finite matching theorem.** For any finite $j_0$ in a characteristic-zero
field, the equations

\[
        (b^2+10b+5)^3-j_0b=0,\qquad
        h^2-(b+11)h-1=0                                  \tag{10}
\]

produce a smooth Lyness moduli value after a field extension of degree at
most $6\cdot2=12$. Indeed, the first polynomial is monic and has nonzero
constant term $125$, so $b\ne0$; the second has root product $-1$, so
$h\ne0$, and $h^2-11h-1=bh\ne0$. Equivalently the monic degree-twelve
polynomial is

\[
             I(h)^3-j_0h^5(h^2-11h-1)=0.                 \tag{11}
\]

The resultant of $I(h)$ and $h^5(h^2-11h-1)$ is $25$, another exact check
that no singular root was silently admitted. The paired parameters are
$h$ and $-1/h$. Forgetting this choice loses a parameter marking. Their
root product is $-1$; it is a field norm when the quadratic is irreducible.
At $b=0$, the same quadratic displays the golden unit $\phi^5$, precisely
on the excluded singular boundary.

This is a theorem about **matching moduli**. A specified elliptic curve
with that $j$ becomes isomorphic to this one over an algebraic closure,
and over some further finite extension. It need not already be isomorphic
over the field containing $h$: a twist and its chosen isomorphism still
matter. No actual point, source puncture, chosen origin, or differential
scale is reconstructed from $j_0$ alone.

The comparison $j_0=2048$ is especially sharp. In Legendre coordinates,

\[
\begin{split}
 &(1-\lambda+\lambda^2)^3-8\lambda^2(1-\lambda)^2\\
 &\qquad=(\lambda^2-\lambda-1)(\lambda^2-3\lambda+1)
                    (\lambda^2+\lambda-1).               \tag{12}
\end{split}
\]

Thus a golden quadratic parameter gives a smooth Legendre curve with this
rational $j$. In contrast, its polynomial (10) is

\[
 b^6+30b^5+315b^4+1300b^3+1575b^2-1298b+125,              \tag{13}
\]

and (11) is

\[
\begin{split}
 h^{12}-36h^{11}+474h^{10}-2700h^9+5775h^8-2984h^7
 +14124h^6\\
 \qquad+2984h^5+5775h^4+2700h^3+474h^2+36h+1.             \tag{14}
\end{split}
\]

Both are irreducible over $\mathbf Q$: modulo three they become
$b^6+b^3+b+2$ and $h^{12}+h^7+2h^5+1$. The JSON contains exact Frobenius
remainders and gcd certificates. For degree $n=6,12$, it verifies
$X^{3^n}=X\pmod f$ and
$\gcd(f,X^{3^{n/\ell}}-X)=1$ for $\ell=2,3$; these are the complete
finite-field irreducibility conditions. Consequently **every** Lyness
$h$ for $j=2048$ has degree exactly twelve over $\mathbf Q$, not merely a
possible degree-twelve upper bound.

This is a concrete cost of the different markings: the Legendre
presentation selects two-torsion data, while the five-cycle selects a
five-torsion translation. It is not an increase in the geometric genus or
in the complexity of the underlying rational $j$ itself.

### A finite Frobenius explanation of the number twelve

There is a separate arithmetic reason to expect this marking cost. The
rational golden model $E_0:Y^2=X^3+X^2-3X-2$ has good reduction at three.
Its five affine points plus infinity give $\#E_0(\mathbf F_3)=6$, hence
Frobenius trace $-2$. Use the classical characteristic equation on prime-to-
characteristic torsion; see
[Sutherland, 18.783 Lecture 8, slide 5](https://math.mit.edu/classes/18.783/2023/LectureSlides8.pdf).
On $E_0[5]$ the polynomial is $t^2+2t+3$, irreducible over $\mathbf F_5$,
so a basis identifies Frobenius with

\[
 M=\begin{pmatrix}0&2\\1&3\end{pmatrix},\qquad
 M^6=3I,\quad M^{12}=-I,\quad M^{24}=I.                   \tag{14a}
\]

**PROVED finite action.** Every nonzero vector has orbit length 24,
every pair $\{v,-v\}$ has orbit length 12, and every projective line has
orbit length 6. Indeed the irreducible algebra $\mathbf F_5[M]$ is
$\mathbf F_{25}$. The displayed element has order 24: none of the proper
divisors of 24 gives the identity (the checker checks every power).
Multiplication on its nonzero elements is therefore a single 24-cycle.
Its first scalar power is 6; its first power equal to either sign of the
identity is 12. The script independently checks all 24 vectors.

There is also an exact characteristic-zero consequence for this fixed
rational model. At the good prime three, its prime-to-three five-torsion
is finite étale, so good-reduction specialization identifies the geometric
torsion and realizes the displayed Frobenius as a Galois action. Its
24-cycle is contained in the global Galois orbit of any nonzero five-point;
there are only 24 such points. Thus each such point has degree **24** over
$\mathbf Q$, and its $X$ coordinate has degree **12**, since the two points
over an $X$ coordinate are exactly the pair $\{P,-P\}$. This uses the
classical good-reduction torsion specialization in addition to (14a); it
does not identify an $h$ coordinate with an $X$ coordinate without a map.

This is a geometric interpretation of the sextic and degree-twelve counts:
a cyclic five-subgroup, and then a five-point modulo sign, retain different
marks. The classical Lyness translation supplies the relevant five-torsion
setting. Identifying these modular marks is **not** used as an alternative
proof that (13)–(14) are irreducible; the direct polynomial certificates
already prove that claim. An actual torsion point on the fixed rational
model can require the further sign choice, consistent with the retained
twist/isomorphism obligation.

## 5. Applying the compiler to the actual row-ten pencil

Retain the source parameter $a$ of THM-4389 and its exact quartic

\[
                    q_a(r)=4A(r)B(r)+c_0^2a^2,
                    \qquad c_0=707284001250,              \tag{15}
\]

where

\[
\begin{split}
 A(r)&=353642000625r^3+6793915500000r^2
             +50541940696500r+7231154026500,\\
 B(r)&=91584545734393856r+631918028977864704.
\end{split}
\]

Compute its ordered quartic invariants $I_a,J_a$ and
$\Delta_a=4I_a^3-J_a^2$. For every smooth source parameter, an exact
matching condition is

\[
 \Delta_a I(h)^3-6912I_a^3h^5(h^2-11h-1)=0.               \tag{16}
\]

Over the source field this is (11) after division by $\Delta_a$.
The API exports the condition and its actual $j(a)$, rather than pretending
to have selected an algebraic root or a curve isomorphism. The canonical
quartic is authenticated from its coefficients; the controls include
$a=0$ and $a=1$, whose $j$-values differ. This does not undo the
non-isotriviality correction in THM-4389.

| Source | Target and map | Preserved | Missing data / cheapest hostile |
|---|---|---|---|
| Literal Lyness point | Five actual steps on (2) | Point, $h$, marked log form | $(-1,2)$ reaches a forbidden denominator |
| Marked cubic point | Equation (7) and its inverse | Actual curve point, differential scale | $x=-1$ requires another chart |
| Smooth canonical quartic | Equation (16) | Geometric $j$ | Root, twist/isomorphism, punctures, differential; $a=0,1$ differ |
| Golden Legendre modulus | Equations (13)–(14) | $j=2048$ | Degree-two parameter does not supply degree-twelve $h$ |

In particular, changing to a five-torsion presentation does not enter the
source-normal chart, solve a higher-row compatibility equation, or produce
a polynomial Keller pair. Those are still separate obligations.

## 6. The differential cannot be flattened by a finite algebraic repair

**Specified-form obstruction.** There are no rational functions $U,V$ with

\[
                    dU\wedge dV=c\Omega\quad(c\ne0),      \tag{17}
\]

where $c$ is a nonzero constant in the characteristic-zero ground field.
This remains impossible after any finite algebraic extension of the
two-variable rational function field over that ground field.

For the rational case, $dU\wedge dV=d(U\,dV)$ is exact. The iterated residue,
first at $x=0$ and then at $y=0$, of $\Omega$ is one. For an arbitrary
rational one-form $A\,dx+B\,dy$, that residue of
$(\partial_xB-\partial_yA)dx\wedge dy$ is zero: an $x$ derivative has no
$x^{-1}$ coefficient, and the remaining coefficient is a $y$ derivative,
which has no $y^{-1}$ coefficient. Rational Laurent expansions justify the
calculation even with mixed denominators. Thus $\Omega$ is not exact.
For a finite extension, trace on Kähler differentials commutes with $d$;
tracing a hypothetical primitive would make the nonzero multiple
$[L:K]\Omega$ exact over the original field, a contradiction.

The fibre version is equally precise. For
$V=xy(H_y\partial_x-H_x\partial_y)$, one has
$\tau_h(V)=1$ on a smooth fibre. There is no rational time coordinate
$t$ with $dt=\tau_h$, even after a finite algebraic cover: the pullback is
a nonzero holomorphic differential on the smooth projective covering
curve. A rational function whose differential is holomorphic has no poles
and is constant. This applies the inherited finite-cover nonexactness
mechanism to an explicit new source; it does not reopen the closed route.

The quantifier in (17) is essential. It excludes flattening **this marked
logarithmic form**. It does not exclude every rational conjugacy of the
finite-order map, or every other invariant form. Analytic logarithms
$U=\log x,V=\log y$ flatten the form locally, but leave the rational and
polynomial categories. Multiplying the form by $xy$ makes it the exact
Cartesian form, but then $L$ fails to preserve it by (3).

The usable remaining JC interface is therefore concrete: a proposed
carrier-to-polynomial construction must supply the actual source map,
the puncture divisor, and the pulled-back differential. A nonzero
holomorphic time form or this nonzero iterated residue rejects the
specified finite algebraic flattening. Passing that test would still not
supply a Keller pair; the full polynomial equations or source-polynomial
dual closure must then be proved. No universal JC conclusion follows.

## 7. Exact reproduction and boundaries

The companion [script](../../04-computation/experiments/jacobian_elliptic_bridge_20261016c.py)
prints deterministic [JSON](jacobian_elliptic_bridge_20261016c.json).
Run from the repository root:

```text
python -B 04-computation/experiments/jacobian_elliptic_bridge_20261016c.py
python -B -O 04-computation/experiments/jacobian_elliptic_bridge_20261016c.py
```

The finite universe uses all 33 distinct rational coordinates $n/d$ with
$-7\le n\le7$, $1\le d\le3$: 934 guarded source points and 155 rejected
full guards. It checks all five literal steps, source fibre, determinant
product, point-level Weierstrass inverse, and both moduli equations.
Symbolic controls check the full polynomial identities, the golden
singularity, differential scale, exact-form residue examples, and the
Legendre factorization. Independent integer polynomial arithmetic modulo
three certifies both irreducibility claims. Actual canonical coefficients,
malformed source/type/parameter controls, and the cancellation hostile
are retained. This is not a census of arbitrary polynomial maps or a
verification of the whole incoming Jacobian theory.

Normal and optimized runs agree on **6,879 checks**. There are 19 distinct
canonical parameter controls. The LF-normalized JSON SHA256 is
`5bcbfae4da591b7b33c555ee7607ea0aca53fc2813d2de7b7b04b2767297da14`.
An independent peer read the complete proof/API and reproduced the optimized
initial output exactly. Its only requested clarification was to make the
constant coefficient in (17) explicit. The final 56-check addition is the
separately stated Frobenius action, including all 24 nonzero vectors; it
does not change the original irreducibility proof or any production API.
