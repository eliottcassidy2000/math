# Joint D moments exclude the closed nonpositive circuit octant

**Status: PROVED + complete FINITE-EXACT certificate + INDEPENDENTLY AUDITED.**
In the anchored degree-five B/C/D model, the three Newton-circuit quotients
cannot all be at most one. The exclusion includes all tie boundaries. It
already follows from two ordinary D-Hankel determinants considered together.
The certificate covers a continuous coefficient region and an unbounded
z half-line; it is not a root sample or a coefficient-height census.

Two exact survivors identify the remaining loss: the joint ordinary D packet
through degree six still admits the words \(--+\) and \(+--\). The latter
survivor even has positive simple beta roots and strict C interlacing. Both
are rejected by the native degree-eight D packet. No full-model point with
two negative circuits is constructed here, and their general exclusion is
still OPEN.

## 1. Inheritance and the exact retained object

The closest proved mechanism is the
[degree-eight moment decoder](continuing4_20260906_moments_packet.md).
It retains the native denominator recurrence and is equivalent to the full
weak two-interlacer geometry. The
[fixed-moment circuit theorem](continuing9_20260907_fixed_moment_circuits.md)
shows that the first two moments alone permit every ternary word. The
[C-only sidecar](continuing9_20260907_interlacer_circuits.md)
preserves that flexibility with strict C interlacing. Thus a circuit
restriction must use the missing D predicate.

The corrected near miss is testing only the fourth D determinant. The exact
hostile in Section 5 has three negative circuits, all strict Newton
inequalities, and positive fourth determinant, but its third determinant
is negative. The least-used operation here is to use the third determinant
as a moving coefficient boundary before controlling the fourth as a
concave quadratic in z.

The concept board is: fixed-moment coefficient coordinates; native Gram
matrices; conditional moment completion; curved coefficient domains; and
exact Bernstein sign certificates. The positive move is to combine the
Gram conditions before discarding coordinates. No theorem transfer through
generic real-rootedness, independent marginal signs, or an original-phase
projection is used.

Retain the original coefficient pencil
\[
 B(t)=t^5-13t^4+55t^3-xt^2+yt-z,
\]
\[
 C(t)=t^4-12t^3+45t^2-\frac23xt+\frac37y,\qquad
 D(t)=t^4-11t^3+36t^2-\frac5{12}xt+\frac17y.
 \tag{1}
\]
At a positive-product model point, the three circuit quotients are
\[
 a=C_2=\frac{831875}{8788x},\qquad
 b=C_3=\frac{13x^3}{166375y},\qquad
 c=C_4=\frac{44y^3}{x^3z}.
 \tag{2}
\]
Here C(t) and the subscripted circuit quotients are different objects.
The beta-root sum is 13 and its square sum is 59. An original Laurent phase
or response is not part of this theorem.

Let \(D(t)/B(t)=\sum_{j\ge0}m_jt^{-j-1}\). Formal division gives
\[
\begin{aligned}
 m_0&=1,&m_1&=2,&m_2&=7,\\
 m_3&=7x/12-19,&
 m_4&=115x/12-632-6y/7,\\
 m_5&=199x/2-7171-92y/7+z,\\
 m_6&=7x^2/12+8969x/12-58463-915y/7+15z.
\end{aligned}
 \tag{3}
\]
Write
\[
 J_3=\det(m_{i+j})_{i,j=0}^2,\qquad
 J_4=\det(m_{i+j})_{i,j=0}^3.
 \tag{4}
\]
Weak D interlacing of a real-rooted B gives nonnegative residues after
cancellation, hence positive-semidefinite moment matrices and
\(J_3,J_4\ge0\). Repeated roots or canceled factors do not invalidate this
necessary direction. We will prove a stronger algebraic exclusion using
only the two determinant inequalities, without assuming B or C geometry.

## 2. The theorem and the bounded coefficient chart

Put
\[
 x_0=\frac{831875}{8788},\quad
 \ell(x)=\frac{13x^3}{166375},\quad
 h(x)=\frac{-343x^2+67788x-3157056}{2592},\quad
 k(x,y)=\frac{44y^3}{x^3}.
 \tag{5}
\]

**Theorem.** There are no real \(x,y,z\) satisfying simultaneously
\[
 x\ge x_0,\qquad y\ge\ell(x),\qquad z\ge k(x,y),
 \qquad J_3\ge0,\quad J_4\ge0.
 \tag{6}
\]
The first three inequalities themselves force \(x,y,z>0\), so their circuit
meaning is exactly \(a\le1,b\le1,c\le1\), with no denominator ambiguity.
Consequently every full B/C/D model point with positive product has at least
one circuit quotient strictly greater than one.

**First reduction.** The exact third determinant is
\[
 J_3=-\frac{343x^2-67788x+2592y+3157056}{1008}
     =-\frac{18}{7}(y-h(x)).
 \tag{7}
\]
Thus \(y\le h(x)\). Combining it with \(y\ge\ell(x)\) gives
\[
 g(x):=343x^2-67788x+3157056+2592\ell(x)\le0.
 \tag{8}
\]
But
\[
 g(99)=537759/125>0,
\]
\[
 g'(99+u)=\frac{101088}{166375}u^2+
          \frac{12195334}{15125}u+\frac{8361378}{1375}>0
 \quad(u\ge0).
 \tag{9}
\]
Therefore \(x<99\). Every point in (6) has
\[
 x_0\le x<99,\qquad \ell(x)\le y\le h(x),\qquad z\ge k(x,y).
 \tag{10}
\]
These are necessary conditions, not a claim that every point of the chart
comes from beta-root geometry.

Use the polynomial map of the closed unit square
\[
 x=x_0+(99-x_0)u,\qquad
 y=\ell(x)+(h(x)-\ell(x))v,\qquad 0\le u,v\le1.
 \tag{11}
\]
It covers every feasible pair in (10). If \(h(x)=\ell(x)\), choose \(v=0\);
no division by the fiber length is needed. The map also includes irrelevant
pairs where \(h(x)<\ell(x)\), but the sign certificate below holds on the
entire square. Enlarging to \(x=99\) is harmless.

## 3. Two complete positive-coefficient identities

Direct calculation from (3) gives
\[
 \frac{\partial^2J_4}{\partial z^2}=-6.
 \tag{12}
\]
Define integer polynomials by clearing **positive** factors:
\[
 P_0(x,y)=-7112448\,x^6J_4(x,y,k(x,y)),
\]
\[
 P_1(x,y)=-1008\,x^3
            \frac{\partial J_4}{\partial z}(x,y,k(x,y)).
 \tag{13}
\]
The source verifies both identities exactly from the full native determinant.
For clarity, the smaller polynomial is
\[
 P_1=-4753x^5+1008x^4y+818748x^4
       -97632x^3y-34727616x^3+266112y^3.
 \tag{14}
\]
The complete 15-term \(P_0\), the moments, and the full \(J_4\) are retained
in the standalone JSON certificate. No terms in the source pencil are
removed when forming these determinants.

After substitution (11), the degrees of \(P_0,P_1\) in \((u,v)\) are
\((18,6)\) and \((9,3)\). If
\(F(u,v)=\sum f_{ij}u^iv^j\) has degree bounds \((d,e)\), its tensor Bernstein
coefficients are exactly
\[
 b_{rs}=\sum_{i\le r,\ j\le s}
    f_{ij}\frac{\binom ri}{\binom di}\frac{\binom sj}{\binom ej}.
 \tag{15}
\]
They give the identity
\[
 F(u,v)=\sum_{r=0}^d\sum_{s=0}^e b_{rs}
 \binom dr u^r(1-u)^{d-r}
 \binom es v^s(1-v)^{e-s}.
 \tag{16}
\]
The basis functions are nonnegative on the closed square and sum to one.

The complete exact arrays have the following positive minima:

| Polynomial after (11) | Bidegree | Number of coefficients | Minimum |
|---|---|---:|---|
| \(P_0\) | \((18,6)\) | 133 | \(38259412467502741144816725/262144\) |
| \(P_1\) | \((9,3)\) | 40 | \(4049307478755/256\) |

All 173 entries are included, not only the minimum or a selected subset.
The compiler also reconstructs each polynomial at a full degree-bounded
tensor interpolation grid, checking the complete identity by a second
finite path. Thus (16) proves \(P_0,P_1>0\) on the whole square, including
every edge and corner.

**Completion of the proof.** Since \(x>0\), equations (13) give
\[
 J_4(x,y,k(x,y))<0,\qquad
 \partial_zJ_4(x,y,k(x,y))<0
 \tag{17}
\]
throughout the feasible coefficient chart. Equation (12) makes the derivative
strictly decrease as z increases. More explicitly, for \(w=z-k(x,y)\ge0\),
\[
 J_4(x,y,z)=J_4(x,y,k)+w\,\partial_zJ_4(x,y,k)-3w^2<0.
 \tag{18}
\]
This contradicts (6). It covers the complete unbounded z half-line and does
not depend on an additional product bound, simple-root assumption, or a
strict determinant hypothesis. The closed nonpositive octant is excluded.
\(\square\)

If a weak model has \(z=0\), the finite quotient \(C_4\) is not defined.
There is no claim to define it artificially. Such a point cannot satisfy
the cleared inequalities (6), since their first two force \(y>0\).

## 4. Why the two determinant conditions cannot be separated

The inherited fourth-determinant hostile is
\[
 x=\frac{4159375}{38207},\quad
 y=\frac{4325100097656250}{42626885365649},\quad
 z=\frac{494717113513946533203125000000}
          {13159366797968048415097695043}.
 \tag{19}
\]
Its circuits are
\[
 \left(\frac{2939}{3380},\frac{29201}{29390},
                    \frac{276701}{292010}\right),
\]
all strictly below one; all four Newton ratios are strictly above one.
Nevertheless \(J_4>0\), while
\[
 J_3=-4492863467194616219/42967900448574192<0.
 \tag{20}
\]
This refutes a universal fourth-determinant sign assertion on the negative
circuit/Newton relaxation. It does not refute the model theorem.

Conversely, set \(r=2047/2048\) and take the exact fixed-moment chart with
all three circuit quotients equal to r:
\[
 x=x_0/r,\quad y=y_0/r^4,\quad z=z_0/r^{10},
\]
where
\[
 y_0=3460080078125/52206766144,\quad
 z_0=791547381622314453125/52414354446428935168.
 \tag{21}
\]
This is the frozen strict C-only negative-word witness from continuing9.
All five C-Hankel leading minors are positive, \(J_3>0\), and \(J_4<0\).
The all-one center has the same D-determinant sign split, so ties do not
escape the argument.

## 5. Precise boundary of a two-negative-pair extension

Joint D positivity through degree six does not already exclude every pair
of nonpositive circuits. Two small exact survivors are:

| \((x,y,z)\) | Strict circuit word | \(J_3\) | \(J_4\) |
|---|---|---|---|
| \((95,68,1)\) | \(--+\) | \(10973/1008\) | \(5415309223/7112448\) |
| \((86,50,9)\) | \(+--\) | \(1571/252\) | \(7483825/444528\) |

The leading D minors begin \(1,3,J_3,J_4\), so the whole ordinary
4-by-4 D-Hankel is positive definite at both points. All four Newton ratios
are strictly above one. Their fifth D determinants are respectively
\[
 -3637138316701405/4182119424,\qquad
 -60271634953/43563744,
 \tag{22}
\]
so the full D predicate rejects them.

The second point is stronger than a coefficient-only hostile: it has five
distinct positive beta roots and strict C interlacing. Its five C leading
minors are
\[
 1,\ 2,\ 173/63,\ 282706/27783,\ 250972/50421,
 \tag{23}
\]
all positive. Thus even actual B/C geometry plus the entire ordinary D
packet through degree six permits \(+--\); the degree-eight D information
is necessary for any exclusion of this word by this route. No claim is
made that the first point retains B/C geometry.

The report also rechecks the four earlier strict full-model controls with
words \(+++,++-,+-+,-++\). These establish nonvacuity of the native joint
positivity domain and preserve its known positive signals. They do not
classify every remaining word. In particular the general exclusion of all
two-negative full-model points remains OPEN.

## 6. Connection contract and reproduction

The source is the exact B/C/D coefficient pencil. The target is a necessary
pair of D moment inequalities. The map is formal division, then a
coefficient-dependent domain reduction, followed by positive-denominator
normalization and the tensor Bernstein basis. It preserves the native
x,y,z coefficients and all z dependence of the fourth determinant.
It discards higher D moments and B/C geometry as a **relaxation**, so proving
infeasibility still applies to the full model. The two-negative survivors
exhibit the information lost in that relaxation. No original phase, Laurent
carry evaluation, or response-sign conclusion is asserted.

Run the [standalone producer](../../04-computation/continuing10_20260907_nonpositive_circuits.py):

~~~text
python 04-computation/continuing10_20260907_nonpositive_circuits.py
python -O 04-computation/continuing10_20260907_nonpositive_circuits.py
~~~

The [complete certificate](continuing10_20260907_nonpositive_circuits_certificate.json)
contains the native moments and determinants, both primitive integer sign
polynomials and their positive normalization factors, every power and
Bernstein coefficient, all domain bounds, and all positive/hostile controls.
No producer is imported. Both
[normal output](continuing10_20260907_nonpositive_circuits.out) and
[optimized output](continuing10_20260907_nonpositive_circuits_optimized.out)
pass **384 always-active exact gates** with identical raw LF stdout.
The source configures newlines explicitly; output bytes were captured
without text normalization.

Frozen pins:

    source d392980046171db2ea1df6cb759c04882fd57305ae811fe3d7f5746df6b6a973
    output 36e5c21f2113d4a3dfbc825af03ddb0d3570a3e359c58ac5a3a78f164db9908b
    JSON   be4762cf527f1add2887ee2b29a98440a83f0c68454571b54d182b7395b58092

The parent owns independent-audit promotion and filing. The general
two-negative circuit exclusion and the original-phase Laurent sign problem
remain open.

The [independent audit](continuing10_20260907_nonpositive_circuits_audit.md) accepts the scoped proof.
Repository filing changes only routing and status prose; frozen source,
output and certificate bytes match the independently reviewed originals.
