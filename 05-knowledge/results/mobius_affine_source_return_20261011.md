# Changing generators: a complete Möbius test and an affine-source return

2026-10-11. **PROVED, scoped and self-audited:** M1–M3 below over any
characteristic-zero field K; M4 adds the fibre-affine stabilizer and its
elliptic interpretation. **FINITE-EXACT:** independent coefficient tests
in the stated bounded universe. **OPEN:** other changing-generator groups,
the remaining planar source problem, and the planar Jacobian conjecture.

## Inheritance, and the coordinate that changes the answer

The closest proved mechanism is F3 of
[source frames](source_frames_eit_20261011.md): no nonconstant polynomial
in A=K[p,y] returns under a nonzero time of the rational generator u=x²t.
The stronger polynomial-carrier exclusions in
[September7](planar_jc48_sep07_universal_genus.md) and
[September8](planar_jc48_sep08_cusp_ideal.md) remain inherited.
Their hypotheses do not exclude arbitrary products of rational flows.

The hostile is a return at an isolated time even when intermediate
coordinate responses have poles. The corrected near miss would be to
replace the affine source family G_H=-u/2+H by its vector space A: a
return of G_H need not be a return of H. The least-used sidecar is the
transformed offset (Phi(u)-u)/2. This is the planar counterpart of retaining
the endpoint/carry rather than only its multiplier in a Collatz receipt.

Board: **two-generator group / polynomial cone / affine offset / exact
return / marked terminal / observer kernel**. The source-to-target map
here is pullback by a birational symplectic transformation; the predicate
is membership in -u/2+A. Passing only to A destroys its affine origin.
The cheap decisive test is the three-flow word below.

Use the inherited coordinates

\[
 p=t+x^2t^2,quad y=xt^2+x^3t^3,quad
 D=p^3-y^2=t^3(1+x^2t)^2,quad u=x^2t=y^2/D.
\]

They give K(p,y)=K(x,t), with x=py/D,
t=D/p². In particular p,y are algebraically independent. Set v=xt, so
p=t+v² and y=vp. All assertions concern these specific embeddings.

## M1. The two flows generate cotangent lifts of SL2

For the bracket {x,t}=1, the flows of t and u are

\[
 T_b(x,t)=(x+b,t),qquad
 L_c(x,t)=\left(\frac{x}{1+cx},t(1+cx)^2\right).
\]

The second is u-time -c. Finite compositions have the form

\[
 \Phi_M(x,t)=\left(\frac{ax+b}{cx+d},t(cx+d)^2\right),
 \qquad M=\begin{pmatrix}a&b\\c&d\end{pmatrix},quad ad-bc=1.       \tag{1}
\]

Direct differentiation gives Jacobian1; direct substitution gives matrix
composition. Elementary upper and lower matrices generate SL2(K), and
M and -M induce the same map. These are rational maps, not asserted to
be polynomial automorphisms of the affine (x,t)-plane.

## M2. Polynomial intersection for the whole group

For every M other than ±I,

\[
                     A\cap\Phi_M^*(A)=K.                         \tag{2}
\]

Here Phi*H means H composed with Phi. The proof has three steps.

**Leading coefficient cone.** A has the basis p^a y^b D^c with
a,c>=0 and b in {0,1}, because y²=p³-D. Its highest t-degree and
corresponding x-degree are

\[
                 L=2a+3b+5c,qquad e=L-c,qquad 4L/5\le e\le L.
\]

For fixed L, distinct c give distinct e; after c is fixed, parity fixes
b and then a. Thus leading terms cannot cancel. Every nonconstant H in A
has highest t-degree L>=2 with coefficient h(x) satisfying

\[
     \operatorname{ord}_0h\ge\lceil4L/5\rceil,qquad \deg h\le L.   \tag{3}
\]

The highest coefficient after (1) is
(cx+d)^(2L) h((ax+b)/(cx+d)). It is a polynomial and must also obey (3)
if Phi*H belongs to A. If c=0 and b!=0, h must vanish at both0 and b/d,
each with multiplicity at least ceil(4L/5), exceeding its degree.
Thus c=0 forces b=0.

If c!=0, let m be the multiplicity of a/c as a root of h. The transformed
coefficient has degree 2L-m. This follows by expansion at infinity:
(ax+b)/(cx+d)=a/c-1/(c(cx+d)). Degree at most L forces m=L,
h=C(x-a/c)^L. Condition (3) forces a=0. Its transformed coefficient is
C b^L(cx+d)^L, whose required positive order at0 forces d=0.
Consequently a nonconstant return can occur only for a diagonal or
anti-diagonal matrix.

**Diagonal case.** Write q=d², so v is fixed and
p maps to q p+(1-q)v². Inside K[p,v], A has the basis p^i v^j,
0<=j<=i. Give p weight2 and v weight1. A homogeneous component of weight
N>0 is v^N f(p/v²), where

\[
 \operatorname{ord}_0 f\ge\lceil N/3\rceil,qquad
 \deg f\le\lfloor N/2\rfloor.
\]

The diagonal map preserves this grading. Return would force both f(z)
and f(qz+1-q) to satisfy these bounds. If q!=1, the distinct roots0 and
1-q require degree at least 2ceil(N/3)>floor(N/2), impossible. If q=1,
M=±I.

**Anti-diagonal case.** Write M=[[0,b],[-1/b,0]], b!=0. In the (t,v)
chart it sends v to -v and t to v²/(b²t). Write
H=sum_{k=0}^K p^k P_k(v), deg P_k<=k, with P_K nonzero. For K>=1,
Phi*H has a pole of order K along t=0, whose coefficient is
(v²/b²)^K P_K(-v). Lower k cannot cancel it. Hence it cannot lie in
K[t,v], which contains A. This proves (2). QED.

The leading-coefficient argument uses only the highest t-degree; it
also applies when Phi*H differs from a member of A by degree at most1
in t. This observation is needed next.

## M3. Restoring the affine origin produces a genuine return

The actual compatibility equation is

\[
 \Phi_M^*G_H=G_K
 \quad\Longleftrightarrow\quad
 \Phi_M^*H-K=\Delta_M,qquad
 \Delta_M=\tfrac12t\big((a^2-1)x^2+2abx+b^2\big).                 \tag{4}
\]

For H,K in A the following list is exhaustive. Matrices are understood
modulo their common sign.

| Matrix | Permitted polynomial corrections |
|---|---|
| I | any H, with K=H |
| [[1,0],[c,1]], c!=0 | H=K is constant |
| [[0,b],[-1/b,0]], b!=0 | H=K=-b²p/2+C, C in K |
| every other matrix | no solution |

**Proof.** If H is constant, (4) and the absence of any nonconstant A
element of highest t-degree1 force K constant and Delta=0. Thus b=0,
a²=1, giving precisely the identity/lower cases with H=K.

If H is nonconstant, Phi*H and K have equal highest t-degree at least2.
Since Delta has degree at most1, the cone argument in M2 forces M
diagonal or anti-diagonal. For a diagonal matrix, the left side of (4)
lies in K[t,v], whereas Delta=(a²-1)v²/(2t). Hence a²=1, leaving the
identity case.

In the anti-diagonal case Delta=(b²t-v²/t)/2 has pole order1 at t=0.
The pole calculation in M2 forces H to have p-degree at most1:
H=C+h_1 p+h_2 y. Matching the pole gives h_2=0, h_1=-b²/2.
Substitution then gives K=H. QED.

In particular, a positive construction is

\[
 T_b L_{-1/b} T_b=\Phi_{\left(\begin{smallmatrix}0&b\\-1/b&0\end{smallmatrix}\right)},
 \qquad
 \Phi_M^*\left(-\frac{u+b^2p}{2}\right)
              =-\frac{u+b^2p}{2}.                                \tag{5}
\]

For b=1 this is (x,t)->(-1/x,x²t), an involution. It swaps u and t
and fixes ut=v², so -(u+p)/2=-(u+t+ut)/2 is visibly fixed. This is
also an independent explanation of the coefficient calculation.

Thus the two-generator group has isolated nonconstant source stabilizers
even though its nonidentity maps have no nonconstant polynomial return.
For H=h p+C with h!=0, the stabilizer is order2 if -2h is a square in K,
and trivial otherwise. For constant H the lower one-parameter group is
the stabilizer. For every other nonconstant H it is trivial. These are
stabilizers of individual sources; the stabilizer of the whole source
family is the identity. No two distinct H,K are exchanged in this group.

## M4. Beyond the two flows: a Klein four group on elliptic fibres

The b=1 example has invariant I=u+p=x²t²+(1+x²)t. Allow a larger
class of rational symplectic maps that are affine in t:

\[
                  (x,t)\longmapsto(R(x),A(x)t+B(x)),\qquad R' A=1.
\]

The **complete stabilizer of I within this class** is the following
four-element group:

\[
\begin{array}{c|cc}
 &x'&t'\\\hline
1&x&t\\
W&-1/x&x^2t\\
V&1/x&-1-x^2(1+t)\\
Z&-x&-1-t-x^{-2}
\end{array}                                                       \tag{6}
\]

All three nonidentity maps square to the identity, commute, and
WV=Z. They preserve G_{-p/2}=-I/2. V and Z retain an additive fibre
offset and therefore are outside the cotangent SL2 family of M1.

**Exhaustiveness proof.** The coefficient of t² in I(R,At+B)=I(x,t)
gives R² A²=x². Since A=1/R', this forces R'/R=±1/x, hence
R=cx or c/x for a constant c!=0. The constant coefficient gives

\[
                     B(R^2B+1+R^2)=0.
\]

For R=cx, A=1/c. If B=0 the coefficient of t forces c=1; if
B=-1-R^{-2} it forces c=-1. For R=c/x, A=-x²/c; the same two cases
force c=-1 and c=1 respectively. These are exactly (6). Direct
substitution checks every map and its Jacobian1. QED.

There is a geometric reason for this four-element structure. On a fibre
I=kappa, use X=-t and Y=xt(t+1). Then

\[
                      Y^2=X(X-1)(X+\kappa).                       \tag{7}
\]

For kappa!=0,-1 this is a smooth genus-one curve: its double cover of
the X-line has four distinct branch points0,1,-kappa,infinity, so
Riemann–Hurwitz gives genus1. The point at infinity chooses its elliptic
group origin. The discriminant of the cubic is kappa²(kappa+1)²; the
two excluded parameters are precisely the singular-fibre boundary.

Translation by a root point (a,0) of this cubic is

\[
 X'=a+\frac{B_a}{X-a},\qquad
 Y'=-\frac{B_aY}{(X-a)^2},\qquad
 B_a=\prod_{e\ne a}(a-e),quad e\in\{0,1,-\kappa\}.              \tag{8}
\]

This follows directly by intersecting the line through (a,0) and (X,Y)
with the cubic and reflecting its third intersection. Formula (8) for
a=1,0,-kappa gives W,V,Z respectively under the chart above. Thus the
three nonidentity maps are exactly translations by the three nonzero
points of order two. This identification is checked as a rational
identity, not inferred from the number four alone.

The single-flow rational-time exclusions and this positive construction
are compatible: an isolated torsion translation made from different
generators is not a nonzero constant-time formula for the excluded
Hamiltonian family. The fibres, flow parameter, and source correction
are different pieces of data and remain explicit.

Over K=C, Z can also be expressed as a fibre shear following
(x,t)->(-x,-t): the shear t->t-1-x^{-2} is time1 of the rational
Hamiltonian x-1/x. The sign map is a cotangent diagonal map and hence
a word in the two elementary flows of M1. Combining it with W produces
V. This gives an explicit further changing-generator construction; the
finite group itself and its classification in (6) are defined over any
characteristic-zero K.

## What this changes, and the next test

The useful reroute is concrete: replace the homogeneous intersection test
by an affine compatibility test with its exact offset. Equation (5)
survives; claiming an all-source repair from it would fail. The involution
has a pole at x=0 and does not solve planar JC.

M4 carries out the next enlargement: adding a fibre offset enlarges the
stabilizer from one involution to a Klein four group. A further family
must either let x' depend on t, change the source, or leave the rational
fibre-affine class. First identify its transformed G0 offset and where
that offset can cancel the leading poles. M3 and M4 are complete
boundaries for their declared classes, not for every symplectic map.

This operation matches the useful part of the
[Vitali recovery](vitali_atom_recovery_20261011.md): preserve a missing
target-bearing coordinate, then test the full compatibility equation.
It does not identify finite observation loss with nonmeasurability.

## Reproduction and scope

[Program](../../04-computation/experiments/mobius_source_anchor_weights_20261011.py),
[output](mobius_source_anchor_weights_20261011.json). The shared checker
has 16,854 always-active gates. Ten exact rational matrices are tested by
independent linear equations in x,t using all 13 basis elements of
highest t-degree <=9. This finite experiment supports, but does not replace,
the proofs for arbitrary degree. Symbolic composition, Jacobian, four
values of b, and the independent (t,v) cancellation are also checked.
The companion weight and historical controls share the same count.
Normal and optimized Python outputs agree byte for byte.

```text
python3 -B 04-computation/experiments/mobius_source_anchor_weights_20261011.py
python3 -B -O 04-computation/experiments/mobius_source_anchor_weights_20261011.py
```
