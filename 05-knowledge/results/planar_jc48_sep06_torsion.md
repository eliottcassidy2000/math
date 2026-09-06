# Vertical component jets and the canonical response connection

**Status: PROVED / independently audited; classical relative-exactness
overlap explicitly retained.** This note gives a
presentation-free version of the vertical response arms in THM-3412, with
an intrinsic connection and exact growth of the pole order of a response.
It is not a new case of JC(2). No literature-priority claim is made.

## 1. Inheritance and classical overlap

Closest proved mechanisms:
[THM-3770 / vertical-principal-part-equalizer-and-log-canonical-dressing-gate](../../01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md)
and [THM-3412 / hamiltonian-principal-part-differential-and-prufer-torsion-arms](../../01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md).
The former describes regularization of one rational mate; the latter
computes a whole response module in a linear-in-one-variable chart.
The hostile is the rationally exact but nonregular primitive of
[THM-3758 / quadratic-radial-carrier-rational-exact-split-fibre-nonentry](../../01-canon/theorems/THM-3758-quadratic-radial-carrier-rational-exact-split-fibre-nonentry.md).
The corrected near miss is generic vanishing of the divergence class
without its integral extension, recorded in MISTAKE-374 and
[THM-3386 / linear-z-canonical-divergence-minimal-polynomial-collision-law](../../01-canon/theorems/THM-3386-linear-z-canonical-divergence-minimal-polynomial-collision-law.md).
The unused operation here is transverse differentiation on the full
response quotient, keeping component labels at coincident target values.

Classical overlap is explicit. Bonnet's
[Relative exactness modulo a polynomial map](https://arxiv.org/html/math/0602223),
Proposition 1.2 and Theorem 1.5, identify relative-exactness torsion and
its vanishing under appropriate fiber hypotheses. The introduction
credits Gavrilov for connected reduced plane fibers and Bonnet--Dimca
for a detailed plane torsion analysis. Those precedents preclude claiming
relative-exactness torsion or the connected-fiber consequence as new.
The direct proof below is self-contained and supplies the component-jet
identification and canonical operator needed for this repo's responses.
The earlier Bonnet--Dimca paper has not yet been compared theorem by theorem.

## 2. Setup and complete torsion statement

Let k be algebraically closed of characteristic zero, R=k[x,y], K=k(x,y),
and let P be a nonconstant polynomial with (P_x,P_y)=R. Put

```
D = P_x partial_y - P_y partial_x,
C_P = R/D(R),               B=k[T], T acts as P.
```

Assume explicitly that ker(D:K->K)=k(P). This holds in the controls below;
no unproved chart-entry assertion is hidden in this hypothesis.
For each c in k, let I_c index the irreducible components of P=c,
let r_c=|I_c|, and put

```
W_c = k^{I_c}/k*(1,...,1),
J_c = u_c^(-1) k[u_c^(-1)],             u_c=T-c.
```

Multiplication by T on J_c is multiplication by c+u_c followed by
removing the nonnegative powers. Equivalently J_c=B[u_c^-1]/B.
Then there is a canonical B-module isomorphism

```
tors_B(C_P)  ~=  direct_sum_c W_c tensor_k J_c.          (1)
```

The summands are finite-support sums. No finiteness of the set of
reducible fibers is needed for (1). In particular

```
dim_k ker((P-c)^n:C_P->C_P) = n*(r_c-1),  n>=1.        (2)
```

Thus every disconnected smooth fiber contributes one full principal-part
arm per component modulo the common diagonal. It contributes neither a
single arm regardless of component count nor a finite module bounded by
one multiplicity.

### 2.1 A polynomial response automatically forbids horizontal poles

Suppose f in K and Df in R. Let g be an irreducible polar factor of f,
and work in the DVR R_(g). Write f=a/g^m with a a unit and m>=1.
If g does not divide Dg, then

```
D(a/g^m) = (g Da - m a Dg)/g^(m+1)
```

has a pole of order m+1: Da is regular and the second numerator term
is a unit. This contradicts Df in R. Hence g divides Dg.
The induced derivation on k(g=0) is nonzero, since otherwise both P_x
and P_y vanish modulo g, contradicting the unit gradient. Its constants
are k: a nonzero derivation of a one-variable function field in
characteristic zero cannot kill a transcendental element. Since D(P)=0,
P has a constant value c on g=0. Thus g divides P-c.
Every pole is vertical. Multiplying f by a suitable polynomial in P
clears all affine poles and gives a polynomial; consequently [Df] is
B-torsion.

### 2.2 Componentwise scalar Laurent coefficients are automatic

A smooth fiber is reduced, and its distinct components are disjoint:
a repeated factor or an intersection point would make the gradient zero.
Thus u=P-c is a uniformizer at each component. If f has pole order m
there, multiply by u^m. The residue alpha of u^m f is killed by the
induced nonzero tangent derivation, since D(u)=0 and Df is polynomial.
Therefore alpha belongs to k. Subtract alpha/u^m and repeat. This yields
unique scalar principal parts

```
pp_{c,i}(f) = sum_{j=1}^{m_i} alpha_{c,i,j} u^(-j).     (3)
```

The recursion is performed separately at every component. It is not
necessary that these coefficients agree across a fiber.

### 2.3 The map and its kernel

For a torsion class [a], choose F(P)a=Dh with F nonzero and h polynomial.
Then f=h/F(P) satisfies Df=a. Associate to [a] the tuple of (3), modulo
the common principal parts at each c. A different polynomial representative
changes f by a polynomial; a different primitive changes f by H(P),
H in k(T), by the constant-field hypothesis. These changes disappear in
W_c tensor J_c, so the map is well-defined and canonical.

If every component tuple is diagonal, subtract the finite sum of those
common principal parts as a rational function of P. The new primitive
has no affine pole, so lies in R, since R is normal. Therefore [a]=0.
This proves injectivity. It also recovers THM-3770's equalizer, now with
no horizontal-pole or scalar-coefficient test left unpaid for an actual
polynomial response.

### 2.4 Surjectivity by simultaneous component jets

Fix c and n. Write P-c=prod_i f_i, absorbing its nonzero scalar in one
factor. Distinct component ideals (f_i) are comaximal, hence so are
(f_i^n). Given arbitrary principal polynomials L_i(u) of orders at most n,
CRT gives h in R with

```
h = u^n L_i(u) mod f_i^n   for every i.                 (4)
```

Set f=h/u^n. Its principal part on component i is exactly L_i.
The ideals (f_i^n) are D-stable because Df_i is divisible by f_i.
The right side of (4) is a polynomial in P and is killed by D; therefore
Dh is divisible by every f_i^n, and Df=Dh/u^n belongs to R.
This realizes every component tuple. Adding such constructions at
finitely many different values realizes the direct sum. This proves (1).
The standard basis u^-1,...,u^-n of the n-torsion of J_c proves (2).

## 3. An intrinsic transverse connection

Choose polynomial A,B_0 with A P_x+B_0 P_y=1, and write
V=A partial_x+B_0 partial_y, m=div(V)=A_x+(B_0)_y.
Then the operator

```
nabla[a] = [V(a)+m a]                                  (5)
```

is well-defined on C_P, is independent of the Bezout choice, and obeys

```
nabla(P[a]) - P nabla[a] = [a].                         (6)
```

Indeed the volume identity or direct differentiation gives
[V,D]=-mD, hence (V+m)Dh=D(Vh). Two Bezout fields differ by hD;
the corresponding difference in (5) is hDa+(Dh)a=D(ha), zero in C_P.
Equation (6) follows from V(P)=1. Thus (5) is an intrinsic algebraic
connection, even when C_P is not finitely generated over k[P].

Under the canonical identification (1), its restriction is exactly

```
nabla(w tensor L(u)) = w tensor dL(u)/du.               (7)
```

To see this, use a rational primitive f. The identity above identifies
the primitive of nabla[Df] with Vf. At a vertical component V(u)=1;
a scalar principal part differentiates ordinarily, while V preserves
the local DVR and cannot create a pole from a regular remainder.
Consequently (3) differentiates coefficientwise, proving (7).
This also proves directly that torsion is preserved by nabla.

For a nonzero c-primary class v of exact order n, the highest pole
coefficient in its W_c-valued principal part is nonzero. Characteristic
zero and (7) imply

```
ord_c(nabla^j v) = n+j,                 j>=0.            (8)
```

The entire derivative sequence is linearly independent over k. In
particular a nonzero finite-dimensional torsion subspace cannot be
nabla-stable. Also no nonzero finitely generated B-submodule of the
torsion can be nabla-stable, since its annihilator would impose a bound
on pole order. This is a precise obstruction to closing a finite jet
packet while discarding its transverse derivative.

## 4. The unit response and the canonical divergence

Put theta=[1], mu=[div V]. Then universally

```
mu=nabla theta,              theta=0 iff mu=0.          (9)
```

Only the converse in the second assertion needs an argument. If m=Dh,
the polynomial 1-form

```
omega=A dy-B_0 dx+h dP
```

is closed because d(A dy-B_0 dx)=m dx wedge dy and
 d(h dP)=-Dh dx wedge dy. The polynomial Poincare lemma gives omega=dQ,
and evaluation on D gives DQ=1. Thus theta=0. The other implication is
immediate from (5), or from a divergence-free Bezout field supplied by Q.

If theta is nonzero torsion with annihilator
prod_c(T-c)^{n_c}, then mu has exact annihilator

```
prod_c(T-c)^{n_c+1}.                                    (10)
```

Both statements concern the fixed P. They do not establish that theta
is torsion or zero for any hypothetical noninvertible Keller component.
For actual Keller P, theta already vanishes; proving its vanishing is
existence of a mate, not polynomial invertibility of that pair.

The converse of 'theta torsion implies mu torsion' is false.
THM-3386's multi-root linear-variable controls have free theta and
nonzero torsion mu. Differentiation can kill a nonzero generic cohomology
class while retaining its vertical integral defect.

## 5. Exact controls and boundary

For P=x+lambda*x^r*y, r>=2 and lambda!=0, use

```
V=(1-lambda*r*x^(r-1)*y) partial_x
  +lambda*r^2*x^(r-2)*y^2 partial_y,
Q0=x^(1-r)/(lambda*(r-1)),               DQ0=1.
```

The zero fiber has two disjoint components x=0 and
1+lambda*x^(r-1)y=0. Q0 has a pole only at the first. Its first scalar
coefficient in u=P is 1/(lambda*(r-1)) at order r-1. Thus (1) gives one
arm, theta has order r-1, and nabla^j theta has order r-1+j.
This recovers the named finite block in THM-3318 while explaining its
unbounded canonical continuation.

For P=x^2+(x^2-1)^2*y, the fiber P=1 has three disjoint components
x=1, x=-1, and 1+(x^2-1)y=0. Its primary torsion therefore has TWO
independent arms, not one arm obtained by merging the two root values.
This agrees with the full THM-3412 calculation. The other fibers are
irreducible: as a linear polynomial in y their coefficient and constant
term are coprime, unless the value is1. A pair of confluent CRT selectors
constructs both arms explicitly in the companion.

For the positive coordinate P=x+y^2, V=partial_x and D=partial_y-2y partial_x;
after the polynomial change (P,y), every polynomial has a D-primitive.
Thus C_P=0 and all fibers have one component, as required.

Smoothness is essential to the proof: for P=xy the zero components meet,
and theta is nonzero (D(y)=y, D(x)=-x, but D has zero constant response),
so there is no polynomial mate. The companion uses the direct fact that
D is diagonal on monomials and 1 is not in its image; the fiber ideals
are not comaximal and the CRT construction fails. No use of (1) is made
for this singular control.

## 6. Connection contract and next question

Source: the polynomial Hamiltonian response complex for a fixed smooth P.
Target: scalar principal parts labelled by actual fiber components,
modulo the common target correction. Map: rational integration, local
principal parts and transgression. Preserved: all torsion responses,
component labels, target value and transverse derivative. Lost: the
free generic de Rham quotient, regular source algebra and global inverse.
Sidecar: generic periods and the actual two-coordinate Keller field.
Cheapest test: equal-value three-component fiber versus one merged arm.

The new usable distinction is between a selected cyclic class and the
full connection module. Finite annihilators of individual responses do
not bound the module closed under differentiation. The next problem is
to connect this transverse operator to the source-normal chart without
replacing its polynomial ring or assuming termination. No such entry
map is claimed here.

## 7. Frozen verification

The [standalone source](../../04-computation/planar_jc48_sep06_torsion.py),
[frozen output](planar_jc48_sep06_torsion.out), and
[independent analytic/source audit](planar_jc48_sep06_torsion_audit.md)
pass all **278** exact gates. Normal and optimized executions reproduce
the same bytes. CRT selectors of orders1..5 supply derivative controls
through order4; four one-arm families use derivative orders0..3.

```
python3 -B 04-computation/planar_jc48_sep06_torsion.py
python3 -B -O 04-computation/planar_jc48_sep06_torsion.py
source 3c01beaeca3abf4567de0a31d2939ae7e1858b0b85746986ccf674f6b336332c
output 4798029319d16097ce37fe1db3e1379255f7fe13e74f414b85bf348a2abfdda3
```
