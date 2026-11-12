# Why this five is a dimension: local isothermal stress classification

**Status: PROVED, self-contained recovery of the classical isothermal-net mechanism.**
The five-parameter statement is not claimed as a discovery. The proofs here
fix the convention and expose its four-plus-one decomposition for transport to
other subjects. The discrete closure and measurement statements below are
applications developed for this session. None proves Chan–Pak.

## 1. Scope, phase and equations

Let U be a connected simply connected open plane region. Work in the
body-force-free, homogeneous isotropic linear-elastic model. We classify
C³ symmetric real stress fields satisfying equilibrium and harmonic trace:

```
σxx,x + σxy,y = 0,      σxy,x + σyy,y = 0,
Δ(σxx+σyy)=0.                                      (1)
```

Harmonic trace is the local elastic compatibility condition in this model.
No boundary traction, material variation, finite deformation, or global
single-valued displacement constraint is silently included. In particular,
this statement is not a theorem about arbitrary prescribed boundary loads.

Put z=x+iy, ∂z=(∂x−i∂y)/2 and ∂bar=(∂x+i∂y)/2. Set

```
p=σxx+σyy,                 q=σxx−σyy−2iσxy.
```

The two real equilibrium equations are exactly

```
∂bar q = −∂z p.                                    (2)
```

Let w=f(z) be a holomorphic isothermal chart with f'(z)≠0 and set
h(z)=f'(z)². The two foliations Re f=constant and Im f=constant are orthogonal.
Their tangent directions are principal stress directions precisely when
q/h is real, at points of distinct principal stresses. This follows because
the principal-axis angle θ satisfies arg q=−2θ and the tangent angle of the
net is −arg f', modulo π/2. We therefore impose q=hR, R real.
A nonvanishing holomorphic h also suffices locally, by taking a square root
and integrating it. Fixing h fixes a normalization as well as the net.

## 2. Exact five-parameter theorem

Every stress satisfying (1) and q/h real has the unique form

```
R(z,bar z)=a|z|²+bz+bar(b)bar(z)+c,
q(z,bar z)=h(z) R(z,bar z),
p(z,bar z)=d−2 Re ∫[h(z)(az+bar(b))] dz,             (3)
```

where a,c,d are real and b is complex. The primitive is normalized at one
fixed base point, so d is the trace there. The solution space has real
dimension five. The entire vector space is included; distinct stress
trajectories are defined only off q=0. On an open patch where R≠0, (3)
gives the prescribed unordered pair of principal foliations.

**Proof.** Harmonic trace says ∂bar∂z p=0. Apply ∂bar to (2):
0=∂bar²q=h∂bar²R, hence ∂bar²R=0. Since R is real, the real and imaginary
parts of this equation are Rxx=Ryy and Rxy=0. The latter makes R locally
F(x)+G(y); the first makes F'' and G'' equal to one constant. Connectedness
propagates the constants, yielding the displayed quadratic R. Equivalently,
its Hessian is the scalar matrix 2aI. Its ∂bar derivative is az+bar(b).
Equation (2) then determines ∂z p=−h(az+bar(b)). A real harmonic p is exactly
the primitive expression in (3), with one free real additive constant.
Conversely (3) directly satisfies (2), harmonic trace, and q/h real.
If a parameter combination gives zero stress, q=0 and h≠0 give R=0,
so a,b,c vanish; then p=d=0. Thus the five parameters are independent. ∎

For b=u+iv, the four amplitude coordinates are
`R=a(x²+y²)+2ux−2vy+c`. This is why a 'set of five' here is not a
classification into five species: it is the dimension of a linear fiber
above an arbitrarily complicated prescribed analytic direction field h.
Infinite geometric complexity is in h, not in these five coefficients.

## 3. Examples and failure boundaries

For the Cartesian net h=1, write two free constants C=(c+d)/2 and
E=(d−c)/2. The complete field is

```
σxx=a y²−2v y+C,
σxy=0,
σyy=−a x²−2u x+E.                                  (4)
```

Thus the five are one shared quadratic curvature, two linear slopes and
two constant stresses. A stress can have arbitrarily large coefficients
without changing these straight-line trajectories.

For the non-Cartesian chart f=z²/2 on a disk avoiding z=0, h=z², and

```
q=z²[a|z|²+bz+bar(b)bar z+c],
p=d−2 Re[a z⁴/4+bar(b)z³/3].                        (5)
```

The same four amplitude modes and pressure anchor survive despite the
curved net. For f=exp(z), h=exp(2z),

```
p=d−2 Re{exp(2z)[a(z/2−1/4)+bar(b)/2]},             (6)
```

up to the prescribed base-point normalization.

Without harmonic trace there is no five-dimensional conclusion: for the
Cartesian net σxx=F(y), σyy=G(x), σxy=0 gives arbitrary functions F,G and
always satisfies equilibrium. Compatibility is the dimension-closing
condition, not the shared direction pattern alone. Removing h≠0 permits
critical points where the chosen foliation chart degenerates. On a
multiply connected domain the primitive can acquire periods; the separate
session global analysis records that extra obstruction.

Adding a constant λI preserves every principal direction and equilibrium;
it shifts p by 2λ and changes no q. This fifth direction is an actual
physical pressure offset, not merely the irrelevant affine gauge of an Airy
potential. Two affine Airy gauges produce identical stresses; pressure
changes the stresses and boundary loads.

## 4. Hermitian carrier, circles and exact measurements

Write ξ=(z,1)^T and

```
H=[[a,bar(b)],[b,c]],       R=ξ* H ξ.
```

The amplitude carrier is Herm(2), a four-real-dimensional vector space.
Writing t=(a+c)/2, s=(a−c)/2, b=u+iv gives

```
det H = t²−s²−u²−v².                               (7)
```

This is a concrete Lorentz-signature quadratic form, not just a shared
word. Congruence H↦M*HM transports the zero conic in homogeneous coordinates,
with det(M*HM)=|det M|²det H. This parameter-space identity does NOT assert
that arbitrary conformal coordinate substitutions preserve the physical
elasticity equations without recomputing (3).

The degeneracy set q=0 is R=0 when h≠0. It is a circle, a point, empty, a
line, or the whole domain, as determined by this real quadratic. Thus the
complexity of principal directions and the simplicity of their isotropic
locus live in different coordinates. Restricting to positive/negative R
selects cones or chambers of the parameter space; it does not create five
combinatorial types.

**Five scalar observations suffice.** Choose four points zj for which the
rows `(xj²+yj²,xj,yj,1)` have rank four. Record R(zj)=q(zj)/h(zj) and one
trace p(z0). These determine a,b,c,d exactly. The rank condition is that
there is no nonzero circle-or-line quadratic of this form vanishing at all
four points. For example `(0,0),(1,0),(0,1),(2,2)` work; four corners of a
square do not. The failure is structural: on a common circle, its defining
quadratic can be added invisibly to all four amplitude readings.

Equivalently one can use the local jet
`(R,Rx,Ry,Rxx,p)` at a single point. Its compatible higher derivatives are
forced: Rxy=0, Ryy=Rxx, every third derivative of R vanishes. A finite jet
rather than a finite list of geometric objects is the actual universal
carrier.

## 5. The real-polyanalytic hierarchy

For any m≥1, a real smooth function satisfying ∂bar^m R=0 on a connected open
complex domain is a polynomial

```
R(z,bar z)=Σ_{0≤i,j<m} H_ij z^i bar(z)^j,
H_ij=conj(H_ji).                                    (8)
```

The real solution dimension is m². One proof uses the local polyanalytic
expansion Σ_{j<m} f_j(z)bar(z)^j. Reality gives ∂z^mR=0; apply successive
∂bar derivatives and descend in j to see each f_j has degree<m. The
Hermitian coefficient symmetry is precisely reality. Local coefficients
agree on overlaps and hence globally on a connected domain. Conversely
every displayed polynomial satisfies the equation.

For m=2, this is the amplitude theorem above. For other m, adding one
formal integration anchor gives m²+1 coordinates, but **no higher-order
elasticity theorem is being asserted**: the actual equilibrium primitive
must still be integrable, and the equation driving it must be specified.
This distinction prevents promoting the numerical sequence2,5,10,17,...
into a universal physical law.

## 6. Exact discrete counterpart: four-mode finite-difference closure

Let R be an array on an integer rectangle `{0,...,N−1}×{0,...,M−1}` with
N,M≥3. Impose, wherever defined,

```
ΔxΔy R=0,                 Δx²R=Δy²R.                (9)
```

Then, uniquely,

```
R(x,y)=A[binom(x,2)+binom(y,2)]+Bx+Cy+D.             (10)
```

Here equality of the second differences is imposed at all common valid
base points `x≤N−3,y≤M−3`. Mixed closure gives R=f(x)+g(y). Equality then
forces every second difference of f and g to the same A. Twice summing
proves (10). For integer R all four parameters are integers, read from
R00, R10, R01, R20. Thus increasing N,M introduces no new free response
coordinate. A missing one of the two constraints destroys this conclusion:
mixed closure alone leaves N+M−1 free coordinates, while radial-curvature
agreement alone retains other modes.

**Exact translation rule.** Under (x,y)↦(x+r,y+s), the coefficient tuple
changes by

```
A'=A,
B'=B+Ar,       C'=C+As,
D'=D+Br+Cs+A[binom(r,2)+binom(s,2)].                 (11)
```

The translations compose because of the binomial addition law. This is a
literal all-size recursive grammar, not extrapolation from small grids.
It gives a testable research program for tournament families: identify two
commuting physical operations, prove their mixed response difference
vanishes and their pure second responses agree, then preserve the absolute
anchor D. A finite difference identity on a scalar count alone does not
prove the operations exist, that the complete root interfaces close, or
that integer target values are hit.

## 7. Source and novelty ledger

- Neményi's survey *Inverse and semi-inverse methods in the mechanics of
  continua*, *Advances in Applied Mechanics*2, pp123–151, states the
  isothermal/isometric five-parameter result and cites the1933 paper
  *Stromlinien und Hauptspannungstrajektorien*, ZAMM13,364–366. The available
  [survey OCR](https://www.scribd.com/document/839714086/Advances-in-Applied-Mechanics-1951-Academic-Press)
  loses displayed formulae. We use it for historical attribution, not as a
  reliable transcription of (3). Bibliographic year differs between the
 1951 volume scan and the2009 reference list's1950; no priority claim here
 depends on settling that discrepancy.
- [Sáenz2009](https://doi.org/10.1007/s00707-008-0118-3) states a related
 2–5 parameter result for more general trajectory nets under isotropic
 elastic equilibrium with conservative body forces. Its publisher page
 flags [an erratum](https://doi.org/10.1007/s00707-010-0370-1). Neither its
 inaccessible full proof nor its erratum is an imported dependency here.
 We do not conflate that result with the explicitly proved body-force-free
 local isothermal theorem.
- Sections2–3 recover the classical mechanism; Sections4–6 give elementary
 consequences, new-to-this-session formulations and a transport template.
 No claim of first historical discovery is made for these standard
 Hermitian/polyanalytic facts.

## 8. A literal five-state flat transport system

There is an exact operational version of the jet description. Let
`Y=(R,Rx,Ry,k,p)` with `k=Rxx=Ryy`, and write h=H+iK. Then

```
∂xY=(Rx, k, 0, 0, −H Rx+K Ry),
∂yY=(Ry, 0, k, 0,  K Rx+H Ry).                      (12)
```

These are linear five-state transport equations with coefficient functions
H,K supplied by the prescribed geometry. Their compatibility is exactly
satisfied by the Cauchy–Riemann identities Hx=Ky and Hy=−Kx. For the last
coordinate the two mixed derivatives are

```
∂y∂x p=−Hy Rx+Ky Ry+Kk,
∂x∂y p= Kx Rx+Hx Ry+Kk,
```

which agree. The other coordinates agree immediately. Thus the local
system is flat and has a five-dimensional solution space; the first four
coordinates form a polynomial jet, while the fifth integrates a geometry-
dependent one-form. This makes the fifth an **extension coordinate**, not
an additional degree of the quadratic amplitude. Global monodromy belongs
to that final coordinate and is handled separately.

This is the precise finite-memory analogy to recursive interfaces:
geometry supplies a changing transport operator, while a closed finite
state carries all subsequent responses. A conserved state dimension is
not a finite set of states, and neither implies that arbitrary prescribed
boundary values lie in the solution image.

**Boundary hostile.** Even for h=1 on a rectangle, requiring zero traction
on its entire boundary kills every field in(4). On a vertical edge the
normal traction gives σxx(y)=0 for all y, forcing a=v=C=0. On a horizontal
edge σyy(x)=0 then forces u=E=0. Thus an exact five-parameter unconstrained
fiber can have only the zero realization after a legitimate target
condition. This is why a tournament response grammar cannot skip the
absolute target/physical compatibility stage.

## 9. Two simple trajectory carriers span every local elastic stress

There is a constructive converse to focusing on one fixed net. Any smooth
body-force-free compatible stress has

```
q=bar(z) A(z)+B(z),       p=d−2 Re ∫A(z)dz,           (13)
```

with holomorphic A,B locally, because ∂bar²q=0. For any real c split

```
q1=A(z)[z+bar(z)+c],       p1=p,
q2=B(z)−(z+c)A(z),        p2=0.                      (14)
```

The first term is exactly the five-parameter construction with h1=A,
R1=z+bar z+c (a=0,b=1); the second is the same construction with h2=q2,
R2=1 and constant trace. Their stresses add to the original stress.
Where h1,h2 and R1 are nonzero, each has an isothermal principal net.
Zeros remain lawful stress values but are degeneracies of the individual
net chart. If A≡0 the original field is already the harmonic component.
At a point where A≠0, c can be chosen to avoid both R1=0 and q2=0 there,
then continuity gives a neighborhood of regular component nets.

This recovers the local algebra behind the superposition attributed to
Wegner1934 in Neményi's survey; it is not asserted as a new historical
result. The powerful operation is **change the carrier before adding**:
a complicated general field is a sum of two fields that separately have
simple trajectory descriptions, but their trajectory nets differ.
This is not a theorem that one fixed five-dimensional space contains every
stress. It also does not transport to tournament Hamiltonian counts until
an actual tournament operation implementing the required addition and
preserving root interfaces is supplied.
