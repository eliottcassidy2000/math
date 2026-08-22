---
id: THM-3333
title: "Gaussian-square Farey/Pythagorean light-cone carrier, Kelvin gauge, triangular in/exradius polarization, and two-scale graph law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Gaussian
  squaring sends primitive projective lattice rays faithfully to canonical
  unreduced Pythagorean null points, where Lorentz pairing is twice squared
  determinant.  Farey adjacency is therefore exactly pairing two on this
  image.  Triangular polarization is the same Lorentz form; on the positive
  right-triangle chamber its four signed null-point potentials are the
  inradius and three signed exradii.  For an unscaled ordered Euclid lift the
  parameter graph pair has rank r; for every positive integral right triple
  the triangle-side pair has rank 2r^2.  Clique filling kills all relative
  homology.  Kelvin-square rescaling turns
  the proved THM-2056 determinant gate into a constant Lorentz cap.  None of
  these identities supplies the missing LRC owner/phase/global-exit data.
source: death-star-2026-08-08-arithmetic-triangular-topology
audit: >
  Independent arithmetic and topology reviews reproduced the Lorentz and
  face-determinant constants on separate full-circle and bounded integer
  sweeps, checked the parity and quotient hostiles, rebuilt both relative
  graph calculations, and verified that flag completion annihilates them.
  Normal and optimized executions byte-match the stored transcript.
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
related:
  - THM-401-pair-sum-sieve-modulus-is-2n-minus-1
  - THM-2000-support-harmonic-abel-dini-figurate-surface
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2975-modular-six-sheet-schreier-graphs-and-farey-partial-cube-boundary
  - THM-3145-bass-serre-two-three-tree-and-tetrahedral-congruence-quotient
  - THM-3157-pointed-resolvent-c3-lift-edge-hexagon-and-relative-tournament
script: 04-computation/gaussian_farey_pythagorean_triangular_thm3333.py
output: 05-knowledge/results/gaussian_farey_pythagorean_triangular_thm3333.out
script_sha256: ddeb800a881fdf3576d3d03d269ced2b486fabe8bc48594416a7c343531441df
output_sha256: 7e8eaabb73fc3f86b3099d3b43aa3c3da925fb24a3081e7611a9faddccebaf79
hash_basis: LF-normalized bytes
---

# THM-3333 -- the Farey graph on the Pythagorean light cone

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem is a repository-new synthesis and proof interface.  No
literature-priority or global-novelty claim is made.

## 1. Inheritance, convention, and theorem

[THM-2056, Kelvin-polar/Farey defect](THM-2056-kelvin-polar-farey-defect-certificate.md)
proves that a primitive coefficient vector in a fixed rank-two LRC relation
plane has exact Gaussian norm denominator, and that all but a finite
uncertified set of primitive parameter directions can be discharged by a
finite acute unimodular fan.  Gate failure there means only uncertified.
[THM-2596, modular Gram-owner covariance](THM-2596-modular-free-factor-farey-gram-owner-cocycle.md)
puts the Berggren reduction grammar on the same rational line and proves
that endpoint defect values alone do not determine a Farey child.  The
faithful inherited state is the Gram-owner pair.

There is one convention change.  THM-2596 writes `0<m<n` and

```text
(a,b,c)=(n^2-m^2,2mn,n^2+m^2).
```

Here the positive chart uses `m>n>0` and

```text
(a,b,c)=(m^2-n^2,2mn,m^2+n^2).                           (1)
```

The parameter swap `J(m,n)=(n,m)` has determinant `-1`: it preserves dot
products and norms, reverses determinant and face orientation, and conjugates
the branch matrices.  Barning matrices from THM-2596 must not be reused here
without that gauge change.

For every integer `k`, put

```text
T(k)=k(k+1)/2.                                             (2)
```

For `u=(m,n)` define the ordered signed Gaussian-square lift

```text
Phi(u)=(A(u),B(u),C(u))
      =(m^2-n^2,2mn,m^2+n^2).                             (3)
```

On triples use the Lorentz pairing

```text
<X,Y>_L=X_C Y_C-X_A Y_A-X_B Y_B.                          (4)
```

For `u=(m,n)` and `v=(p,q)`, set

```text
Delta=det(u,v)=mq-np,             h=u dot v=mp+nq.        (5)
```

Write `X=Phi(u)=(A,B,C)` and `Y=Phi(v)=(A',B',C')`.  Then

```text
C C'=h^2+Delta^2,                                         (6)
<Phi(u),Phi(v)>_L=2 Delta^2.                              (7)
```

If `F` is the graph on primitive projective rays
`Z^2_prim/{+-1}` with edge predicate `|det(u,v)|=1`, and `P_Phi` is the
graph on the image of (3) with edge predicate

```text
<Phi(u),Phi(v)>_L=2,                                      (8)
```

then

```text
F isomorphic to P_Phi.                                   (9)
```

Thus the light-cone graph is exactly a Pythagorean-coordinate copy of the
Farey graph, not a new graph with unknown topology.  The qualification
"on the image of Phi" is load-bearing: arbitrary integral null triples with
pairing two need not carry compatible raw Gaussian lifts.

The same carrier unifies the one-vertex and edge laws.  Define

```text
R(X)=T(X_A)+T(X_B)-T(X_C).                                (10)
```

Then, for all integer triples `X,Y`,

```text
R(X)+R(Y)-R(X+Y)=<X,Y>_L.                                (11)
```

On every positive integral right triple `(a,b,c)` with inradius `r`,

```text
R(a,b,c)=T(a)+T(b)-T(c)=r.                               (12)
```

The radius at a null vertex and the squared-intersection index on an edge
are therefore the value and polarization of one triangular potential.

## 2. Gaussian squaring is the rank-one symmetric-square cone

Write `z=m+in`.  Then

```text
A(u)+iB(u)=z^2,                  C(u)=z conjugate(z),      (13)
```

so `A(u)^2+B(u)^2=C(u)^2`.  There is also a real structural model.  For
`X=(A,B,C)`, put

```text
M(X)=[C+A   B ].                                          (14)
     [ B   C-A]
```

Then

```text
det M(X)=C^2-A^2-B^2,
M(Phi(u))=2u u^T.                                        (15)
```

Thus the Pythagorean points are rank-one positive semidefinite integral
forms on the Lorentz null cone.  For `g in GL_2(Z)`,

```text
M(Phi(gu))=g M(Phi(u)) g^T.                              (16)
```

This symmetric-square covariance is the mechanism behind the modular and
Lorentz appearances; it is not an active symmetry of a fixed Euclidean LRC
owner predicate.

The map `Phi` is injective on primitive projective rays.  Away from `m=0`,

```text
n/m=B/(C+A),                    C+A=2m^2,                 (17)
```

and the endpoint `m=0` is unique.  Equivalently, equality of two Gaussian
squares forces their square roots to differ by sign.

Primitivity has one exact parity boundary.  If `gcd(m,n)=1`, then

```text
g(u):=gcd(A,B,C)=2  iff m,n are both odd,
g(u)=1              iff m,n have opposite parity.        (18)
```

No odd common prime is possible.  In the odd/odd case `B` and `C` have
2-adic valuation exactly one, so the common content is exactly two.

Every Farey face `{u,v,u+v}` reduces modulo two to the three nonzero vectors
of `F_2^2`.  Consequently its three raw Pythagorean contents are

```text
sort(g(u),g(v),g(u+v))=(1,1,2).                          (19)
```

After primitive-normalizing each triple, the sorted triple of Lorentz edge
weights becomes `(1,1,2)` instead of `(2,2,2)`.  This is the light-cone scale version
of the three Farey parity channels in THM-2632.

For primitive `u`, both coordinates of the Kelvin inverse

```text
I(u)=u/C(u)                                               (20)
```

have exact common denominator `C(u)`, since
`gcd(m,C(u))=gcd(n,C(u))=1`.  This recovers the represented
sum-of-two-squares denominator used by THM-2056.

## 3. Triangular polarization, torus intersection, and the exact graph

Triangular numbers are the integral quadratic refinement of multiplication:

```text
P_T(x,y):=T(x+y)-T(x)-T(y)=xy                            (21)
```

for all integers `x,y`.  Applying (21) to (10) proves (11), and (7) then
gives

```text
R(Phi(u))+R(Phi(v))-R(Phi(u)+Phi(v))=2 Delta^2.           (22)
```

Notice that the last argument is addition in the ambient triple lattice.
It is not the spinor mediant:

```text
Phi(2,1)+Phi(3,2)=(8,16,18),
Phi(5,3)=(16,30,34).                                     (23)
```

Equation (22) polarizes a coordinate potential; it does not say that `Phi`
is additive.

There is a literal topological reading.  Primitive projective vectors are
unoriented primitive slope classes of essential simple closed curves on the
torus.  Their minimal geometric intersection number is

```text
i([u],[v])=|det(u,v)|.                                   (24)
```

Straight representatives prove attainability, while algebraic intersection
gives the lower bound.  Therefore

```text
<Phi(u),Phi(v)>_L=2 i([u],[v])^2.                        (25)
```

The Farey graph is precisely the minimal positive intersection shell.  The
Lorentz scalar remembers geometric intersection squared and forgets the
oriented algebraic-intersection sign.

The converse in (9) is scoped to canonical raw images.  For example,

```text
X=(4,3,5),                 Y=(6,8,10),
<X,Y>_L=2,                                                (26)
```

but `X` is not `Phi(u)` for an integral spinor, because `C+A=9` is not twice
an integer square.  More generally, compatible scaled lifts satisfy

```text
<k Phi(u), l Phi(v)>_L=2kl det(u,v)^2.                    (27)
```

## 4. An oriented Farey face and its scalar exchange law

Choose ordered signed representatives `(u,v)` before using face signs.
Replacing one spinor by its negative leaves its `Phi` image fixed, reverses
`Delta`, and swaps `[u+v]` with `[u-v]`.  Thus an oriented face has covariant
orientation data, but a bare projective edge does not intrinsically choose a
pre-existing spinor determinant sign.

Put `w_+=u+v` and `w_-=u-v`.  Quadraticity gives the vector diamond

```text
Phi(w_+)+Phi(w_-)=2Phi(u)+2Phi(v).                        (28)
```

With triples written as columns,

```text
det[Phi(u),Phi(v),Phi(u+v)]=-4 Delta^3.                  (29)
```

The sign records the chosen oriented face.  The two incident hypotenuses
satisfy

```text
C(w_+)=C(u)+C(v)+2h,
C(w_-)=C(u)+C(v)-2h.                                    (30)
```

An unlabeled face pair therefore recovers `|h|`; labeled plus/minus faces
recover signed `h`, and an acute flank fixes `h>=0`.

There is a useful scalar exchange form.  On a unimodular edge write

```text
x=C(u),       y=C(v),       z_+=C(u+v),       z_-=C(u-v).
```

Then

```text
h^2=xy-1,
z_++z_-=2(x+y),
z_+ z_-=(x-y)^2+4,                                      (31)
```

and either incident face norm `z` lies on

```text
x^2+y^2+z^2-2xy-2xz-2yz=-4.                             (32)
```

Thus a labeled edge has an exact Vieta-type face mutation
`z -> 2(x+y)-z`.  It is a local operation law, not a global classification
by norm values: Section 6 gives the representation collision.

On an acute unimodular flank,

```text
h=sqrt(C(u)C(v)-1),                                      (33)
```

so the two endpoint norms determine the full Gram matrix.  This strictly
repairs THM-2596's endpoint-defect hostile on a fixed labeled owner cone.
Indeed, for

```text
F_w(d)=||d||^2-91 w dot d,                               (34)
```

the endpoint pairs `(C,F_w)` give

```text
w dot u=(C(u)-F_w(u))/91,
w dot v=(C(v)-F_w(v))/91,                               (35)
```

and the mediant update is

```text
C(u+v)=C(u)+C(v)+2 sqrt(C(u)C(v)-1),
F_w(u+v)=F_w(u)+F_w(v)+2 sqrt(C(u)C(v)-1).               (36)
```

Hence `(C,F_w,owner-id)` is an exact scalar Farey state while all descendants
stay in one acute owner cone.  Its entries are unbounded, and an owner switch
or tie still needs the signed owner sidecar.

Projectivizing a nonzero image to `(A/C,B/C)` gives

```text
(m,n)=rho(cos theta,sin theta)
       -> (A/C,B/C)=(cos 2theta,sin 2theta).              (37)
```

Because `theta` is taken modulo `pi`, this is a homeomorphism
`RP^1 -> S^1`.  Drawing hyperbolic geodesics between adjacent boundary
points gives the ordinary ideal Farey triangulation of the open disk and its
trivalent dual tree.  The cusp set is infinite-valence and non-locally-finite;
it is not a closed triangulated-disk boundary.

The graph itself contains Farey triangles, so it is neither bipartite nor a
partial cube.  Its edge relation is symmetric, so a tournament orientation
would be an added gauge.  This is distinct from THM-3145's subdivided-`K4`
incidence partial cube and from THM-3157's finite, triangle-rich Johnson
octahedron.  THM-2975 supplies the general warning that a lossy finite shadow
may be a partial cube when the literal carrier is not.

Finally, Pythagorean vertices do not make this the Berggren/Barning tree.
The `C_0` branch in the present convention sends `(2,1)` to `(4,1)`, but

```text
|det((2,1),(4,1))|=2,
<Phi(2,1),Phi(4,1)>_L=8.                                (38)
```

It is a length-two Farey prefix, not an edge of `P_Phi`.

## 5. Two radial gauges: constant graph edge versus constant LRC gate

Extend (3) to a real vector `x=(x_1,x_2)` by

```text
H(x)=(x_1^2-x_2^2,2x_1x_2,x_1^2+x_2^2).                 (39)
```

The proof of (7) works over the reals:

```text
<H(x),H(y)>_L=2 det(x,y)^2.                              (40)
```

For a primitive integer direction `u`, put `q_u=C(u)` and define

```text
P_u=H(u)=Phi(u),
Q_u=H(I(u))=P_u/q_u^2,             (Q_u)_C=1/q_u.         (41)
```

Let the labelled rank-two relation columns in THM-2056 be `c_i`, and put

```text
S_i=H(c_i),
D(u)=max_i |det(u,c_i)|.                                 (42)
```

Then

```text
max_i <Q_u,S_i>_L=2D(u)^2/q_u^2,                         (43)
```

and the proved THM-2056 sufficient gate is exactly

```text
D(u)<=q_u/91
iff max_i <Q_u,S_i>_L<=2/91^2.                           (44)
```

The integer gauge `P_u` makes Farey adjacency the constant value `2`; the
Kelvin-square gauge `Q_u` makes the LRC determinant gate the constant cap
`2/91^2`.  Projectively they are the same rational Pythagorean-circle point,
but `q_u` is the load-bearing radial scale relating the two simple laws.

Equation (44) is an exact reformulation of a sufficient gate, not an LRC
criterion.  It still lacks the signed hull owner, non-hull runner labels,
phase height, pair-sum clock, endpoint owner word, and global first-exit data.

## 6. Why represented sums of two squares and leg order matter

The norm value `C` does not retain its Gaussian representation.  A minimal
exact hostile is

```text
u=(8,1),       Phi(u)=(63,16,65),
u'=(7,4),      Phi(u')=(33,56,65),
v=(1,0),       Phi(v)=(1,0,1).                            (45)
```

Although the first two norms are both `65`,

```text
|det(u,v)|=1,        <Phi(u),Phi(v)>_L=2,
|det(u',v)|=4,       <Phi(u'),Phi(v)>_L=32.               (46)
```

The graph therefore lives on represented sums of two squares.  Its ordered
signed Pythagorean legs are the representation sidecar.

Primitive normalization changes the raw scale.  The positive Farey edge

```text
u=(3,1),       Phi(u)=(8,6,10),
v=(2,1),       Phi(v)=(3,4,5)                             (47)
```

has raw pairing `2`; replacing the first point by `(4,3,5)` changes it to
`1`.  An ordered normalized triple still determines whether its raw content
was one or two: in the retained `(A,B,C)` order, `g(u)=1` iff `B` is even
and `g(u)=2` iff `B` is odd.  Thus `g(u)` in (18) is reconstructible if
ordered legs are retained.

Forgetting leg order is genuinely destructive.  Both normalized triples in
(47) become the same unordered geometric triple `{3,4,5}`, collapsing a
Farey edge to one vertex.  Even without normalization, swapping one leg in
the adjacent pair `(3,4,5),(5,12,13)` changes the Lorentz pairing from `2`
to `9`.  The exact graph is on compatible raw spinor lifts, not on classical
unordered primitive triples.

## 7. Right-triangle radii and the polygonal hierarchy

For arbitrary integers `a,b,c`, expansion gives

```text
T(a)+T(b)-T(c)-(a+b-c)/2=(a^2+b^2-c^2)/2.                (48)
```

Hence

```text
T(a)+T(b)-T(c)=(a+b-c)/2
iff a^2+b^2=c^2.                                         (49)
```

For a positive right triangle,

```text
(a+b-c)(a+b+c)=2ab,                                      (50)
```

so

```text
r=(a+b-c)/2=ab/(a+b+c)=area/semiperimeter.               (51)
```

For the unscaled ordered lift (1),

```text
r=n(m-n)=T(m)-T(n)-T(m-n).                               (52)
```

For every positive integer `k`, the radius of `k Phi(m,n)` is `k n(m-n)`;
formula (52) without the factor `k` is not a statement about every
nonprimitive triple.

The four THM-2606 barycentric sign representatives expose more structure.
For `epsilon=(epsilon_A,epsilon_B,epsilon_C) in {+-1}^3`, define

```text
R_epsilon(X)=T(epsilon_A X_A)+T(epsilon_B X_B)
                                      -T(epsilon_C X_C). (53)
```

Since `T(epsilon k)=(k^2+epsilon k)/2`, for all integer triples `X,Y`,

```text
R_epsilon(X)+R_epsilon(Y)-R_epsilon(X+Y)=<X,Y>_L,         (54)
```

independently of `epsilon`.  On a null triple,

```text
R_epsilon(a,b,c)
 =(epsilon_A a+epsilon_B b-epsilon_C c)/2,
R_(-epsilon)=-R_epsilon.                                 (55)
```

Thus the magnitude descends to the four sign classes
`{+-1}^3/<(-1,-1,-1)>`, which form `V_4`.  With the THM-2606 representatives

```text
I=(+++),       I_A=(-++),       I_B=(+-+),       I_C=(++-),
```

one obtains exactly

```text
(R_I,R_(I_A),R_(I_B),R_(I_C))=(r,-r_a,-r_b,r_c).         (56)
```

For the Euclid lift,

```text
(r,r_a,r_b,r_c)
= (n(m-n),m(m-n),n(m+n),m(m+n)).                         (57)
```

This is also a literal graph toggle.  For `k>=0`, put
`Gamma_k=K_(k+1)^(1)`.

```text
T(k)=|E(Gamma_k)|,
T(-k)=T(k-1)=beta_1(Gamma_k).                            (58)
```

Writing `E_k=|E(Gamma_k)|` and `beta_k=beta_1(Gamma_k)`, the four entries in
(56) are respectively

```text
 E_a+E_b-E_c,          beta_a+E_b-E_c,
 E_a+beta_b-E_c,       E_a+E_b-beta_c.                   (59)
```

Sign flips toggle an edge count to a cycle rank.  The common polarization
(54) cannot recover which in/excenter origin was chosen.  Feuerbach supplies
the all-positive origin in THM-2606; this radius identity does not create a
`V_4` motion of circles or a quartic-resolvent transfer.

The phenomenon sits in a complete polygonal hierarchy.  For `s>=3`, let

```text
P_s(k)=((s-2)k^2-(s-4)k)/2.                              (60)
```

Then

```text
P_s(x+y)-P_s(x)-P_s(y)=(s-2)xy.                         (61)
```

If `G_s(X)=P_s(X_A)+P_s(X_B)-P_s(X_C)`, then

```text
G_s(X)+G_s(Y)-G_s(X+Y)=(s-2)<X,Y>_L,                    (62)
G_s(a,b,c)=(4-s)r                  on a right triple.    (63)
```

Thus a Farey edge has polygonal polarization `2(s-2)`.  Triangular numbers
`s=3` are uniquely the usual polygonal family with both unit Lorentz
polarization and positive vertex value `r`; squares `s=4` retain the edge
polarization but have zero null-vertex defect.

## 8. A two-clique surgery lemma at parameter and triangle scales

The graph topology follows from one reusable lemma.  Let finite sets
`V_1,V_2` cover `V`, with nonempty overlap, and put

```text
X=K_V^(1),
Y=K_(V_1)^(1) union K_(V_2)^(1),
p=|V_1-V_2|,                   q=|V_2-V_1|.              (64)
```

The missing edges are exactly `K_(p,q)`.  Because `X,Y` have the same
vertices, `C_0(X,Y)=0`; because they are graphs, `C_2(X,Y)=0`.  Hence

```text
H_1(X,Y;Z) isomorphic to Z^(pq).                          (65)
```

There is no relative-kernel ambiguity: every missing edge is a free relative
one-cycle.  But after clique filling,

```text
Cl(X)=Delta^(|V|-1),
Cl(Y)=Delta^(|V_1|-1) union_(Delta^(|V_1 intersect V_2|-1))
                                      Delta^(|V_2|-1).   (66)
```

Both complexes are contractible, so

```text
H_j(Cl(X),Cl(Y);Z)=0                   for every j.       (67)
```

At the Euclid-parameter scale, for `m>n>=1` take `|V|=m+1`,
`|V_1|=n+1`, `|V_2|=m-n+1`, and make the overlap one vertex.  Then

```text
p=n,             q=m-n,             pq=r.                (68)
```

This realizes (52) as an actual relative graph group:

```text
H_1(X_parameter,Y_parameter;Z) isomorphic to Z^r.         (69)
```

Both the edge and cycle invoices are `r`:

```text
T(m)-T(n)-T(m-n)=r,
T(m-1)-T(n-1)-T(m-n-1)=r.                               (70)
```

At the triangle-side scale, take `|V|=c+1`,
`|V_1|=a+1`, and `|V_2|=b+1`.  The forced overlap is

```text
|V_1 intersect V_2|=a+b-c+1=2r+1.                       (71)
```

The exclusive sizes are `c-b` and `c-a`.  Every positive right triple obeys

```text
(c-a)(c-b)=2r^2,                                         (72)
```

and for the unscaled Euclid lift these are

```text
c-b=(m-n)^2,                  c-a=2n^2.                  (73)
```

Therefore

```text
H_1(X_triangle,Y_triangle;Z) isomorphic to Z^(2r^2).      (74)
```

The separately named complete graphs also give the edge-cycle dual pair

```text
|E(Gamma_a)|+|E(Gamma_b)|-|E(Gamma_c)|=r,
beta_1(Gamma_c)-beta_1(Gamma_a)-beta_1(Gamma_b)=r.        (75)
```

For the included two-clique graph,

```text
beta_1(Y_triangle)
 =T(a-1)+T(b-1)-T(2r-1),
beta_1(X_triangle)-beta_1(Y_triangle)=2r^2.               (76)
```

Equations (69) and (74) are two graph scales of the same inradius: the
parameter cut has rank `r`, while the side/hypotenuse cut has rank `2r^2`.
Both vanish after flag filling.  They are one-skeleton cycle invoices, not a
Cech/Hodge obstruction or physical LRC current.  The octahedral current
discussion routed through HYP-2887/HYP-3061 remains hypothetical and is not a
dependency here.

## 9. Transfer contract and frontier boundary

| Source | Target / map | Preserved predicate | Destroyed information | Repair / cheapest decisive test |
|---|---|---|---|---|
| primitive torus/Gaussian ray `[u]` | raw `Phi(u)` | slope, exact norm, squared intersection, Farey adjacency | only the already-quotiented spin sign | (17) reconstructs the ray |
| raw image | ordered primitive triple | classical primitive triangle | numerical raw scale | reconstruct `g(u)` from retained `B`-parity before using threshold `2` |
| ordered triple | unordered triangle | side multiset | Euclid-leg order and graph injectivity | edge-collapse hostile (47) |
| represented norm | integer `C` | sum-of-two-squares value | representation and adjacency | retain ordered legs; norm-`65` hostile (45)--(46) |
| oriented Farey face | Lorentz edge scalar | determinant magnitude | spin orientation and signed Gram cross term | retain a spin lift and labeled face; otherwise only magnitudes survive |
| integer gauge `P_u` | Kelvin gauge `Q_u` | projective rational-circle point | radial denominator if scale is discarded | retain `q_u`; compare (8) with (44) |
| fixed owner flank | endpoint state `(C,F_w,owner)` | exact Gram-owner mediant recursion | owner changes, phase, clocks, global exit | replay the finite THM-2056 fan and deliberately test every switch/tie |
| graph pair `(X,Y)` | clique-complex pair | vertices and cliques | all relative graph cycles | `(3,4,5)` has graph ranks `1`/`2`, filled ranks `0` |

The strongest LRC-facing gain is now exact: (44) turns the determinant gate
into a Lorentz cap, while (36) compresses a fixed-owner acute Farey recursion
to scalar endpoint pairs `(C,F_w)`.  The next cheap experiment is to replay
every THM-2056 finite residual-fan decision using only this state plus an
owner identifier and to stop at the first owner switch or tie.  This is an
**OPEN experiment**, not a proved LRC terminal.

Even a successful fixed-owner replay would still omit non-hull runner labels,
phase height, pair-sum clocks, endpoint owner words, and global first exit.
The opposite face of an edge may leave the owner cone.  Consequently this
theorem excludes no LRC row and does not prove `LRC(14)`.

MISTAKE-222 blocks treating common triangular syntax as a cross-problem
mechanism; here the source, target, map, and preserved bilinear predicate are
all explicit.  MISTAKE-229 blocks calling the Gaussian gate a Heegner-form
criterion.  MISTAKE-352 blocks identifying a finite parity torsor or the
Berggren ternary grammar with the faithful modular/Farey carrier.

## 10. Scope, controls, and reproducibility

Equations (3), (6)--(7), (11), (15)--(16), (21)--(23), (28)--(32),
(39)--(43), (54), and (61)--(62) hold for all integer vectors/triples in the
displayed domains.  The graph equivalence (9) uses primitive projective rays
and canonical raw images.  Equations (48)--(51), (56), (63), (72), and (75)
hold for every positive integral right triple, primitive or not.  Equations
(52), (57), and (73) use the unscaled ordered Euclid lift; multiplying the
triple by a positive integer `k` multiplies all four radii and the exclusive
sizes by `k`.
Equation (58) is stated for `k>=0`.
The graph pairs are unique only up to relabeling; their isomorphism types and
ranks depend only on their displayed cardinalities.

Run

```bash
python3 04-computation/gaussian_farey_pythagorean_triangular_thm3333.py
python3 -O 04-computation/gaussian_farey_pythagorean_triangular_thm3333.py
```

Both modes reproduce the stored `1,568`-byte transcript exactly.  The main
universe is `713` primitive first-octant rays with `0<=n<=m<=48`, all
`253,828` unordered pairs, `1,423` Farey edges and face diamonds, all `435`
ordered Euclid rows with `1<=n<m<=30`, and five independently constructed
parameter/triangle graph and clique-complex controls.  The audit checks four
signed triangular potentials and polygonal orders `3..8` on every ray pair,
polygonal orders `3..12` and all in/exradii on every Euclid row, and five
positive/hostile Kelvin-cap controls.  It also checks symmetric-square
covariance, the Vieta face equations, and `5,692` fixed-owner scalar updates.

Positive graph controls are `(3,4,5)`, `(5,12,13)`, `(15,8,17)`,
`(7,24,25)`, and `(21,20,29)`.  Hostiles cover the norm-`65` representation
collision, positive odd/odd primitive normalization, unordered-leg edge
collapse, orientation reversal, determinant two, ambient-versus-spinor
addition, leg swap, arbitrary null pairing two, Berggren-versus-Farey edges,
and graph-versus-clique filling.  The code computes `Phi` both directly and
by Gaussian multiplication, and computes cycle ranks and relative clique
boundaries from explicit edge and triangle sets rather than trusting the
closed forms.

The computation corroborates the consequence objects.  The proof is the
integer algebra, symmetric-square model, and relative chain-complex argument
above.
