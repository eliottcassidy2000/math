---
id: THM-3336
title: "Primitive Gaussian multiplication, charge cohomology, content-curved Farey faces, and Boolean groupoid"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Primitive
  Gaussian multiplication is a raw Lorentz/norm similitude whose primitive
  reduction carries an exact content cocycle, split-prime charge H^1 classes,
  and content-curved Farey-face reconstruction.  Fixed-hypotenuse Boolean
  directions have canonical folded weights, but multiplication acts only
  through a section-dependent source groupoid.  Per-column contents can
  reverse the THM-2056 sufficient determinant gate on lawful saturated decks,
  while comparing different planes; LRC(14), owner/phase transport, a
  tournament orientation, and JC flux remain open.
source: codex-2026-08-12-primitive-gaussian-content-curvature
audit: >
  The promotion audit rederived the identities and four 13-column controls,
  swept 1,008 positive-determinant primitive-entry matrices over 233,856
  Farey faces and 1,968,128 ordered matrix cocycles, and reproduced the primary
  transcript.  A second independent hostile audit covered 1,600 Smith
  multipliers, all content patterns through norm 96, 984,960 primitive-matrix
  Farey faces, and 8,632 represented norms at least 91; it exposed the
  composition-domain issue recorded in MISTAKE-370, accepted the repair, and
  verified the universal gate construction.  A separate charge/Boolean audit
  repaired the C8 torsion, conjugation, section, shell, and H^1 typing and then
  accepted the result.  Both companions match in normal and optimized modes,
  and all four recorded hashes match.
depends_on:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
related:
  - THM-2053-rank-two-parameter-plane-geodesic-terminal
  - THM-2055-determinant-gate-normal-fan-and-tangent-sector-reduction
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors
script: 04-computation/primitive_gaussian_content_curvature_thm3336.py
output: 05-knowledge/results/primitive_gaussian_content_curvature_thm3336.out
script_sha256: a63d60669e11cd22ee0ef7afd619e54885308ce55d452a95e7131d19263e881f
output_sha256: 74cad2e39423f3f558d03f3be0b91e1d6102cf3106b6e2686c6d6f6f0a0a5709
secondary_script: 04-computation/gaussian_content_curved_farey_thm3336.py
secondary_output: 05-knowledge/results/gaussian_content_curved_farey_thm3336.out
secondary_script_sha256: 70d9d105ba30eae75262d564fa40319432d515be54817701656110260d35f0cc
secondary_output_sha256: 29247ef153043b0298414557e190d443ce20e4403ac87c8e2adaf823f88428e2
hash_basis: working-tree bytes (LF)
---

# THM-3336 -- primitive Gaussian multiplication curves the Farey labels

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

No literature-priority or global-novelty claim is made.  The sum-of-two-squares
and Smith-normal-form ingredients are elementary; the intended new payload is
their typed operation-level assembly with primitive Pythagorean composition,
split-prime charge cohomology, Farey-face index labels, fixed-hypotenuse
Boolean fibres, signed radii, and the current LRC determinant gate.

## 1. Inheritance, objects, and the raw similitude

[THM-3333, the Gaussian-square Farey light cone](THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md)
uses primitive projective spinors `u=(m,n)` and

```text
Phi(u)=(m^2-n^2,2mn,m^2+n^2).                            (1)
```

It proves

```text
<Phi(u),Phi(v)>_L=2 det(u,v)^2,                           (2)
```

so Farey adjacency is the Lorentz shell `2` on the raw labelled image.
[THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md)
classifies the equal-hypotenuse factor-choice fibre, while
[THM-2632](THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar.md)
shows that parity is a genuine retained channel.  Neither theorem gives a
binary operation on two primitive spinors.  That operation is the object here.
[THM-3341](THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md)
proves a special diagonal self-square transplant on the positive Berggren
middle ray.  Here multiplier and operand are independent and both Gaussian
and parity contents are retained; no Berggren-depth transplant is claimed.

Identify `(m,n)` with `m+in`.  Fix a primitive nonzero Gaussian multiplier

```text
s=(a,b),            gcd(a,b)=1,
G_s=[a -b],         N=N(s)=a^2+b^2.                      (3)
    [b  a]
```

For a primitive vector `u`, define

```text
d_s(u)=gcd((G_s u)_1,(G_s u)_2),
mu_s(u)=G_s u/d_s(u).                                    (4)
```

Absolute values are understood inside `gcd`.  The map is on projective rays,
so changing `u` to `-u` changes `mu_s(u)` only by sign.  Before primitive
normalization, Gaussian multiplication is an exact similarity:

```text
N(G_su)=N N(u),
(G_su) dot (G_sv)=N(u dot v),
det(G_su,G_sv)=N det(u,v).                               (5)
```

Write `Phi(s)=(A_s,B_s,C_s)` and `Phi(u)=(A,B,C)`.  On
triples define the Brahmagupta product

```text
(A,B,C) star (A',B',C')
  =(AA'-BB', AB'+BA', CC').                              (6)
```

Then Gaussian multiplication and (1) intertwine exactly:

```text
Phi(s) star Phi(u)=Phi(su),
<Phi(su),Phi(sv)>_L=N^2 <Phi(u),Phi(v)>_L.                (7)
```

Thus the raw operation is a Lorentz similitude, not a new Pythagorean
parameterization.  All curvature below enters when the two coordinates of
`su` are divided by their common content.

## 2. Primitive reduction is a multiplicative content cocycle

For every primitive `u`, equations (4)--(5) give

```text
N(mu_s(u))=N N(u)/d_s(u)^2,                              (8)
det(mu_s(u),mu_s(v))
  =N det(u,v)/(d_s(u)d_s(v)).                            (9)
```

Projectively, `mu_s` is just the rational boundary action of `G_s`; it is a
bijection, with inverse induced by the conjugate multiplier.  The primitive
normalization is nevertheless not dispensable.  For arbitrary nonzero
Gaussian integers `s,t` and primitive `u`, extend (4) in the evident way.
Then

```text
d_(st)(u)=d_t(u) d_s(mu_t(u)),
mu_s(mu_t(u))=mu_(st)(u).                               (10)
```

Indeed, `tu=d_t(u)mu_t(u)`, and the gcd of an integral vector scales by the
same positive integer.  Hence primitive Gaussian multiplication is a genuine
commutative projective action equipped with a multiplicative content cocycle.
The bare projective action forgets precisely the integers in (10).

The first hostile is already decisive:

```text
(2+i)(4+3i)=5+10i=5(1+2i),                              (11)
```

although both input spinors are primitive.  On raw triples,

```text
(3,4,5) star (7,24,25)=(-75,100,125)
                           =25(-3,4,5).                 (12)
```

Primitive Pythagorean triples are therefore not closed under the unscaled
Brahmagupta product.  The naive primitive product does not preserve Berggren
depth, and no Farey-distance additivity follows.

### 2A. The signed primitive group and its exact charge cohomology

The fixed-multiplier notation hides a symmetric group law.  For nonzero
signed primitive Gaussian integers define

```text
d(z,w)=d_z(w),                  z star_p w=mu_z(w).       (12a)
```

The subscript distinguishes this parameter operation from the triple product
in (6).  Positive content is removed without an additional unit choice.
Equations (10) imply literally

```text
(z star_p w) star_p t=z star_p (w star_p t),
d(z,w)d(z star_p w,t)=d(w,t)d(z,w star_p t).             (12b)
```

The operation is commutative, has identity `1`, and has inverse `conj(z)`.
It is therefore already an abelian group on signed primitive pairs.  Its
finite torsion is `C8`, generated by `1+i`:

```text
(1+i)^star_p2=i,                    (1+i)^star_p8=1.      (12c)
```

Let

```text
tau(z)=v_(1+i)(z) in {0,1};                                 (12d)
```

equivalently, `tau(z)=1` exactly when both coordinates are odd.  The
ramified-two content and primitive parity obey

```text
v_2(d(z,w))=tau(z)tau(w),
tau(z star_p w)=tau(z) XOR tau(w).                       (12e)
```

Define the odd content

```text
k_odd(z,w)=d(z,w)/2^(tau(z)tau(w)).                      (12f)
```

The binary carry identity makes it a multiplicative two-cocycle:

```text
k_odd(z,w)k_odd(z star_p w,t)
 =k_odd(w,t)k_odd(z,w star_p t).                         (12g)
```

This cocycle gives the clean primitive-triple composition law.  With

```text
T(z)=Phi(z)/2^tau(z),                                    (12h)
```

one has, for the signed triple product (6),

```text
T(z) star T(w)=k_odd(z,w)^2 T(z star_p w).               (12i)
```

Absolute values or canonical leg sorting cannot be inserted into (12i)
without recording the resulting unit/leg gauge.  If

```text
H(z)=N(z)/2^tau(z)                                       (12j)
```

is primitive hypotenuse grade, then

```text
H(z star_p w)=H(z)H(w)/k_odd(z,w)^2,
k_odd(z,w)^2=H(z)H(w)/H(z star_p w).                     (12k)
```

Thus the square content is intrinsically the multiplicative coboundary of
the grade.

The faithful signed coordinates are even sharper.  For every rational split
prime `p=1 mod 4`, choose one Gaussian prime `pi_p` above it and set

```text
lambda_p(z)=v_(pi_p)(z)-v_(conj(pi_p))(z).               (12l)
```

Modulo the four Gaussian units, the primitive group is

```text
G = C2 direct-sum (direct sum over p=1 mod4 of Z),        (12m)
```

with coordinates `(tau,(lambda_p)_p)`.  Every finitely supported coordinate
occurs, and

```text
lambda_p(z star_p w)=lambda_p(z)+lambda_p(w).             (12n)
```

With trivial coefficient actions, these are precise group-cohomology classes

```text
tau in H^1(G;F_2),                 lambda_p in H^1(G;Z). (12o)
```

They are not classes on the Farey tessellation, Berggren tree, Boolean graph,
LRC base, or a JC complement.

Odd content is coordinate cancellation:

```text
v_p(k_odd(z,w))
 =(|lambda_p(z)|+|lambda_p(w)|
   -|lambda_p(z)+lambda_p(w)|)/2.                        (12p)
```

After weighting coordinate `p` by `log p`, `log k_odd` is half the `l1`
triangle defect.  It is itself an exact coboundary after choosing the prime
orientations.  Put

```text
F_-(z)=product_p p^max(0,-lambda_p(z)).                  (12q)
```

Then

```text
k_odd(z,w)=F_-(z)F_-(w)/F_-(z star_p w).                (12r)
```

Equation (12r) is orientation-dependent; (12k) is the intrinsic statement.
The `H^1` classes are the signed charges, while nonnegative content is their
lossy cancellation shadow.

### 2B. Binary Farey content equals the endpoint norm gcd

If primitive `u,v` satisfy `|det(u,v)|=1`, then

```text
d(u,v)=gcd(N(u),N(v)),                                   (12s)
```

and this integer is odd.  At an odd shared split prime, the two norm-zero
vectors modulo `p` must occupy opposite isotropic Gaussian lines; the same
line would make their determinant vanish.  Opposite orientations contribute
exactly the minimum endpoint valuation to both sides of (12s).  Inert primes
cannot divide a primitive Gaussian norm, and two odd/odd endpoints would have
even determinant.

Writing `h=u dot v`, one obtains

```text
N(u)N(v)=h^2+1,
N(u star_p v)=(h^2+1)/gcd(N(u),N(v))^2,
h^2=-1 mod gcd(N(u),N(v))^2.                             (12t)
```

Adjacency is load-bearing: `d((2,1),(2,1))=1`, not `5`.  The sharp positive
control is

```text
u=(2,1),       v=(7,4),       det(u,v)=1,
uv=10+15i=5(2+3i),                                     (12u)
```

so the parameter content is `5`, the primitive output norm is `13`, and the
triple content in (12i) is `25`.

## 3. A root of minus one linearizes every content

Choose Bezout coefficients `r,t` with

```text
ar+bt=1,
h=at-br  (mod N).                                       (13)
```

Then

```text
L G_s=[1 h],       L=[ r t] in SL_2(Z),                 (14)
      [0 N]          [-b a]
1+h^2=N(r^2+t^2),       so h^2=-1 (mod N).              (15)
```

Different Bezout choices change `h` by a multiple of `N`.  Since unimodular
row operations preserve coordinate gcd, for every primitive `u=(m,n)` one
has the single-linear-form law

```text
d_s(m,n)=gcd(m+hn,N).                                   (16)
```

To see the last equality directly, (14) first gives
`gcd(m+hn,Nn)`.  The first entry is coprime to `n`, because a common prime
would divide both `m` and `n`; hence only the factor `N` remains.  Consequences
include:

1. `d_s(u)` always divides `N`.
2. Every divisor `d|N` occurs, already at `u=(d-h,1)`.
3. With `B=[1 -h;0 1]`, one has

   ```text
   L G_s B=diag(1,N).                                   (17)
   ```

Thus the primitive-content geometry is the cyclic Smith embedding of index
`N`, together with the fixed-coordinate root `h` selecting its absorbing
cusp.  The sum-of-two-squares structure is exactly what strengthens the
generic Smith parameter to `h^2=-1`.

The root is also the operation-level sidecar of THM-3334's representation
fibre.  Equations (13)--(15) imply

```text
a h=-b (mod N),              b h=a (mod N).              (18)
```

There is an explicit factorization

```text
h+i=(a+ib)(t+ir),
det((a,-b),(t,r))=ar+bt=1.                              (18a)
```

The two Gaussian factors in the second determinant are coprime, so the
Gaussian gcd of `N=(a+ib)(a-ib)` and `h+i` recovers `a+ib` up to a unit.
Conjugation sends `h` to `-h`.  For coprime multiplier norms, the other
multiplier is invertible modulo either norm, so the product root restricts to
the original root in both factors and is their Chinese-remainder gluing.
Shared represented primes can instead cancel into the content: in (11), the
roots of the two norm-five representations are opposite modulo `5`.

The primitive hypothesis on `s` is load-bearing.  For `s=(2,0)`, the axis
Farey face has contents `(2,2,2)` although `N=4`; all pairwise-coprimality and
index formulas in the next section would fail.

## 4. Farey edges become exact Lorentz shells

Let primitive `u,v` satisfy `det(u,v)=delta` with `delta=+-1`, and abbreviate

```text
d_u=d_s(u),     d_v=d_s(v),
U=mu_s(u),      V=mu_s(v).                               (19)
```

Then

```text
d_u|N,       d_v|N,       gcd(d_u,d_v)=1.               (20)
```

For if a prime divided both contents, `G_s` would kill both members of a
basis modulo that prime, forcing it to divide all four entries of `G_s`,
contrary to `gcd(a,b)=1`.  Therefore

```text
lambda_s(u,v)=N/(d_ud_v) is a positive integer,          (21)
det(U,V)=delta lambda_s(u,v),
<Phi(U),Phi(V)>_L=2 lambda_s(u,v)^2.                     (22)
```

An individual Farey edge survives as a standard Farey edge exactly when its
two endpoint contents exhaust `N`:

```text
|det(U,V)|=1  iff  d_ud_v=N.                             (23)
```

For `s=1+i`, the edge `((1,0),(1,1))` has endpoint contents `(1,2)` and
survives.  The edge `((1,0),(2,1))` has contents `(1,1)` and moves to
determinant shell `2`, or Lorentz shell `8`.  Hence no nonunit Gaussian
multiplier preserves the whole standard Farey graph.  In particular the
axis edge already moves to determinant `N`.

## 5. One transformed face reconstructs the lost norm

Retain the hypotheses of Section 4 and put

```text
w=u+v,      W=mu_s(w),      d_w=d_s(w),
D=d_ud_vd_w,                kappa=N/D.                  (24)
```

Every two vertices of `(u,v,w)` form a unimodular pair.  Repeating (20)
therefore proves that `d_u,d_v,d_w` are pairwise coprime.  Since all three
divide `N`, `D|N` and `kappa` is a positive integer.  Primitive reduction
turns the ordinary mediant and its three determinants into the weighted laws

```text
d_w W=d_u U+d_v V,                                      (25)
det(U,V)= delta kappa d_w,
det(V,W)=-delta kappa d_u,
det(W,U)=-delta kappa d_v.                              (26)
```

At `p=2`, the same coprimality can also be seen directly: a Farey face has
the three nonzero residues of `F_2^2`, hence exactly one odd/odd vertex, so
two divides at most one endpoint content.

Equations (22) and (26) use raw `Phi` lifts.  With the primitive normalized
triple `T(x)=Phi(x)/2^tau(x)`, every image edge instead obeys

```text
<T(x_i),T(x_j)>_L
 =2 det(x_i,x_j)^2/2^(tau(x_i)+tau(x_j)).                (26a)
```

Projectively the image is an ideal triangulation under `[G_s] in PGL_2(Q)`,
but its weighted edges are generally not determinant-one edges of the
standard Farey graph.

Label each edge by the absolute determinant opposite a vertex:

```text
ell_u=|det(V,W)|,   ell_v=|det(W,U)|,
ell_w=|det(U,V)|.                                       (27)
```

The primitive image face itself then reconstructs every scalar that
normalization hid:

```text
kappa=gcd(ell_u,ell_v,ell_w)
      =[Z^2:<U,V,W>],                                   (28)
d_i=ell_i/kappa,
N=lcm(ell_u,ell_v,ell_w)
 =ell_u ell_v ell_w/kappa^2.                            (29)
```

The lattice index in (28) is the gcd of the three `2 x 2` minors.  Thus the
unlabelled transported triangulation is abstractly just another Farey
triangulation and loses `N`; the embedding relative to the standard lattice,
or equivalently the three edge labels, recovers it exactly.

There is also a Pythagorean volume invoice.  For all integer vectors `x,y,z`,

```text
det[Phi(x),Phi(y),Phi(z)]
 =-4 det(x,y)det(y,z)det(z,x).                           (30)
```

Applying (26) gives

```text
det[Phi(U),Phi(V),Phi(W)]=-4 delta kappa^2 N.            (31)
```

No nonunit can send a whole Farey face to a standard Farey face: the product
of its three absolute edge labels is `kappa^2 N>1`.

### Exact face-content range

Let `ell_h(m,n)=m+hn`.  If

```text
alpha=ell_h(u),      beta=ell_h(v),
```

then `gcd(alpha,beta)=1` and

```text
(d_u,d_v,d_w)
 =(gcd(alpha,N),gcd(beta,N),gcd(alpha+beta,N)).           (32)
```

Conversely, for odd `N`, an ordered triple of positive divisors of `N` occurs
in (32) exactly when the three divisors are pairwise coprime.  If `N` is even,
primitivity of `s` forces `v_2(N)=1`, and occurrence is equivalent to pairwise
coprimality plus exactly one even divisor.

The converse is prime-local.  At an odd prime power, assign the desired
divisibility to `alpha`, `beta`, or `alpha+beta` by the residue pairs

```text
(p^j,1),       (1,p^j),       (1,-1+p^j),                (33)
```

and use `(1,1)` when that prime stays in `kappa`.  At `2`, the three nonzero
row states assign the unique factor `2` to exactly one of those entries.
Chinese remaindering and the elementary lift of a unimodular residue row to
`SL_2(Z)` produce the required Farey face.

Two hostiles delimit the allocation language.  With `s=1+2i`, `N=5`, the
axis face has contents `(1,1,1)` and `kappa=5`: an odd prime may remain wholly
unabsorbed.  With

```text
s=100+71i,  N=15041=13^2*89,  u=(100,59),
d_s(u)=gcd(5811,13000)=13,                               (34)
```

one has `gcd(13,N/13)=13`.  A content need not be a unitary divisor; the
correct statement is valuation allocation, not all-or-nothing prime-power
allocation.

## 6. Weighted Vieta exchange and the radius matrix

Across a Farey edge, put `w_+=u+v`, `w_-=u-v`, with primitive images `W_+,
W_-` and contents `d_+,d_-`.  Linear and quadratic parallelogram identities
give

```text
d_+ W_+ + d_- W_-=2d_u U,
d_+ W_+ - d_- W_-=2d_v V,                               (35)
d_+^2 Phi(W_+)+d_-^2 Phi(W_-)
 =2d_u^2 Phi(U)+2d_v^2 Phi(V).                          (36)
```

Equation (36) is the content-weighted Vieta exchange.  It is not a Berggren
parent-child law.

The four in/exradii from THM-3333 admit a compact operation law.  For a
triple `X=(A,B,C)`, define the all-positive-radius coordinate matrix

```text
H(X)=1/2 [ A+B-C   A-B+C]
          [-A+B+C  A+B+C].                              (37)
```

On a positive ordered right triangle its entries are `(r,r_a;r_b,r_c)`.
Outside that chamber they are signed algebraic coordinates.  Directly,

```text
det H(X)=(A^2+B^2-C^2)/2=-<X,X>_L/2.                    (38)
```

Thus a nonzero triple is null exactly when `H(X)` has rank one, and

```text
det(H(X)+H(Y))-det H(X)-det H(Y)=-<X,Y>_L.              (39)
```

For `u=(m,n)`,

```text
H(Phi(u))=[m-n][n m]
          [m+n]
 =[(n(m-n))  (m(m-n))].                                 (40)
  [(n(m+n))  (m(m+n))]
```

Gaussian multiplication has the unusual but exact two-sided covariance

```text
H(Phi(G_su))=G_s H(Phi(u)) G_s,
H(Phi(mu_s(u)))=G_s H(Phi(u)) G_s/d_s(u)^2.             (41)
```

The right factor in (41) is `G_s`, not `G_s^T`, and this is not a congruence
action.  If `L_0=[1 -1;1 1]` and `P=[0 1;1 0]`, then
`H(Phi(u))=(L_0u)(Pu)^T`, while `L_0G_s=G_sL_0` and
`PG_s=G_s^TP`; these two identities prove (41).

For two null points `Phi(U),Phi(V)`, (39) becomes

```text
|det(U,V)|^2=-det(H(Phi(U))+H(Phi(V)))/2.               (42)
```

Consequently the three pairwise determinant polarizations of the four-radius
matrices recover the labels (27), and hence `kappa` and `N` by (28)--(29).
The matrix `H` is an injective linear repackaging of the triple, not a new
quotient invariant; the content reconstruction is the substantive statement.

## 7. Primitive Pythagorean composition needs two contents

For primitive `x`, put

```text
epsilon(x)=gcd(Phi(x)) in {1,2},
P(x)=Phi(x)/epsilon(x).                                  (43)
```

Here `epsilon(x)=2` exactly for an odd/odd spinor.  Combining (7) with (4)
gives the exact operation on primitive ordered triples:

```text
P(s) star P(u)
 =[d_s(u)^2 epsilon(mu_s(u))/(epsilon(s)epsilon(u))]
    P(mu_s(u)).                                          (44)
```

The bracket is an integer, but it is not generally `d_s(u)^2`.  For the
positive hostile `s=u=(3,1)`, both input parity contents equal `2`, while
`d_s(u)=2` and `epsilon(mu_s(u))=1`; accordingly

```text
(4,3,5) star (4,3,5)=(7,24,25),                         (45)
```

with coefficient one, not four.  Primitive normalization therefore requires
both the Gaussian product content `d_s(u)` and THM-3333's parity content
`epsilon`; unordered legs erase the spinor gauge needed to evaluate them.

In the symmetric notation of Section 2A,

```text
epsilon(x)=2^tau(x),
d_s(u)^2 epsilon(mu_s(u))/(epsilon(s)epsilon(u))
 =k_odd(s,u)^2.                                         (45a)
```

Thus (44) and (12i) are the same exact law.  The two-content presentation in
(44) is convenient for fixed multipliers; the charge-cancellation presentation
in (12i) is convenient for composition and Boolean fibres.

## 7A. Fixed-hypotenuse fibres are weighted, not globally acted on

Let an admissible odd hypotenuse be

```text
c=product_(j=1)^r p_j^e_j,                 p_j=1 mod 4, (45b)
```

and choose `pi_j` above each `p_j`.  The raw Gaussian allocation cube is

```text
z_x=product_j pi_j^(e_j(1-x_j)) conj(pi_j)^(e_j x_j),
x in F_2^r.                                               (45c)
```

Global conjugation sends `x` to `x+1`.  After choosing the prime orientations
and one affine chart, THM-3334's parent torsor is therefore

```text
X_c isomorphic to F_2^r/<1>.                             (45d)
```

For raw lifts, (12p) specializes to

```text
k_odd(z_x,z_y)=product_(j:x_j!=y_j) p_j^e_j.             (45e)
```

Conjugating one lift replaces the displayed divisor by its complement in
`c`.  Hence a quotient direction `s=[x-y]`, represented by a subset `S`, has
the intrinsic folded weight

```text
K_c(s)={P_S,c/P_S},
P_S=product_(j in S)p_j^e_j.                             (45f)
```

Unique factorization makes (45f) injective on nonzero directions of (45d):
equality of unordered divisor pairs means equality or complementation of
subsets, exactly the quotient relation.  Every nonzero direction of an
affine binary space is one perfect matching, so (45f) canonically weights the
Cayley one-factorization of the complete parent graph.

At `c=65`, its one matching has folded weight

```text
{5,13}.                                                   (45g)
```

At THM-3334's first affine `V4` fibre

```text
c=1105=5*13*17,                                         (45h)
```

the three prime-XOR matchings have exact weights

```text
{5,221},                    {13,85},                    {17,65}. (45i)
```

These weights distinguish all three directions but are invariant under every
translation of the four parent vertices.  They select no owner, ancestry
path, or tournament orientation.

### Multiplication does not descend through conjugation

The operation `star_p` acts on signed primitive lifts or associate classes.
It does not directly act on the conjugation quotient (45d).  The minimal
hostile is

```text
2+i ~ 2-i,                     w=2+i,
(2+i) star_p w=3+4i with H=25,
(2-i) star_p w=1 with H=1.                              (45j)
```

In charge coordinates, the raw cube (45c) is an antipodal orthotope

```text
lambda_(p_j)(z_x) in {+e_j,-e_j}.                        (45k)
```

Multiplication by a fixed lift translates charges by a constant `mu`; a bit
flip reflects one coordinate.  No translation reflects both signs.  More
intrinsically, representative independence would require

```text
[mu+lambda]=[mu-lambda]                                  (45l)
```

for every antipodal class.  For nonzero `lambda`, equality forces `mu=0`.
Only zero split-charge gauges descend canonically, and they contribute no
odd-prime Boolean motion.

A chosen orientation section can make a nontrivial fixed multiplier look
like a quotient translation.  At `c=65`, take

```text
pi_5=2+i,       pi_13=3+2i,
z_0=4+7i,       z_1=8-i,
alpha=(2-i)^2=3-4i.                                     (45m)
```

Then

```text
alpha star_p z_0=8+i ~ z_1,
alpha star_p z_1=4-7i ~ z_0.                            (45n)
```

On both conjugate representatives, however, the output grade is `1625`
rather than `65`.  The action belongs to the chosen section and does not
descend.

### The lawful Boolean motion is a Gaussian groupoid

Let `g_(j,x)` denote the entire selected Gaussian prime-power factor at
coordinate `j` of `z_x`.  For a subset `S`, define the source-dependent
multiplier

```text
A_(S,x)=product_(j in S) conj(g_(j,x))^2.                (45o)
```

Then

```text
d(A_(S,x),z_x)=P_S,
A_(S,x) star_p z_x=z_(x+S).                             (45p)
```

These arrows compose as a groupoid over the raw allocation cube.  Their
dependence on `x` is not optional: it records which prime orientation must be
cancelled before its conjugate is installed.

The smallest one-sided hostile uses

```text
pi=2+i,                  alpha=pi^2=3+4i.                (45q)
```

Then `alpha conj(pi)=5pi` reduces back to grade `5`, whereas
`alpha pi=pi^3=2+11i` has grade `125`.  A fixed multiplier flips one
orientation and exits the grade on the other.

Finally, suppose a primitive multiplier preserves `H(z)` for every primitive
`z`.  Taking `z=1` in (12k) forces `H(alpha)=1`, so `alpha` is a Gaussian unit
or an associate of `1+i`.  Conversely those eight lifts have
`k_odd(alpha,z)=1` and preserve every grade.  Up to signed unit/leg gauge,
multiplication by `1+i` on an odd/odd pair is

```text
(m,n)->((m-n)/2,(m+n)/2),                                (45r)
```

the normalizer behind THM-3339's third Fibonacci ray.  It supplies the
ramified-two correction and no split-prime Boolean translation.

### Cohomology and ancestry scope

The cohomology in (12o) is group cohomology of the associate-class
multiplication group `G`.  It makes no statement about graph cohomology.  In
particular the external `K4` at `c=1105`, viewed as a graph, has first Betti
number `6-4+1=3`; only its clique filling is a tetrahedron with vanishing
positive-degree cohomology.

After choosing one parent as origin, XOR displacements reconstruct all other
vertices along paths in the fibre graph.  This is an affine calibration, not
an intrinsic owner.  Equal-hypotenuse parents form a directed ancestry
antichain.  The ambient Berggren tree has a unique undirected path between
them, but that path exits the fixed-hypotenuse fibre and supplies no canonical
fibre origin.

## 8. The determinant gate is content-weighted, not invariant

For a labelled primitive coefficient deck `c_1,...,c_13` and primitive
parameter direction `u`, THM-2056 uses

```text
D(u)=max_i |det(u,c_i)|,       q(u)=N(u),
D(u)<=q(u)/91.                                         (46)
```

Apply `mu_s` separately to the direction and every column, and write
`d_i=d_s(c_i)`.  Equations (8)--(9) give the exact ratio law

```text
|det(mu_s(u),mu_s(c_i))|/N(mu_s(u))
 =(d_s(u)/d_i) |det(u,c_i)|/N(u),                        (47)
```

and hence

```text
D'(mu_s(u))/q'(mu_s(u))
 =max_i (d_s(u)/d_i)|det(u,c_i)|/q(u).                  (48)
```

There is no uniform primitive-normalized gate transport: the per-column
contents are load-bearing.

This failure persists on lawful saturated positive `13`-column planes and in
both directions.  Take `s=1+i`.

```text
pass -> fail:
u=(9,7),  c_k=(5,4)+k(9,7),  0<=k<=11;
mu_s(u)=(1,8),  mu_s(c_k)=(2k+1,16k+9).

fail -> pass:
u=(8,1),  c_k=(7,1)+2k(8,1),  0<=k<=11;
mu_s(u)=(7,9),  mu_s(c_k)=(7k+3,9k+4).                  (49)
```

In each line the deck is `{u,c_0,...,c_11}`.  Both the source and target
column-minor gcds are one.  Primitive positive covectors are respectively

```text
(-38,49), (-50,7)     and     (-6,49), (-50,39),        (50)
```

giving positive, nonzero, distinct speed rows, so all four coefficient planes
meet the standing rank-two hypotheses.  At the displayed parameter rays,

```text
(D,q)=(1,130) -> (1,65),
(D,q)=(1,65)  -> (1,130).                               (51)
```

The first certificate changes from pass to fail and the second from fail to
pass.

The distortion is unbounded.  For `s=10+i`, `N=101`, the deck

```text
d=(10,-1),       c_k=(1,0)+k d                           (52)
```

has contents `(101,1,...,1)` and moves `(D,q)=(1,101)` to `(1,1)`.
Conversely,

```text
d=(1,0),        c_k=(10,-1)+101k d                       (53)
```

has contents `(1,101,...,101)` and moves `(1,1)` to `(1,101)`.  Taking
`s=a+i` gives both displayed directions for the unbounded sequence
`N=a^2+1`.  More generally, the pass-to-fail construction works for every
primitive Gaussian norm `N=a^2+b^2>=91`.  Choose signs with `a,b>0`, solve
`ar+bt=1`, and put

```text
d=(a,-b),       c=(t,r),       h=at-br,
c_j=c+(K+j)d,                                      0<=j<=11. (53a)
```

Then

```text
G_s d=(N,0),            G_s c_j=(h+(K+j)N,1).            (53b)
```

For `K` large enough that `h+KN>N`, the source and target decks
`{d,c_0,...,c_11}` are saturated, lie in open positive half-planes, and admit
generic primitive covectors with positive, nonzero, distinct speeds.  Their
content patterns are `(N,1,...,1)`, and the displayed direction changes
`(D,q)=(1,N)` to `(1,1)`.  Thus the sufficient gate passes and then fails.
No corresponding all-representations claim is made for the reverse direction
at the fixed threshold `91`.

These are certificate flips, not LRC safety flips.  Gate failure means only
uncertified.  More importantly, raw left multiplication of every coefficient
column by `G_s` is merely a row-basis change of one rational plane, whose
parameter transforms contragrediently by `G_s^{-T}`.  Dividing different
columns by different `d_i`, as in (47)--(53), generally produces a different
rational plane, and the displayed row lattice need not remain saturated.
The four examples are independently saturated and compare two lawful
primitive-column decks; they do not prove or disprove invariance of the LRC
maximum on one fixed plane.

## 9. Topological meaning and maximal generalization

Primitive projective vectors are unoriented essential torus curves, with
minimal geometric intersection `|det(u,v)|`.  The integer matrix `G_s` induces
an oriented degree-`N` torus covering.  After each image curve is divided by
its covering multiplicity `d_s(u)`, equation (9) says

```text
i(mu_s(u),mu_s(v))
 =N i(u,v)/(d_s(u)d_s(v)).                               (54)
```

Thus the edge labels and face index are literal intersection and covering
invoices.  On the rational boundary, `[G_s]` sends the ideal Farey
triangulation to a commensurator translate.  Only Gaussian units preserve the
standard triangulation; a nonunit produces a labelled overlay.

The face mechanism has a useful maximal extension, but its operation class
must be typed carefully.  For every nonsingular integral `2 x 2` matrix `M`
and primitive `u`, define

```text
d_M(u)=gcd(Mu),             mu_M(u)=Mu/d_M(u).            (55)
```

This full nonsingular class is closed under composition and obeys

```text
d_(ML)(u)=d_L(u)d_M(mu_L(u)),
mu_M(mu_L(u))=mu_(ML)(u).                                (55a)
```

If `g=cont(M)` is the gcd of the four entries and `M_0=M/g`, then

```text
mu_M=mu_(M_0),       d_M=g d_(M_0),
Delta(M)=|det M|/g^2=|det M_0|.                          (55b)
```

Thus `Delta`, not `|det M|`, is the effective primitive degree after scalar
content is discarded.  In particular, the content-one subclass is **not**
closed under composition.  The minimal Gaussian witness is

```text
A=B=[1 -1;1 1],       AB=[0 -2;2 0],                    (55c)
```

where `A,B` have entry-content one and determinant two, while `AB` has
entry-content two.  Correspondingly

```text
Delta(A)=Delta(B)=2,                 Delta(AB)=1.         (55c')
```

Thus `Delta` is an objectwise effective degree, not a multiplicative grading
of the primitive-normalized action.  The content cocycle is exactly the
sidecar that records this loss.

Now let `A` have positive determinant `N` and entry-content one.  Its Smith
normal form is `diag(1,N)`.  Replacing `G_s,N` by `A,N` preserves the
determinant identity (9), the edge and weighted-face identities (19)--(31) and
(35)--(36), and the intersection law (54).  The content formula (32) also
persists after replacing `ell_h` by the primitive linear row supplied by Smith
normal form; a general matrix need not have a fixed-coordinate root of minus
one.  For odd `N`, the exact range is again all pairwise-coprime divisor
triples.  If `e=v_2(N)>0`, the contents remain pairwise coprime and exactly one
is even, but its valuation may be any `1<=j<=e`; the residual factor
`2^(e-j)` remains in `kappa`.  For example,

```text
A=diag(1,4),       u=(2,1),       v=(1,1)                (55d)
```

gives the Farey-face contents `(2,1,1)` and `kappa=2`.
It does **not** preserve the Euclidean norm formula (8): for example,
`A=diag(1,2)` and `u=(1,1)` give `N(Au)=5`, not
`|det A|N(u)=4`.  Gaussian structure contributes conformality, the extra
root `h^2=-1`, the Brahmagupta product, the sum-of-two-squares norm
interpretation, and the single-factor radius law (41).  This boundary
explains why content curvature is a primitive cyclic lattice phenomenon,
while the Pythagorean/radius and determinant-gate compilers are genuinely
Gaussian.

Each hypothesis and boundary matters independently:

```text
A=2I                         breaks content coprimality;
A=I, u=2(1,0)                breaks d_A(u)|det A;
A=diag(1,2), u=(0,1),v=(2,1) breaks pairwise contents
                              when det(u,v) is not +-1.  (56)
```

## 10. Exact consequence boundary and reproduction

What is proved here is:

1. raw Gaussian multiplication is a Lorentz/norm similitude;
2. primitive reduction is governed by the exact content cocycle (10) and
   the root-of-minus-one linear form (16);
3. every transformed Farey face carries enough index data to reconstruct
   all three endpoint contents and the Gaussian norm;
4. Brahmagupta composition of primitive triples requires both Gaussian and
   parity contents;
5. the THM-2056 determinant gate acquires the per-column weights (47), which
   can reverse its sufficient verdict on lawful positive saturated decks;
6. signed primitive multiplication has explicit `C8`/split-prime charge
   coordinates, with `tau,lambda_p` the exact associate-group `H^1` classes
   and odd content an orientation-dependent coboundary;
7. fixed-hypotenuse Boolean directions have canonical folded content weights,
   while actual flips are source-dependent groupoid arrows and multiplication
   does not descend through conjugation without a chosen section.
8. the primitive-matrix extension lives in the full nonsingular integral
   composition class, with scalar content removed by the effective degree
   `Delta`; its larger `2`-adic range is not the Gaussian range.

It does **not** prove an LRC(14) row, a safety equivalence, a global Berggren
endomorphism, a canonical phase/owner, a tournament orientation, a JC flux
class, or a global exit from the finite uncertified fan.  The surviving LRC
opportunity is more precise: content is a new exact coordinate for comparing
primitive deck operations, and any use of Gaussian composition must transport
saturation, the conjugation section, owner, phase, and clocks in addition to
the determinant ratio.

Reproduce the exact audit with

```bash
python3 04-computation/primitive_gaussian_content_curvature_thm3336.py
python3 -O 04-computation/primitive_gaussian_content_curvature_thm3336.py
python3 04-computation/gaussian_content_curved_farey_thm3336.py
python3 -O 04-computation/gaussian_content_curved_farey_thm3336.py
```

For each companion, normal and optimized modes must match its stored
transcript after LF normalization (the repository evidence convention).  The
primary companion checks finite Gaussian and Smith-form censuses, face
reconstruction, Vieta/radius, parity, determinant-gate examples, and the three
explicit non-Gaussian hostiles.  The arbitrary primitive-matrix extension and
the universal `N>=91` construction are proved above rather than exhaustively
enumerated by that companion.  The secondary companion independently checks
the signed group/content cocycles,
charge coboundary, binary Farey gcd, raw and normalized shells, Boolean
weights/groupoid, conjugation and section hostiles, and universal grade
preservers.
