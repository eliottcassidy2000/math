---
id: THM-3336
title: "Primitive Gaussian multiplication, content cocycle, and the weighted Farey triangulation"
status: "RESERVED / UNPROVED PROVISIONAL CANDIDATE UNDER AUDIT"
source: codex-2026-08-12-primitive-gaussian-content-curvature
depends_on: []
related: []
script: 04-computation/primitive_gaussian_content_curvature_thm3336.py
output: 05-knowledge/results/primitive_gaussian_content_curvature_thm3336.out
---

# THM-3336 -- primitive Gaussian multiplication curves the Farey labels

**RESERVED / UNPROVED PROVISIONAL CANDIDATE UNDER AUDIT.**

This proof candidate is not yet part of proved canon.  It may not be cited as
a result until the exact companion, hostile audit, hashes, and dependencies
have been independently checked and the status is explicitly promoted.

No literature-priority or global-novelty claim is made.  The sum-of-two-squares
and Smith-normal-form ingredients are elementary; the intended new payload is
their typed operation-level assembly with primitive Pythagorean composition,
Farey-face index labels, signed radii, and the current LRC determinant gate.

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
Brahmagupta product.  Berggren depth and Farey distance cannot be additive
under this operation.

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

Consequently the Gaussian gcd of `N` and `h+i` recovers `a+ib` up to a unit,
and conjugation sends `h` to `-h`.  For coprime norms, multiplication glues
the two roots by the Chinese remainder theorem.  Shared represented primes
can instead cancel into the content: in (11), the roots of the two norm-five
representations are opposite modulo `5`.

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
The matrix `H` is an invertible linear repackaging of the triple, not a new
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
primitive Gaussian norm `N>=91`: use the content-`N` direction `(a,-b)`, a
Bezout neighbor, and a harmless positive shift.  No corresponding
all-representations claim is made for the reverse direction at the fixed
threshold `91`.

These are certificate flips, not LRC safety flips.  Gate failure means only
uncertified.  More importantly, raw left multiplication of every coefficient
column by `G_s` is merely a row-basis change of one rational plane, whose
parameter transforms contragrediently by `G_s^{-T}`.  Dividing different
columns by different `d_i`, as in (47)--(53), generally produces a different
saturated plane.  The examples compare two lawful primitive-column decks;
they do not prove or disprove invariance of the LRC maximum on one fixed plane.

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

The face mechanism has a useful maximal extension.  Let `A` be any integral
`2 x 2` matrix with positive determinant and gcd of its four entries equal to
one.  For primitive `u`, define

```text
d_A(u)=gcd(Au),             mu_A(u)=Au/d_A(u).            (55)
```

Smith normal form is `diag(1,|det A|)`.  Replacing `G_s,N` by
`A,|det A|` preserves the ordered content cocycle
`mu_A(mu_B(u))=mu_(AB)(u)`, the determinant identity (9), the edge and
weighted-face identities (19)--(31) and (35)--(36), and the intersection law
(54).  The content-range statement (32)--(33) also persists after replacing
`ell_h` by the primitive linear row supplied by Smith normal form; a general
matrix need not have a fixed-coordinate root of minus one.
It does **not** preserve the Euclidean norm formula (8): for example,
`A=diag(1,2)` and `u=(1,1)` give `N(Au)=5`, not
`|det A|N(u)=4`.  Gaussian structure contributes conformality, the extra
root `h^2=-1`, the Brahmagupta product, the sum-of-two-squares norm
interpretation, and the single-factor radius law (41).  This boundary
explains why content curvature is a primitive cyclic lattice phenomenon,
while the Pythagorean/radius and determinant-gate compilers are genuinely
Gaussian.

Each hypothesis matters independently:

```text
A=2I                         breaks content coprimality;
A=I, u=2(1,0)                breaks d_A(u)|det A;
A=diag(1,2), u=(0,1),v=(2,1) breaks pairwise contents
                              when det(u,v) is not +-1.  (56)
```

## 10. Exact consequence boundary and reproduction

What is proved here, after audit and promotion, is:

1. raw Gaussian multiplication is a Lorentz/norm similitude;
2. primitive reduction is governed by the exact content cocycle (10) and
   the root-of-minus-one linear form (16);
3. every transformed Farey face carries enough index data to reconstruct
   all three endpoint contents and the Gaussian norm;
4. Brahmagupta composition of primitive triples requires both Gaussian and
   parity contents;
5. the THM-2056 determinant gate acquires the per-column weights (47), which
   can reverse its sufficient verdict on lawful positive saturated decks.

It does **not** prove an LRC(14) row, a safety equivalence, a global Berggren
endomorphism, a canonical phase/owner, or a global exit from the finite
uncertified fan.  The surviving LRC opportunity is more precise: content is a
new exact coordinate for comparing primitive deck operations, and any use of
Gaussian composition must transport saturation, owner, phase, and clocks in
addition to the determinant ratio.

Reproduce the exact audit with

```bash
python3 04-computation/primitive_gaussian_content_curvature_thm3336.py
python3 -O 04-computation/primitive_gaussian_content_curvature_thm3336.py
```

Both modes must byte-match
`05-knowledge/results/primitive_gaussian_content_curvature_thm3336.out`.
