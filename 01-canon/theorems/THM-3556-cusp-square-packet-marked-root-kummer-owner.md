---
id: THM-3556
title: "Cusp-square packet as a marked-root/Kummer inverse-cubic owner"
status: >
  PROVED + VERIFIED-EXACT + FINITE-EXACT POSITIVE PACKET.  For every
  polynomial U(v,y), the cusp-square packet (L,T,U,S) is the coefficient
  packet of an inverse cubic L X^3+T X+2U that factors into the marked root
  X=-1/v and a quadratic Kummer pair.  Its discriminant is -4 L S^2, and
  the marked root escapes projectively at v=0.  None of the six natural
  two-coordinate projections is Keller.  One explicit four-coordinate
  packet nevertheless has differential rank two everywhere and has an exact
  pair of distinct geometric source points with equal packet value.  No
  polynomial target projection in any degree has nonzero constant Jacobian;
  more strongly, no arbitrary descended coefficient six-tuple combines the
  packet tangent minors to a nonzero constant.
  In integral coordinates its image is exactly a prime (3,4) complete
  intersection; adjoining the marked root gives its smooth normalization.
  The parameterization is birational but nonfinite, and its conductor has a
  genus-two unordered-pair model on the main open set.  This closes the
  explicit packet, not other packets or the planar Jacobian conjecture.
source: kps-s188
depends_on: []
related:
  - THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor
  - THM-3535-fixed-keller-full-wreath-and-all-level-linear-primitivity
  - THM-3548-planar-keller-conductance-shadow-gates
  - THM-3554-punctured-kummer-collision-surface-normal-form
  - THM-3555-catalan-thickening-universal-cubic-root-cover
companion: 04-computation/cusp_square_marked_root_kummer_owner_kps_s188.py
output: 05-knowledge/results/cusp_square_marked_root_kummer_owner_kps_s188.out
script_sha256: a3dc5f8c58c0fd18a81cdc77e96fe1f9717b76d40777ef0db81633f14257152b
output_sha256: e11cc70931087e21cff412cb244262eb64f1c5466eb9f1bc182871c8181072dd
quadratic_no_go_script: 04-computation/cusp_square_quadratic_projection_no_go_thm3556.py
quadratic_no_go_output: 05-knowledge/results/cusp_square_quadratic_projection_no_go_thm3556.out
quadratic_no_go_script_sha256: 2db8317c35df10a882a7f87837b7b09e48887b17fe6d67e67af1650681b36b7f
quadratic_no_go_output_sha256: a21b0cc00f077bb67a9eda80ee39628c3988fe66a5d534a070ff0c09865bd975
cubic_no_go_script: 04-computation/cusp_square_cubic_projection_no_go_thm3556.py
cubic_no_go_output: 05-knowledge/results/cusp_square_cubic_projection_no_go_thm3556.out
cubic_no_go_script_sha256: 0b0a7f129ebbed95cc85e4ed8dad2f8996ab3dc0d8c66b468e1452d7bcaeb281
cubic_no_go_output_sha256: 41841c39426d6a42872fe1ea09bf9cc9169c3efe8473097e55bd84871770b933
image_ideal_script: 04-computation/cusp_square_packet_image_ideal_thm3556.py
image_ideal_output: 05-knowledge/results/cusp_square_packet_image_ideal_thm3556.out
image_ideal_script_sha256: 392b4eb2b3cae0cdfe985a1ca4e15dc36121a1a7d2bf1e66e23824a502728755
image_ideal_output_sha256: 494c04e2a63b88d70e9b94076002d92d6d3f4d7d8747213ade0e3b6bcf7d850a
all_degree_no_go_script: 04-computation/cusp_square_all_degree_descended_bivector_no_go_thm3556.py
all_degree_no_go_output: 05-knowledge/results/cusp_square_all_degree_descended_bivector_no_go_thm3556.out
all_degree_no_go_script_sha256: 78de3e52d2095c7cedfd47b1d5533a8a37f2c0c7aab053a73409062fae428818
all_degree_no_go_output_sha256: 323a884e710e4a2698f1267fc2596b69547da7b5cbc4e76ec6979a32fbfa23dc
hash_basis: LF-normalized bytes
---

# THM-3556 -- cusp-square packet as a marked-root/Kummer inverse-cubic owner

**PROVED + VERIFIED-EXACT + FINITE-EXACT POSITIVE PACKET.**  The
four-coordinate cusp packet is not merely a discriminant identity.  It is an
explicit `1+2` inverse-cubic owner: one marked root and one quadratic Kummer
pair.  The marked root reaches infinity exactly where the leading cubic
coefficient vanishes.

The field has characteristic zero.

## 1. The packet and its inverse cubic

Let `U=U(v,y)` be any polynomial and define

```text
T=y^2-6vU,
S=y^3-9vUy,
L=v^2(8vU-y^2).                                       (1)
```

Then the polynomial identity

```text
S^2=T^3+27LU^2                                        (2)
```

holds.  Associate to `(1)` the cubic

```text
E(X)=LX^3+TX+2U.                                      (3)
```

It has the exact factorization

```text
E(X)=(vX+1)
     [v(8vU-y^2)X^2+(y^2-8vU)X+2U].                  (4)
```

For `v!=0`, the first factor gives the marked root

```text
X_0=-1/v.                                              (5)
```

The factorization is better read projectively.  The homogenized linear
factor is `vX+Z`, so its root is `[-1:v]`.  At `v=0` it becomes `[1:0]`, the
point at infinity.  Thus the marked root genuinely escapes when the cubic
leading coefficient `L`, which contains `v^2`, vanishes.

## 2. Kummer pair and discriminant-square law

The discriminant of the quadratic factor in `(4)` is

```text
Delta_2=y^2(y^2-8vU)=-y^2 L/v^2             in k[v,v^-1,y]. (6)
```

After removing the displayed square `(y/v)^2`, its quadratic square class is
`-L`.  The resultant of the linear and quadratic factors is

```text
Res=-2v(y^2-9vU).                                     (7)
```

Since `S=y(y^2-9vU)`, the product-discriminant formula gives

```text
disc_X(E)=Delta_2 Res^2=-4LS^2.                        (8)
```

The same result follows directly:

```text
disc(LX^3+TX+2U)
  =-4LT^3-108L^2U^2
  =-4L(T^3+27LU^2)
  =-4LS^2.                                             (9)
```

Thus `(2)` is exactly the discriminant-square owner equation.  The sole odd
discriminant factor is `L`; `S^2` records projection collisions between
finite roots.  This is the same typed split used by the fixed Keller inverse
cubic in THM-3535, but `(1)` gives a flexible two-source-parameter packet
rather than asserting that its coefficients already descend to a planar
Keller target.

### 2.1 The dual visible marked-root cubic

The source coordinate `y` is itself a marked root of a second cubic, now
written directly in the visible packet coordinates `(T,S)`:

```text
Y^3-3TY+2S
 =(Y-y)[Y^2+yY+2(9vU-y^2)].                           (9a)
```

Its discriminant is

```text
disc_Y(Y^3-3TY+2S)
 =108(T^3-S^2)=-(54U)^2L,                             (9b)
```

and the residual quadratic has discriminant

```text
9(y^2-8vU)=-9L/v^2               in k[v,v^-1,y].      (9c)
```

Thus both marked-root presentations have the same Kummer square class
`-L`.  The first presentation `(3)` makes the escaping root `-1/v`
explicit; the dual presentation `(9a)` makes the finite source coordinate
`y` a root of a cubic whose coefficients are already among the packet
observables.  Any nonlinear projection that hides the marked sheet must
preserve this common resolvent sidecar, not merely the scalar cusp identity
`(2)`.

## 3. Every natural two-coordinate projection is obstructed

Write `Jac` for `Jac_(v,y)`.  Exact differentiation gives

```text
Jac(T,S)=54vU(U+vU_v),                                 (10)
v | Jac(L,T), Jac(L,U), Jac(L,S),                      (11)
Jac(T,U)=-2(yU_v+3UU_y),                              (12)
Jac(U,S)|_(v=0)
 =3y[y U_v(0,y)+3U(0,y)(d/dy)U(0,y)].                 (13)
```

Equations `(10)`, `(11)`, and `(13)` visibly exclude a nonzero constant.
For `(12)`, suppose `U=sum_(i=0)^d a_i(y)v^i`.  If `d>=1`, let `j` be the
largest index with `a_j'!=0`.  The highest `v`-degree involving a derivative
in `UU_y` has coefficient `a_d a_j'`, with no competing term; hence it must
vanish.  Descending gives `a_i'=0` for every `i`.  Equation `(12)` then
reduces to `-2yU_v`, which cannot be a nonzero constant.  If `d=0`, a
nonconstant `U(y)` makes `UU_y` have degree `2 deg(U)-1>0`, while constant
`U` gives zero.  This handles the final pair.

Therefore none of

```text
(L,T), (L,U), (L,S), (T,U), (T,S), (U,S)               (14)
```

is a planar Keller map, for any polynomial `U`.

## 4. An everywhere differential-rank-two positive packet

The obstruction in `(14)` is not a rank defect of the full packet.  Take

```text
U_*=1+y-y^2/2-(3/2)vy(y-3)                             (15)
```

and define `Z=(L,T,U_*,S):A^2->A^4` by `(1)`.  Let `M_ij` be its six
two-by-two Jacobian minors.  Exact Groebner reduction over `Q` gives

```text
(M_12,M_13,M_14,M_23,M_24,M_34)=(1) in Q[v,y].         (16)
```

Hence the differential of `Z` has rank two at every geometric affine point:
`Z` is differential-immersive.  This does not assert that `Z` is injective or
an algebraic immersion/embedding.  A second exact linear solve gives

```text
sum_(i<j) c_ij M_ij != 1                               (17)
```

for every six-tuple of constants `c_ij`.  Since the Jacobian of any constant-
linear projection `A^4->A^2` is such a constant combination, no
constant-linear projection of this packet is Keller.

The companion checks `(1)`--`(17)` with exact symbolic arithmetic and keeps
all truth gates active under optimized execution.  The Groebner and linear
systems are finite exact calculations, not numerical evidence.

### 4.1 An exact double packet fibre

The warning after `(16)` is active: the packet is not geometrically
injective.  In `Q[v,y]` consider the reduced degree-two scheme

```text
I=(27v-20y+22, y^2-3y+9/20).                         (17a)
```

The second generator has discriminant `36/5`, so over `Q(sqrt(5))` its two
distinct points are

```text
p_+=(8/27+4sqrt(5)/9, 3/2+3sqrt(5)/5),
p_-=(8/27-4sqrt(5)/9, 3/2-3sqrt(5)/5).               (17b)
```

Direct substitution, and independently Groebner reduction modulo `I`, give

```text
Z(p_+)=Z(p_-)=(-6724/3645,57/20,27/40,27/40)         (17c)
```

in `(L,T,U_*,S)` order.  Thus every target projection of `Z` retains this
collision.  This is a geometric collision over a quadratic extension, not a
claim that the displayed points are rational over `Q`.

There is no first-jet incompatibility between the two sheets.  With

```text
B_0=(9627/304384)T-(205687/1369728)U_*,
```

one has

```text
Jac(L,B_0)=1 on V(I).                                 (17d)
```

Hence the bounded no-go below is genuinely global; it is not obtained by
forcing incompatible tangent planes at the selected double fibre.

## 5. The nonlinear projection equation

The unit ideal `(16)` says an arbitrary source-polynomial combination of the
six minors can equal one.  It does **not** yet produce a legal projection.
For polynomials `A,B in k[Z_1,Z_2,Z_3,Z_4]`, the chain rule is

```text
Jac(A(Z),B(Z))
 =sum_(i<j) [A_iB_j-A_jB_i](Z) M_ij.                  (18)
```

The six coefficient functions in `(18)` are highly constrained:

1. they must be pullbacks of ambient polynomials `C_ij in k[Z_1,...,Z_4]`,
   not merely arbitrary elements of the source ring;
2. before pullback, the ambient tuple is the Pluecker coordinate tuple of
   the decomposable two-form `dA wedge dB`, so it obeys

   ```text
   c_12 c_34-c_13 c_24+c_14 c_23=0;                  (19)
   ```

3. more strongly, all six coefficients must arise simultaneously as
   `C_ij=A_iB_j-A_jB_i` for polynomial potentials `A,B`; descent, the
   Pluecker relation, and de Rham closedness/exactness alone do not imply this
   polynomial-potential realization; and
4. for this explicit packet the nontrivial fibre is already automatic by
   `(17c)`: every pair `(A(Z),B(Z))` identifies `p_+` and `p_-`.

This first converts the nonlinear-projection search into a typed syzygy
problem: find a descending, polynomial-potential-realizable decomposable
Bezout certificate for the minor ideal.  Constant coefficients fail by
`(17)`; arbitrary source coefficients exist by `(16)`; the gap between them
is exactly stated here and then closed for this packet in Section 5.4.
Conditions checked only after quotienting by the image ideal are necessary
but need not lift to such ambient `A,B`.

### 5.1 Total target degree at most two is empty

Constants in `A` and `B` do not affect their Jacobian.  The pullbacks of the
fourteen nonconstant target monomials of total degree at most two form the
exactly `14`-dimensional space

```text
W_2=<L,T,U_*,S,L^2,LT,LU,LS,T^2,TU_*,TS,U_*^2,U_*S,S^2>.
                                                               (19a)
```

Order this basis as displayed.  For each of its `91=binom(14,2)` unordered
pairs `i<j`, form the column

```text
Jac(f_i,f_j)=(f_i)_v(f_j)_y-(f_i)_y(f_j)_v.           (19b)
```

Rows are the coefficients of the `139` occurring source monomials `v^a y^b`,
including the constant row.  If `C_2` is this `139 by 91` rational matrix and
`C_2^+` is obtained by deleting the constant row, exact reduction gives

```text
rank_Q(C_2)=rank_Q(C_2^+)=67.                         (19c)
```

Thus the constant row lies in the span of the nonconstant rows.  Any bivector
coefficient vector that kills every nonconstant coefficient also kills the
constant coefficient.  This already allows the `91` wedge coefficients to
vary independently, forgetting their Pluecker decomposability and
polynomial-potential origin.

There is a second, broader descent relaxation.  Allow an arbitrary
coefficient `q(Z)` multiplying each original minor `M_ij`, where `q` ranges
over the fifteen target monomials of total degree at most two, including
`1`.  Its `139 by 90` matrix has exact full/deleted-constant ranks

```text
80,80.                                                (19d)
```

For quadratic `A,B`, every coefficient `A_iB_j-A_jB_i` in `(18)` lies in
this coefficient space.  Therefore, over every characteristic-zero field,

```text
deg_target(A),deg_target(B)<=2
   ==> Jac(A(Z),B(Z)) is not a nonzero constant.      (19e)
```

Both rank pairs were rederived modulo two large primes; an independent audit
reproduced `(19c)` over six further primes and checked the collision by both
radical substitution and the ideal `(17a)`.

### 5.2 Total target degree three is also empty

For an integral calculation put `W=2U_*` and `R=2S`.  The `34` nominal
nonconstant monomials of total target degree at most three have pullback
coefficient matrix of shape `111 by 34` and rank `33`.  Its unique relation
is the rescaled cusp equation

```text
R^2-4T^3-27LW^2=0.                                   (19f)
```

Retain all `34` functions despite this relation and form all
`561=binom(34,2)` pairwise brackets.  Their source-coefficient matrix has
shape `336 by 561`.  Exact integer reduction gives

```text
rank(C_3)=rank(C_3^+)=187,                            (19g)
```

where again `+` means deletion of the constant row.  Four independent odd-
prime reductions give the same pair.  This is an arbitrary-bivector
relaxation, stronger than imposing Pluecker decomposability.  Consequently

```text
deg_target(A),deg_target(B)<=3
   ==> Jac(A(Z),B(Z)) is not a nonzero constant.      (19h)
```

This finite exact statement did not by itself close total target degree four.
The all-degree singular-minor obstruction in Section 5.4 now supersedes that
degree boundary for this fixed packet, but not for any other packet.

### 5.3 Full image ideal, normalization, and conductor

The finite shadow after `(19h)` admits an exact completion.  Use the integral
target coordinates

```text
(L,T,W,R)=(L,T,2U_*,2S)
```

and let `phi:k[L,T,W,R]->k[v,y]` be their pullback.  In addition to

```text
F=R^2-4T^3-27LW^2,                                  (19i)
```

put

```text
G= 27LR^2-243LRT-81LRW+729LR-243LTW^2+486LTW
   -243LW^3+1215LW^2-1458LW
   +3R^2T+7R^2W-18R^2-15RT^2-12RTW+48RT
   +3RW^2-16RW+36R
   +21T^2W^2-138T^2W+192T^2
   +6TW^3-36TW^2+48TW+W^4-6W^3+12W^2-8W.           (19j)
```

Then, over every characteristic-zero field,

```text
ker(phi)=(F,G).                                      (19k)
```

In particular, the packet image is an absolutely integral `(3,4)` complete
intersection.  There are no hidden components supported on `W=0`, and the
Hilbert-function agreement previously seen only through degree seven is now
a theorem in every degree.  Indeed the top total-degree forms are

```text
F_top=-4T^3-27LW^2,
G_top=-W^2(243LT+243LW-21T^2-6TW-W^2),
```

and are coprime.  They are therefore a homogeneous regular sequence, so the
filtered complete-intersection lemma identifies the associated graded
Hilbert series with that of degrees `(3,4)`.

Here is a proof that also records the saturation boundary.  Let

```text
C(Y)=Y^3-3TY+R,
Q(Y)=(W+2T)Y^2-(2W+6T+R)Y+(W^2-2W+3R),              (19l)
H=Res_Y(C,Q),       A=partial G/partial L.
```

Exact elimination gives

```text
27W^2G+AF=27H.                                      (19m)
```

The polynomial `H` is monic of degree four in `R`.  At `(T,W)=(-3,-3)` it
specializes to

```text
R^4+63R^3+891R^2+4320R+216000,
```

whose reduction modulo `7` is the irreducible polynomial
`R^4+2R^2+R+1`.  Monicity therefore proves that `H` is irreducible over
`Q[T,W,R]`.  It is absolutely irreducible: the rational point
`(T,W,R)=(0,2,0)` of `H=0` has `H_R=64`, while a rational point on two
Galois-conjugate geometric components would be singular.

Let `J=(F,G)`.  The primitive polynomial `F`, linear in `L`, is prime, and
`G` is nonzero modulo `F`, so `F,G` is a regular sequence and `J` is
unmixed of height two.  After inverting `W`, `(19m)` eliminates `L` and
identifies `J_W` with the prime hypersurface ideal `(H)`.  A height-two
minimal prime of `J` containing `W` would have to be

```text
(W,R^2-4T^3),
```

but `F=0` and `G=270` at `(L,T,W,R)=(0,1,0,2)`.  Thus every minimal prime
avoids `W`; contraction from the unique prime after localization gives a
unique minimal prime `P`.  Unmixedness makes `J` `P`-primary.  If `f in P`,
then `W^n f in J` for some `n`, and `W notin P` forces `f in J`; hence
`J=P`.  Finally `Jac(T,W)(0,0)=-12`, so the pullback image has
transcendence degree two.  Its prime kernel also has height two and must
equal `J`.  Absolute irreducibility of `H` makes the same proof valid after
every characteristic-zero field extension.

#### 5.3.1 The smooth normalization is obtained by adjoining the marked root

The linear subresultant of `C,Q` is

```text
alpha Y+beta,

alpha=R^2+6RT+RW-12T^3-12T^2W+36T^2
      -5TW^2+28TW-W^3+6W^2,

beta=-3R^2+4RT^2+4RTW-18RT-4RW
     -6TW^2+12TW-2W^3+4W^2.                         (19n)
```

On the packet,

```text
beta(Z)=-y alpha(Z),       alpha(Z)!=0.               (19o)
```

Consequently `y=-beta/alpha` in the image function field, and then

```text
v=(y^2-T)/(3W)                                    on W alpha!=0. (19p)
```

This proves that `Z` has generic degree one.  Since `C(y)=0`, the ring

```text
B=k[L,T,W,R]/(F,G),       Bbar=B[y] subset k(v,y)     (19q)
```

is finite and birational over `B`.  Put `p=y^2-T=3vW`.  Exact lexicographic
elimination of `v` from

```text
p-3vW,
W-(2+2y-y^2)-3vy(3-y),
L-v^2(4vW-y^2)                                      (19r)
```

gives four generators for the kernel in `k[L,p,W,y]`.  For those four
generators, the ideal generated by them and all `2 by 2` minors of their
`4 by 4` Jacobian matrix is the unit ideal.  This is an exact rational
Groebner calculation with its universe and term order fixed in the
companion.  Hence `Bbar` is a smooth surface, in particular normal, and
`(19q)` is the normalization of `B`.  The presentation embeds in
`Q[v,y]` and has the rational smooth point `(L,p,W,y)=(0,0,2,0)`, so the
same Galois-component argument makes it geometrically integral; the
normalization statement therefore survives every characteristic-zero base
extension.

This normalization statement does **not** make the original affine plane
the normalization scheme.  The inclusion `Bbar subset k[v,y]` induces a
birational but nonfinite map `A^2->Spec(Bbar)`.  Indeed, normalizing the cusp
`F|_(W=0)=0` by `T=t^2,R=2t^3` gives

```text
G=6t^3(2t+3)
  [9L(t-3)^2+t^4-4t^3+8t+4].                        (19s)
```

Thus the `W=0` boundary consists set-theoretically of

```text
ell_0={(lambda,0,0,0)},
ell_1={(lambda,9/4,0,-27/4)},
Gamma_0={(-(t^4-4t^3+8t+4)/(9(t-3)^2),t^2,0,2t^3)}. (19t)
```

Both vertical lines are genuine limits from source infinity.  Set
`v=t^-1`, take either `y_b=tu` or `y_b=3+tu`, and impose `L=lambda`.  With

```text
E=W(t^-1,y_b)-(t y_b^2+t^3 lambda)/4,
```

the equation `E=0` has a unique formal solution near

```text
u(0)=-2/9, E_u(0)=9,       or       u(0)=-1/9, E_u(0)=-9.
```

On either branch one has exactly

```text
L=lambda,
W=(t y_b^2+t^3 lambda)/4,
T=(y_b^2-3t^2 lambda)/4,
R=-(y_b^3+9t^2 lambda y_b)/4.                        (19u)
```

The limits are respectively `ell_0` and `ell_1`, while `v=t^-1` has no
affine limit.  The valuative criterion therefore excludes properness and
finiteness of the plane parameterization.  This is the precise failure of
the tempting but false inference “birational parameterization equals
normalization.”

#### 5.3.2 Singular and multiple-fibre locus

On `W!=0`, the complete intersection is isomorphic to the resultant
hypersurface `H=0`.  Exact ideal containments

```text
alpha^2,beta^2 in (H,H_T,H_W,H_R),
beta^2,H_T^2,H_W^2,H_R^2 in (H,alpha)                (19v)
```

give

```text
Sing(X) intersect D(W)=V_X(alpha) intersect D(W).    (19w)
```

This is also the multiple-marked-root locus: there `C` and `Q` have at
least two common roots.  The boundary in `(19t)` is the exact failure
boundary for this simple `alpha` description.  Direct Jacobian reduction
gives

```text
Sing(X) intersect V(W)
 =ell_0 union {(-169/2916,9/4,0,-27/4)}.             (19x)
```

The isolated displayed point is the intersection of `ell_1` with
`Gamma_0` at `t=-3/2`; the rest of `ell_1` is smooth even though the
projected cubics `C,Q` alone have two common roots there.  The coordinate
`L`, lost when one divides by `W`, removes that spurious resultant signal.

There is a compact model for the main conductor.  For two common roots let
`s` and `p` be their sum and product.  Then

```text
T=(s^2-p)/3,       R=sp,

E_1=-3Ws+6W+5ps-6p-2s^3+6s^2=0,
E_2=3W^2-3Wp-6W+2p^2-2ps^2+9ps=0.                  (19y)
```

Eliminating `p` gives `-3P(s,W)`, where

```text
P=-16W^2s^2+36W^2s-24W^2
  +12Ws^4-79Ws^3+206Ws^2-228Ws+72W
  +4s^6-42s^5+126s^4-108s^3.                       (19z)
```

Its discriminant as a quadratic in `W` is

```text
disc_W(P)=(5s-6)^2 D_6(s),
D_6=16s^6-168s^5+601s^4-1000s^3+1048s^2-672s+144,
disc(D_6)=-26132324317743978381312!=0.               (19aa)
```

Thus the smooth projective model of this unordered-pair curve is
hyperelliptic of genus two.  Away from `5s-6=0`,

```text
p=(3Ws-6W+2s^3-6s^2)/(5s-6),
delta=s^2-4p
     =-3(4Ws-8W+s^3-6s^2)/(5s-6),                   (19ab)
```

and adjoining `sqrt(delta)` orders the two source sheets:

```text
y_+/-=(s+/-sqrt(delta))/2,
v_+/-=(y_+/-^2-T)/(3W).                              (19ac)
```

At `s=3,W=27/20`, one has `p=9/20` and `delta=36/5`, recovering `(17b)`.
The denominator in `(19ab)` is a genuine chart boundary: at `s=6/5`, the
equations `(19y)` give `W=-54/25` and
`p=-18/5 +/- (18/25)sqrt(-1)`.  There is also one triple fibre on `W!=0`:

```text
(L,T,W,R)=(1/54,1/2,-1,-1),                         (19ad)
```

where `Q` vanishes identically and the cubic `C` has three distinct roots.
Accordingly, `(19z)` is a generic unordered-pair model with explicitly
recorded exceptional charts, not a claim that every fibre is exactly
double.

#### 5.3.3 Exact normal-multiplier reformulation of the projection search

Use columns `(L,T,W,R)` and rows `(dF,dG,dA,dB)`.  For every target pair
`A,B`, the six oriented Hodge-normal identities combine to

```text
det d(F,G,A,B)(Z)=9 alpha(Z) Jac_(v,y)(A(Z),B(Z)).   (19ae)
```

The sign and factor `9` are fixed by direct exact substitution in the
`(L,W)` component.  Both gradient rows annihilate the two packet tangent
rows, and that component is nonzero, so two-dimensional normal-space linear
algebra gives the other five identities over `k(v,y)`; polynomial equality
extends them everywhere.  Since `(19k)` is the full kernel and
`alpha(Z)!=0`, a scalar
`c in k*` is the Jacobian of a target projection exactly when

```text
det d(F,G,A,B) = 9c alpha          modulo (F,G).      (19af)
```

Thus `alpha` is simultaneously the rational-inverse denominator, the main
conductor owner, and the required ambient normal multiplier.  This is a
strict reformulation, not merely a necessary quotient test; Pluecker and
polynomial-potential realizability remain encoded by the fact that the last
two rows really are `dA,dB`.

The exact companion verifies `(19i)`--`(19af)`, including the elimination,
smoothness, radical, boundary, and formal-chart gates, in normal and `-O`
execution.

### 5.4 All-degree descended-bivector obstruction

For `i<j`, let

```text
N_ij=det(dF,dG,e_i,e_j),                              (19ag)
```

with columns ordered `(L,T,W,R)`.  Direct exact pullback gives all six
oriented identities

```text
N_ij(Z)=9 alpha(Z) M_ij,                              (19ah)
```

where `M_ij=Jac_(v,y)(Z_i,Z_j)`.  These are the coefficientwise form of
`(19ae)`.

Now form the ambient ideal

```text
J_sing=(F,G,N_12,N_13,N_14,N_23,N_24,N_34).          (19ai)
```

An exact rational grevlex Groebner basis has thirteen elements, nine of
degree four and four of degree three.  Reduction gives a nonzero eighteen-
term normal form for `alpha`.  A compact quotient-dual witness is

```text
Phi(h)=[coefficient of L R^2 in NF_(J_sing)(h)],
Phi(J_sing)=0,                  Phi(alpha)=81/50.      (19aj)
```

Hence `alpha notin J_sing`.  Reductions modulo `1009`, `10007`, and
`1000003` independently reproduce the thirteen-element degree profile and a
nonzero eighteen-term remainder; the characteristic-zero conclusion is the
exact rational calculation.

Suppose, even in the relaxed universe, that arbitrary ambient polynomials
`C_ij` satisfied

```text
sum_(i<j) C_ij(Z) M_ij=c in k*.                       (19ak)
```

Multiplying by `9 alpha(Z)`, using `(19ah)`, and then using the full-kernel
identity `(19k)` would give

```text
sum C_ij N_ij-9c alpha in (F,G),
```

contradicting `(19aj)`.  Thus `(19ak)` has no solution in any target degree.
For an actual target pair `A,B`, its coefficients
`C_ij=A_iB_j-A_jB_i` are a special descended six-tuple, so in particular

```text
Jac_(v,y)(A(Z),B(Z)) notin k*                         (19al)
```

for every `A,B in k[L,T,W,R]`.  This is an all-degree theorem for the fixed
packet and is stronger than the bounded arbitrary-bivector ranks in Sections
5.1--5.2.  It does not constrain a different packet whose normal multiplier
does lie in its image singular-minor ideal.

### 5.5 Cited degree-cap transfer -- superseded for this packet

The source total degrees of `(L,T,U_*,S)` are `(6,4,3,5)`.  Hence target
polynomials of total degree at most `E` pull back to source polynomials of
degree at most `6E`.  If such a pair had nonzero constant Jacobian, `(17c)`
would make it a noninjective planar Keller pair.  The **CITED** reduced-height
bound of
[Guccione--Guccione--Horruitiner--Valqui](https://arxiv.org/abs/2204.14178),
as routed in `05-knowledge/reference/CORE-PAPERS.md`, therefore gives

```text
E>=18.                                               (19am)
```

At `E=18`, the same cited sub-`125` classification forces the target-reduced
source degree pair `(72,108)`, up to order.  In the raw packet-weight ledger,
the unique cap-`18` monomial of source weight `108` is `L^18`.  Since

```text
L_top=-12v^4y^2=-12(v^2y)^2,
```

the common-leading-base law forces the degree-`72` source top form to be a
scalar multiple of `(v^2y)^24`, equivalently `(L_top)^12`.  The direct
packet-leading-monomial ledger has `L^12` as the unique monomial with leading
exponent `v^48y^24`.  Before `(19al)` was known, this identified cap `18` as
the first numerical search boundary and required quotienting higher-weight
cancellations before treating the coefficient of `L^12` as invariant.  The
all-degree obstruction now proves that no cap survives; this ledger is
retained as an independent degree/leading-face hostile control, not as an
open packet lane.

## 6. Counterexample architecture and scope

The packet combines the two boundary models from THM-3554 and THM-3555:

```text
marked factor vX+1       -> one root, escaping at v=0;
quadratic factor         -> Kummer pair with square class -L;
disc=-4LS^2              -> odd L is the prospective infinity owner;
four-coordinate packet  -> no differential rank defect.              (20)
```

The selected collision `(17c)` means that any unit-Jacobian projection would
have been a genuine planar counterexample.  The exact nonmembership
`(19aj)` proves that no such projection exists, even after relaxing actual
polynomial potentials to arbitrary descended minor coefficients.  This does
not exclude other packets or planar counterexamples.  What is proved is the
owner factorization, all natural projection no-gos, the differential-rank-two
but geometrically noninjective packet, the prime `(3,4)` image ideal, its
smooth marked-root normalization and nonfinite boundary, the determinant
reformulation `(19af)`, and the all-degree obstruction `(19al)`.  No `JC(2)`
conclusion beyond closing this explicit construction is claimed.
