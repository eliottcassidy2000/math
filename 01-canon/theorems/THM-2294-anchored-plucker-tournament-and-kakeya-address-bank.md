---
id: THM-2294
title: "Anchored Plucker tournaments and the Kakeya address bank"
status: >
  PROVED + INDEPENDENTLY AUDITED + VERIFIED-EXACT + CITED KNOT APPLICATION.
  A rank-three real
  relation packet with a nonzero anchored minor has an intrinsic alternating
  pair observable: its anchored Plucker two-form. Its sign shadow, when
  tie-free, is a locally transitive tournament; if the packet annihilates a
  strictly positive vector, that tournament is also strong. Exact anchored
  minors together with one affine offset per coordinate reconstruct the full
  rank-three plane. For the THM-2284 packet, either one coordinate column is
  zero modulo thirteen, or at every depth n there are at least
  72*13^(2n-2) anchored relation addresses whose nine coefficients are all
  thirteen-units, of height at most 22143(13^n-1)/2. The raw line-union
  comparison gives 65*13^(2n-2); rank three sharpens 65 to 72, and seven
  parallel lines plus one transversal show that 72 is sharp for abstract
  rank-three packets. Quadratic-character orientation works at primes
  congruent to three modulo four but is symmetric, not a tournament, at
  thirteen. Every universal four-vertex Pfaffian tournament gauge has exactly
  one cyclic triangle and is therefore incompatible with an anchored
  chirotope tournament. On knots, directional catalytic imbalance is an
  honest pair observable but the first Brittenham--Hermiller example is a
  complete tie with positive symmetric bypass, so its tournament shadow
  cannot classify unknotting interaction. No LRC(14) profile is excluded.
source: codex-2026-07-25-anchored-plucker-kakeya-bank
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2284-thirteen-adic-anchored-rank-three-plucker-lift
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
related:
  - THM-2267-static-owner-coverage-is-flag-and-transition-holonomy-is-a-cut-kernel
  - THM-2278-two-shallow-proper-root-spectrum-and-gap-ancestry-activation
  - THM-2282-thirteen-adic-saturation-and-unit-anchored-minor
script: 04-computation/anchored_plucker_kakeya_bank_thm2294.py
output: 05-knowledge/results/anchored_plucker_kakeya_bank_thm2294.out
script_sha256: b11e393b1a18e0090d6539c97bc45fc4da84c22ed893f36472eef1c442f03900
output_sha256: 942beae70378c74ece3a7d03530abe33040333ebd780ed2a22693c2e06956e20
hash_basis: working-tree bytes (LF)
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number is not additive under connected sum, arXiv:2506.24088v2."
---

# THM-2294 -- the anchored Plucker object is a valued edge field

THM-2284 supplies three bounded integer relations and a three-by-three
minor, containing the shallow blocker `c_1`, which is nonzero modulo
thirteen. There are two tempting but inequivalent ways to read that minor:

```text
real sign:             an oriented-matroid / order-type shadow;
thirteen-adic unit:    an arithmetic address and Smith certificate.       (1)
```

Neither coordinate determines the other. The exact pair object retains both
and, after choosing a pivot, becomes an arrangement of affine lines. This
gives a lawful tournament shadow, an exact account of what it forgets, and a
finite-field line-avoidance theorem for the whole THM-2284 address bank.

The comparison with THM-2290 is deliberately hostile. Its Pfaffian
tournament orients matching signs; the tournament here orients projected
coefficient columns. The two constructions already give disjoint
four-vertex tournament classes.

## 1. Anchored contraction and exact reconstruction

Let `R` be a real `3 by N` matrix of row rank three, with labelled columns

```text
v_i in R^3,                         i in I.             (2)
```

For ordered distinct labels, write

```text
Delta_(ijk)=det(v_i,v_j,v_k).                          (3)
```

Fix an anchor `a` and an ordered pivot `(a,b,c)` with

```text
Delta_(abc)!=0.                                       (4)
```

The **anchored Plucker edge field** is

```text
h_(ij)=Delta_(aij)/Delta_(abc),        i,j!=a.         (5)
```

It is alternating, has `h_(bc)=1`, and is independent of the chosen row
basis. A different ordered pivot multiplies every `h_(ij)` by one common
nonzero scalar; over the reals its sign shadow therefore changes at most by
global reversal.

Define three normalized coordinates for every column:

```text
Z_i=Delta_(ibc)/Delta_(abc),
X_i=Delta_(aic)/Delta_(abc),
Y_i=Delta_(abi)/Delta_(abc).                          (6)
```

> **Anchored reconstruction theorem.** After the unique row-basis change
> sending `(v_a,v_b,v_c)` to the standard basis, the `i`th column is
>
> ```text
> (Z_i,X_i,Y_i)^T.                                    (7)
> ```
>
> Moreover
>
> ```text
> h_(ij)=X_iY_j-Y_iX_j,                               (8)
> ```
>
> and hence every four labels satisfy
>
> ```text
> h_(ij)h_(kl)-h_(ik)h_(jl)+h_(il)h_(jk)=0.          (9)
> ```

### Proof

Let `B=(v_a v_b v_c)`. Cramer's rule gives

```text
B^(-1)v_i
 =(
   Delta_(ibc)/Delta_(abc),
   Delta_(aic)/Delta_(abc),
   Delta_(abi)/Delta_(abc)
  )^T,                                                (10)
```

which is (7). Taking the determinant of the last two coordinates gives
(8), and the decomposable two-form in (8) satisfies (9). QED.

The edge field already recovers its two-dimensional contraction:

```text
X_i=h_(ic),                  Y_i=h_(bi).              (11)
```

It does not recover `Z_i`. The smallest literal witness has pivot columns
`e_1,e_2,e_3` and a fourth column equal respectively to

```text
(0,1,1)^T             or             (1,1,1)^T.      (12)
```

The complete anchored edge fields agree, but the two `Z` coordinates and
the two row planes differ. Thus the exact `h` field plus the node offset
ledger `(Z_i)` reconstructs the plane; the tournament of signs alone is two
quotients further away.

There is a useful Kakeya-style dual description. In the affine chart of row
combinations whose anchor coefficient is one, write the combination as
`(1,s,t)` in the normalized basis. Its coefficient at coordinate `i` is

```text
L_i(s,t)=Z_i+sX_i+tY_i.                              (13)
```

Consequently coefficient zero is an affine line in the `(s,t)` plane, with
normal direction `(X_i,Y_i)` and offset `Z_i`. Equation (8) detects
parallelism and cyclic order of the normals. It forgets the offsets which
place the lines. In this precise sense the tournament is a direction-only
Kakeya shadow and `(Z_i)` is its translation sidecar.

## 2. The exact real tournament shadow

Put

```text
u_i=(X_i,Y_i),                    i!=a.              (14)
```

If every `h_(ij)` on a vertex set `S` is nonzero, orient its pairs by

```text
i -> j             iff             h_(ij)>0.        (15)
```

This is intrinsic after choosing the orientation of the pivot; reversing
the pivot reverses every arc. When some `h_(ij)=0`, it remains an honest
tie: `u_i,u_j` are parallel, or one is zero.

> **Tournament theorem.**
>
> 1. Every tie-free tournament obtained from (15) is locally transitive:
>    the in-neighborhood and out-neighborhood of every vertex induce
>    transitive tournaments.
> 2. Suppose additionally that there is a strictly positive vector
>
> ```text
> w=(w_i)_(i in I),             w_i>0,
> Rw=0.                                               (16)
> ```
>
>    If the full anchored contraction on `I\{a}` is tie-free, its tournament
>    is strongly connected.
> 3. Without the tie-free hypothesis, (16) forces a positive circuit in the
>    contraction: either a zero vector, an antipodal pair, or a triple whose
>    origin lies in its positive triangle.

### Local transitivity

Choose polar angles `theta_i` for the nonzero vectors `u_i`. Then

```text
h_(ij)>0
 iff theta_j-theta_i mod 2pi lies in (0,pi).          (17)
```

For fixed `i`, its out-neighbors lie in one open semicircle. In increasing
angular order every arc among them points forward, because the angular
difference remains strictly less than `pi`. They are transitive. The same
argument applies to the opposite semicircle of in-neighbors.

### Positive dependence forces strength

Projecting (16) to the last two normalized coordinates gives

```text
sum_(i!=a) w_i u_i=0.                                (18)
```

If the tournament were not strong, its strongly connected components would
give a nontrivial dominance cut

```text
A -> B,                 A disjoint_union B=I\{a}.    (19)
```

Put `U_A=sum_(i in A)w_i u_i` and similarly `U_B`. By (18),
`U_B=-U_A`, and therefore

```text
0=det(U_A,U_B)
 =sum_(i in A,j in B) w_iw_j det(u_i,u_j)>0,         (20)
```

because every cross determinant has the positive orientation in (19).
This contradiction proves strength.

Finally, (18) puts zero in the convex hull of the `u_i`. A support-minimal
positive dependence in the plane has size at most three. Size one is a zero
vector, size two is an antipodal pair, and size three is a positive
triangle. In the tie-free case only the third alternative remains and its
three signs form a directed cycle. Hence the LRC packet has the exact
alternative

```text
anchored Plucker tie,
or a strong locally transitive tournament with a directed triangle.     (21)
```

For the THM-2284 relation matrix, (16) is not an extra hypothesis: its rows
annihilate the strictly positive scalar profile vector `w_*` of that theorem.

## 3. Why `chi_7` can orient and `chi_13` cannot

Let `p` be an odd prime and let `chi_p` be its quadratic character. For a
nonzero finite-field anchored minor,

```text
chi_p(h_(ji))=chi_p(-1)chi_p(h_(ij)).                 (22)
```

Thus the rule

```text
i -> j iff chi_p(h_(ij))=+1                          (23)
```

is antisymmetric exactly when

```text
chi_p(-1)=-1,
equivalently p=3 mod 4.                              (24)
```

This explains the sharp `chi_7`/`chi_13` boundary. At seven, a nonzero
quadratic-character edge is an orientation, up to global reversal under a
change of ordered pivot by a nonsquare factor. Row-basis change cancels in
the normalized ratio (5). At thirteen, `-1` is a square, so

```text
chi_13(h_(ji))=chi_13(h_(ij)).                       (25)
```

The result is an undirected square/nonsquare edge coloring, not a
tournament. More generally every multiplicative homomorphism
`F_13^* -> {+1,-1}` sends `-1` to `+1`, so no multiplicative binary
character repairs (25). A quartic character retains four phases and can
change sign under negation, but converting those phases to two arcs requires
an additional noncanonical sector selector.

The exact residue-valued edge field remains lawful. A tournament quotient
does not preserve whether a Plucker coordinate is a thirteen-unit: multiplying
one row of an integral packet by thirteen leaves every real determinant sign
unchanged while multiplying every maximal minor by thirteen and destroying
mod-thirteen row rank.

## 4. The all-depth Kakeya address bank

Apply the preceding algebra to THM-2284's three scalar relation rows

```text
rho_1,rho_2,rho_3 in Z^9                              (26)
```

with height bounds

```text
||rho_1||_infinity<=9841,
||rho_2||_infinity<=4921,
||rho_3||_infinity<=7381.                            (27)
```

Let `a=c_1`. THM-2284 gives a three-column minor containing `a` which is a
unit modulo thirteen. For every coordinate `i`, let the integral column and
its reduction be

```text
Vtilde_i=((rho_1)_i,(rho_2)_i,(rho_3)_i) in Z^3,
V_i=Vtilde_i mod 13 in F_13^3.                       (28)
```

Call `i` **dark** if `V_i=0`.

For `n>=1`, let

```text
A_n={
 lambda in (Z/13^n Z)^3:
 lambda.Vtilde_a=1 mod 13^n
}.                                                   (29)
```

The anchored minor makes `V_a` primitive, so

```text
|A_n|=13^(2n).                                       (30)
```

> **All-depth dark-column/address-bank dichotomy.**
>
> - If some coordinate is dark, every relation in the three-row packet has
>   coefficient divisible by thirteen there, so no packet word is
>   all-unit.
> - If no coordinate is dark, then for every `n>=1` at least
>
> ```text
> 72*13^(2n-2)                                       (31)
> ```
>
>   members of `A_n` give relation words whose nine coefficients are all
>   units modulo thirteen. Each has a centered integral representative of
>   height at most
>
> ```text
> 22143(13^n-1)/2.                                   (32)
> ```

### The base affine-plane count

On the `13^2`-point affine plane `A_1`, the zero condition at a nonanchor
coordinate is one of three things:

```text
V_i=0:                         the whole plane;
V_i in F_13^* V_a:            the empty set;
otherwise:                    one affine line of 13 points.              (33)
```

The first case is absent in the no-dark branch. There are at most eight
forbidden lines. The raw union bound covers at most

```text
8*13=104
```

points and leaves the direction-free comparator (still within the no-dark
branch)

```text
169-104=65.                                         (34)
```

Rank three gives seven more points. The quotient columns of the two other
pivot labels are independent modulo the anchor, so the forbidden family
contains at least two nonparallel line directions.

For completeness, first deduplicate coincident forbidden lines. Suppose the
resulting `m<=8` affine lines use at least two directions. Choose a maximal
parallel class of size `r`, with `1<=r<=m-1`. Its union has `rp` points.
Every remaining line is nonparallel to this class and meets its `r` distinct
lines in `r` distinct points, so it adds at most `p-r` new points. At
`p=13`,

```text
|union lines|
 <=rp+(m-r)(p-r)
 =mp-r(m-r)
 <=mp-(m-1)
 <=8*13-7
 =97.                                                 (35)
```

Therefore at least

```text
169-97=72                                            (36)
```

base addresses avoid every coefficient-zero line.

The constant is sharp for abstract rank-three packets. Take `V_a=e_1`,
seven further columns

```text
(-r,1,0)^T,                    r=0,...,6,
```

and the last column `(0,0,1)^T`. On `A_1`, their zero lines are

```text
s=0,...,6                     and                  t=0.                  (37)
```

Seven parallel lines plus one transversal cover exactly

```text
7*13+(13-7)=97
```

points. This sharpness witness is for the affine line lemma; no claim is
made that it also satisfies the positive-kernel constraints of an LRC row.

### Self-similar lifting and height

Reduction `A_n -> A_1` is surjective and every base point has exactly

```text
13^(2n-2)                                           (38)
```

lifts. A coefficient nonzero modulo thirteen remains a unit at every lift.
Equations (36) and (38) prove (31). The bank is therefore a union of at
least seventy-two full thirteen-adic cylinders, of constant relative density
`72/169` at every depth.

Choose centered integer representatives for the three combination
coefficients:

```text
|lambda_j|<=(13^n-1)/2.                             (39)
```

Equations (27) and (39) give

```text
||sum_j lambda_j rho_j||_infinity
 <=(9841+4921+7381)(13^n-1)/2
 =22143(13^n-1)/2,                                  (40)
```

which is (32). Every such sum is an exact integer relation, not merely a
modular word.

This theorem does not decide which branch occurs for a given LRC packet. A
dark coordinate is a precise stopping object, not a contradiction. In the
other branch, an all-unit coefficient bank still does not identify one
Fourier atom or one owner-ancestry return.

## 5. Pfaffian signs and chirotope signs are different tournaments

THM-2290 identifies the only universal tournament-gauge range for turning a
hafnian into a Pfaffian. At order four, write its six signs as

```text
s_(ij) in {+1,-1},                   i<j.            (41)
```

Universality up to one global sign says that, for some `epsilon`,

```text
s_12 s_34=-s_13 s_24=s_14 s_23=epsilon.             (42)
```

For a triangle `i<j<k`, its cyclic indicator is

```text
C_(ijk)
 =(1+s_(ij)s_(jk))(1-s_(ij)s_(ik))/4.               (43)
```

Put

```text
x=s_12,        y=s_13,        z=s_14.
```

Equation (42) gives

```text
s_34=epsilon*x,
s_24=-epsilon*y,
s_23=epsilon*z.                                     (44)
```

Substitution in the four instances of (43) gives the exact identity

```text
C_123+C_124+C_134+C_234=1.                          (45)
```

Thus every universal four-vertex Pfaffian gauge has exactly one cyclic
triangle. The fourth vertex must beat all three vertices of that triangle
or lose to all three: any mixed adjacency creates a second cyclic triangle.
Equivalently the gauge tournament is, up to converse,

```text
K_1 -> C_3.                                         (46)
```

It is not strong, and one of its in- or out-neighborhoods is cyclic, so it
is not locally transitive.

Conversely, (42) has sixteen sign solutions, while there are exactly

```text
4 choices of cyclic triple
*2 cyclic orientations
*2 choices of source/sink fourth vertex
=16                                                  (47)
```

labelled tournaments with exactly one cyclic triangle. Hence (45) is also
an exact tournament classification of the up-to-sign universal gauges.

An anchored real chirotope tournament is locally transitive by Section 2,
so it cannot be a universal `K_4` Pfaffian gauge. At larger even orders
THM-2290 proves that no universal Pfaffian tournament exists at all. The two
uses of Plucker signs are therefore distinct:

```text
Krenn/Pfaffian:   signs compensate the parity of matching monomials;
LRC/chirotope:    signs orient projected coefficient columns.            (48)
```

A graph-specific Pfaffian orientation may survive on restricted support,
but it must retain that support/cycle sidecar.

## 6. Knot interactions are not tournament edges

In the metric-monoid notation of THM-2176, define the directional
translation saving

```text
C_y(x)=ell(x)-d(x+y,y)>=0                            (49)
```

and the symmetric connected-sum defect

```text
sigma(x,y)=ell(x)+ell(y)-ell(x+y)>=0.                (50)
```

The antisymmetric comparison

```text
tau(x,y)=sign(C_y(x)-C_x(y))                         (51)
```

is an intrinsic, tie-aware tournament shadow. It preserves only which
summand benefits more from the other as a translator. It loses a common
positive saving and, more importantly for the known examples, the symmetric
geodesic bypass.

The first Brittenham--Hermiller pair is an exact hostile witness. Put

```text
K=T(2,7),                    Kbar=mirror(K),
X=K#Kbar.                                           (52)
```

THM-2176 records

```text
u(K)=u(Kbar)=3,
d_G(X,K)=d_G(X,Kbar)=3,
u(X)<=5.                                            (53)
```

Consequently

```text
C_Kbar(K)=C_K(Kbar)=0,
sigma(K,Kbar)=6-u(X)>=1.                            (54)
```

Edges from the unknot `U` also have both directional savings zero and
symmetric defect zero. Thus on the three vertices

```text
{U,K,Kbar}
```

the entire catalytic-imbalance tournament is tied, while the symmetric
interaction graph has the nonzero edge `{K,Kbar}`. No orientation of that
edge can be extracted from (53)--(54), or from this catalytic/symmetric
ledger, without an additional asymmetric sidecar.

There is also an all-pairs symmetry no-go. Let `mu` be an isometric monoid
involution, and let a real pair score `A` be both antisymmetric and
simultaneously `mu`-invariant:

```text
A(y,x)=-A(x,y),
A(mu(x),mu(y))=A(x,y).                              (54a)
```

Then every involution pair is forced to tie, because in one line

```text
A(x,mu(x))=A(mu(x),x)=-A(x,mu(x)),
```

so `A(x,mu(x))=0`. Knot mirroring is such an isometric connected-sum
involution. Therefore every **mirror-blind** antisymmetric score ties
`K` and `Kbar`; a strict orientation must retain a mirror-odd chirality
sidecar or choose an external gauge. This does not claim that every knot
score is mirror-blind.

The exact first pair ledger is instead

```text
(C_y(x), C_x(y), sigma(x,y)).                        (55)
```

Together with the two root lengths it reconstructs

```text
d(x+y,y),            d(x+y,x),            ell(x+y), (56)
```

and the directional bypasses are `sigma-C_y(x)` and `sigma-C_x(y)`.
Even (55) is only a pair compression. The operation-complete knot object
remains THM-2176's min-plus continuation kernel `P_x(a,b)`.

## 7. Connection contract and exact scope

For the LRC packet, the connection is:

```text
source:
  three bounded relation rows with a c_1-anchored unit minor;

target:
  the anchored alternating edge field h, its sign/tie oriented-matroid
  shadow, and the affine line arrangement L_i(s,t)=0;

map:
  columns v_i |-> Delta_(aij), then normalize at one pivot;

preserved by the tournament:
  the sign of every real anchored minor, anchored basis-extension versus
  tie, cyclic order of projected columns, local transitivity, and the
  positive-kernel strongness obstruction;

destroyed by the tournament:
  determinant magnitude, thirteen-adic valuation and angular component,
  affine offset Z_i, integral lattice content, exact Fourier frequency,
  root history, phase, and owner ancestry;

minimal reconstruction sidecars:
  the exact normalized edge values h_(ij), the offsets Z_i, the integral
  Plucker content/valuation ledger, and the original LRC labels;

cheapest hostile tests:
  the two columns in (12), row scaling by thirteen, chi_13(-1)=+1,
  the one-cycle K_4 Pfaffian gauge, and the dark-coordinate branch.       (57)
```

The verification companion checks the Cramer/Plucker identities on an exact
integer bank, the sharp `65 -> 72` affine-line counts and all-depth cylinder
numbers, a positive eight-vector strong locally transitive model, the
`chi_7/chi_13` antisymmetry boundary, and all `64` labelled four-vertex
tournament gauges.

Reproduce with

```bash
python3 04-computation/anchored_plucker_kakeya_bank_thm2294.py
python3 -O 04-computation/anchored_plucker_kakeya_bank_thm2294.py
```

The two transcripts are byte-identical to the stored output. The theorem
provides a new representation and a positive-density relation-address bank.
It does not choose the no-dark branch, land an exact Fourier atom, exclude
one of the `120` interior profiles, handle the `45` boundary/repeated rows,
or prove LRC(14). QED.
