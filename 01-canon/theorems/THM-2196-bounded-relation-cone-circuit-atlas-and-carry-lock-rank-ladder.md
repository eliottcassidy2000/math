---
id: THM-2196
title: "Bounded relation-cone circuit atlas and the carry-lock rank ladder"
status: >
  PROVED (integral circuit atlas and explicit rank-six/rank-seven
  specializations) +
  PROVED RELATIVE TO CITED LOWER-DIMENSIONAL LRC through THM-609
  and THM-2053 (seven-label divergence, the carry-lock rank ladder, and
  projective finiteness of zero-Haar rows). Every bounded rank-r relation
  matrix with a strictly positive kernel vector has a finite semilinear atlas
  whose generators are primitive nonnegative circuits of support at most
  r+1 and explicitly bounded height. Applying THM-2190/2193 gives finite
  seven-carry/rank-six and six-carry/rank-seven atlases. THM-2053 turns every
  stabilized carry subspace of dimension at least two into a positive-Haar
  torus unless a new fixed independent relation survives. Consequently every
  infinite zero-Haar sequence has a common-dilation subsequence, and the
  primitive zero-Haar locus is finite. This does not prove LRC(14).
source: codex-2026-07-24-relation-carry-spectrum
depends_on:
  - THM-2190-basis-safe-floor-and-height-500-rank-six-harvest
  - THM-2193-uniform-rank-six-safe-torus-floor
  - THM-2053-rank-two-parameter-plane-geodesic-terminal
  - THM-609-base-good-region-floor
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - LRCUpTo13
related:
  - THM-2187-height-preserving-saturation-and-universal-radix-carrier
  - THM-2188-finite-phase-bank-and-pairwise-overlap-no-go
  - THM-2174-endpoint-phase-scale-obstruction
  - THM-2191-catalytic-localization-of-the-gordian-metric
  - THM-763-strict-finite-height-for-tight-lrc-instances
---

# THM-2196 -- bounded relation-cone circuit atlas and the carry-lock rank ladder

The bounded relation harvest does more than supply a finite list of linear
constraints. Its nonnegative kernel is a rational cone, so the speed row
itself has a finite relation-and-carry normal form. The missing coordinate is
then visible: it is the unbounded carry vector, not another bounded relation.

Throughout, for a finite speed set `A`, write `G_A` for its `1/14`-safe set
and `mu(G_A)` for Haar measure in `R/Z`.

## 1. Integral circuit atlas

Let

```text
B in Z^(r x 13),          rank(B)=r,          ||B||_infinity<=H,
C_B={x in R_(>=0)^13:B x=0},                 d=13-r,              (1)
```

and suppose `ker(B)` contains a strictly positive vector. Define

```text
U_r(H)=max_(0<=k<=r) (sqrt(k) H)^k,                         (2)
```

where the `k=0` term is `1`.

> **Integral circuit-atlas theorem.** There is a finite collection of
> simplicial charts
>
> ```text
> sigma=cone(u_1,...,u_d)                                  (3)
> ```
>
> covering `C_B`, where the `u_j` are linearly independent primitive
> nonnegative integral extreme-ray generators. For every chart there is a
> finite residue set
>
> ```text
> R_sigma subset C_B intersection Z^13                     (4)
> ```
>
> such that every integral point `x` of the chart has an exact expression
>
> ```text
> x=b+sum_(j=1)^d n_j u_j,
> b in R_sigma,                  n_j in Z_(>=0).             (5)
> ```
>
> Every generator and residue obeys
>
> ```text
> |supp(u_j)|<=r+1,             ||u_j||_infinity<=U_r(H),
> ||b||_infinity<d U_r(H).                                  (6)
> ```

The cover need not be disjoint. A fixed pulling order and half-open face
convention makes the chart choice deterministic if a normal form is wanted.

### Extreme rays are bounded circuits

Fix an extreme ray `R=R_(>=0)u` of `C_B` and normalize `u` so that its
nonzero coordinates are positive. Put

```text
S=supp(u),                  s=|S|,                  k=rank(B_S).   (7)
```

Because `u_S in ker(B_S)`, one has `k<=s-1`. If `k<=s-2`, choose
`z in ker(B_S)` not proportional to `u_S` and extend it by zero off `S`.
For sufficiently small positive `epsilon`, both

```text
u+epsilon z,             u-epsilon z                         (8)
```

are nonnegative kernel vectors and are not on the same ray. Their midpoint
is `u`, contradicting extremality. Hence

```text
k=s-1,                     s=k+1<=r+1.                       (9)
```

Choose `k` independent rows of `B_S`; call the resulting
`k x (k+1)` matrix `A`. These rows span the row space of `B_S`, so their
kernel is one-dimensional. For `k>=1`, the signed maximal-cofactor vector

```text
w_j=(-1)^j det(A with column j deleted)                       (10)
```

is a nonzero integral kernel generator. All its coordinates are nonzero
because the kernel vector has support exactly `S`; one global sign makes
them positive. Dividing by their gcd gives the primitive generator `u`.
Each cofactor matrix has `k` rows of Euclidean norm at most `sqrt(k)H`.
Hadamard's inequality therefore gives

```text
|w_j|<=(sqrt(k)H)^k<=U_r(H).                                 (11)
```

Division by the gcd only decreases the coordinates. When `k=0`, `S` is a
singleton and the primitive ray is a coordinate vector. This proves every
claim in (6) concerning a generator.

If the relation rows have individual height bounds
`h_1<=...<=h_r`, the same proof retains their anisotropy. A `k`-minor uses
some `k` rows, so

```text
||u||_infinity
 <=max_(0<=k<=r) k^(k/2) product_(j=r-k+1)^r h_j,            (11a)
```

with the `k=0` term again equal to one. This profile bound is often much
sharper than replacing every `h_j` by `H`.

### Existing-ray triangulation

The cone is pointed because it lies in the nonnegative orthant. A strictly
positive kernel vector has a relative neighborhood in `ker(B)` that remains
positive, so

```text
dim(C_B)=dim(ker(B))=d.                                      (12)
```

Intersect `C_B` with `sum_i x_i=1`. The result is a compact rational
polytope of dimension `d-1`; its vertices are precisely the normalized
extreme rays. A pulling triangulation introduces no new vertices. Every
maximal simplex has `d` affinely independent normalized rays, equivalently
`d` linearly independent primitive ray generators. Coning these simplices
from zero gives the finite cover (3).

### Integral parallelepiped residues

In one chart every `x` has a unique expansion

```text
x=sum_j alpha_j u_j,              alpha_j>=0.                (13)
```

For integral `x`, set

```text
n_j=floor(alpha_j),
b=x-sum_j n_j u_j=sum_j {alpha_j}u_j.                        (14)
```

The first formula makes `b` integral; the second puts it in `C_B` and in the
half-open fundamental parallelepiped

```text
Pi_sigma={sum_j theta_j u_j:0<=theta_j<1}.                   (15)
```

Coordinatewise,

```text
0<=b_i<sum_j (u_j)_i<=d U_r(H).                              (16)
```

Thus `R_sigma=Pi_sigma intersection Z^13` is finite, and (5)--(6) follow.
The strict middle inequality is harmless even in a coordinate absent from
some rays: because the chart spans `ker(B)` and that kernel contains a
strictly positive vector, every coordinate occurs in at least one chart ray.

## 2. The explicit rank-six LRC atlas

Let `v in Z_(>0)^13` have pairwise distinct coordinates and
`mu(G_v)=0`. THM-2190 supplies six independent relation rows with respective
height bounds

```text
(105,105,178,204,262,450).                                  (17)
```

Use them as the rows of `B`. Then `rank(B)=6`, `Bv=0`, and `v` is strictly
positive, so Section 1 applies with `d=7`. Consequently `v` belongs to one
of finitely many seven-carry charts, over the finite collection of possible
matrices `B`,

```text
v=b+sum_(j=1)^7 n_j u_j,                                    (18)
```

with

```text
|supp(u_j)|<=7,
||u_j||_infinity
 <=6^3 105^2 178 204 262 450
  =10,195,213,482,720,000,
||b||_infinity
 <71,366,494,379,040,000.                                   (19)
```

The maximum in (11a) is its `k=6` term for this profile. The coarser
`U_6(450)=1,793,613,375,000,000,000` remains valid.

THM-2193 supplies the explicit seventh-relation height

```text
H_7=78 7^21=43,566,577,398,496,152,546.                     (19a)
```

Appending that row gives a finite rank-seven six-carry atlas. Its primitive
rays have support at most eight. A convenient wholly integral wrapper for
the anisotropic cofactor height is

```text
V_7=7^4 105^2 178 204 262 450 H_7
   =4,937,284,759,496,105,853,936,801,653,089,320,000,       (19b)
```

and every rank-seven chart residue has height `<6V_7`. The sharper real
Hadamard factor in (11a) is `7^(7/2)` in place of `7^4`.

## 3. Seven-label divergence in every stabilized atlas

> **Stabilized divergence lemma.** Assume the lower-dimensional
> `LRCUpTo13` citation node used by THM-609. Suppose an infinite
> pairwise-distinct sequence of positive, repetition-free thirteen-speed rows
> has a fixed representation
>
> ```text
> v^(N)=c+sum_(j=1)^k n_j^(N)u_j,                           (20)
> ```
>
> where `c,u_j in Z_(>=0)^13`, and every carry is either constant or tends to
> infinity. If
>
> ```text
> mu(G_(v^(N)))=0,
> D=union_(n_j^(N) -> infinity) supp(u_j),                  (21)
> ```
>
> then `|D|>=7`.

At least one carry diverges, since otherwise (20) would make all rows equal.
If `i` lies outside `D`, then `v_i^(N)` is fixed. If `i` lies in `D`,
nonnegativity gives

```text
v_i^(N) -> infinity.                                       (22)
```

Suppose `q=|D|<=6`. The complementary coordinate list `E` is fixed,
positive, and repetition-free, with `7<=|E|<=12`. THM-609 gives

```text
m_E=mu(G_E)>0.                                             (23)
```

Let `r_E` be the finite number of positive-length components of `G_E`.
The simultaneous tail `F_N={v_i^(N):i in D}` satisfies

```text
sum_(w in F_N) 1/w -> 0
  <(7-q)m_E/(sqrt(2)r_E)                                  (24)
```

for all sufficiently large `N`. THM-735 applies to the fixed body and all
`q<=6` tails at once, with no separation hypothesis among the tails, and
gives `mu(G_(v^(N)))>0`, a contradiction.

For the rank-six application, order every row and select one of THM-2190's
independent relation six-tuples. Only finitely many matrices have the
profile (17), and each has finitely many charts and residues. Pigeonhole
extraction fixes `B`, a chart, and its residue. Successive extraction makes
each carry constant or divergent, so the lemma applies. The same argument
applies to any finite bounded relation-matrix family, in particular the
explicit rank-seven family supplied by THM-2193.

## 4. The two-plane terminal forces a carry-lock rank ladder

For a rational subspace `S subset Q^13`, put

```text
Lambda_S=S^perp intersection Z^13,
K_S={x in T^13:a.x=0 for every a in Lambda_S}.               (25)
```

This is the connected subtorus whose character relations are exactly
`Lambda_S`.

> **Positive-subspace terminal.** Assume `LRCUpTo13`. If `S` has dimension
> at least two and contains a positive integer row, then
>
> ```text
> Haar_(K_S)(K_S intersection [1/14,13/14]^13)>0.            (26)
> ```

Choose a primitive positive row `v in S intersection Z^13` and an integral
`z in S` independent of it. The rational plane

```text
U=span_Q{v,z}                                               (27)
```

contains a primitive positive row. THM-2053 therefore gives a point of
`K_U` whose thirteen coordinate clearances are at least `1/13`. Since

```text
U subset S              implies              K_U subset K_S, (28)
```

this point lies strictly inside the `1/14`-safe cube in `K_S`. Its relative
neighborhood has positive Haar measure, proving (26).

### Carry lock or dilation

Fix a relation matrix `B`, stabilize one of Section 1's charts and residues,
and absorb every constant carry into `c`. The remaining sequence has the
form

```text
v^(N)=c+sum_(j=1)^k n_j^(N)u_j,       n_j^(N)->infinity.     (29)
```

Put

```text
S=span_Q{c,u_1,...,u_k}.                                   (30)
```

If `dim(S)=1`, all rows in (29) are positive integer multiples of the unique
primitive positive generator of `S intersection Z^13`. This is the common
dilation alternative.

Suppose `dim(S)>=2`. If every `a in Z^13\Lambda_S` obeyed
`a.v^(N)!=0` for all sufficiently large `N`, then the Haar measures on the
Kronecker lines generated by `v^(N)` would converge weakly to Haar measure on
`K_S`. Indeed, their Fourier coefficients are

```text
1_(a.v^(N)=0),                                             (31)
```

and the proposed limit has coefficient `1_(a in Lambda_S)`. Pointwise
convergence on characters implies weak convergence. Every coordinate
character is nonzero on `K_S`, because `S` contains the positive row
`v^(N)`. Hence the safe-cube boundary has limiting Haar measure zero.
Equation (26) would force

```text
lim_N mu(G_(v^(N)))>0,                                     (32)
```

contrary to a zero-Haar sequence.

Therefore some fixed `a notin Lambda_S` satisfies `a.v^(N)=0` infinitely
often. Pass to that subsequence. Every row in the span of `B` annihilates
`c,u_1,...,u_k` separately, so it belongs to `Lambda_S`. Thus `a` is
independent of `B`: it is a genuine rank lift.

### Projective finiteness

Start with THM-2193. Its seventh-relation height is uniform, so an infinite
zero-Haar sequence has a subsequence with one fixed rank-seven relation
matrix. Stabilize its finite circuit atlas and apply the preceding
dichotomy. A sequence of primitive rows cannot take the dilation branch,
because a rational ray has only one primitive positive lattice generator.
It therefore acquires a fixed eighth relation.

Repeat with the augmented matrix. Each rank lift is fixed after passing to a
subsequence, so Section 1 applies again even though no uniform height is
known for the new row. The rank rises successively until at most twelve. A
rank-twelve positive kernel is one-dimensional and again contains only one
primitive positive row. This contradicts an infinite pairwise-distinct
primitive sequence. Consequently:

> **Projective zero-Haar finiteness.** Relative to `LRCUpTo13`, there are
> only finitely many primitive positive, repetition-free thirteen-speed rows
> with zero Haar-safe mass. Equivalently, every infinite zero-Haar sequence
> has an infinite subsequence consisting of positive integer dilations of one
> primitive row.

This is qualitative. THM-763 remains the effective finite-height theorem;
the rank-ladder proof neither improves its ceiling nor enumerates the finite
locus.

### Explicit proper-support one-carry terminal

The rank-two mechanism also gives a sharp constructive carry gate. Let

```text
c in Z^13,       u in Z_(>=0)^13\{0},
D=supp(u) proper subset {1,...,13},
c_l>0 for l outside D,             C=max_l |c_l|.            (33)
```

Then `c+nu` has positive Haar-safe mass for every integer

```text
n>91C.                                                       (34)
```

Choose `j outside D`, put `e=c_j`, and choose `i in D` minimizing
`c_i/u_i`. Restrict the two-phase forms `c_l t+u_l beta` to

```text
(t,beta)=s(u_i,e-c_i).                                      (35)
```

The induced positive integer speeds are

```text
a_l=c_lu_i+u_l(e-c_i)
   =u_i u_l(c_l/u_l-c_i/u_i)+eu_l       for l in D,          (36)
```

and `a_l=c_lu_i>0` off `D`. Moreover `a_i=eu_i=a_j`, so there
are at most twelve distinct speeds. `LRCUpTo13` supplies `s` with every
clearance at least `1/13`. Put

```text
t_0=u_i s,                 beta_0=(e-c_i)s.                  (37)
```

Among the `n` solutions of `nt=beta_0` in `T`, choose `t_n` nearest `t_0`.
Then `||t_n-t_0||<=1/(2n)`, and

```text
||(c_l+nu_l)t_n||
 >=||c_l t_0+u_l beta_0||-|c_l| ||t_n-t_0||
 >=1/13-C/(2n)>1/14.                                       (38)
```

Thus `t_n` is a strict safe time. The threshold is independent of the size
of `u`.

## 5. Boundary and exact loss

The atlas is lossless only while the carry vector is retained. Forgetting
the carries leaves a finite carrier but destroys endpoint current, phase
scale, and safe measure, exactly as THM-2174 and THM-2188 warn.

The theorem handles arbitrary changing cores by subsequence stabilization
and removes every infinitary projective escape: only common dilation can
survive indefinitely. It does not:

1. enumerate or discharge the finite primitive zero-Haar locus;
2. improve THM-763's effective height ceiling;
3. bound the new relations produced by the carry-lock argument;
4. preserve weak safe-set nonemptiness after forgetting the carries; or
5. prove LRC(14).

The surviving object is therefore a genuinely graded coefficient monoid
with a finite residue sidecar. Independent phase growth forces positivity;
failure of phase growth is an exact additional relation. This is the precise
relation-and-carry spectrum promised by the bounded-rank harvest. QED.
