---
id: THM-2193
title: "Uniform rank-six safe-torus floor"
status: >
  PROVED + VERIFIED-EXACT. Every saturated rank-six relation lattice
  L subset Z^13 with no coordinate vector e_i has
  mu_(K_L)(K_L intersection J^13) >= 7^(-21), where
  J=[1/14,13/14]. The proof resolves the critical cancellation in
  THM-2190 by disintegrating each missing facet exactly and treating three
  arithmetic regimes. Consequently every zero-safe distinct positive
  thirteen-speed row has seven independent integer relations of height at
  most 78*7^21=43566577398496152546. In combination with THM-2190's first
  six relations, every twelve-speed deletion of such a zero-safe row has six
  independent relations of height at most
  900*(78*7^21)=39209919658646537291400. This does not prove LRC(14).
source: codex-2026-07-24-rank-six-safe-torus-floor
depends_on:
  - THM-2190-basis-safe-floor-and-height-500-rank-six-harvest
related:
  - THM-2185-rank-two-safe-cube-floor-and-height-500-rank-three-harvest
  - THM-2164-relative-packet-rank-harvesting
script: 04-computation/lrc14_uniform_rank_six_safe_torus_floor_thm2193.py
output: 05-knowledge/results/lrc14_uniform_rank_six_safe_torus_floor_thm2193.out
script_sha256: 8fd4d9ec70711e90c4b0c513465c12a6f212f9e6fafb3d4f9e69737112dfcbd5
output_sha256: 7a8d5d3418bfb744d71c18d99f6b16db5ca6940db2488db4143af701e6db0cf1
hash_basis: working-tree bytes (LF)
---

# THM-2193 -- uniform rank-six safe-torus floor

Put

```text
I={x in R/Z:||x||<1/14},
J=(R/Z)\I=[1/14,13/14],
p=mu(I)=1/7,
q=mu(J)=6/7.                                         (1)
```

For a saturated lattice `L subset Z^13`, let

```text
K_L={x in (R/Z)^13:l.x=0 mod 1 for every l in L}.     (2)
```

The new point is that the rank-six positivity in THM-2190 is uniform.
More precisely:

> **Uniform rank-six safe-torus theorem.** If
>
> ```text
> rank(L)=6,        e_i notin L for 1<=i<=13,         (3)
> ```
>
> then
>
> ```text
> mu_(K_L)(K_L intersection J^13)
>     >=1/343^7
>      =1/7^21.                                       (4)
> ```

No height or coefficient bound on `L` is assumed.

## 1. Coordinate cover and the critical deficit identity

Saturation makes `X=Z^13/L` a free character lattice of rank seven.
The thirteen coordinate classes

```text
c_i=e_i+L in X                                       (5)
```

are all nonzero by (3), and they generate `X`. Choose seven of them that
form a rational basis. Exactly as in THM-2190, their finite-index span
gives an integer `D>=1` and a finite surjective cover

```text
rho:(R/Z)^7 -> K_L                                   (6)
```

on which the selected coordinates and the six remaining coordinates pull
back to

```text
c_j(rho(y))=D y_j                      (1<=j<=7),
c_i(rho(y))=a_i.y                      (8<=i<=13),    (7)
```

where every `a_i` is a nonzero vector in `Z^7`. The cover preserves Haar
probability.

Write

```text
B={y:D y_k in J for every 1<=k<=7},
A_i={y:a_i.y in I}.                                  (8)
```

Then

```text
rho^(-1)(K_L intersection J^13)=B\union_(i=8)^13 A_i,
mu(B)=q^7.                                           (9)
```

Fix one extra vector `a`, choose `j` with `a_j!=0`, and put

```text
C_j={y:D y_k in J for every k!=j}.                   (10)
```

For fixed values of the other six coordinates, the nonzero character
`y_j -> a_j y_j` is Haar-uniform. Hence

```text
mu({a.y in I} intersection C_j)=p q^6.               (11)
```

Up to null endpoints, the part of (11) missing from `A intersection B`
is

```text
Delta(D,a,j)
 =mu{y:
       a.y in I,
       D y_j in I,
       D y_k in J for k!=j}.                          (12)
```

Thus the exact facet-deficit identity is

```text
mu(A intersection B)=p q^6-Delta(D,a,j).             (13)
```

The task is to bound (12) uniformly in `D` and `a`.

## 2. Exact disintegration of one facet

Set

```text
g=gcd(D,a_1,...,a_7),
h=D/g,
b=a/g.                                               (14)
```

Then

```text
gcd(h,b_1,...,b_7)=1,        b_j!=0.                  (15)
```

Disintegrate the multiplication-by-`D` map coordinatewise. For
`z=D y mod 1`, its fibre is represented by

```text
y=(m+z)/D,       m in (Z/DZ)^7,                      (16)
```

and

```text
a.y=(b.m+b.z)/h mod 1.                               (17)
```

The homomorphism

```text
(Z/DZ)^7 -> Z/hZ,       m -> b.m mod h               (18)
```

is surjective by (15), so all its fibres have the same size. Define

```text
N_h(u)=#{r in Z/hZ:||(r+u)/h||<1/14}.                (19)
```

Let `R_j` be the product box in `z`-space with its `j`-th coordinate in
`I` and every other coordinate in `J`. Its volume is

```text
mu(R_j)=p q^6.                                       (20)
```

Equations (16)--(19) give the exact formula

```text
Delta(D,a,j)
   =integral_(R_j) N_h(b.z)/h dz.                    (21)
```

All subsequent endpoint statements are essential, or equivalently may be
read after replacing `I,J` by their interiors. Endpoints have zero Haar
mass.

## 3. Three arithmetic regimes

### 3.1 Many residue phases: `h>=7`

The `h` residue points in (19) have unit spacing on a circle of circumference
`h`. The danger arc has length `h/7`. Except when an endpoint hits a residue,

```text
N_h(u)>=floor(h/7).                                  (22)
```

If `h=7s+r` with `0<=r<=6`, then `s>=1` and

```text
h<=7s+6<=13s.
```

Consequently

```text
floor(h/7)/h>=1/13.                                  (23)
```

Equations (20)--(23) imply

```text
Delta(D,a,j)>=p q^6/13.                              (24)
```

### 3.2 Few phases but a large coefficient:
`1<=h<=6` and `max_k |b_k|>=7`

For `h<=6`, the arc in (19) has length less than one residue spacing.
Therefore, almost everywhere,

```text
N_h(u)=1_{||u||<h/14}.                               (25)
```

Choose `ell` with

```text
B_0=|b_ell|=max_k |b_k|>=7.                          (26)
```

Fix the other six `z`-coordinates. As `z_ell` crosses its interval of
length

```text
lambda_ell=p if ell=j,
lambda_ell=q if ell!=j,                              (27)
```

the real variable `b_ell z_ell+constant` crosses an interval of length
`B_0 lambda_ell`. The periodic set

```text
{w:||w||<h/14}
```

has measure `h/7` in every interval of length one. Partitioning from the
left endpoint into complete unit intervals shows that the target preimage
in the `z_ell` interval has length at least

```text
h floor(B_0 lambda_ell)/(7B_0).                      (28)
```

If `ell=j`, write `s=floor(B_0/7)`. Since
`B_0<=7s+6<=13s`,

```text
floor(B_0/7)/B_0>=1/13.                              (29)
```

Using (21), (25), and the other-coordinate volume `q^6` gives

```text
Delta(D,a,j)
 >=q^6 floor(B_0/7)/(7B_0)
 >=q^6/91
 =p q^6/13.                                          (30)
```

If `ell!=j`, then

```text
floor(6B_0/7)/B_0>=6/13                              (31)
```

because `floor(6B_0/7)>=6B_0/7-1>=6B_0/13` for
`B_0>=7`. The other-coordinate volume is now `p q^5`, so

```text
Delta(D,a,j)
 >=p q^5 floor(6B_0/7)/(7B_0)
 >=p q^5(6/91)
 =p q^6/13.                                          (32)
```

Thus the same strong bound (24) holds throughout the large-coefficient
regime.

### 3.3 The finite local core:
`1<=h<=6` and `max_k |b_k|<=6`

Use the centered real lifts

```text
z_j in (-1/14,1/14),
z_k in (1/14,13/14) for k!=j.                        (33)
```

Their centre and coordinate radii are

```text
c_j=0,          r_j=1/14,
c_k=1/2,        r_k=3/7 for k!=j.                    (34)
```

Also

```text
sum_k |b_k|<=42.                                     (35)
```

At the centre,

```text
b.c=(1/2)sum_(k!=j)b_k.                              (36)
```

If the sum in (36) is even, `b.c` is an integer. The open
`l_infinity` cube of radius `1/686` around `c` remains in the box (33), and
the linear-form variation in that cube is less than

```text
42/686=3/49<1/14<=h/14.                              (37)
```

It is therefore contained in the target strip in (25).

Suppose instead that the sum in (36) is odd. It has a half-integer residue.
Put

```text
R=|b_j|/14+(3/7)sum_(k!=j)|b_k|.                     (38)
```

Here `b_j!=0`, and oddness forces at least one nonzero coefficient away from
`j`. Hence

```text
R>=1/14+3/7=1/2.                                    (39)
```

Move from `c` in the coordinatewise support direction that reduces the
half-integer residue, through the fraction

```text
theta=(1/2-h/98)/R.                                  (40)
```

This is legal because

```text
0<theta<=1-h/49<1.                                   (41)
```

Call the resulting point `z_*`. Then

```text
||b.z_*||=h/98,                                      (42)
```

its distance to every face of (33) is at least

```text
(1-theta)/14>=h/686,                                 (43)
```

and its remaining target-strip margin is

```text
h/14-h/98=3h/49.                                     (44)
```

For `||u||_infinity<1/686`, (35) gives

```text
|b.u|<42/686=3/49<=3h/49.                            (45)
```

Equations (43)--(45) show that the open cube of radius `1/686` around
`z_*` lies in both the box and the target strip. This is the same cube as in
the even case. Its side is `1/343`, so in either parity

```text
mu{z in R_j:||b.z||<h/14}
 >=1/343^7.                                          (46)
```

By (21) and (25), and because `h<=6`,

```text
Delta(D,a,j)
 >=1/(6*343^7).                                      (47)
```

The deliberately small cube in (46) is what makes the estimate independent
of all lattice coefficients.

## 4. Summing the six deficits

Apply (47), or the stronger (24), to each of the six extra characters in
(7), choosing a nonzero basis coordinate separately for each. The union
bound, (9), and (13) give

```text
mu_(K_L)(K_L intersection J^13)
 >=q^7-sum_(i=8)^13 mu(A_i intersection B)
 =q^7-6p q^6+sum_(i=8)^13 Delta_i.                   (48)
```

At the rank-six threshold,

```text
q^7=6p q^6.                                          (49)
```

The two critical terms cancel exactly, but the six facet deficits do not:

```text
mu_(K_L)(K_L intersection J^13)
 >=6/(6*343^7)
 =1/343^7
 =1/7^21.                                            (50)
```

This proves (4). Notice that overlaps among the six danger events can only
improve the union bound; no disjointness assumption occurs.

The source-target map is now explicit:

```text
rank-six torus
  -> seven coordinate clocks plus six extra characters
  -> six one-facet conditional deficits
  -> residue count / maximal-coefficient slice / bounded local cube
  -> uniform safe mass.                              (51)
```

The cover forgets the original coefficient heights, while `(h,b)` retains
exactly the arithmetic needed to decide which of the three estimates applies.

## 5. Why a degenerating family cannot evade the floor

There is also a qualitative compactness explanation. Give sublattices of
`Z^13` the pointwise membership topology inside `{0,1}^{Z^13}`. Any sequence
of saturated rank-six lattices has a subsequence whose membership indicators
converge. Its limit `L_infinity` is saturated and has rank at most six.
Moreover, `e_i notin L_infinity` if no lattice in the sequence contains
`e_i`.

The Haar Fourier transform of `K_L` is the membership indicator of `L`.
Thus the associated Haar measures converge weakly to Haar measure on
`K_(L_infinity)`. Every coordinate character remains nontrivial, so the
boundary of `J^13` has limiting Haar measure zero and safe mass is continuous
along the subsequence. Ranks at most five have THM-2190's positive basis
floor; rank six has (50). Hence no sequence of rank-six tori can have safe
mass tending to zero. This compactness argument explains the obstruction,
while (50) supplies the effective constant that compactness alone would not.

## 6. An explicit seventh-relation height

The uniform floor makes THM-2190's existential height effective. We record a
simple self-contained Jackson bound.

For `N>=2`, let

```text
P_N(t)=sum_(k=0)^(N-1) exp(2 pi i k t),
C_N=integral_(R/Z)|P_N(t)|^4 dt
   =N(2N^2+1)/3,
Q_N(t)=|P_N(t)|^4/C_N.                               (52)
```

Then `Q_N` is a nonnegative normalized kernel of degree `2N-2`. Put

```text
f=1_J,
f_N=Q_N*f,
eta_N=||f_N-f||_1.                                   (53)
```

For `x=||t|| in [0,1/2]`,

```text
|P_N(t)|<=min(N,1/(2x)),
C_N>=(2/3)N^3,
Q_N(t)<=min(3N/2,3/(32N^3 x^4)).                    (54)
```

Splitting at `x=1/(2N)` gives

```text
integral_(R/Z) ||t|| Q_N(t) dt
 <=3/(4N)-3/(8N^3)
 <3/(4N).                                            (55)
```

Since translating an interval indicator by `t` changes its `L1` norm by at
most `2||t||`,

```text
eta_N
 <=integral Q_N(t)||f(. -t)-f||_1 dt
 <3/(2N).                                            (56)
```

Choose

```text
N_*=39*7^21+1
   =21783288699248076274,
H_*=2N_*-2
   =78*7^21
   =43566577398496152546.                            (57)
```

Then

```text
26 eta_(N_*)<39/N_*<1/7^21.                          (58)
```

Let `v` be a distinct positive thirteen-speed row with zero safe measure,
and suppose

```text
dim_Q W_(H_*)(v)<=6.                                 (59)
```

Take the full saturated bounded-relation lattice

```text
L=W_(H_*)(v) intersection Z^13.                      (60)
```

Every Fourier resonance of the degree-`H_*` tensor `product_i f_(N_*)(x_i)`
lies in `L`, and conversely `L` is the annihilator of `K_L`. Hence its line
average and `K_L` average agree exactly. Product telescoping on the line and
on `K_L` gives

```text
mu_(K_L)(K_L intersection J^13)<=26 eta_(N_*).       (61)
```

If `rank(L)=6`, (50) contradicts (58)--(61). If `rank(L)<=5`,
THM-2190's basis-safe floor is larger still and gives the same contradiction.
Therefore

> **Explicit height-`H_*` rank-seven harvest.** Every zero-safe distinct
> positive thirteen-speed row satisfies
>
> ```text
> dim_Q W_(43566577398496152546)(v)>=7.               (62)
> ```

The constant is intentionally crude; it is the first completely explicit
height obtained from the uniform rank-six argument, not an optimized height.

Finally, combine a seventh relation of height `H_*` with THM-2190's six
independent relations of heights at most

```text
(105,105,178,204,262,450).                            (63)
```

The kernel of evaluation at any fixed coordinate has dimension at least six.
Using the first relation nonzero in that coordinate as a pivot produces six
independent deletion relations, each of height at most

```text
2*450*H_*
 =900H_*
 =39209919658646537291400.                           (64)
```

## 7. Scope

The theorem closes the uniform rank-six safe-torus question left open in
THM-2190 and makes its seventh-relation height effective. It does not:

1. optimize the floor `7^(-21)` or the height in (57);
2. provide a comparable safe-torus floor in rank seven;
3. turn the seven global relations into the full terminal certificate needed
   by the current LRC(14) proof graph; or
4. prove LRC(14).

Reproduction:

```bash
python3 04-computation/lrc14_uniform_rank_six_safe_torus_floor_thm2193.py
python3 -O 04-computation/lrc14_uniform_rank_six_safe_torus_floor_thm2193.py
```
