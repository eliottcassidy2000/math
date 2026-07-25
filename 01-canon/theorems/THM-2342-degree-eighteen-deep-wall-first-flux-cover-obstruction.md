---
id: THM-2342
title: "Degree-eighteen deep-wall first-flux cover obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. On each of
  THM-2341's two genus-zero
  deep-wall spectral curves, the recovered first-flux function Z has a
  pole of order three at all three points at infinity. Hence the
  connected double cover T^2=Z has at least four branch points and
  positive genus. A rational Keller trajectory would lift to a
  nonconstant map from P^1 into this cover, which is impossible. Together
  with THM-2338 and THM-2341 this empties the complete deep common-root
  wall 126D=25B^2, 21W=-20BC. The off-wall degree-eighteen cone and
  JC(2) remain open.
source: codex-2026-07-25-deep-wall-first-flux-cover
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2338-degree-eighteen-deep-common-root-wall-hurwitz-quartet
  - THM-2341-degree-eighteen-deep-wall-local-genus-split
related:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
script: 04-computation/jc2_degree18_deep_wall_first_flux_cover_thm2342.py
output: 05-knowledge/results/jc2_degree18_deep_wall_first_flux_cover_thm2342.out
script_sha256: c33ff230410dcb7a470bf52a120ea7b56055bd1d06cf99489b917712bdc5e9e8
output_sha256: f33e9d884e40d8e652641cf2aae97c547070dec7fd286a4cea6ded764ab33612
hash_basis: working-tree bytes (LF)
---

# THM-2342 -- the first-flux cover empties the deep wall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2338 finds exactly three low-square-class points on the simultaneous
wall

```text
126D=25B^2,                 21W=-20BC.              (1)
```

THM-2341 eliminates the rational ratio by genus and proves that the two
quadratic ratios `t_+,t_-` have genus-zero spectral normalization. Those
two curves nevertheless cannot carry a Keller trajectory.

The coordinate forgotten by the spectral cubic is the first-flux square

```text
Z=T^2.                                              (2)
```

On each rational normalization, `Z` has three odd poles. The resulting
quadratic cover has positive genus, so a rational trajectory cannot lift
through (2). Consequently the entire wall (1) is empty in the genuine
degree-eighteen Keller branch.

## 1. The two survivors are nodal plane cubics

Normalize `B=1`. On (1), THM-2341 blows up `v=yz` and obtains

```text
Phi(z,y)
 =z^3
  +(16/964467)(1890+245y^2)z
  +(64/703096443)(11340y+183708C+539y^3)
 =0.                                               (3)
```

This is a plane cubic. At each of `t_+,t_-`, it is absolutely
irreducible and has the split node identified in THM-2341.

The homogeneous degree-three part of (3) is

```text
Phi_3(z,y)
 =z^3+(3920/964467)y^2z+(34496/703096443)y^3.      (4)
```

There is no infinity point with `y=0`, because (4) would then give
`z^3=0`. Put `h=z/y` on the infinity line. Its three points are the roots
of

```text
B_infinity(h)
 =h^3+(80/19683)h+704/14348907.                    (5)
```

They are distinct:

```text
Disc(B_infinity)
 =-94208/282429536481 !=0.                         (6)
```

Let `(r,z_0)` be the node and intersect (3) with a line through it:

```text
y=r+s,                    z=z_0+h s.               (7)
```

Because `(r,z_0)` is a double point and `Phi` is cubic,

```text
Phi(z_0+hs,r+s)
 =s^2(A(h)+B_infinity(h)s).                        (8)
```

Thus the third intersection is

```text
s=-A(h)/B_infinity(h),                             (9)
```

which is the standard line-through-node parameterization of the
normalization.

If `A` and `B_infinity` had a common root, the corresponding line would
make (8) vanish identically and would be a component of `Phi`, contrary
to THM-2332's absolute irreducibility. Therefore `A` is nonzero at all
three roots of (5). Equations (6) and (9) show that `s`, hence `y` and
`z`, has a simple pole at each of the three infinity parameters.

## 2. Recover the first-flux function

The inverse depressed-coordinate translation on (1) is

```text
u
 =yz-(49601160+1607445y^2)/(3*(-26040609))

 =(35y^2+1701yz+1080)/1701.                        (10)
```

At `alpha=0`, `B=1`, and `D=25/126`, THM-2262's first-flux numerator is

```text
N_2
 =45927u^2-58320u-5670uy^2+2160y^2
  +93312(25/126)-15552Cy+35y^4.                    (11)
```

Substitution of (10) gives the unexpectedly linear factor

```text
N_2
 =-(y/9)(
    139968C+560y^3+34020y^2z
    -413343yz^2+12960y).                           (12)
```

Hence the exact first-flux equation

```text
Z=-2N_2/(5103y)
```

becomes the regular function on the affine plane cubic

```text
Z
 =(2/45927)(
    139968C+560y^3+34020y^2z
    -413343yz^2+12960y).                           (13)
```

The cancellation of `y` in (12) is legitimate on the function field;
THM-2262 already excludes the separate constant wall `y=0`.

## 3. All three infinity valuations are odd

Along (7), with `s` tending to infinity and fixed direction `h`, the
degree-three part of (13) is

```text
(2/45927)s^3 L(h),

L(h)=560+34020h-413343h^2.                         (14)
```

There is no leading cancellation at any infinity point, because

```text
Res_h(B_infinity,L)
 =1024192512
 =2^12 3^6 7^3
 !=0.                                              (15)
```

At every root of `B_infinity`, equation (9) gives a simple pole for `s`,
and (14)--(15) give

```text
ord(Z)=-3.                                         (16)
```

Thus `Z` has an odd pole at each of the three distinct infinity points
of the rational normalization.

## 4. The full radical tower has positive genus

Let `C_tilde` be either genus-zero normalization from THM-2341 and let
`D_tilde` be the normalization of the quadratic cover

```text
D_tilde -> C_tilde,                 T^2=Z.          (17)
```

An odd valuation of `Z` is a branch point of (17). Equation (16) gives
three. It also proves that `Z` is not a square, so the cover is
connected.

The number `r` of odd valuations of a rational function on `P^1` is
even: reducing its principal divisor modulo two leaves total degree
zero. Therefore

```text
r>=4.                                              (18)
```

Riemann--Hurwitz for the connected double cover gives

```text
2g(D_tilde)-2=2(-2)+r,

g(D_tilde)=(r-2)/2>=1.                             (19)
```

Now suppose a degree-eighteen Keller trajectory survived at `t_+` or
`t_-`. Its nonconstant pair `(u,y) in C(x)^2` gives a nonconstant map

```text
P^1_x -> C_tilde.
```

The retained rational function `T in C(x)` and equation `T^2=Z` lift it
to

```text
P^1_x -> D_tilde.                                  (20)
```

Properness extends (20) across poles, and it is nonconstant because its
projection is nonconstant. But there is no nonconstant map from `P^1`
to a curve of genus at least one. This contradiction eliminates both
quadratic ratios.

Notice why "`Z` is nonsquare on `C_tilde`" alone would not suffice: a
nonsquare rational function can become a square after pullback along a
rational cover. The positive genus of the full cover (17), not
nonsquareness by itself, is the invariant obstruction.

## 5. Frontier effect

The exact deep-wall ledger is now closed:

```text
H=1 cyclic square class:               absent       THM-2335;

H_2 square class on the wall:          absent       THM-2338;

rational H_4 ratio t_0:                genus one    THM-2341;

quadratic H_4 ratios t_+,t_-:          flux-cover
                                      genus >=1     this theorem.      (21)
```

By THM-2332, every Keller survivor must enter one of the low-square-class
strata. Therefore no survivor lies on (1).

The structural lesson is a cover tower:

```text
P^1_x  --->  D_tilde  --degree 2-->  C_tilde
                                      --degree 3--> P^1_y.             (22)
```

Genus must be audited after restoring every radical sidecar, not only on
the first spectral quotient. The trigonal quotient can be rational while
the retained flux cover is not. This mechanism is reusable on other
Newton edges and at other terminal degrees.

The mixed and four-simple square-class strata away from the wall, other
degree-eighteen parameter walls, `JC(2)`, and `DC(2)` remain open.

## 6. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_deep_wall_first_flux_cover_thm2342.py
python3 -O 04-computation/jc2_degree18_deep_wall_first_flux_cover_thm2342.py
```

Both transcripts are byte-identical to the stored output. The companion
verifies the plane-cubic degree, the complete infinity polynomial
(5)--(6), the inverse depression (10), the first-flux factorization
(12)--(13), the leading flux direction (14), the exact resultant (15),
and the double-cover genus lower bound. The line-through-node,
irreducibility, valuation, lifting, and Riemann--Hurwitz arguments are
the mathematical proof above. No executable check uses Python `assert`.

The independent hostile audit rebuilt the infinity cubic discriminant,
inverse depressed-coordinate translation, first-flux factorization, and
resultant (15) using custom rational and Sylvester arithmetic. It also
checked that the three infinity intersections are transverse, that both
`y` and `z` have order `-1` there, and that a common zero of `A` and
`B_infinity` would force a line component. Finally, it audited the
function-field lifting argument even for a non-birational map
`P^1 -> C_tilde`, the parity/Riemann--Hurwitz step, the complete wall
ledger, both execution modes, and both hashes. QED.
