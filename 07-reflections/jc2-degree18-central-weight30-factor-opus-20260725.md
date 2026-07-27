---
source: codex-2026-07-25-jc2-central-weight30-factor
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED; SUPERSEDED AS
  A CLOSURE ROUTE BY THM-2345. On
  D=25B^2/126, remove the
  universal y^2 from the degree-eighteen spectral branch discriminant.
  The repeated-root resultant of the residual degree-ten polynomial is
  exactly a nonzero scalar times R^6*S^3*T. The primitive weight-thirty
  factor T is reconstructed explicitly with 21 integer monomials. This
  makes the last central T=0 wall concrete. The later independently
  audited `jc2-degree18-central-twall-closure-opus-20260725` packet now
  excludes that wall and the full central ratio; this factorization packet
  alone does not prove JC(2).
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
related:
  - THM-2345-degree-eighteen-common-root-wall-saturation
  - jc2-degree18-central-wall-s0-closure-opus-20260725
  - jc2-degree18-central-wall-r0-closure-opus-20260725
  - jc2-degree18-central-twall-closure-opus-20260725
---

# The last degree-eighteen central wall has an explicit weight-thirty equation

> **CANONICAL SUPERSESSION.** THM-2345 closes the full common-root wall.
> This packet remains the exact `R^6 S^3 T` factorization and
> weight-thirty equation, not the current closure authority.

## 1. Exact central resultant

Use THM-2297's normalized spectral cubic

```text
G(u,y)
 =-5878656Wy-26040609u^3+49601160Bu^2+1607445u^2y^2
  -20995200B^2u-2857680Buy^2-52907904Du-138915uy^4
  +777600B^2y^2+33592320BD-5598720BCy+78120By^4
  +1959552Dy^2-435456Cy^3+1127y^6.                 (1)
```

The weights are

```text
wt(y)=1,  wt(u)=2,  wt(B,C,D,W)=(2,3,4,5),
```

and `G` has weight six. On the central ratio

```text
D=25B^2/126,                                       (2)
```

put

```text
Delta(y)=Disc_u G(u,y).
```

Direct cubic-discriminant arithmetic gives

```text
y^2 divides Delta,

delta_bar=Delta/y^2 has degree 10 and weight 10.   (3)
```

Consequently its repeated-root resultant has weight `10*9=90`.
The exact factorization is

```text
Res_y(delta_bar,d delta_bar/dy)
 =K R^6 S^3 T,                       K in Q*,       (4)
```

where

```text
R=20BC+21W,                                         (5)

S=2888B^5+108864B^2C^2+571536BCW+750141W^2,        (6)
```

and `T` is the primitive weight-thirty polynomial below. This proves the
multiplicities which the earlier central-wall routing used only through
its radical:

```text
R exponent 6,             S exponent 3,             T exponent 1. (7)
```

## 2. The 21-term primitive factor

Write a row `(i,j,k):q` for the term `q B^i C^j W^k`. Then

```text
(0,0,6):   1857584148891020160556640625
(0,5,3):   -720893165611681749052800000
(0,10,0):  187775146320922282243129344
(1,1,5):   6184255303140290720859375000
(1,6,2):   -1511892036737051357606400000
(2,2,4):   8180469698416791719343750000
(2,7,1):   -1038920491424993478174720000
(3,3,3):   5592912523479745046100000000
(3,8,0):   -224453394609952422813696000
(4,4,2):   2217787332403775132700000000
(5,0,4):   12754325524369432500000000
(5,5,1):   550162720634152586284800000
(6,1,3):   31485168419060160000000000
(6,6,0):   73571138739313427520000000
(7,2,2):   25958526361650288000000000
(8,3,1):   7389030553435238400000000
(9,4,0):   206537825241216000000000
(10,0,2):  23073416987040000000000
(11,1,1):  40577760518400000000000
(12,2,0):  17795598912000000000000
(15,0,0):  256000000000000000.                     (8)
```

Every row satisfies `2i+3j+5k=30`, and the gcd of all 21 coefficients
is one. Thus (8) fixes `T` without a residual scalar ambiguity.

## 3. Why the reconstruction is an identity

The exact companion uses no symbolic-algebra package.

First set `B=1`. A weight-ninety polynomial then has

```text
degree_C<=30,                  degree_W<=18.         (9)
```

For every point of the tensor grid

```text
0<=C<=30,                     0<=W<=18,             (10)
```

the companion:

1. forms the cubic discriminant in (3) by the five-term cubic formula;
2. removes the checked universal `y^2`;
3. computes the degree-ten resultant by exact rational Euclidean
   recursion; and
4. interpolates first in `W` and then in `C`.

All `31*19` values are exact fractions. Coefficient comparison on the
whole rectangle proves (4) on `B=1`, including every zero coefficient.
Weighted homogeneity then lifts the identity globally: for any surviving
pair `(j,k)`, the exponent

```text
i=(90-3j-5k)/2
```

is uniquely determined. Therefore setting `B=1` does not identify two
different weight-ninety monomials.

The archived modular scout independently confirms at three large primes
the profile

```text
degree_W 18 = R^6 * S^3 * squarefree residual of degree_W 6,
```

and its optional reconstruction routine, manually replayed at those same
three primes, matches the normalized coefficients of (8). The stored
modular transcript does not certify a larger CRT bank. The rational tensor
identity, rather than this modular scout, is the proof over characteristic
zero.

## 4. A smaller quotient model for the geometric attack

On `B!=0`, introduce the invariant coordinates

```text
x=C^2/B^3,             q=W^2/B^5,             z=CW/B^4,

z^2=xq.                                             (11)
```

Because every term in (8) has `j` and `k` of the same parity, define

```text
A(x,q)
 =sum_(j,k even) t_(i,j,k) x^(j/2)q^(k/2),

B_1(x,q)
 =sum_(j,k odd) t_(i,j,k) x^((j-1)/2)q^((k-1)/2).  (12)
```

The first sum has 13 terms and the second has eight. Then

```text
T/B^15=A(x,q)+z B_1(x,q).                          (13)
```

Thus away from `B_1=0`, eliminating the selected sign lift makes the last
wall birational to its exact plane image

```text
A(x,q)^2-xq B_1(x,q)^2=0.                          (14)
```

Equation (14) is the cheapest next geometric object. It preserves the
weighted-projective orbit while removing the simultaneous sign
`(C,W)->(-C,-W)`. It does not preserve which lift of that sign quotient
the Keller trajectory uses; the first-flux square and Keller one-form
remain necessary sidecars.

## 5. Frontier effect and scope

The two companion closures already prove that the central faces `S=0`
and `R=0` are empty in the genuine nonsplit polynomial branch. Equations
(4)--(8) therefore reduce the entire remaining central problem to

```text
B!=0,                 R!=0,                 S!=0,                 T=0.
                                                                  (15)
```

At the stage of this identity, this was a major reduction in explicitness,
not a closure. The later singular-atlas packet
`jc2-degree18-central-twall-closure-opus-20260725` supplies the missing
normalization/genus argument and closes (15). The factorization (4) by
itself does not:

```text
make the spectral normalization rational;

retain the square root from the first flux;

exclude a positive-genus or singular component;

or prove the planar Jacobian or two-variable Dixmier conjecture.   (16)
```

That decisive test has now been executed: the parameter curve has exactly
16 nodal and six cuspidal singular points, and every corresponding spectral
normalization has positive genus. The next target has moved to the other
weighted components of the degree-eighteen repeated-branch discriminant.

## 6. Exact reproduction

Characteristic-zero identity:

```text
04-computation/jc2_degree18_central_t_factor_exact_probe.py

05-knowledge/results/jc2_degree18_central_t_factor_exact_probe.out
```

Normal, optimized (`-O`), and stored transcripts match. LF-byte SHA-256:

```text
script  5dc5ebae4e09a930ff81a977989383b7f5461bbb786f60c4cfe5b26175d4a3c3
output  8fb7687d0a0f185a698de3fca607fc06ed261fa6ab9987b2c87af0123c8b82b7
```

Independent modular scout:

```text
04-computation/jc2_degree18_central_t0_modular_scout.py

05-knowledge/results/jc2_degree18_central_t0_modular_scout.out
```
