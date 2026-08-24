---
id: THM-3959
title: "Centered degree-five root maps are scalar or fold a ramification arm"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the centered
  trace-zero linear-color binary-cubic grammar with deg A=3 and repeated-root
  degree five, THM-3941 leaves seven collision-free non-root-regular pole
  carriers. Exact color division makes three carriers empty and leaves five
  explicit polynomial-color families. One endpoint is a literal scalar value
  in an actual normal maximal cubic order. Every nonscalar survivor has an
  irreducible one-place implicit component occurring to exponent one in the
  order discriminant and a forced two- or three-address target fibre.
  Rank-three finite flatness coalesces those addresses at one source
  ramification point, so the source arm is non-unibranch and THM-3920 forbids
  an affine-plane Keller open. This closes the stated centered degree-five
  grammar, not arbitrary root gauges, nonlinear color planes, degree at least
  six, or JC(2).
source: jc-cohn3709 + jc-zero-debt-lift / post-THM-3941 degree-five color closure, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc_zero_debt_lift and
  jc_extra_debt_local, 2026-08-24). Both audits independently reconstructed
  all seven trace forms; every color resultant and exceptional seam; the
  saturated J, M, N, and P systems; the denominator-complement L/O duality;
  all five survivor incidences, coefficient ideals, implicit discriminants,
  and address fibres; the scalar maximal-order R1+S2 gate; and the
  rank-three source non-unibranch conclusion. The P row was explicitly
  checked in the centered gauge: omitting its constant shift gives nonzero
  trace and is invalid. A hostile replay corrected the address-discriminant
  scalar from 1/27 to 1/9 and replaced a multivariate-gcd heuristic on J by
  a genuine saturation certificate. A final completeness audit restored the
  independent triple-pole modulus k in row M; the full pk!=0 packet has
  saturated color ideal <1>, so the row remains empty. All repairs preceded
  promotion. Normal and optimized 105-gate runs byte-match the frozen output;
  hashes and documentation checks pass.
depends_on:
  - THM-3801-cubic-etale-normalization-nonmonogenic-and-companion-sheet-gate
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3929-root-regular-one-place-linear-color-cubic-is-monogenic
  - THM-3941-all-degree-centered-cubic-pole-carrier-routing
related:
  - THM-3933-centered-degree-three-root-map-pole-partition-octic-nonentry
  - THM-3938-centered-degree-four-root-map-nonentry
  - THM-3952-minimal-mobius-internal-split-carriers-are-four-critical-colors-and-nonentry
script: 04-computation/jc2_centered_degree_five_root_map_thm3959.py
output: 05-knowledge/results/jc2_centered_degree_five_root_map_thm3959.out
script_sha256: 6cb3e9f43114124d416a80de2e5f7de628b2944b699b586c637069fdb521b8ac
output_sha256: 69c474ed97dcfd1a080ed9ed9b0820f864e20a2795c9290561c83dfa1ddeedc3
semantic_sha256: 2a8a5aab0520c4edeab61b1717d2a6b9a4bf7261ef64dc9a4fb0a6ac253b8808
hash_basis: raw LF bytes
---

# THM-3959 -- degree five has seven doors, and five fold the same arm

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Consider

```text
Phi=a(A)U^3+C U^2V+c(A)UV^2+d(A)V^3.                     (1)
```

Let an irreducible one-place discriminant component have normalization
`A1_u`, and assume

```text
A=A(u), C=C(u) in k[u],       deg A=3,
deg(t:P1_u -> P1_[U:V])=5,    t=U/V,                      (2)
```

with centered primitive incidence

```text
a(A)t^3-c(A)t-2d(A)=0.                                   (3)
```

Thus `Tr_(k(u)/k(A))(t)=0`. We also impose the coefficient and normalization
gates

```text
(a,C,c,d)=k[A,C],                 k(A,C)=k(u).             (4)
```

The conclusion is that no such component occurs in the finite normalization
attached to a planar Keller map. The centered root gauge and linear `C`
coefficient plane are hypotheses: translating a root generally mixes the
coefficient slots of `(1)`.

## 1. The seven complete trace forms

THM-3941 proves that degree five has exactly seven collision-free
non-root-regular pole carriers. Solving the lower Laurent-jet trace equations,
and using diagonal root/color scaling to normalize the infinity jet when it
is present, gives

```text
row  carrier                 A(u)          t(u)

J    m=0; C3(5)             u^3           p/u+q/u^2+r/u^4+1/u^5

K    m=0; C2(5)             u^3+u^2       p/u-3q/(2u^2)+q/u^3
                                            -5/(2u^4)+1/u^5

L    m=1; C3(4)             u^3           u+p/u+q/u^2+r/u^4

M    m=1; C2xC2(1,3)        u^3-3u        u+p/(u-1)+q/(u+1)
                                            +k/[2(u+1)^2]+k/(u+1)^3

N    m=2; C2(3)             u^3+u^2       u^2+pu+(p-1)/3+q/u
                                            -3k/(2u^2)+k/u^3

O    m=4; C3(1)             u^3           u^4+pu^2+qu+r/u

P    m=4; C2(1)             u^3+u^2       u^4+4u^3/3+pu^2+qu
                                            +(1/9-p/3+q/3)+r/u.
                                                                  (5)
```

The exact-pole conditions used below are

```text
L,O: r!=0;       M:pk!=0;      N:k!=0;       P:r!=0.     (6)
```

For example, if `v=1/u` over `A=u^3+u^2`, then

```text
Tr(h+z1v+z2v^2+z3v^3+z4v^4+z5v^5)
 = [3A^2h+A(2z2+3z3)+2z4+5z5]/A^2.                      (7)
```

This forces both cancellation rows in K. The same trace gives the constants
and principal-part ratios in N and P. In particular, after the surviving P
parameters are inserted, its constant is `(9q+2)/27`. Omitting that term
would give trace `-(9q+2)/9`; it is not a centered representative.

The polynomial row `m=5` is the separate root-regular exit and is already
closed by THM-3929.

## 2. Exact color resultants and seam exhaustion

Eliminate `u` from `A-A(u)` and the denominator-cleared equation `T-t(u)`.
The trace-zero resultant supplies `a,c,d` in `(3)`, and the repeated-root
derivative reconstructs the color:

```text
C(u)=-(3a(A(u))t(u)^2+c(A(u)))/(2t(u)).                  (8)
```

The exact resultants of the numerator of `t` and the numerator in `(8)` are

```text
J:  3^4 p^2(p-qr)^6,
K:  2^21 p^2(2p+5q)^6(2p+9q+27)^3,
L:  3^5 r^12(pq-r)^6,
M: -8549802417586176 k^9p^3 E_M^6,
N:  1184595334580404224 k^9 E_N1^3 E_N2^6,
O:  3^5 r^6(pq-r)^6,
P: -109418989131512359209 r^3 E_P1^3 E_P2^6,             (9)
```

where

```text
E_M =kp+kq+6k+6p^2+8pq+24p+2q^2+16q+24,
E_N1=243k+12p+54q-4,
E_N2=9k+4p^2+6pq-8p-6q+4,
E_P1=18p-54q-243r-14,
E_P2=9p^2-9pq-6p+3q+9r+1.                               (10)
```

Resultant vanishing is only necessary, so every factor in `(9)` is divided
exactly.

### 2.1 J: two families

Write the numerator of `t` as

```text
N=uX+Y,       X=pu^3+r,       Y=qu^3+1.                 (11)
```

The color is polynomial exactly when `N` divides `XY`. Since `N(0)=1`, `u`
is a unit modulo `N`, and `Y=-uX mod N`; equivalently `N|X^2`. Exact
division gives only

```text
(p,r)=(0,0), q arbitrary;        or
(p,q)=(r^4,r^3), r!=0.                                  (12)
```

For completeness, on `p=0,q!=0` the remainder is

```text
3r^2[(6qr^2-r^5)u^2+(4qr-5r^4)u+q-4r^3]/q^5.           (13)
```

Saturating its three coefficients by `r!=0` gives the unit ideal; at `q=0`
the remainder is `3/r^10`. On the other seam `p=qr`, with `qr!=0`, the final
remainder is `-3(q-r^3)/r^13`. Thus `(12)` is exhaustive.

### 2.2 K: one isolated point

The `p=0` seam has a three-coefficient univariate ideal of gcd one; its
degree-drop endpoint `p=q=0` also fails. The seam `2p+9q+27=0` has a
two-coefficient ideal of gcd one, and its exceptional endpoint `p=0,q=-3`
fails directly. On `2p+5q=0`, division leaves

```text
-2^16(28q+125)/5^13.                                    (14)
```

Hence the unique K point is

```text
p=625/56,                   q=-125/28.                   (15)
```

### 2.3 L and O: a denominator-complement duality

Both mixed C3 rows have the same numerator

```text
u^5+pu^3+qu^2+r,                                        (16)
```

but denominators `u^4` and `u`. Their only live resultant seam is `r=pq`,
where `(16)=(u^3+q)(u^2+p)`. Cancellation leaves the common denominator
`2(u^2+p)` and the respective remainders

```text
L: -3p^6(pu-q),               O: 3p^3(pu-q).             (17)
```

Since `r=pq!=0`, neither can vanish. Thus L and O are uniformly empty. The
two rows are not unrelated failed cases: exchanging which complementary
power of `u` is placed at infinity preserves their numerator and the final
linear obstruction.

### 2.4 M: the conic is empty after the color obligations

The infinity coefficient has already been normalized to one; the triple-pole
coefficient `k` is an independent modulus and cannot also be normalized away.
The only live seam is `E_M=0`, since `pk!=0`. Reducing the color numerator by
its denominator gives five coefficient obligations. The exact saturated ideal
is

```text
<E_M, five remainder coefficients, zpk-1> = <1>
    in k[z,p,q,k].                                        (18)
```

Thus no affine point of the conic passes color division. This is an exact
Groebner certificate, not a numerical conic search.

### 2.5 N and P: the same C2 branch factor at two infinity orders

On either N seam `E_N1=0` or `E_N2=0`, saturating the exact remainder
coefficients by `k!=0` puts `(3p-4)^2` or `(3p-4)^3` in the ideal. Conversely
substitution divides exactly. Hence

```text
N: p=4/3,             k=-2(9q+2)/81,       q!=-2/9.      (19)
```

Likewise the two P seams, saturated by `r!=0`, contain `(9p-1)^2` and
`(9p-1)^3`; exact substitution gives

```text
P: p=1/9,             r=-2(9q+2)/81,       q!=-2/9.      (20)
```

The same factor `9A+9q+2` will control both address packets. This is the C2
analogue of the L/O denominator-complement pairing, except that here both
complementary infinity orders survive.

## 3. The five polynomial-color survivors

Put `J_i=Res_u(A-A(u),C-C(u))` for each displayed parametrization.

### 3.1 J1: the free C3 family

Let `L_1=qA+1`. Then

```text
t=(qu^3+1)/u^5,
C=-3u^10(qu^3+1)/2,
a=A^5,                  c=0,                 d=L_1^3/2,
Disc(Phi)=-2L_1^3 J_1.                                  (21)
```

For `q!=0`, `J_1` has `C`-degree three and total degree thirteen. All three
roots of `u^3=-1/q` map to `(-1/q,0)` and are distinct because their
discriminant is `-27/q^2`.

At `q=0`, however, `d=1/2` is a literal scalar value. This is the unique
scalar seam and is treated in Section 4.

### 3.2 J2: the factored C3 family

For `r!=0`, put

```text
L_2=Ar^3+1,             Q_2=r^2u^2-ru+1.
```

Then

```text
t=(ru+1)^2Q_2/u^5,
C=-3u^10Q_2(r^2u^2+3ru+1)/2,
a=A^5,       c=3A^2rL_2^2,       d=L_2^4/2,
Disc(Phi)=-2L_2^4 J_2.                                  (22)
```

Here `J_2` has `C`-degree three and total degree fourteen. The identity

```text
(ru+1)Q_2=r^3u^3+1                                      (23)
```

puts both roots of `Q_2` over `(-r^-3,0)`; they are distinct because
`Disc(Q_2)=-3r^2`. The third point of that cubic `A`-fibre is `u=-1/r`,
where `Q_2=3` and `C=9/(2r^10)!=0`; hence the packet has exactly two
addresses.

### 3.3 K: the isolated C2 family

Set

```text
L_K=125A-28,       Q_K=25u^2+35u+14,
R_K=75u^7+385u^6+602u^5+280u^4-140u^3-140u^2+16.
```

The survivor `(15)` is

```text
t=(5u-2)^2Q_K/(56u^5),             C=-u^5Q_KR_K/14,
a=8A^5,
c=L_K^2(70A^2-35A+4)/392,
d=L_K^4/43904,
Disc(Phi)=-L_K^4 J_K/10976.                              (24)
```

The implicit prime has `C`-degree three and total degree fourteen. Moreover

```text
(5u-2)Q_K=125(u^3+u^2)-28,       Disc(Q_K)=-175,          (25)
```

so its address packet has exactly two points over `(28/125,0)`.
Indeed, the third point in that cubic `A`-fibre is `u=2/5`, where
`Q_K=32`, `R_K=-1024/3125`, and therefore `C!=0`.

### 3.4 N and P: the paired C2 family

Put

```text
L=9A+9q+2.                                               (26)
```

For N, equations `(19)` give

```text
t=(3u-1)(3u+2)L(A(u))/(81u^3),
C=-4u^3(3u+2)L(A(u))(9u^4+30u^3+24u^2-4)/3,
a=216A^3,
c=8(6A-1)(27A-4)L^2/243,
d=4(27A-4)^2L^3/19683,
Disc(Phi)=-16(27A-4)^2L^3 J_N/19683.                    (27)
```

The prime `J_N` has `C`-degree three and total degree eleven.

For P, divide the incidence and color by the harmless common scalar `729`.
Equations `(20)` become

```text
t=(3u-1)(3u+2)L(A(u))/(81u),
C=-9u(3u+2)(9u^2+6u-4)L(A(u))/2,
a=729A,
c=-(27A-4)L^2/9,
d=(27A-4)^2L^3/1458,
Disc(Phi)=-2(27A-4)^2L^3 J_P/729.                       (28)
```

Here `J_P` has `C`-degree three and total degree seven.

Both address fibres are the roots of

```text
u^3+u^2+q+2/9=0.                                        (29)
```

Its discriminant is

```text
-(9q+2)(27q+10)/9.                                      (30)
```

The first exceptional value `q=-2/9` is precisely `k=r=0`, which deletes the
required pole. At the other value,

```text
q=-10/27:       (29)=(3u-1)(3u+2)^2/27,                 (31)
```

so there are two distinct addresses `u=1/3,-2/3`. At every other admissible
parameter there are three. Thus neither N nor P has a one-address escape.

## 4. Coefficient ideals and the scalar endpoint

Every nonscalar survivor passes the coefficient-ideal gate `(4)`:

```text
J1: gcd(A,qA+1)=1,
J2: gcd(A,Ar^3+1)=1,
K:  L_K(0)=-28,
N:  d(0)=64(9q+2)^3/19683,
P:  d(0)=8(9q+2)^3/729.                                 (32)
```

The last two are nonzero exactly away from the pole-deletion value.

At the J1 endpoint `q=0`, exact specialization gives

```text
J_0=(27A^10+8C^3)/8,               Disc(Phi_0)=-2J_0.    (33)
```

The binomial `J_0` is irreducible because `gcd(10,3)=1`. The generic binary
cubic `A^5T^3+CT^2+1/2` is primitive and linear in `C`: any factorization
over `k(A)[C,T]` would have a `C`-independent factor dividing both `T^2` and
`A^5T^3+1/2`, whose gcd is one. Thus it is irreducible over `k(A,C)` and the
finite-free cubic order is a domain. Its only height-one discriminant prime
occurs to exponent one. An
order-to-maximal-order index changes discriminant valuation by an even
integer, hence the order is already maximal at `J_0`; it is etale and normal
at every other height-one prime. The binary-cubic order is finite free of
rank three over the regular ring `k[A,C]`, hence Cohen--Macaulay and `S2`.
Thus `R1+S2` makes the whole order normal. Therefore `d=1/2` is a scalar
value in the actual maximal normalization, not in an auxiliary nonmaximal order.
THM-3801 excludes a nontrivial Keller plane atlas.

## 5. Every nonscalar survivor is a genuine folded source arm

Each `J_i` in `(21)--(28)` occurs to exponent exactly one in the displayed
order discriminant. The other factors depend only on `A`, so none contains
the `C`-degree-three prime `J_i`. Since an order/maximal-order correction is a
square, maximalization cannot remove an odd valuation. Thus `J_i` supports a
genuine reduced tame residue-degree-one `(2,1)` ramification prime `E_i` in
the maximal cubic normalization.

The implicit relations are irreducible and the displayed parametrizations
are birational. First, `d!=0`, and a factor of the cubic in
`k(A)[C,T]`, which is linear in `C`, would have a `C`-independent factor
dividing both `T^2` and `aT^3+cT+d`; this is impossible. Second, the repeated
root is not in `k(A)`: at a finite pole its order modulo the inertia index is

```text
J1,J2: 5 mod 3;       K:5 mod 2;       N:3 mod 2;       P:1 mod 2. (34)
```

Since `[k(u):k(A)]=3` is prime, `k(A,t)=k(u)`. At the generic exact-double
root, Euclid's algorithm applied to the cubic and its derivative recovers
`t` from `A,C`; hence `k(A,C)=k(u)`. Consequently `E_i^nu=A1_u`, and `J_i`
is not a parametrization power.

Now let `B` be the maximal normalization of the target plane in this cubic
field. It is normal and finite over the regular two-dimensional target,
hence Cohen--Macaulay and finite flat of rank three. Every point of `E_i` is
non-etale and consumes geometric fibre length at least two. A rank-three
fibre cannot contain two distinct points of `E_i`.

Sections 3.1--3.4 exhibit two or three distinct points of `E_i^nu` over one
target point. They must therefore coalesce to a single source point of
`E_i`. The source ramification prime has two or three normalization branches
at that point. THM-3920 says every irreducible boundary curve in a normal
affine surface containing an `A2` open is rational and unibranch. Since a
Keller open must delete every ramification prime, no nonscalar survivor can
belong to a Keller affine-plane completion.

Combining this with the scalar endpoint and THM-3929's root-regular row closes
the full centered degree-five grammar.

## 6. Scope and structural handoff

Degrees three, four, and five are now closed in this centered linear-color
root-map lane. Degree five adds two reusable structural observations:

1. complementary placements of one numerator between finite and infinity
   poles can have the same final color obstruction, as in L/O; and
2. the surviving families have a two-type address dichotomy. In J1 and N/P,
   `C(u)` retains a nonconstant base factor `L(A(u))`, so a whole `A`-fibre
   collapses. In J2 and K, a base factor splits as `L(A(u))=ell(u)Q(u)` and
   `C(u)` retains the complementary quadratic `Q`, forcing exactly its two
   addresses. The only constant-base-factor endpoint is J1 at `q=0`, and it
   is exactly the scalar seam.

The theorem does not close degree at least six, an uncentered root gauge,
another coefficient plane, a nonlinear color slot, an intrinsically
nonmonogenic cubic completion, or JC(2). **QED.**

## Reproduction

```bash
python3 04-computation/jc2_centered_degree_five_root_map_thm3959.py
python3 -O 04-computation/jc2_centered_degree_five_root_map_thm3959.py
```

The frozen transcript is
`05-knowledge/results/jc2_centered_degree_five_root_map_thm3959.out`.
