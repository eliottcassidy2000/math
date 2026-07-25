---
id: THM-2262
title: "Degree-eighteen trigonal spectral discriminant reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In the polynomial exact-square-prefix, genuine
  nonsplit terminal quartic branch, every reduced degree-eighteen Keller
  mate lies on the explicit singular-discriminant locus of a trigonal
  spectral curve G(u,y)=0. Off y=0 the first two Faber fluxes eliminate
  Z=T^2 and give a cubic in u=dT whose branch discriminant has exact degree
  twelve. A squarefree discriminant gives a smooth irreducible genus-four
  cover with three distinct unramified points at infinity, so no
  nonconstant rational Keller trajectory exists. The exceptional wall y=0
  is empty: its Keller one-form forces a monomial cusp, while the exact
  degree-eighteen Faber sidecar has an uncancellable -21T^3/1024 pole.
  The generic third flux and Keller one-form are retained for the singular
  locus. This is a strict degree-eighteen reduction, not a closure of that
  locus or a proof of the planar Jacobian conjecture.
source: codex-2026-07-25-degree18-trigonal-spectral
depends_on:
  - THM-2071-quadratic-fiber-square-parity-gate
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2247-nonsplit-terminal-quartic-degree-fourteen-closure
related:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2217-square-prefix-pole-alternative-and-odd-leading-degree-terminal-wall
  - THM-2245-degree-fourteen-spectral-quartic-discriminant-reduction
script: 04-computation/jc2_degree18_trigonal_spectral_discriminant_thm2262.py
output: 05-knowledge/results/jc2_degree18_trigonal_spectral_discriminant_thm2262.out
script_sha256: 5c78356494dfd9973d7e33bc5df3d8cc4850ca746842e389b345b2acdae1ce22
output_sha256: 3aed261096fdc3ab811ee95fb2c1c61e3cd596ba1f1661054a503ea7bf665804
hash_basis: working-tree bytes (LF)
---

# THM-2262 -- degree eighteen lands on a singular trigonal spectrum

THM-2247 closes the reduced degree-fourteen branch left by THM-2245. The
next nonsplit Faber degree is eighteen. At this degree the first two
constant fluxes no longer complete to a hyperelliptic equation: after the
first elimination they give a cubic cover of a centered coefficient line.
The new point is that the generic cover has genus four, while the apparent
exceptional center is killed by information that the flux quotient forgets:
the polynomial Faber sidecar.

## 1. Inherited terminal coordinates

Work over

```text
R=C[x],                         K=C(x),
```

where `C` is algebraically closed of characteristic zero. Use the
polynomial exact-square prefix and genuine nonsplit terminal branch of
THM-2214 and THM-2247:

```text
P=H^2+L,

H=V z^2+Bz+C_0 in R[z],        V!=0,
L=A z+E in R[z],               A!=0,                 (1)
```

with `V` nonsquare in `K`. In the quadratic extension used to depress the
quartic, the intrinsic coordinates are

```text
P=w^4+2d w^2+q w+(d^2-s),      q^2=T,
T=A^2/V.                                               (2)
```

The deck sends `(w,q)` to `(-w,-q)` and fixes `d,s,T`. It is genuine, so

```text
q!=0,                          T!=0.                  (3)
```

After THM-2230's exact target-shear quotient and full Faber gauge, a
reduced degree-eighteen mate has the intrinsic form

```text
Q=E_18+alpha E_14+beta E_10+gamma E_6+delta E_2,     (4)
```

where all four displayed coefficients lie in `C`. THM-2129's Hamiltonian
identity and nonsplit parity give

```text
Phi_Q=0,             Psi_Q=Psi0 in C,
R_Q'=kappa/U,        kappa in C*,                    (5)
```

where `U^2=V`, and a prime holds `w` fixed. As before, if

```text
F=R_Q/q,
```

then `q=A/U` and `q^2=T` convert the last equation in (5) into the retained
Keller one-form

```text
A(2T F'+F T')=2kappa T.                              (6)
```

## 2. Exact degree-eighteen Laurent row

The new Laurent observables are

```text
Phi_18/q
 =63(2T^2d^2+T^2s-20Td s^2+10s^4)/128,

Psi_18
 =-63(-T^3d+20T^2d^2s+5T^2s^2
       -40Td s^3+4s^5)/256,

R_18/q
 =-3(-3T^3+84T^2d^3+210T^2ds
       -840Td^2s^2-280Ts^3+420ds^4)/512.            (7)
```

The exact companion obtains (7) in two independent ways: the
fibre-derivative Laurent recurrence and a finite multinomial expansion. It
rechecks the rows of degrees `2,6,10,14` at the same time.

Put

```text
u=dT,                  Z=T^2,                  y=9s-2alpha.       (8)
```

Before centering, the first flux divided by `q` is

```text
Phi_Q/q
 =[126u^2+560u alpha s-160u beta-1260u s^2
   -14Z alpha+63Zs-560alpha s^3+480beta s^2
   -384gamma s+256delta+630s^4]/128.                (9)
```

The coefficient of `Z` in (9) is `7y/128`. Consequently the analysis
splits honestly into the rational-function branches `y!=0` and `y=0`.

## 3. The trigonal spectral curve on `y!=0`

Define

```text
N_2
 =45927u^2+22680u alpha^2-58320u beta-5670u y^2
  -1680alpha^4-2240alpha^3y+8640alpha^2beta
  -840alpha^2y^2+8640alpha beta y-31104alpha gamma
  +2160beta y^2+93312delta-15552gamma y+35y^4.       (10)
```

The first equation in (5) is exactly

```text
Z=-2N_2/(5103y).                                     (11)
```

Substitute (11) in the second equation in (5). No squaring or root choice
is used. The result is the cubic

```text
G(u,y)=0,                                            (12)
```

where

```text
G
 =-5878656Psi0 y
  -26040609u^3
  -19289340u^2alpha^2+49601160u^2beta+1607445u^2y^2
  -2222640u alpha^4+11430720u alpha^2beta
  +1111320u alpha^2y^2+17635968u alpha gamma
  -20995200u beta^2-2857680u beta y^2
  -52907904u delta-138915u y^4
  +235200alpha^6+326144alpha^5y-1814400alpha^4beta
  +82320alpha^4y^2-2096640alpha^3beta y
  +4354560alpha^3gamma-62720alpha^3y^3
  +3110400alpha^2beta^2-423360alpha^2beta y^2
  -13063680alpha^2delta+2612736alpha^2gamma y
  -30380alpha^2y^4+3110400alpha beta^2y
  -11197440alpha beta gamma+241920alpha beta y^3
  -2612736alpha delta y-653184alpha gamma y^2
  +777600beta^2y^2+33592320beta delta
  -5598720beta gamma y+78120beta y^4
  +1959552delta y^2-435456gamma y^3+1127y^6.         (13)
```

Let

```text
Delta(y)=Disc_u G(u,y).                              (14)
```

The leading coefficient of `G` in `u` is the nonzero constant
`-26040609`, and the exact branch discriminant satisfies

```text
deg_y Delta=12,
lc_y Delta=-153384762202971019112448!=0.             (15)
```

The degree-twelve polynomial is frozen by its canonical symbolic digest

```text
05bf2aef4d6209c73e4d16ef40cd9a26cd9b3ee409a4af621ad0de9c5bd0ebee.
                                                               (16)
```

The condition below is not an identically vanishing artifact. At

```text
(alpha,beta,gamma,delta,Psi0)=(0,1,1,1,1),
```

the exact Euclidean algorithm gives

```text
gcd(Delta,Delta')=1.                                 (17)
```

### 3.1 The squarefree cover has genus four

Use the weighted coordinate

```text
v=u/y^2
```

at infinity. Dividing (13) by `y^6` gives the limiting cubic

```text
L_infinity(v)
 =1127-138915v+1607445v^2-26040609v^3.              (18)
```

Its discriminant is

```text
Disc_v L_infinity
 =-153384762202971019112448!=0.                      (19)
```

Thus the projective completion has three distinct smooth branches over
infinity. In the local coordinate `r=1/y`, the implicit-function theorem
applied to

```text
L_infinity(v)+O(r)=0
```

shows that all three are unramified for the degree-three map to the
`y`-line.

Suppose `Delta` is squarefree. A finite singular point of (12) would make
both `G_u` and `G_y` vanish and hence give a multiple zero of the branch
discriminant. Therefore the affine curve is smooth, and every one of the
twelve finite branch values in (15) is simple.

It is also irreducible. Indeed, because the leading `u` coefficient is a
unit, every nontrivial factor has positive `u`-degree. Distinct factors
split the three distinct weighted leading roots in (18). Their cross
resultant therefore has positive `y`-degree and nonzero leading
coefficient. The product formula

```text
Disc(fg)=Disc(f)Disc(g)Res(f,g)^2                   (20)
```

would put a nonconstant square factor in `Delta`, contrary to
squarefreeness.

The smooth irreducible degree-three cover therefore has twelve simple
finite branch points and no branch at infinity. Riemann--Hurwitz gives

```text
2g-2=3(-2)+12=6,                    g=4.             (21)
```

### 3.2 The Keller trajectory forces singular parameters

The pair `(u,y) in C(x)^2` defines a rational map from `P^1` to the smooth
projective curve in Section 3.1. Properness extends it across every pole.
A nonconstant map from `P^1` to a genus-four curve is impossible by
Riemann--Hurwitz, so `u,y` are constant. Equation (11) makes `Z=T^2`
constant. Since the constant field is `C`, first `T` and then `q` are
constant. But the genuine deck fixes `C` and sends the nonzero `q` to
`-q`, a contradiction.

Hence every degree-eighteen survivor with `y!=0` satisfies

```text
gcd(Delta,Delta')!=1,

equivalently                    Disc_y(Delta)=0.     (22)
```

This is an explicit algebraic hypersurface in the five constant
parameters `(alpha,beta,gamma,delta,Psi0)`.

## 4. The exceptional center `y=0` is empty

The division by `y` in (11) loses no survivor. On `y=0`, the first flux is

```text
[15309u^2+7560u alpha^2-19440u beta
 -560alpha^4+2880alpha^2beta
 -10368alpha gamma+31104delta]/15552=0.              (23)
```

Thus `u=dT` is algebraic over `C`, hence constant. The coefficient of `Z`
in the second flux is

```text
(567u+140alpha^2-360beta)/2304.                      (24)
```

If it were nonzero, the second flux would make `Z` constant and the deck
contradiction of Section 3.2 would apply. Every exceptional survivor must
therefore satisfy

```text
u=u_0=-20(7alpha^2-18beta)/567,                      (25)

delta
 =(490alpha^4-2520alpha^2beta+3402alpha gamma
    +2025beta^2)/10206,                              (26)

Psi0
 =2(14alpha^2-45beta)
    (56alpha^3-225alpha beta+486gamma)/45927.        (27)
```

Under (25)--(27), the third flux collapses to

```text
F=R_Q/q=T(2187T^2-K)/124416,                         (28)

K=4480alpha^3-17280alpha beta+31104gamma.            (29)
```

Substitution into the Keller one-form (6), followed by cancellation of
the nonzero `T`, gives

```text
A(1701T^3-KT)'=82944kappa.                           (30)
```

### 4.1 The one-form forces a monomial cusp

Use the rational-primitive lemma of THM-2071/2214. If

```text
A in C[x]\{0},       S in C(x),       AS' in C*,
```

then either `A` is constant and `S` is affine-linear, or, after translating
`X=x-xi`,

```text
A=a_0X^m,        m>=2,
S=c_0+c_1X^(1-m),        c_1!=0.                    (31)
```

The constant case is impossible for `S=1701T^3-KT`: a finite pole of `T`
would give a pole of `S`, so `T` would be a polynomial, but a nonconstant
cubic in a polynomial cannot have degree one.

In the second case put `t=X^-1`. The same pole argument makes `T=R(t)` a
nonconstant polynomial, and differentiation gives

```text
(5103R(t)^2-K)R'(t)=c t^(m-2),       c!=0.           (32)
```

If `K!=0`, the quadratic factor in (32) has two distinct roots. The
nonconstant polynomial `R` assumes both values at finite points, whereas
the right side can vanish only at `t=0` (and never vanishes when `m=2`).
Both values cannot equal `R(0)`. Hence

```text
K=0.                                                 (33)
```

Equation (31) now says that a cube is a binomial:

```text
1701R(t)^3=c_0+c_1t^(m-1).                          (34)
```

If both terms on the right were nonzero, all of its nonzero roots would
be simple, impossible for a cube. Therefore `c_0=0`; for some `k>=1`,

```text
m=3k+1,                 T=t_0X^-k,       t_0!=0.     (35)
```

For reference, (2) then gives

```text
V=A^2/T=v_0X^(7k+2).                                 (36)
```

The genuine nonsplit condition would force `k` odd, since a monomial over
the algebraically closed field is a rational square exactly when its
exponent is even. The stronger pole contradiction below does not need
this parity.

### 4.2 The polynomial Faber sidecar forbids the cusp

The flux quotient alone cannot see polynomiality of the reduced mate. To
restore it, for `j>=1` define the exact binomial truncation

```text
R_j(P,H)
 =sum_(i=0)^(j-1) binom(j-1/2,i)
    H^(2j-1-2i)L^i.                                  (37)
```

Both `P,H,L` are polynomial, so every `R_j(P,H)` belongs to `R[z]`.
Independent whole-polynomial expansion of the Faber seeds gives

```text
E_6-R_2
 =3T/8,

E_10-R_3
 =-5T(-L+2s)/16,

E_14-R_4
 =35T(L^2-2Ls-Td+3s^2)/128,

E_18-R_5
 =-21T[-12L^3+24L^2s+L(12Td-36s^2)
        +T^2-48Tds+48s^3]/1024.                     (38)
```

On the exceptional branch, `s=2alpha/9` and `Td=u_0` are constants.
Therefore the sidecar of the whole combination (4) has the form

```text
Q-[R_5+alpha R_4+beta R_3+gamma R_2+delta R_1]

 =-21T^3/1024+T M(x,z),              M in C[x,z].    (39)
```

At `X=0` in (35), the first term of (39) has valuation `-3k`, while the
second has valuation at least `-k`. It cannot cancel. Thus the left side
has a pole. But `Q` and every truncation in brackets are polynomials, a
contradiction. This eliminates the entire `y=0` wall.

The sidecar is the decisive missing coordinate: the spectral fluxes alone
leave the monomial cusp (35), while polynomial descent kills it
immediately.

## 5. Third flux retained on the singular locus

For the remaining branch (22), the third observable should not be
recomputed from the discriminant. Before substituting (11), it is

```text
F=R_Q/q=-N_3/(373248T),                              (40)
```

where

```text
N_3
 =183708u^3+90720u^2alpha^2-233280u^2beta-22680u^2y^2
  +51030uZy-6720u alpha^4-8960u alpha^3y
  +34560u alpha^2beta-3360u alpha^2y^2
  +34560u alpha beta y-124416u alpha gamma
  +8640u beta y^2+373248u delta-62208u gamma y+140u y^4
  -6561Z^2+13440Z alpha^3+10080Z alpha^2y
  -51840Z alpha beta-25920Z beta y+93312Z gamma-840Zy^3.
                                                               (41)
```

Equations (6), (11), (12), (14), and (40)--(41) are the complete retained
continuation state. The next task is not another generic genus argument:
it is to normalize the repeated-branch locus (22), pull back the one-form
(6), and simultaneously retain the polynomial sidecar (38).

## 6. Scope and reproduction

This theorem concerns only the polynomial exact-square-prefix, genuine
nonsplit terminal quartic branch in reduced degree eighteen. It does not
close the singular locus (22), the split deck, even-leading descent, other
short Newton edges, `JC(2)`, or `DC(2)`.

Run

```bash
python3 04-computation/jc2_degree18_trigonal_spectral_discriminant_thm2262.py
python3 -O 04-computation/jc2_degree18_trigonal_spectral_discriminant_thm2262.py
```

Both runs are byte-identical to the stored output, including an independent
replay and hostile proof audit. The companion checks
the five Faber rows by two coefficient engines; the `Z` elimination and
the full cubic (13); the degree and leading coefficient of (14); the
squarefree specialization (17); the infinity cubic and its discriminant;
the retained third flux; all exceptional relations; and the
whole-polynomial sidecar (38). The genus, irreducibility,
rational-primitive, and valuation arguments are the mathematical proof
above rather than delegated computer conclusions.
