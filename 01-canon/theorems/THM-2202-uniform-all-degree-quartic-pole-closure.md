---
id: THM-2202
title: "Uniform all-degree quartic pole closure"
status: >
  PROVED. Let a planar Keller component be
  P=V^2z^4+Vbz^3+gamma z^2+delta z+epsilon after THM-2180, and suppose
  the reduced mate has any twice-odd fibre degree n=4r-2. Then
  V divides 4gamma-b^2, without any square or nonsquare assumption on V.
  At a hypothetical bad place, the boundary/flux triple first forces the
  exact square cage of THM-2189. For every even Faber seed m_s=4s-2,
  the expansion (K^2+B)^(s-1/2) has a universal two-chamber filtration.
  Below the cusp r nu(Lambda)=(2r-1)nu(c), the top first and second
  fluxes both have order r nu(Lambda); an odd seed which cancels the first
  flux has its second flux exactly nu(Lambda) orders shallower. Above the
  cusp, the top boundary has the unique c^(2r-1) tooth; an odd seed which
  cancels it has a uniquely deepest first-flux pole. Lower even seeds are
  strictly shallower in both chambers. This proves polynomiality of the
  canonical quadratic approximate root in every remaining reduced degree
  and removes the quartic finite-pole obstruction completely. The explicit
  degree-ten bank is retained as an exact referee. The terminal nonmonic
  square-prefix/quadratic-member descent remains open, so this is not a
  proof of JC(2) or DC(2).
source: codex-2026-07-24-split-all-degree-quartic-pole
depends_on:
  - THM-2129
  - THM-2158
  - THM-2180
  - THM-2189
  - THM-2194
related:
  - THM-2181
script: 04-computation/jc2_quartic_degree10_phase_bank_codex_20260724.py
output: 05-knowledge/results/jc2_quartic_degree10_phase_bank_codex_20260724.out
script_sha256: ce91b7916fd9262bb23d7875ca6b8b577ca1b9580fb762a66809420990352863
output_sha256: 02760315d6d552c1e218f46609d0134b7a166a093342a0c3f8dc1d4079454a0d
hash_basis: working-tree bytes (LF)
---

# THM-2202 -- every twice-odd quartic pole closes

Let

```text
P=V^2z^4+beta z^3+gamma z^2+delta z+epsilon
```

belong to a planar polynomial Keller pair over `C`. Reduce the mate by
polynomial target shears and suppose its remaining fibre degree is

```text
n=4r-2,                         r>=1.                 (1)
```

THM-2180 gives `beta=Vb`, with `b in C[x]`. Put

```text
D=4gamma-b^2.                                         (2)
```

The case `r=1`, namely reduced degree two, is the uniform degree-two
addendum of THM-2189. We may therefore assume `r>=2` below; lower even seeds
still include `m_1=2`.

Then, with no assumption on whether `V` is a square in `C(x)`,

```text
V divides D.                                          (3)
```

By THM-2158, this is exactly the polynomiality of the canonical quadratic
approximate root

```text
H_0=Vz^2+(b/2)z+D/(8V).                               (4)
```

## 1. The inherited square cage

Fix an irreducible `pi|V` and suppose for contradiction that, for an
extended `pi`-valuation `nu`,

```text
e=nu(V)>0,                 d=nu(D)<e.                 (5)
```

The initial negative-scale argument in THM-2189, Section 2, does not use
deck parity. The top degree `n` is uniquely deeper than every lower Faber
degree in the polynomial boundary and both fluxes. Regularity and constancy
therefore force the same three consecutive boundary coefficients to vanish.
THM-2129 gives the exact square face.

Writing `f=nu(b)` and `ell_0=nu(delta)`, one obtains

```text
2f<d<e,             nu(gamma)=2f,
ell_0+e>3f.                                          (6)
```

Choose `U^2=V`, on either the field or split deck, and put

```text
a=b/(2U),                 A=-nu(a)=e/2-f>0,
C=d-2f>0,                 H=e-d>0.                   (7)
```

Then

```text
A=(C+H)/2,                         0<C<2A.            (8)
```

In the canonically scaled reverse variable `X=au`, THM-2189 gives

```text
T
 =1+2X+(1+c)X^2+lX^3+mX^4
 =K^2+B,                                                (9)

K=1+X+(c/2)X^2,
B=Lambda X^3+Omega X^4,
c=D/b^2,
Lambda=l-c,
Omega=m-c^2/4.                                      (10)
```

The exact valuation data are

```text
nu(c)=C,                 nu(Lambda)>0 or Lambda=0,
mu:=Omega+c^2/4=m,       nu(mu)>=4A,
nu(Omega)=2C.                                           (11)
```

The last equality follows from `4A>2C`. Equivalently,

```text
T=(1+X)^2+cX^2(1+X)+Lambda X^3+mu X^4.               (12)
```

The identities (10)--(12), especially `nu(mu)>2C`, retain the coordinate
which is lost if `Omega` is treated as an arbitrary quartic coefficient.

## 2. The all-degree Faber observables

For every reduced degree `j`, define

```text
B_j=[X^j]T^(j/4),
F_j=[X^(j+1)]T^(j/4),
G_j=[X^(j+2)]T^(j/4)
       +(1/2)[X^(j+1)]T^(j/4).                       (13)
```

The exact contribution of `E_j` to the polynomial boundary, first flux
divided by four, and second flux divided by four is

```text
a^j B_j,                  a^(j+1)F_j,
a^(j+2)G_j.                                           (14)
```

For `s>=1`, put

```text
m_s=4s-2.                                             (15)
```

The even Faber series has the exact expansion

```text
T^(s-1/2)
 =sum_(q>=0) binom(s-1/2,q)
      K^(2s-1-2q) B^q.                               (16)
```

For `q<s`, the summand is a polynomial of degree at most

```text
2(2s-1-2q)+4q=4s-2=m_s,                              (17)
```

so it contributes neither `F_(m_s)` nor `G_(m_s)`.

### Lemma 1: the Lambda-dominated flux pair

Let `L=nu(Lambda)<2C`. Then

```text
nu(F_(m_s))=nu(G_(m_s))=sL.                          (18)
```

Indeed, the unique least-order part comes from `q=s` in (16), taking
`Lambda X^3` from every factor of `B`. Modulo the positive-valuation ideal,
`K=1+X`, so the relevant coefficients of `1/K` are

```text
[X^(s-1)](1+X)^-1=(-1)^(s-1),

[X^s](1+X)^-1
 +(1/2)[X^(s-1)](1+X)^-1
 =(1/2)(-1)^s.                                      (19)
```

Both are nonzero. Every use of `Omega`, every positive-order correction
from `K`, and every `q>s` term has strictly larger valuation.

### Lemma 2: the c-dominated first flux

If `L>2C`, then

```text
nu(F_(m_s))=L+2(s-1)C.                               (20)
```

The leading term again has `q=s`; it uses one `Lambda X^3` and
`s-1` copies of `Omega X^4`, so its `X`-degree is exactly `4s-1`.
Its coefficient is a nonzero multiple of

```text
s Lambda Omega^(s-1).                                (21)
```

For a term with `q=s+h`, `h>=1`, to reach degree at most `4s-1` it must
use at least `4h+1` copies of `Lambda X^3`. If it uses `k` such copies,
its valuation exceeds (20) by

```text
(k-1)(L-2C)+2hC>0.                                   (22)
```

At the boundary `L=2C`, the same argument gives the sufficient lower bound

```text
nu(F_(m_s))>=2sC.                                    (23)
```

If `Lambda=0`, every even first flux vanishes. This is also immediate in
the centered coordinate: an even quartic raised to an even Faber degree
has no `Z^-1` Laurent coefficient.

### Lemma 3: the boundary cusp

If

```text
sL>(2s-1)C
```

or `Lambda=0`, then

```text
nu(B_(m_s))=(2s-1)C.                                 (24)
```

To prove this without a genericity assumption, expand (12). A monomial
obtained by choosing `p` copies of `cX^2(1+X)`, `q` copies of
`Lambda X^3`, and `h` copies of `mu X^4` has the form

```text
c^p Lambda^q mu^h X^(2p+3q+4h)
 (1+X)^(2s-1-p-2q-2h).                              (25)
```

Give `c,Lambda,mu` the cusp weights

```text
1,              (2s-1)/s,              2.            (26)
```

If the final exponent of `(1+X)` in (25) is nonnegative, contribution to
`[X^(4s-2)]` forces

```text
p+q+2h>=2s-1.                                        (27)
```

Its cusp weight is therefore at least `2s-1`. If that exponent is
negative, then `p+2q+2h>=2s`; for `q<=s` this again gives weight at least
`2s-1`, while for `q>s` the `Lambda` weight alone is already larger.

Equality in this weight bound is possible only for

```text
Lambda^s
```

or for monomials `c^p mu^h` with

```text
p+2h=2s-1.                                           (28)
```

The actual inequality `nu(mu)>=4A>2C` makes every equality candidate with
`h>0` strictly higher than `(2s-1)C`. The strict hypothesis in Lemma 3
makes `Lambda^s` higher as well. The surviving coefficient is

```text
[c^(2s-1)]B_(m_s)=binom(s-1/2,2s-1)!=0.              (29)
```

This proves (24), including `s=1`, where it reads `B_2=c/2`.

### Lemma 4: the odd bank

For every odd `j`,

```text
nu(B_j)=nu(F_j)=0,                                   (30)
```

because at `c=Lambda=Omega=0`,

```text
(B_j,F_j)
 =(binom(j/2,j), binom(j/2,j+1)),                    (31)
```

and both generalized binomial coefficients are nonzero.

Moreover,

```text
G_j=Lambda U_j,                       nu(U_j)=0.      (32)
```

When `Lambda=0`, the quartic is even in the translated Laurent coordinate
whose `Z^-2` observable is encoded by the combination defining `G_j` in
(13). Its odd Faber Laurent series contains only odd powers, so that
`Z^-2` coefficient, equivalently `G_j`, vanishes. Thus `Lambda` divides
`G_j`. For `j>=3`, differentiating at the exact square gives

```text
U_j mod m
 =-j/(8(j-1)) binom(j/2-2,j-2)!=0,                  (33)
```

while `U_1 mod m=1/4`. Hence `U_j` is a unit.

These four lemmas are the complete bank needed below. They are identities
and valuation consequences, not finite evidence.

## 3. Below and on the cusp

Normalize the nonzero coefficient of `E_n=E_(m_r)` to one. The reduced
Faber normal form is

```text
Q=J(P)+E_(m_r)
   +sum_(1<=s<r) u_s E_(m_s)
   +sum_(odd j<n) v_j E_j,                           (34)
```

with constant coefficients. Multiples of four have been absorbed into the
target shear `J(P)`.

Suppose first that `Lambda!=0` and

```text
rL<=(2r-1)C.                                         (35)
```

Then `L<2C<4A`. Lemma 1 gives top first- and second-flux coefficient
orders `rL`, and order `sL` for every lower even degree `m_s`.
The top first-flux valuation is strictly below every lower even one:

```text
[rL-(4r-1)A]-[sL-(4s-1)A]
 =(r-s)(L-4A)<0.                                     (36)
```

Odd first-flux terms have the pairwise distinct valuations

```text
-(j+1)A.                                             (37)
```

Any active odd term strictly deeper than the top is therefore unique and
impossible. After eliminating all such terms, first-flux constancy can hold
only if one odd degree `j` ties the top:

```text
rL-(4r-1)A=-(j+1)A,
rL=(n-j)A.                                           (38)
```

If no such `j` exists, the top first-flux pole is unique.

Under (38), Lemma 4 gives the matching odd second-flux valuation

```text
L-(j+2)A=(r+1)L-4rA,                                 (39)
```

whereas the top second-flux valuation is

```text
rL-4rA.                                              (40)
```

Thus the top is strictly deeper by `L>0`. It is also strictly below every
lower even second flux by the same difference (36). Lower odd degrees are
at least two powers of `a` shallower, and higher odd coefficients already
vanished. Therefore (40) is a unique negative pole, contradicting constancy
of the second flux.

This includes the cusp equality in (35). No leading-residue genericity or
next-tooth analysis is needed: Lemma 1 makes the first and second fluxes
simultaneously exact there.

## 4. Above the cusp

Suppose now that

```text
Lambda=0
  or
rL>(2r-1)C.                                          (41)
```

For every `s<r`,

```text
sL>(2s-1)C,                                          (42)
```

because `(2r-1)/r>(2s-1)/s`. Lemma 3 gives boundary
coefficient orders `(2r-1)C` for the top and `(2s-1)C` below it.
The top boundary valuation is strictly below every lower even one:

```text
[(2r-1)C-(4r-2)A]
 -[(2s-1)C-(4s-2)A]
 =2(r-s)(C-2A)<0.                                    (43)
```

The top boundary is polar because `(2r-1)C<(4r-2)A`. As before, every
active odd boundary term strictly deeper than the top is unique and must
vanish. A cancellation of the top boundary can therefore occur only with
one odd degree `j`, necessarily satisfying

```text
(2r-1)C-(4r-2)A=-jA,
(2r-1)C=(n-j)A.                                      (44)
```

The corresponding odd first-flux term has valuation

```text
-(j+1)A=(2r-1)C-(4r-1)A<0.                           (45)
```

We show that it is unique.

The top first-flux coefficient order is strictly greater than
`(2r-1)C`:

```text
rL>(2r-1)C                         if L<2C,

nu(F_(m_r))>=2rC>(2r-1)C           if L=2C,

L+2(r-1)C>(2r-1)C                  if L>2C,

nu(F_(m_r))=infinity               if Lambda=0.       (46)
```

For a lower even degree `m_s`, subtracting (45) from its first-flux
valuation gives, respectively,

```text
sL-(2r-1)C+4(r-s)A
 >(r-s)(4A-(2r-1)C/r)>0            if L<2C,           (47)

2sC-(2r-1)C+4(r-s)A
 =C+2(r-s)(2A-C)>0                 if L=2C,           (48)

L+2(s-1)C-(2r-1)C+4(r-s)A
 =(L-C)+2(r-s)(2A-C)>0             if L>2C.           (49)
```

If `Lambda=0`, every even first flux vanishes. All lower odd degrees are
shallower than (45), while all higher odd coefficients were already forced
to vanish by the boundary. Thus (45) is a unique negative first-flux pole,
a contradiction.

Together with the cited degree-two case `r=1`, Sections 3 and 4 exhaust
`Lambda`, both cusp boundaries, every `r>=1`, and every support in (34).
The bad-place assumption (5) is impossible, proving (3).

## 5. The explicit degree-ten bank

For `r=3`, the top degree is ten. The exact recurrence

```text
T(T^(j/4))'=(j/4)T'T^(j/4)                           (50)
```

gives

```text
B_10=(-10Lambda^3+30Lambda^2Omega+30Omega^2c
        +10Omega c^3+c^5)/32,                        (51)

F_10=-(5Lambda/32)
       (Lambda^2c-2Lambda^2+6Lambda Omega-6Omega^2),
                                                               (52)

G_10=(5/128)(
       -Lambda^4+6Lambda^3c-4Lambda^3
       -12Lambda^2Omega c+12Lambda^2Omega
       -12Lambda Omega^2+8Omega^3).                  (53)
```

The lower even rows are

```text
B_2=c/2,
F_2=Lambda/2,
G_2=-(Lambda-2Omega)/4;                              (54)

B_6=(3Lambda^2+6Omega c+c^3)/8,
F_6=-3Lambda(Lambda-2Omega)/8,

G_6=3((Lambda-2Omega)^2
       +Lambda^2(1-2c))/32.                          (55)
```

For the odd degrees `1,3,5,7,9`, the exact square residues are

```text
 j       B_j mod m       F_j mod m       (G_j/Lambda) mod m
 1          1/2            -1/8                    1/4
 3         -1/16            3/128                  3/32
 5          3/256          -5/1024                -5/512
 7         -5/2048         35/32768                7/4096
 9         35/65536       -63/262144             -45/131072.
                                                               (56)
```

This is the `r=3` instance of Lemmas 1--4, not an extra hypothesis.

## 6. Exact consequence and remaining debt

The argument proves (3) at every divisor of `V` on either deck. Hence the
canonical approximate root (4) is polynomial throughout the remaining
twice-odd quartic branch. Combined with THM-2129's odd-degree and square-face
classification, THM-2181's monic depressed closure, and the preceding
reductions, there is no remaining quartic finite-pole survivor.

This does not yet close the quartic source-fibre stratum. The polynomial
approximate root (4) has leading coefficient `V`, not one. The exact square
prefix

```text
P=H_0^2+(linear remainder)
```

still needs a terminal nonmonic monicization/quadratic-member theorem.
Neither the present pole closure nor THM-2181 supplies that descent. Thus
`JC(2)` and `DC(2)` remain open.

The all-degree mechanism can be summarized as

```text
first rational tooth
  -> cusp rL=(2r-1)C
  -> Lambda chamber: equal top F/G orders
  -> c chamber: unique c^(2r-1) boundary tooth
  -> odd cancellation loses one flux order
  -> lower even seeds remain strictly shallower.                 (57)
```

## 7. Exact referee

Run

```bash
python3 04-computation/jc2_quartic_degree10_phase_bank_codex_20260724.py
```

The first path independently derives every degree
`1,2,3,5,6,7,9,10` boundary/two-flux row from (50), verifies
(51)--(56), and checks the exact `Omega=mu-c^2/4` identities. The second
path uses an exact `Fraction` Laurent solver on

```text
1<=A<=6,              1<=C<2A,             1<=E<=4A,

T=1+2X+(1+t^C)X^2+lambda t^E X^3,

lambda in {-2,-1,-1/2,0,1/2,1,2}.                    (58)
```

All `4508` degree-ten hostile rows reject simultaneous boundary and
two-flux regularity. The `Lambda=0` row
`A=2,C=E=3,lambda=1` is a named negative control. The exact square
`T=(1+X)^2` is a positive algebraic control. Normal and optimized Python
execute the same explicit exceptions and have byte-identical output.

The output labels this sweep `FINITE_REFEREE_ONLY`. It checks the first
nontrivial all-degree instance beyond THM-2194, but no finite bound supplies
a quantifier in Sections 2--4. QED.
