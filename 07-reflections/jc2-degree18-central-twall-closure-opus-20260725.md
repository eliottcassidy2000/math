---
source: codex-2026-07-25-jc2-central-twall-closure
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED; SUPERSEDED AS
  A CLOSURE STATEMENT BY THM-2345. In the
  genuine nonsplit polynomial
  exact-square-prefix degree-eighteen branch of THM-2262/2297, the whole
  remaining central sublocus D=25B^2/126, B*R*S!=0, T=0 is empty. On
  B=1, Disc_W(T)=K*P8(C^2)^2*P3(C^2)^3. The 16 P8 points are ordinary
  nodes and the six P3 points are ordinary cusps. At smooth, nodal, and
  cuspidal points the irreducible trigonal spectral normalization has
  positive genus, so a rational Keller trajectory is constant and the
  retained square/deck or y=0 polynomial sidecar gives a contradiction.
  The complementary R*S*T!=0 locus has ten simple residual branch values
  and genus four. Together with the separately audited R=0 and S=0
  closures, this closes the entire central ratio inside the THM-2262/2297
  branch. It does not close the other degree-eighteen discriminant
  components or prove JC(2) or DC(2).
depends_on:
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - jc2-degree18-central-weight30-factor-opus-20260725
  - jc2-degree18-central-wall-r0-closure-opus-20260725
  - jc2-degree18-central-wall-s0-closure-opus-20260725
related:
  - THM-2311-degree-eighteen-two-sparse-weighted-ratio-bank
  - THM-2313-degree-eighteen-bd-linear-ratio-closure
  - THM-2345-degree-eighteen-common-root-wall-saturation
---

# The complete central degree-eighteen wall has positive-genus normalization

> **CANONICAL SUPERSESSION.** THM-2345 now closes the entire common-root
> wall `126D=25B^2`, not merely this `T=0` central chart.  This file is
> retained for its independent smooth/node/cusp spectral atlas and
> positive-genus mechanism; it is no longer the closure authority.

## 1. Target and inheritance

After THM-2297's legal target translation, use the spectral cubic

```text
G(u,y)
 =-5878656Wy-26040609u^3+49601160Bu^2+1607445u^2y^2
  -20995200B^2u-2857680Buy^2-52907904Du-138915uy^4
  +777600B^2y^2+33592320BD-5598720BCy+78120By^4
  +1959552Dy^2-435456Cy^3+1127y^6.                 (1)
```

On the central ratio

```text
D=25B^2/126,                                       (2)
```

the exact weight-thirty-factor packet proves

```text
Disc_y(Disc_u G/y^2)=K R^6 S^3 T,                 K in Q*,          (3)

R=20BC+21W,

S=2888B^5+108864B^2C^2+571536BCW+750141W^2,       (4)
```

where `T` is the explicit primitive 21-term polynomial recorded in
`jc2-degree18-central-weight30-factor-opus-20260725.md`. The separately
audited central packets close `R=0` and `S=0`. It remains to close

```text
B!=0,                   R!=0,                   S!=0,                   T=0.
                                                                    (5)
```

The operation used here is not another coefficient domination. We normalize
the discriminant curve `T=0`, classify every place where that parameter
curve itself is singular, and count ramification on the *spectral*
normalization over each parameter point.

Weighted covariance permits `B=1`. All statements below then lift back to
every `B!=0` orbit.

## 2. The complement of the three factors has genus four

Before entering `T=0`, suppose

```text
R*S*T!=0.                                                           (5a)
```

Equation (3) says that the degree-ten polynomial `delta_bar` is
squarefree, so it has ten distinct finite nonzero roots. The central fibre
`y=0` is smooth and totally ramified by (27)--(28) below, contributing
two. The ten simple residual roots contribute ten, and infinity is
unramified by (29). The polynomial-root argument of Section 6 also proves
irreducibility here: a rational root would force `C=W=0`, hence `R=0`,
contrary to (5a). Riemann--Hurwitz therefore gives

```text
Ram=2+10=12,                 2g-2=-6+12,             g=4.           (5b)
```

The rational Keller trajectory is constant; the first-flux square and
genuine deck give the same contradiction as in Section 8. Thus the
factorization complement (5a) is empty as well. The rest of this packet
closes the only still-unhandled factor `T=0` with `R*S!=0`.

## 3. The singular atlas of the parameter curve

Exact reconstruction of the discriminant of `T(1,C,W)` in `W` gives

```text
Disc_W T(1,C,W)=K_0 P8(C^2)^2 P3(C^2)^3,          K_0 in Q*,          (6)
```

Here `T` has weight 30, degree six in the weight-five variable `W`, so
its `W`-discriminant has weight

```text
(6-1)(2*30-6*5)=150.                                  (6a)
```

On `B=1`, its `C`-degree is consequently at most `150/3=50`. Thus 51
exact consecutive `C`-evaluations determine the whole discriminant in
(6); this is the load-bearing interpolation bound used by the companion.

where the coefficients are written from constant term upward:

```text
P3(X)=
  72250000
 +26254935000 X
 +1530102783825 X^2
 +24450881901768 X^3,                              (7)
```

and

```text
P8(X)=
 -1087948019845507812500000000
 +92930586076313704687500000000 X
 +8954570361875562292640625000000 X^2
 -313103637508400358291807187500000 X^3
 +13873912636573086466704217312500000 X^4
 -588322473659000442440395871400000000 X^5
 +8764467368105916596083808036836642500 X^6
 -68598223636291030733181729239483693550 X^7
 +57091662724422970717844096005443494577 X^8.      (8)
```

The exact finite-field certificates are

```text
P3(C^2) irreducible modulo 11,                 degree 6;

P8(C^2) irreducible modulo 173,                degree 16.             (9)
```

In particular both characteristic-zero polynomials are irreducible,
squarefree, and coprime. Put

```text
K_i=Q[C]/(P_i(C^2)),                  i in {3,8}.                     (10)
```

Euclidean arithmetic in `K_i[W]` gives

```text
gcd(T,T_W)=W+g_i(C),                                                (11)

T_C(C,-g_i(C))=0.                                                   (12)
```

For the cubic factor one may take

```text
g_3(C)
 =965C/714-104733C^3/3400-85562001C^5/53125.                        (13)
```

The leading `W^6` coefficient of `T(1,C,W)` is a nonzero constant.
Therefore (6), (9), and (11)--(12) find every affine solution of

```text
T=T_C=T_W=0                                                        (14)
```

without a lost leading-coefficient point. There are exactly

```text
16+6=22                                                           (15)
```

singular parameter points. Direct reduction in both quotient fields also
shows that `R` and `S` are nonzero there. Thus none belongs to one of the
already closed central faces.

## 4. Sixteen nodes and six cusps

At a point of (14), exact quotient arithmetic gives `T_WW!=0`. The
implicit equation `T_W(C,w(C))=0` therefore defines the critical section
`W=w(C)`. Put

```text
w'=-T_CW/T_WW,

Phi2=T_CC-T_CW^2/T_WW,                                             (16)

Phi3=T_CCC+3T_CCW w'+3T_CWW(w')^2+T_WWW(w')^3.                     (17)
```

The two irreducible strata have the exact profiles

```text
P8(C^2)=0:       T_WW!=0,       Phi2!=0;

P3(C^2)=0:       T_WW!=0,       Phi2=0,       Phi3!=0.              (18)
```

For example, nonzero reduced coordinates on the `P3` stratum include

```text
[T_WW]_(C^0)=611393420358000000000000/23,

[Phi3]_(C^1)=5660096667257733120000000000/8993.                    (19)
```

The first line of (18) is the nondegenerate Hessian test, so all 16
`P8` points are ordinary nodes. In the second line the Hessian has rank
one and the first nonzero derivative in its null direction has order
three. Hence all six `P3` points are ordinary cusps. This is a complete
analytic singularity atlas of the affine `T`-curve, not merely a count of
zeros of its discriminant.

## 5. The branch profiles on the two singular strata

Let

```text
Delta(y)=Disc_u G(u,y),                  delta_bar=Delta(y)/y^2.     (20)
```

On (2), `delta_bar` has degree ten. Reducing its coefficients at the
critical section (11), exact gcds in `K_i[y]` give

```text
stratum     deg gcd(delta_bar,delta_bar')   repeated-gcd defect
P8          2                              squarefree
P3          2                              one repeated root.        (21)
```

Equivalently, the finite nonzero branch-value profiles are

```text
P8:       6 simple roots + 2 double roots;

P3:       7 simple roots + 1 triple root.                            (22)
```

At a smooth point of `T=0` with `R*S!=0`, the universal discriminant
hypersurface is smooth. Indeed,

```text
gradient_(C,W)(K R^6 S^3 T)!=0.
```

By the chain rule, the universal degree-ten discriminant gradient is
nonzero. Its singular locus consists exactly of polynomials with two
double roots or one root of multiplicity at least three. A degree-ten
polynomial at our smooth point therefore has exactly one double root and
all remaining roots simple. Thus the smooth profile is

```text
8 simple roots + 1 double root.                                    (23)
```

Equations (21)--(23) are deliberately stated as root-multiplicity
profiles. A double zero of the discriminant need not itself be counted as
one simple ramification point; the genus lower bound below uses only the
simple zeros and the independently controlled central fibre.

## 6. Irreducibility of the trigonal spectral curve

For every point of `T=0`, the cubic (1) is irreducible over `C(y)`.
Because its leading coefficient in `u` is a nonzero constant, a rational
root would be integral over `C[y]`, hence polynomial. Weighted degree at
infinity restricts it to

```text
u=a y^2+b y+c.                                                       (24)
```

The squarefree infinity cubic

```text
1127-138915a+1607445a^2-26040609a^3=0                              (25)
```

and the next coefficient of (1) force `b=0`. The odd `y^3` coefficient
then forces `C=0`, and the odd `y` coefficient forces `W=0`. But

```text
T(1,0,0)=256000000000000000!=0.                                    (26)
```

Thus (24) cannot occur on the wall. Since a cubic is reducible exactly
when it has a rational root, the spectral normalization is connected of
degree three over the `y`-line.

## 7. Uniform ramification at zero and infinity

At `y=0`, the cubic specializes exactly to

```text
G(u,0)=-26040609(u-40/63)^3.                                       (27)
```

Moreover

```text
G_y(40/63,0)=-279936 R!=0.                                         (28)
```

The local spectral curve is smooth there and the map to `y` is totally
ramified. This fibre contributes two to the ramification divisor. It is
disjoint from the roots of `delta_bar`, because with
`a=-26040609`,

```text
delta_bar(0)=-27 a^2(279936R)^2!=0.                   (28a)
```

At infinity, the cubic (25) has nonzero discriminant

```text
-153384762202971019112448,                                         (29)
```

so the three branches are distinct, smooth, and unramified.

Every simple root of `delta_bar` supplies one further simple ramification
point on the normalization. We may ignore the multiple roots and still
obtain the lower bounds

```text
smooth T point:       Ram>=2+8=10,        so g>=3;

P8 node:              Ram>=2+6=8,         so g>=2;

P3 cusp:              Ram>=2+7=9.                              (30)
```

For a connected degree-three cover, Riemann--Hurwitz reads

```text
2g-2=-6+Ram.                                                        (31)
```

The left side is even, so the last lower bound improves automatically to
`Ram>=10`, hence `g>=3`. In all three cases the spectral normalization has
positive genus.

## 8. Keller contradiction and full central closure

A rational Keller trajectory gives a rational map from `P^1` to the
projective spectral normalization. Positive genus forces this map to be
constant.

If its constant `y` is nonzero, THM-2262's first-flux equation makes
`Z=T_flux^2` constant, and then makes `T_flux` and `q` constant. The
genuine deck fixes the base invariants and changes the sign of the nonzero
square root, which is impossible over the algebraically closed constant
field. If `y=0`, THM-2262's retained third-flux, Keller one-form,
rational-primitive lemma, and whole-polynomial Faber sidecar exclude the
monomial cusp. Therefore (5) is empty.

Combining this result with the independently audited `R=0` and `S=0`
packets and the exact factorization (3) proves:

```text
Within the genuine nonsplit polynomial exact-square-prefix branch of
THM-2262/2297, no degree-eighteen survivor lies on

                 D/B^2=25/126.                                    (32)
```

The quantifiers in (32) matter. This closes one complete weighted
discriminant ratio. It does not close the other components of
`Disc_y(Disc_u G)=0`, the split/excluded branches outside THM-2262's scope,
the planar Jacobian conjecture, or the two-variable Dixmier conjecture.

## 9. Exact reproduction and information ledger

The standard-library companion is

```text
04-computation/jc2_degree18_central_t_singular_atlas_probe.py

05-knowledge/results/jc2_degree18_central_t_singular_atlas_probe.out. (33)
```

It reconstructs (6), checks the modular irreducibility certificates,
performs the quotient-field Euclidean algorithms in (11)--(22), and
checks (26)--(29). Normal and optimized Python executions must reproduce
the stored transcript byte for byte before promotion. The independently
replayed LF-normalized hashes are

```text
script:
c0408580c03cfdccfc593d48e982bb0e1b9a5351513930f63fd52c375839c04b

output:
44bbc875e444064724810aeac0c7880ad891a080c445e975bef72f2643058628. (33a)
```

```text
source:
  the explicit central factor T and the retained trigonal spectral
  cubic/flux/Faber package;

map:
  singular parameter point -> quotient-field critical section ->
  residual branch profile -> normalized trigonal cover;

preserved:
  every affine T-singularity, R/S face membership, branch multiplicity,
  the central total ramification, and the genuine deck;

destroyed:
  the original source-coordinate presentation and all information about
  the other discriminant ratios;

needed sidecar:
  first-flux square/deck for y!=0 and the third-flux/Keller/Faber
  polynomial obstruction for y=0;

new sharp target:
  audit the remaining rational B--D ratio D/B^2=4075/85176 from
  THM-2311/2313 on the quotient x=y^2 before pulling back to y. Its
  repeated-root resultant mixes the base-change factor delta(0) with
  Disc_x(delta)^2, so multiplicity in the y-resultant is not by itself a
  singular-normalization certificate. Factor delta(x), its derivative,
  and the singular ideal of H(u,x), then count ramification on the
  x-normalization before the quadratic base change.                 (34)
```
