# Complete sign on the original shared-root singular fibres

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
Every singular fibre of either shared-positive-root chart has strictly negative
original Laurent response whenever it meets the anchored model. This closes
whole coefficient intervals, including an explicitly nonempty high-phase
interval. General sign away from these singular fibres remains **OPEN**.

The proof needs only the C interlacer, including weak or repeated-root cases.
The prescribed D polynomial remains in the response and in the D shared-root
chart, but its interlacing is not required for the sign argument. Thus the
conclusion applies in particular to the full two-interlacer model. No converse
integer-support realization, general two-return theorem, or JC/LRC implication
is asserted.

[Independent proof and exact audit](continuing7_20260906_singular_fibres_audit.md) passes.

## 1. Inheritance, map, and statement

The closest proved mechanism is
[the original shared-wall chart and its admissible fibre](continuing6_20260906_shared_wall.md).
The positive-response hostile in
[the zero-beta sign theorem](continuing6_20260906_zero_boundary.md)
shows why both the original phase and the geometric moment packet must survive.
The corrected near miss is division by a vanishing coefficient of x: the
shared-wall note proves that this discards genuine model intervals. The
previously unused sidecar here is the ordinary four-by-four C moment
determinant, requiring moments through degree six.

The concept board is: the original Laurent carrier; singular chart algebra;
the low moment packet; the next Hankel determinant; full root admissibility;
and the quadratic response on a whole fibre. The map eliminates y,z by a shared
root r, then retains every simultaneous coefficient/constant zero of the
original phase. It preserves the phase and response exactly. Projecting to low
moments loses the actual interlacing predicate; the degree-six determinant is
the necessary repair on the one problematic fibre.

Put

```text
B(v)=v^5−13v^4+55v^3−xv^2+yv−z,
C(v)=v^4−12v^3+45v^2−(2x/3)v+3y/7,
D(v)=v^4−11v^3+36v^2−(5x/12)v+y/7.
```

Assume B has five nonnegative roots, counted with multiplicity, and C weakly
interlaces B. Their sum is13 and their square sum is59. Let r>0 be a common
root of B with A, where A is either C or D. Let s>0 satisfy the **original**
phase equation

`g(s)=zs^4−(12/7)ys^3+xs^2−10s+1/11=0`.

Define u=rs and

`A_C(u)=5u^2−24u+9`, `A_D(u)=23u^2−60u+12`.

**Theorem.** If the corresponding A_C(u) or A_D(u) vanishes, then Q(−s)<0,
where Q is the original carried second response defined below. The theorem
quantifies over every x on each such fibre. It concerns its singular original
phase s; it does not claim to handle the other, nonsingular phases of that
same coefficient shape.

## 2. Original coefficients and the necessary moment packet

Use the Laurent rows

```text
beta=t^−1+13+55t+xt^2+yt^3+zt^4,
c=t^−1+12+45t+(2x/3)t^2+(3y/7)t^3,
d=t^−1+11+36t+(5x/12)t^2+(y/7)t^3.
```

For O(t)=sum_j binom(14,2j+1)t^j and E(t)=sum_j binom(14,2j)t^j,
the first row is P=O star beta, with coefficientwise product star, and

`Q=(O^2+t^−1 E^2) star (beta^2+2tcd)`.

The source reconstructs both convolutions and independently checks the closed
weight identity `[t^j](O^2+t^−1 E^2)=binom(28,2j+2)`.
In particular, P(−s)=2002g(s), Q has full support −1 through8, and q_−1=28.
No zero of a transformed first polynomial is substituted.

Weak interlacing of monic C and B gives a positive probability measure on the
distinct B roots with Stieltjes transform C(v)/B(v). This remains valid at
repeated roots: weak interlacing cancels all but simple poles, and the remaining
residues are nonnegative; equivalently take a limit of strict interlacing
arrays. Hence every ordinary and shifted Hankel matrix of its moments is PSD.
Polynomial division at infinity gives

```text
m0=1, m1=1, m2=3,
m3=x/3−16,
m4=16x/3−373−4y/7,
m5=54x−59y/7+z−3969,
m6=x^2/3+1178x/3−568y/7+14z−31082.
```

The shifted two-by-two and ordinary three-by-three determinants are

```text
det(m_(i+j+1))_(i,j=0..1)=(x−75)/3,
det(m_(i+j))_(i,j=0..2)=(x−75)(135−x)/9−8y/7.
```

Together with nonnegative elementary coefficients, these imply

`x>=75`, `0<=y<=cap(x):=(7/72)(x−75)(135−x)`, `z>=0`.

The extra determinant used below is
`F4(x,y,z)=det(m_(i+j))_(i,j=0..3)>=0`.
The certificate contains its complete polynomial and the seven moments, all
rebuilt by formal division rather than imported from an earlier producer.

Finally every beta root M satisfies
`M^2+(13−M)^2/4<=59`, by Cauchy on the other four roots. Thus M<71/10;
the left side is already greater than59 at71/10 and is increasing thereafter.

## 3. Complete singular-pair enumeration

The shared-root equations give exactly

```text
Y_C=(14/9)xr−(7/3)r^2(r^2−12r+45),
Z_C=r^2[5x/9−r(4r^2−45r+150)/3],
Y_D=(35/12)xr−7r^2(r^2−11r+36),
Z_D=r^2[23x/12−r(6r^2−64r+197)].
```

After s=u/r, the original phase is
`x u^2 A_C(u)/(9r^2)−H_C(r,u)/(33r)` or
`x u^2 A_D(u)/(12r^2)−H_D(r,u)/(11r)`, respectively, where

```text
H_C=44r^2u^4−132r^2u^3−495ru^4+1584ru^3−3r
    +1650u^4−5940u^3+330u,
H_D=66r^2u^4−132r^2u^3−704ru^4+1452ru^3−r
    +2167u^4−4752u^3+110u.
```

Thus a singular original phase requires A=H=0. Their complete eliminants are

```text
P_C=705672r^4−15079284r^3+115450055r^2−375696750r+441317250,
P_D=120434688r^4−2310536448r^3+15142183583r^2
    −39453584808r+35022385200.
```

Their rational u reconstructions are

```text
U_C=(11286r^2−109085r+245025)/(27126r^2−261360r+586025),
U_D=(446688r^2−4023865r+7744176)/(1978416r^2−17739216r+34339426).
```

The source verifies the resultants up to nonzero rational constants, every
reconstruction identity modulo P, and coprimality of all reconstruction
denominators with P. Each quartic has exactly four positive roots. Ordered
by r, the eight pairs have the following coarse exact isolators (the two
listed numerators are divided by10^9):

| Pair | Lower numerator | Upper numerator | Outcome |
|---|---:|---:|---|
| C1 |3528453826|3528453827|negative on whole fibre|
| C2 |3550824804|3550824805|negative using degree-six determinant|
| C3 |6081020695|6081020696|negative on whole fibre|
| C4 |8208387541|8208387542|r exceeds beta maximum|
| D1 |2049574992|2049574993|negative on whole fibre|
| D2 |2831681376|2831681377|negative on whole fibre|
| D3 |6130326664|6130326665|already empty under low C packet and z>=0|
| D4 |8173391713|8173391714|r exceeds beta maximum|

Exact Sturm counts establish completeness. For coefficient sign evaluation,
the source refines each root to a verified decimal interval of width10^−28.
Every U is positive there. No x is divided out: at these pairs the original
phase holds identically for every x on the corresponding chart.
There is only one common u for each listed r. Indeed A has two distinct
roots, and its resultant with H is, up to a nonzero constant, the product
of H evaluated at those two roots. Two common u values at the same r
would make that resultant have a multiple r root, whereas each quartic
has four distinct roots. Thus the rational reconstruction misses no
second phase above a listed r.

## 4. Whole-fibre sign certificate

For each chart substitute y=Y_A, z=Z_A, s=U_A/r into the **complete** Q(−s).
Reduce each coefficient modulo P_A, with every denominator certified invertible.
The result is

`q_A(x;r)=a_A(r)x^2+b_A(r)x+c_A(r)`,

where a_A,b_A,c_A are rational polynomials of degree at most3 in r.
All six polynomials are saved explicitly in the certificate. Quotient
identities verify them exactly; root approximation is used only inside
rational sign enclosures.

First, for C1,C3,D1,D2 and D3, put f(x)=cap(x)−Y_A(x,r).
Its second derivative is −7/36. The exact certificate shows f(L)<0 and
f'(L)<0 at the following rational caps:

| Pair | L | Interval containing f(L) | Interval containing f'(L) |
|---|---:|---:|---:|
| C1 |94|(−2,−1)|(−4,−3)|
| C3 |89|(−3,−2)|(−7,−6)|
| D1 |102|(−4,−3)|(−6,−5)|
| D2 |98|(−5,−4)|(−7,−6)|
| D3 |95|(−4,−3)|(−16,−15)|

Consequently no feasible x is at least L. For D3, Z_D(95) lies in
(−102,−101) and its x slope lies in(72,73); therefore Z_D(x)<0 for every
x<95 as well. D3 has no feasible x.

For C2 use F4>=0. After replacing x by88+X and reducing its coefficients
modulo P_C, the five coefficients, in descending powers of X, lie in

`{−2/81}, (−2,−1), (−39,−38), (−304,−303), (−171,−170)`.

They are all strictly negative. Hence F4(88+X,Y_C,Z_C)<0 whenever X>=0,
so every feasible point on C2 has x<88. This is the entire additional use
of the degree-six packet; no D determinant is required.

The response certificate is now small. For each of the five remaining
pairs it proves a_A<0, q_A'(L)>0, and q_A(L)<0:

| Pair | L | Exact upper bound for q_A(L) |
|---|---:|---:|
| C1 |94|−4,139,408|
| C2 |88|−84,680,595,317|
| C3 |89|−5,979,289,215|
| D1 |102|−3,289,168|
| D2 |98|−6,096,143,057|

Since q_A''=2a_A<0, its derivative at every x<=L is at least its positive
derivative at L. Thus q_A is increasing throughout that whole ray, and
q_A(x)<=q_A(L)<0. This proves the theorem, including weak endpoints of
the geometric domain. The rational caps are deliberately sufficient bounds,
not claims to compute the exact admissible interval endpoints.

## 5. The failed projection and an actual full-model interval

Let r be the C2 algebraic root, with u=U_C(r), s=u/r, and set x=92,
y=Y_C(92,r), z=Z_C(92,r). These are exact algebraic data. The shared-root
and original phase equations hold identically. They satisfy x>75, y>0,
z>0 and cap(x)−y in(4,5). In fact both C and D ordinary three-by-three
and shifted three-by-three moment matrices are **positive definite**:
every leading principal minor is certified strictly positive. These are
the two degree-five Stieltjes packets.

Nevertheless the same original response satisfies

`37,071,436,949<Q(−s)<37,071,436,950`,

while `−2117<F4<−2116`. This is a hostile to dropping the degree-six
constraint, not a counterexample in the anchored model. It does not claim
to pass every possible degree-five upper-support localizing test, nor to
establish degree-six minimality among all conceivable encodings.

The repaired fibre is genuinely populated by full-model points. Keep the
same r,u and allow

`86999/1000 <= x <= 87001/1000`.

Uniform rational endpoint signs give the following root brackets; listed
integers are divided by1000:

```text
B/(v−r): (80,100), (560,610), (2440,2470), (6320,6340),
C/(v−r): (400,410), (1910,1930), (6120,6140),
D:       (180,200), (1470,1500), (3380,3400), (5930,5950).
```

The number of disjoint sign brackets equals each quotient's degree. Thus
all these roots are simple, real, positive, and exhaust the polynomial.
Their ordering, with r between3.55 and3.56, proves C weakly interlaces B
with its one prescribed shared root, while D strictly interlaces B.
The source also verifies y,z>0 uniformly. The anchors follow from B's
coefficients. Every point of this explicit interval is therefore in the
full model, at the same singular original phase s approximately1.23632539,
and the whole-fibre theorem applies.

## 6. Consequence and stopping boundary

The existing full-model simplicity theorem excludes repeated beta roots;
the zero-boundary theorem closes beta zero. Combined with the incoming
fixed-original-phase extremum reduction, this theorem removes all singular
denominator fibres from the remaining shared-root boundary search. If a
nonnegative response exists anywhere in the full model at some phase s,
there must be a **possibly different** shared-root model point at that same
phase with nonnegative response. Such a point must lie in a chart where
the corresponding A_C(rs) or A_D(rs) is nonzero. This is an existence/test
reduction, not an assertion that every possible counterexample is itself
on the boundary.

The general shared-root chart off these finite singular fibres is still
OPEN. The theorem supplies no uniform sign on that two-parameter remainder.

## 7. Reproduction and frozen artifacts

- [Standalone source](../../04-computation/continuing7_20260906_singular_fibres.py)
- [Normal output](continuing7_20260906_singular_fibres.out)
- [Complete rational certificate](continuing7_20260906_singular_fibres_certificate.json)

```text
python -B 04-computation/continuing7_20260906_singular_fibres.py
python -B -O 04-computation/continuing7_20260906_singular_fibres.py
```

Both runs pass167 always-active exact gates. Python's stdout newline is
explicitly configured to LF; the frozen outputs were captured as raw
subprocess bytes, compared byte-for-byte, and checked to contain no CR.
SymPy supplies exact polynomial division, resultants, Sturm counts and
rational root refinement. All sign conclusions use standard-library
Fraction interval arithmetic, not numerical root signs. The exhaustive
universe is the four positive singular pairs in each of two charts, not
a coefficient or geometric sample census. The source writes its JSON
alongside itself outside the repository, or under05-knowledge/results
when filed under04-computation.

SHA256:

```text
source 13897869eb93ca25e2af4e9125cd7610c1734ae6228edeb6a3e22e8b6f75cef5
output 12e7d6b17369eb452edc99a27124a7573a6e261b7dcd0eb991fdf88a2c43938a
JSON   34f178da766757e31968103e47543d71ab8e975100a8ac0d1f4dd1d69c027ee9
```

The independent audit above passed. The session manifest pins the filed
source, proof, output and complete certificate.
