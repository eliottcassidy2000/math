---
id: THM-2636
title: "Degree-twenty-two B-C-D triple spectral-square Kummer closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the inherited genuine nonsplit degree-twenty-two branch on the open
  first-flux chart, the support-three stratum B,C,D nonzero and E=W=0 has a
  44-term weighted-scale eliminant, monic of degree five in v.  Its t=0
  section is the irreducible squarefree quintic L_5.  The five roots and ten
  unordered root pairs are one degree-five and one degree-ten field.  The
  unique Hensel lifts have five order-eleven equations supported on
  c^3,c*d^2,c*d,c; in both fields the resulting 5-by-4 coefficient matrix
  has rank four uniformly, excluding every line and quadratic factor.  The
  retained base-field square T^2=Z=rho^3*zeta/t^3 is odd at all five smooth
  fixed points, so its connected double cover has at least six branch places
  and genus at least two.  The y=0 boundary reduces to a nonzero quartic.
  Hence the entire B-C-D stratum is empty.  No other support-three stratum,
  JC(2), or DC(2) follows.
source: jc-degree22-2026-07-28-bcd-hensel-kummer
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
related:
  - THM-2480-degree-twenty-two-BC-plane-hensel-ramification-closure
  - THM-2617-degree-twenty-two-BDW-triple-fixed-section-and-last-quadratic-type
script: 04-computation/jc2_degree22_bcd_weighted_hensel_kummer_thm2636.py
output: 05-knowledge/results/jc2_degree22_bcd_weighted_hensel_kummer_thm2636.out
script_sha256: 0866a29f665aedc6d2c226f35943852e56907ff821e705a0dbca2651e71fa15c
output_sha256: 52901047964006429d5b8b243e455b9025f4186c56a29777c97d9e090d50ad1d
hash_basis: working-tree bytes (LF)
---

# THM-2636 -- the B-C-D triple closes by its retained spectral square

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The support-two coefficient planes are closed, and THM-2617 closes the first
support-three torus `B,D,W != 0`.  The next coprime-weight triple is

```text
B,C,D != 0,                         E=W=0.             (1)
```

The odd weight of `C` destroys THM-2617's coordinate `p=B/y^2`: the
eliminant sees a genuine signed scale `t`, not only `t^2`.  Nevertheless the
spectral coordinate retains a different square class.  Once the signed-scale
eliminant is proved irreducible, its five fixed-section points are all odd
places of `Z=T^2`.  This supplies the positive-genus lift without any
two-parameter discriminant classification.

## 1. Inherited coordinates and weighted scale

Use THM-2411's target-translated coordinates

```text
y=11s,                    u=d_0 T,                    Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).              (2)
```

The genuine nonsplit branch has `T=q^2 != 0`; `T` belongs to the base
rational function field, while `q` carries the quadratic deck character.
Work on the inherited open chart

```text
mathcal A=616B-1089u+63y^2 != 0.                     (3)
```

First suppose `y != 0` as a function.  Choose the constant `rho in C*` with
`rho^2=B` and put

```text
c=C/rho^3,       d=D/rho^4,       t=rho/y,
v=u/y^2,         zeta=Z/y^3.                         (4)
```

Both `c,d` are nonzero constants.  Dividing the two inherited fluxes by
`y^5,y^6` gives

```text
F_1=
 819896t^2zeta-1449459v zeta+83853zeta
 -2981440t^2v+24640t^2
 +9370240c t^3v-232320c t^3+2044416d t^4
 +3689532v^2-101640v+252,                            (5)

F_2=
 15944049zeta^2+65591680t^2zeta-206145280c t^3zeta
 -162339408v zeta+2236080zeta
 +1443016960t^2v^2-71554560t^2v+98560t^2
 +449771520c t^3v-1239040c t^3
 -1978994688d t^4v+16355328d t^4
 -1190488992v^3+147581280v^2-1219680v+672.           (6)
```

The open chart becomes

```text
616t^2-1089v+63 != 0,                                (7)
```

so (5) reconstructs `zeta` in the function field.

## 2. Universal eliminant and the fixed quintic

For a hostile normalization check, the companion first retains all five
dimensionless coefficients

```text
C=c rho^3, D=d rho^4, E=e rho^5, W=w rho^6.          (8)
```

The primitive universal resultant `R_univ=prim Res_zeta(F_1,F_2)` has

```text
integer content before primitivization = 28344976,
number of terms                         = 60,
degree in (t,v)                         = (10,5),
t-support                               = {0,2,3,4,5,6,7,8,9,10}.          (9)
```

After `e=w=0`, write `R=R_BCD`.  It has 44 terms, degree nine in `t`, degree
five in `v`, and Newton polygon

```text
conv{(0,0),(9,0),(8,1),(0,5)}.                       (10)
```

Its leading `v` coefficient is the nonzero constant

```text
[v^5]R=-88239118492602.                              (11)
```

Thus `P=R/[v^5]R` is monic in `v`.  At `t=0`, equivalently before monic
normalization,

```text
R(0,v)=-567L_5(v),                                  (12)

L_5
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                        (13)
```

Exact factorization proves `L_5` irreducible over `Q`; it is squarefree.
There is no order-one term in `t`.

## 3. Every possible factor has a fixed root or root pair

Suppose `P=QS` is a nontrivial factorization over `C`.  Since `P` is monic
in `v` with constant leading coefficient, both factors can be normalized
monic in `v`; neither can have `v`-degree zero.  Their `v`-degrees add to
five and remain unchanged at `t=0`.  After exchanging `Q,S`, the fixed
factor therefore has degree one or two.

For degree one, put

```text
A_1=Q[gamma]/(L_5(gamma)),       q_0=v-gamma.         (14)
```

This is one irreducible degree-five field whose embeddings are the five
roots.  For degree two, divide `L_5` by

```text
q_0=v^2+alpha v+beta.                                 (15)
```

The two remainder coefficients have a triangular Groebner basis: `beta` is
a polynomial in `alpha`, and `alpha` satisfies an irreducible squarefree
degree-ten polynomial.  Hence

```text
A_2=Q[alpha]/(P_10(alpha))                            (16)
```

is one reduced degree-ten field whose embeddings are precisely the ten
unordered pairs of distinct roots of `L_5`.  The complementary fixed
cofactors have degrees four and three.  Thus (14)--(16) omit no factor
component or factor/cofactor swap.

## 4. One order-eleven obstruction closes both fields

Normalize (12) to its monic form and expand

```text
P(t,v)=sum_(n=0)^10 r_n(v)t^n,
Q(t,v)=sum_(n>=0)q_n(v)t^n,
S(t,v)=sum_(n>=0)s_n(v)t^n.                          (17)
```

In either field `A_k`, `q_0s_0=r_0` and squarefreeness makes `q_0,s_0`
coprime.  Hensel uniqueness gives the exact recursion

```text
f_n=r_n-sum_(i=1)^(n-1)q_i s_(n-i),

q_n = rem_(q_0)(f_n (s_0 mod q_0)^(-1)),
s_n = (f_n-q_n s_0)/q_0.                             (18)
```

The companion performs (18) in the sparse ring `A_k[c,d][v]`, checks every
division, and directly multiplies the reconstructed factor and cofactor at
all eleven orders.  At order eleven the two factor coefficients together
with the complementary cofactor coefficients give five equations.  In both
fields their complete parameter support is

```text
{c^3, c d^2, c d, c}.                                (19)
```

Since `c != 0`, divide by `c` and write the equations as

```text
M_k (c^2,d^2,d,1)^T=0,              M_k in Mat_(5x4)(A_k).                (20)
```

Exact maximal-minor reduction gives

```text
                    A_1 root field       A_2 pair field
field degree                 5                    10
nonzero 4x4 minors           2                     2
first numerator degree       4                     9
first numerator terms        5                    10
gcd(numerator,modulus)        1                     1.                    (21)
```

Thus `rank M_k=4` at every embedding of either field.  Equation (20) would
force `(c^2,d^2,d,1)=0`, impossible.  A polynomial factor of `P`, whose
`t`-degree is nine, has zero factor and cofactor coefficients by order
eleven; hence (20) is a necessary condition for every actual factorization.
There are no line or quadratic factors.  Section 3 exhausts the possible
smaller side, so

```text
P_(c,d)(t,v) is absolutely irreducible for every c != 0.                 (22)
```

In particular this is uniform on the physical locus `c,d != 0`.

## 5. Five fixed points are odd places of the spectral square

Let `C_(c,d)` be the smooth projective normalization of (22).  The five
points above `t=0` are distinct and smooth by squarefreeness of `L_5`.
Since `L_5'` is nonzero at each root, the implicit-function criterion alone
makes `t` a uniformizer at all five.  The vanishing order-one coefficient in
(17) is a separately verified stronger tangency control, not the reason that
`t` is a parameter.

At `t=0`, the first flux is

```text
(83853-1449459v)zeta
 +(3689532v^2-101640v+252)=0.                         (23)
```

Both displayed coefficient polynomials have degree below five and are
coprime to the irreducible `L_5`.  Therefore `zeta` is a unit at every one
of the five points.

Now use the correctly typed retained square.  From (4),

```text
rho^3 zeta/t^3=Z=T^2,                                (24)
```

where `T` is in the base rational function field.  Put

```text
H=rho^3 zeta/t^3 in C(C_(c,d)).                       (25)
```

At each fixed point, `ord(H)=-3`.  Hence `H` is not a square, the quadratic
Kummer cover

```text
D_(c,d): X^2=H                                      (26)
```

is connected, and it ramifies at all five points.  A quadratic branch
divisor has even degree, so (26) has at least six branch places.  Without
assuming any genus for `C_(c,d)`, Riemann--Hurwitz gives

```text
2g(D)-2=2(2g(C)-2)+R >= -4+6=2,

g(D)>=2.                                             (27)
```

This is the missing uniform genus gate.  It uses the spectral square
`T^2=Z`, not the fourth power involving the deck-character coordinate `q`.

## 6. The nonzero-y trajectory is impossible

A rational Keller trajectory satisfies (5)--(7), (22), and (24).  If `t`
is nonconstant, absolute irreducibility embeds the function field of
`C_(c,d)` in `C(x)`, and the physical element `T` realizes the connected
cover (26) inside `C(x)`.  This gives a nonconstant rational map from
`P^1` to the genus-at-least-two normalization `D_(c,d)`, impossible.

If `t` is constant, then `y` is constant.  Equation (22) makes `v`, hence
`u`, algebraic over the algebraically closed constant field and therefore
constant.  The open first flux reconstructs constant `zeta`; then (24)
makes `Z,T,q` constant, contradicting the genuine nonsplit deck.  Thus the
branch `y != 0` is empty.

## 7. The identically-zero boundary

The case `y=0` is treated before defining `t`.  A constant weighted scaling
normalizes `B=1,C=c,D=d`.  Put

```text
A=616-1089u != 0.                                    (28)
```

The first two fluxes become

```text
N_1=1331 A Z+9370240c u,

N_2=15944049Z^2-206145280cZ+1443016960u^2
     -1978994688d u-1190488992u^3.                   (29)
```

Eliminating `Z` through `N_1=0`, multiplying by `(1331A)^2`, and taking
primitive part gives exactly

```text
-u Q_4(u),                                           (30)

Q_4=
 2264031u^4-5305608u^3+(3829056+3763584d)u^2
 +(1267200c^2-4257792d-878080)u
 -1433600c^2+1204224d.                               (31)
```

The leading coefficient `2264031` is constant and nonzero.  If `u=0`,
(29) gives `Z=0`, contrary to `T!=0`.  Otherwise (31) makes `u` algebraic
over `C`, hence constant; the open equation reconstructs constant `Z`, and
then `T,q` are constant, again contradicting the deck.  The boundary is
empty.

Combining Sections 6--7 proves that (1) has no genuine degree-twenty-two
trajectory on the inherited branch.

## 8. Exact reproduction and scope

Run

```bash
python3 04-computation/jc2_degree22_bcd_weighted_hensel_kummer_thm2636.py
python3 -O 04-computation/jc2_degree22_bcd_weighted_hensel_kummer_thm2636.py
```

The companion reconstructs (5)--(13), the fixed root and unordered-pair
fields, the shared Hensel recursion through order eleven, all direct product
controls, every maximal-minor coprimality in (21), the exact squared
first-flux wall, both unit checks in (23), the Kummer branch/genus invoice,
and the boundary factor (30)--(31).  All truth-bearing checks use explicit
exceptions and remain active under optimized Python.

An independent hostile audit checked the three logical gates separately.
First, monicity and the constant `v`-leading coefficient forbid a
nonconstant `v`-degree-zero factor and preserve fixed-section degree, so the
root/pair atlas is exhaustive.  Second, degree in `t` is additive under a
polynomial factorization; since the BCD specialization has degree nine,
both factor series vanish by order eleven and the rank-four equations are
necessary.  Third, nonconstant `t` makes the function-field map injective;
the five odd valuations make `H` nonsquare, and the base-field element `T`
embeds the connected quadratic cover in `C(x)`.  The audit also separated
the squarefree uniformizer argument from the stronger missing-order-one
tangency check.

This theorem closes only the physical support-three stratum (1) inside the
inherited genuine nonsplit degree-twenty-two reduction.  It does not close
another support-three pattern, a split/even short-edge branch, integral
order raising, `JC(2)`, or `DC(2)`.

**QED.**
