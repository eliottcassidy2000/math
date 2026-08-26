---
id: THM-4233
title: "Pair-specific primitive-observable oscillation Haar charts"
status: >
  PROVED ANALYTIC REDUCTION + FINITE-EXACT + INDEPENDENTLY AUDITED FIXED-PAIR
  CERTIFICATE. For the displayed thirty-label pool, every nine-body is
  Haar-safe after adjoining every common dilation of the coprime pair
  (3713,5149). A second analytic chart proves the same conclusion for every
  dilation of (5k+2,7k+3), k>=180785. The adjacent coprime family has limiting
  primitive oscillation 36/2401 and is a rigorous hostile control. Arbitrary
  finite pair entry and LRC(14) remain OPEN.
source: codex-lrc14-pair-specific-oscillation-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
related:
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
primary_script: 04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.cpp
primary_output: 05-knowledge/results/lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.out
independent_audit_script: 04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_independent_audit_thm4233.cpp
independent_audit_output: 05-knowledge/results/lrc14_pair_specific_primitive_observable_gcd_one_chart_independent_audit_thm4233.out
primary_script_sha256: 9c02b42562d3b911e9df2a553844c6e1acb1ab6d64dee176d9d76c861f54dcc9
primary_output_sha256: 03e0c222ee7ba9e2384de76c33217c05887990c64ec24c7635ca59bc05b9df07
independent_audit_script_sha256: b6691337c81199e7459fbabfeec10e49d0899ee40368d7b1d729deabd7a0ad20
independent_audit_output_sha256: e43dc469c08d89a5e5ca20e0c734bf56f44d35aa84d5d4cac8ef304464487253
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary integer common-grid sweep and an independent
  reduced-rational safe-teeth intersection/integration audit reproduce the
  fixed-pair mass, primitive extrema, oscillation, threshold margins,
  collision ticks, and hostile controls. Their independently serialized
  semantic ledgers are 4725493ac81fe903 and 635148be446ddc28. O0/O3 replays
  match; the independent path also passes ASan/UBSan.
---

# THM-4233 -- pair-specific primitive-observable oscillation Haar charts

**PROVED ANALYTIC REDUCTION + FINITE-EXACT + INDEPENDENTLY AUDITED FIXED-PAIR
CERTIFICATE; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Retain the thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive label set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

> **Fixed primitive-ratio ray.** For every integer `g>=1` and every
> `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {3713g,5149g}))>=4/63.                  (3)
> ```

In particular, `g=1` is a genuinely coprime comparable pair. Both labels are
below the [THM-4231](THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift.md)
quadrant, their gcd is below the
[THM-4228](THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray.md)
gate, and the pair is outside the
[THM-4227](THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge.md)
scale-separated wedge.

There is also an analytic infinite chart. For `k>=1`, define

```text
u_k=5k+2,                 v_k=7k+3.                      (4)
```

> **Canceled-resonance ray family.** For every `k>=180785`, every integer
> `g>=1`, and every `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {g u_k,g v_k}))>=4/63.                  (5)
> ```

The conclusion `(5)` is numerically subsumed by THM-4231, since already
`u_k,v_k>17548`. Its new content is the cyclotomic-zero mechanism that makes
the primitive-observable oscillation `O(1/k)` without a common divisor.

For every body in `(3)` or `(5)`,
[THM-4150](THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer.md)
gives, for every positive integer `c` and every two distinct positive odd
integers `a,b`, an `x in R/Z` such that

```text
min_(w in 2c(B union {gu,gv}) union {a,b})||wx||>=1/14, (6)
```

where `(u,v)` is respectively `(3713,5149)` or `(u_k,v_k)`. These are
infinite thirteen-speed LRC(14) families, not a proof of LRC(14).

The closest proved mechanism is THM-4228's periodic-observable discrepancy
lemma and its exact depth-eight deck `E_3467`. The canonical hostile is the
adjacent coprime family `(n,n+1)`, whose primitive oscillation does not tend
to zero. The corrected near miss is to interpret failure of a sufficient
oscillation gate as literal Haar danger. The least-used sidecar is the
centered primitive of the **joint** two-comb observable, rather than either
marginal comb. The live concept board is

```text
depth-eight repair deck | pair density | primitive oscillation
endpoint collisions | slow-fast averaging | cyclotomic zeros.             (7)
```

## 2. Pair-specific transfer through the THM-4228 deck

Let

```text
j(t)=1_(||t||>=1/14),
A_(u,v)=G_u intersect G_v,
beta_(u,v)=mu(A_(u,v)),                                  (8)
```

where changing the values of `j` at its two endpoints is immaterial. Define
the one-periodic centered primitive

```text
H_(u,v)(t)=integral_0^t(1_(A_(u,v))(s)-beta_(u,v))ds,
omega_(u,v)=max H_(u,v)-min H_(u,v).                    (9)
```

If `U` is a union of `c_U` positive-length circular intervals, the oriented-
lift proof of THM-4228, equation `(24)`, gives for every integer `g>=1`

```text
mu(U intersect G_(gu) intersect G_(gv))
 >=beta_(u,v)mu(U)-c_U omega_(u,v)/g.                   (10)
```

There is no independence assumption in `(10)`: it is the identity
`G_(gu) intersect G_(gv)=m_g^(-1)(A_(u,v))`, followed component by component
by the endpoint cocycle of `H_(u,v)`.

Put

```text
beta_0=66/91,
T=(1650/8281)/3467=1650/28710227.                       (11)
```

For `R in binom(P,8)`, write

```text
U_R=G_(P\R),       M_R=mu(U_R),       c_R=#components(U_R). (12)
```

THM-4228 constructs a finite deck `E_3467` such that

```text
tau(E_3467)>9,
beta_0 M_R-c_R T>=4/63              for every R in E_3467. (13)
```

Consequently, any distinct positive pair `(u,v)` satisfying

```text
beta_(u,v)>=beta_0,             omega_(u,v)<=T           (14)
```

supplies every dilation `g>=1`. Indeed, for any `B in binom(P,9)`, `(13)`
provides `R in E_3467` disjoint from `B`. Equations `(10)`, `(13)`, and `(14)`
give

```text
mu(G_((P\R) union {gu,gv}))>=4/63.                     (15)
```

Since `B union {gu,gv}` is contained in the label set in `(15)`, safe-set
monotonicity proves the corresponding conclusion for `B`. More generally,
the same proof works when `omega_(u,v)/g<=T`; `(14)` is the convenient gate
that includes the primitive dilation `g=1`.

## 3. Exact coprime certificate for `(3713,5149)`

Set

```text
u=3713,                 v=5149,                 gcd(u,v)=1,
N=14uv=267655318.                                         (16)
```

Both `u` and `v` exceed `max(P)=290`. On the integer grid `y=t/N`, the
circular endpoint multisets of `G_u` and `G_v` are

```text
v(14i+1), v(14i-1) mod N,        0<=i<u,
u(14j+1), u(14j-1) mod N,        0<=j<v.                  (17)
```

There are `17724=2(u+v)` raw circular entries. Because
`u+v=8862=14(633)`, precisely two endpoint pairs coincide, at ticks

```text
95591185,                 172064133.                      (18)
```

These are null collisions, not positive-length cells. After adjoining the
linear sentinels `0,N`, the sweep has `17726` raw points, `17724` unique
points, and `17723` open cells. It finds

```text
safe cells=7595,          unsafe cells=10128,
safe components=7595,
S=safe ticks=196644720.                                  (19)
```

Thus

```text
beta_(u,v)=S/N=98322360/133827659,
beta_(u,v)-66/91=16387038/1739759567>0.                 (20)
```

For completeness, the exact primitive computation can be stated without any
floating-point convention. Let

```text
0=e_0<e_1<...<e_m=N
```

be the unique sorted linear points from `(17)` and the two sentinels, put
`d_i=e_(i+1)-e_i`, and let `a_i in {0,1}` be the joint-safe value on that
open cell. Starting with `C_0=0`, set

```text
C_(i+1)=C_i+d_i(Na_i-S).                                (21)
```

Then `H_(u,v)(e_i/N)=C_i/N^2`; since `H_(u,v)` is affine on each cell, its
extrema occur among these endpoint values. The exhaustive recurrence `(21)`
gives

```text
min C_i=-2057535171948 at t=207775767=55959N/72086,
max C_i= 2057535171948 at t= 59879551=16127N/72086,
N^2=71639369253681124,                                  (22)

min H=-138535899/4823550313337,
max H= 138535899/4823550313337,
omega_(u,v)=277071798/4823550313337.                    (23)
```

The decisive strict margin is

```text
T-omega_(u,v)
 =82934716896/2826229070241355051>0.                    (24)
```

Equations `(20)`, `(23)`, and `(24)` verify `(14)`. Section 2 proves `(3)`
for every `g>=1`, and THM-4150 proves `(6)`. This is the fixed-pair theorem.

The exact artifacts retain the hostile and positive controls

```text
(1,13):       beta=66/91, omega=990/8281,
(1,2):        beta=11/14, omega=11/98,
(2584,4181):  omega=42967355/553336008694>T.             (25)
```

The last line is a certificate-gate failure, not an unsafe-set witness. The
primary semantic ledger is `4725493ac81fe903`; the independent ledger is
`635148be446ddc28`. The primary path is the integer-grid recurrence
`(17)--(23)`. The audit path independently intersects reduced rational safe
teeth and integrates the resulting rational primitive. Neither path uses
floating point, imports the other's geometry, or treats the pair-selection
process as a proof of minimality.

## 4. The adjacent coprime family is a hostile limit

Let

```text
A_n=G_n intersect G_(n+1),
beta_n=mu(A_n),
H_n(t)=integral_0^t(1_(A_n)(s)-beta_n)ds.               (26)
```

Then

```text
lim_(n->infinity) osc(H_n)=36/2401.                     (27)
```

Here is a direct slow-fast proof. Define the autocorrelation

```text
C(z)=integral_(R/Z) j(x)j(x+z)dx.                       (28)
```

The danger arc has length `1/7`, so, with `||z||` the circle distance,

```text
C(z)=5/7+(1/7-||z||)_+,
integral_(R/Z)C(z)dz=36/49.                             (29)
```

On a cell `[a/n,(a+1)/n]`, substitute `t=(a+x)/n`. Apart from null walls,

```text
j(nt)j((n+1)t)=j(x)j(x+t).                              (30)
```

Replacing `t=(a+x)/n` in the second factor by `a/n` costs at most `4/n` in
the `x`-integral: if the two values differ, `x+a/n` lies within `1/n` of one
of the two walls of `j`. After the outer factor `1/n`, each complete cell
therefore costs at most `4/n^2`. Summing complete
cells, treating one final partial cell identically, and using the uniform
Riemann convergence of the continuous piecewise-linear function `C` proves

```text
integral_0^t j(ns)j((n+1)s)ds
 -> integral_0^t C(s)ds                                 (31)
```

uniformly for `0<=t<=1`. At `t=1`, `(31)` also gives `beta_n->36/49`.
Therefore `H_n` converges uniformly to

```text
K(t)=integral_0^t(C(s)-36/49)ds.                        (32)
```

On the first side of the triangular peak,
`C(s)-36/49=6/49-s`; its positive area on `[0,6/49]` is
`18/2401`. By reflection symmetry around zero, `K` has maximum `18/2401`,
minimum `-18/2401`, and oscillation `36/2401`. Oscillation is continuous
under uniform convergence, proving `(27)`.

In particular,

```text
36/2401-T=21012378/1406801123>0.                        (33)
```

Thus primitive comparable pairs do not acquire the tiny THM-4228 oscillation
gate merely by moving to large adjacent labels. This identifies why a pair-
specific cancellation, not size alone, is needed below the common-gcd or
cofinal-quadrant routes.

## 5. A cyclotomic-zero family with vanishing oscillation

### 5.1. Exact zero fast mean

For the pair `(4)`, note first that

```text
7u_k-5v_k=-1,                                           (34)
```

so every pair is coprime. Put

```text
F(x,y)=j(5x+2y)j(7x+3y)-36/49.                          (35)
```

The Fourier coefficients of `j` are

```text
c_0=6/7,
c_m=-sin(pi m/7)/(pi m)       for m!=0.                 (36)
```

In the `x`-average of the product in `(35)`, a Fourier pair `(a,b)` can
survive only when `5a+7b=0`, hence `(a,b)=(7h,-5h)`. If `h!=0`, its first
coefficient is `c_(7h)=0`; if `h=0`, the contribution is `c_0^2=36/49`.
Fejer approximation justifies the calculation for the step function. Thus

```text
integral_0^1 F(x,y)dx=0                 for every y.    (37)
```

This is the load-bearing cyclotomic cancellation.

### 5.2. Exact density quasipolynomial

Let

```text
L_k=14u_kv_k,
k=7q+r,                  0<=r<=6,                       (38)
d_0,...,d_6=(16,0,-16,-18,-6,6,18).                    (39)
```

Then the joint density is exactly

```text
beta_(u_k,v_k)=36/49+d_r/(7L_k).                        (40)
```

To prove `(40)`, let `D_m=(R/Z)\G_m` up to null endpoints and write
`gamma_k=mu(D_(u_k) intersect D_(v_k))`. The exact two-comb overlap formula
from THM-4228 is

```text
L_k gamma_k=W_k,
W_k=w(0)+2sum_(t>=1)w(t),
w(t)=[min(2u_k,u_k+v_k-14t)]_+.                         (41)
```

For `k=7q+r`, the plateau ends at

```text
floor((v_k-u_k)/14)=q,                                  (42)
```

and the positive support ends at

```text
h=6q+b_r,
(b_0,...,b_6)=(0,1,2,2,3,4,5).                         (43)
```

Consequently,

```text
W_k=2u_k+2[q(2u_k)+sum_(t=q+1)^h(u_k+v_k-14t)].         (44)
```

Substitution of `k=7q+r` into the finite arithmetic sum gives, residue by
residue,

```text
7W_k-2u_kv_k=d_r.                                       (45)
```

Finally, inclusion-exclusion and `L_k=14u_kv_k` give

```text
beta_(u_k,v_k)-36/49
 =5/7+W_k/L_k-36/49
 =(7W_k-2u_kv_k)/(7L_k),                                (46)
```

which is `(40)`. Since `|d_r|<=18`, `L_k>=980` for `k>=1`, and

```text
36/49-66/91=6/637>9/3430>=18/(7L_k),                   (47)
```

every member of the family satisfies `beta_(u_k,v_k)>beta_0`.

### 5.3. Primitive bound and the explicit tail

By `(37)`, define the one-periodic fast primitive

```text
P(x,y)=integral_0^x F(s,y)ds.                           (48)
```

For fixed `y`, the uncentered indicator in `(35)` has mean `36/49`.
THM-4228's sharp primitive estimate therefore gives

```text
osc_x P(.,y)<=(36/49)(13/49)=468/2401.                 (49)
```

Moreover, translating the first safe comb by `2 delta` changes it on measure
at most `4|delta|`, and translating the second by `3 delta` changes it on
measure at most `6|delta|`. Hence

```text
|P(x,y)-P(x,y')|<=10|y-y'|.                             (50)
```

The function `P` is Lipschitz and piecewise affine; hence its absolutely
continuous restriction to the line `(kt,t)` obeys the chain rule away from
finitely many wall crossings, and integration yields

```text
integral_0^t F(ks,s)ds
 =P(kt,t)/k-(1/k)integral_0^t partial_y P(ks,s)ds.      (51)
```

Because `P(0,y)=0`, `(49)` implies `|P(x,y)|<=468/2401`. The first term on
the right of `(51)` therefore has oscillation at most `936/(2401k)`, while
the primitive of `partial_y P` has oscillation at most `10`. Centering from
`36/49` to the exact density `(40)` adds at most
`|d_r|/(7L_k)`. Thus

```text
omega_(u_k,v_k)
 <=24946/(2401k)+18/(7L_k)
 <=24946/(2401k)+9/(1715k^2).                           (52)
```

The right side decreases with `k`. At `k=180785`, it is

```text
22549313113/392362010781125,                            (53)
```

and the exact threshold margin is

```text
T-22549313113/392362010781125
 =28971847951/229893926442909103375>0.                  (54)
```

Equations `(47)`, `(52)`, and `(54)` verify the pair-specific gate `(14)`
for every `k>=180785`. Section 2 proves `(5)` for every common dilation
`g>=1`, and THM-4150 gives the odd-tail consequence `(6)`.

### 5.4. General unimodular resonance-zero template

The cancellation above is not isolated. Let `a,c` be positive integers and
`b,d` be integers such that

```text
ad-bc=+-1,                  and                  7|a or 7|c. (55)
```

For every integer `k>=1` for which

```text
u=ak+b,                     v=ck+d                         (56)
```

are distinct and positive, they are coprime, because every common divisor
divides `cu-av=cb-ad=+-1`. Define

```text
F_(a,b,c,d)(x,y)=j(ax+by)j(cx+dy)-36/49.                (57)
```

The `x`-resonances are `(ch,-ah)`. Condition `(55)` makes one of their two
nonzero Fourier coefficients a multiple-of-seven coefficient of `j`, hence

```text
integral_0^1 F_(a,b,c,d)(x,y)dx=0.                       (58)
```

The same primitive proof as `(48)--(52)` has fast oscillation at most
`468/2401` and slow Lipschitz constant at most `2(|b|+|d|)`. On the other
hand, Fourier orthogonality in the two danger combs forces the modes
`(vh,-uh)`; product-to-sum and
`sum_(h>=1)cos(2pi hz)/h^2=pi^2 B_2({z})` give

```text
beta_(u,v)-36/49
 =[B_2({(v-u)/14})-B_2({(v+u)/14})]/(uv),               (59)
```

where braces denote fractional part and `B_2(z)=z^2-z+1/6` on `[0,1)`.
Since
`-1/12<=B_2<=1/6`, equations `(58)--(59)` prove the uniform compiler bounds

```text
|beta_(u,v)-36/49|<=1/(4uv),
omega_(u,v)
 <=[936/2401+2(|b|+|d|)]/k+1/(4uv).                    (60)
```

Whenever the lower density and upper oscillation in `(60)` pass `(14)`, the
THM-4228 deck supplies every common dilation. The family `(4)` is the
instance `(a,b,c,d)=(5,2,7,3)`; its exact quasipolynomial `(40)` improves the
last term of `(60)` and yields the displayed cutoff.

## 6. Exact audit and reproduction

The fixed-pair computation has the explicit finite universe `(17)--(21)`:
all joint endpoints, all resulting cells, and every primitive prefix. Its
inherited filters are exactly the two safe-comb predicates; there is no
search truncation or floating comparison. The positive control `(1,2)`, the
sharp-density hostile `(1,13)`, and the failed-gate Fibonacci predecessor in
`(25)` test the density, primitive normalization, and inequality direction.

The primary implementation uses the common integer grid and midpoint cell
classification. The independent implementation builds reduced-rational safe
teeth, intersects intervals directly, and integrates its own rational
primitive. Reproduce from the repository root with

```bash
clang++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_thm4233.cpp \
  -o /tmp/thm4233-primary
/tmp/thm4233-primary

clang++ -std=c++20 -O3 -DNDEBUG \
  04-computation/lrc14_pair_specific_primitive_observable_gcd_one_chart_independent_audit_thm4233.cpp \
  -o /tmp/thm4233-independent
/tmp/thm4233-independent
```

Both output streams must byte-match their frozen files.
Repeating both compilations with `-O0` gives the same two output streams; the
independent executable also passes combined AddressSanitizer/UBSan replay.

## 7. Connection contract, scope, and open boundary

```text
source:       one primitive pair observable A_(u,v) and the THM-4228 deck
target:       a lawful fixed-pair or canceled-family depth-eight repair
map:          A_(u,v) -> (beta,centered primitive H,omega), then pull back
              by the common dilation and intersect each base component
preserved:    literal pair ratio, exact density, endpoint phase, base mass,
              component count, common dilation, and target threshold
destroyed:    endpoint addresses after reduction to (beta,omega)
sidecar:      exact primitive extremizers or an analytic fast primitive
decisive test: beta>=66/91 and omega/g<=1650/28710227.   (61)
```

The proof and failure boundaries are:

1. The exact computation certifies `(3713,5149)` and its dilations; it makes
   no claim that this is the first, smallest, or unique coprime pair passing
   the gate.
2. Equation `(27)` is an obstruction to one proof mechanism, not a claim
   that adjacent pairs are unsafe. Common dilation still reduces their
   discrepancy, and THM-4231 handles their sufficiently large coordinates.
3. The canceled family proves an `O(1/k)` mechanism and the conservative
   analytic cutoff `180785`. Exploratory data suggesting a permanent exact
   residue-class transition at `k=748` remain **OPEN** and are not used.
4. The deck `E_3467` is sufficient, not necessary. Missing `(14)` says only
   that this transfer does not fire.
5. Arbitrary finite pair entry, replacement of the fixed pool, and full
   LRC(14) remain open.

**QED.**
