---
id: THM-2823
title: "Degree-twenty-six triple-pole Faber valuation obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In a nonsplit
  polynomial exact-square-prefix degree-twenty-six chart, any response
  point with one finite local pair
  (ord V,ord M)=(3,1) is impossible.  The exact top-row identity
  K=-(d/2)(Phi/q)+T H_26 reduces the first-flux locus to K=T H_26.
  A complete valuation split gives regular valuation mismatches, two
  coprime cubic initial forms, or one exceptional two-scale unit ideal.
  Both THM-2817 sextic carriers have such a triple pole, so neither enters
  this first open degree.  Other degrees/charts, JC(2), and DC(2) remain
  open.
source: root/sextic-degree26-triple-pole-obstruction-2026-07-28
audit: >
  thm2823-hostile-audit-2026-07-28 independently derived the Faber
  recurrence and translation normalization, checked the K decomposition,
  audited every valuation lane for cancellation and exhaustiveness,
  reproduced both resultants and the exceptional unit ideal, verified both
  carrier multiplicities, and replayed normal, optimized, and stored
  transcripts with matching LF hashes: ACCEPT.
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2784-nonsplit-response-square-potential-divisor-and-infinity-classification
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2817-sextic-e3-maximal-pole-power-chebyshev-accessory-classification
related:
  - THM-2760-exact-prefix-even-faber-flux-gcd-and-smooth-boundary-exclusion
  - THM-2822-sextic-response-centered-lift-mod-three-faber-obstruction
script: 04-computation/jc_degree26_triple_pole_faber_valuation_obstruction_thm2823.py
output: 05-knowledge/results/jc_degree26_triple_pole_faber_valuation_obstruction_thm2823.out
script_sha256: 338e106b2326be8ab59e784aca197c85c05f89df88c62db20eb106070b51fbcc
output_sha256: 9fdd094358adae1c61808d78d3d26ef6f6c554e33377425aef745df34aeaa812
independent_script: 04-computation/jc_degree26_triple_pole_faber_valuation_independent_audit_thm2823.py
independent_output: 05-knowledge/results/jc_degree26_triple_pole_faber_valuation_independent_audit_thm2823.out
independent_script_sha256: bde1e5d90db29c8e19b12832400fc9f8dde43b742ab21e76385d0a2e31b212c3
independent_output_sha256: 60b92d4ee89b4663b55b4315bb398f6ec056c30935312856251c0cb175faad84
hash_basis: LF-normalized bytes
---

# THM-2823 -- one triple pole excludes degree twenty six

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2822 excludes the canonical centered lifts of THM-2817's two sextic
response carriers in every Faber degree.  The first open degree admits a
stronger result: no choice of the remaining exact-square source data can
lift either carrier there.  The obstruction is local at one multiplicity-
three zero of the response potential.

## 1. The normalized degree-twenty-six bank

In the polynomial exact-square-prefix chart write

```text
H=Vz^2+B_src z+C_0,                L=A_src z+E_0,
P=H^2+L,                           U^2=V,              (1)
```

and use the standard quartic invariants

```text
d=C_0-B_src^2/(4V),
s=A_src B_src/(2V)-E_0,
T=A_src^2/V=q^2.                                       (2)
```

Start with a reduced degree-twenty-six full Faber representative.  Under a
legal constant translation `P_c=P+c`, the exact binomial law gives

```text
E_26(P)=E_26(P_c)-(13/2)c E_22(P_c)+lower rows.        (3)
```

Thus `c=2a_22/13` kills the original `E_22` coefficient.  This translation
preserves the Keller bracket and only changes the named invariant `s` by a
constant.  In the normalized coordinate,

```text
Q=E_26+c_18E_18+c_14E_14+c_10E_10+c_6E_6+c_2E_2.     (4)
```

The nonsplit transfer of THM-2784 and the Hamiltonian flux identities give

```text
Phi_Q=0,                         Psi_Q in C,
K=R_Q/q,                        A_src K=lambda M,      (5)
```

where `lambda` is a nonzero constant and `M` is the response carrier of
THM-2796.

## 2. Exact top-row decomposition

For a row of reduced degree `m`, let

```text
(1+2dt^2+qt^3+(d^2-s)t^4)^(m/4)=sum_(n>=0)c_n t^n,

Phi_m=4c_(m+1),                 Psi_m=4c_(m+2),
R_m=4c_(m+3)+2d c_(m+1).                            (6)
```

Extracting `(6)` for the six rows in `(4)`, dividing the odd expressions
by `q`, and replacing `q^2` by `T` gives the exact identity

```text
K=-(d/2)(Phi_Q/q)+T H_26,                             (7)
```

where

```text
4096 H_26 =
    72c_18T^2 -4032c_18Tds +6720c_18s^3
  +896c_14Td -4480c_14s^2
  +2560c_10s -1024c_6
  -143T^3d +5148T^2d^2s +1287T^2s^2
  -24024Tds^3 +12012s^5.                             (8)
```

Consequently the first flux equation makes the response quotient

```text
K=T H_26.                                             (9)
```

The companion derives `(7)--(8)` twice: from the Faber differential
recurrence and independently by direct generalized multinomial extraction.

## 3. Local valuation lemma

Let `beta` be a finite point and write `v=ord_beta`.  Assume

```text
v(V)=3,                          v(M)=1.               (10)
```

Put

```text
a=v(A_src),                     b=v(B_src),            (11)
```

where `b=infinity` is allowed when `B_src=0`.  The response identity in
`(5)` with `M!=0` makes `A_src` nonzero, so `a` is finite.  Since all four
source coefficients in `(1)` are polynomials, `a,b` are nonnegative and
`C_0,E_0` are regular at `beta`.  Equations `(2)` and `(5)` give

```text
v(T)=2a-3,                      v(K)=1-a,
v(d)>=min(0,2b-3),             v(s)>=min(0,a+b-3).    (12)
```

Whenever the minimum on the right of `(12)` is negative, equality holds:
a regular term cannot cancel a pole.

We prove that `(5)` and `(9)` are incompatible for every pair `(a,b)`.

## 4. Three regular valuation regions

### 4.1 `a>=3`

Here

```text
v(T)>=3,                        v(d)>=-3,
v(s)>=0.                                               (13)
```

Every monomial in `(8)` has nonnegative valuation, so
`v(H_26)>=0`.  Equation `(9)` gives `v(K)>=3`, whereas `(12)` gives
`v(K)=1-a<=-2`.

### 4.2 `a=2,b>=1`

Now

```text
v(T)=1,                         v(d)>=-1,
v(s)>=0.                                               (14)
```

Again `(8)` gives `v(H_26)>=0`; hence `(9)` gives `v(K)>=1`, while
`(12)` gives `v(K)=-1`.

### 4.3 `a<=1,b>=2`

In this region `v(d)>=0`.  The only subcase with negative `s` is
`(a,b)=(0,2)`, where `v(s)=-1`.  Directly in the exact polynomial
`Phi_Q/q`, the unique least-valuation monomial is

```text
(143/16384)T^4.                                      (15)
```

For `a=0` it has valuation `-12`; for `a=1` it has valuation `-4`.
All other monomials have strictly larger valuation, including in the
`(0,2)` boundary.  Thus `Phi_Q=0` is impossible.

These arguments include `b=infinity`.

## 5. The generic polar lanes

The cases

```text
(a,b)=(0,0),(1,0),(2,0),(1,1)                       (16)
```

have

```text
v(Td)=2v(s)<0.                                       (17)
```

Let

```text
r=lc(Td)/lc(s^2).                                    (18)
```

The least terms of the two fluxes are nonzero constants times

```text
s^6 f(r),                         s^7 g(r),            (19)

f(r)=r^3-21r^2+35r-7,
g(r)=7r^3-35r^2+21r-1.
```

The first expression must vanish because `Phi_Q=0`; the second must vanish
because `Psi_Q` is constant while its displayed valuation is negative.
But

```text
Res_r(f,g)=-2^21.                                    (20)
```

Hence no lane in `(16)` survives.

There is also a sharper source-typed shortcut.  In every lane in `(16)`,
the polar summands in `(2)` dominate, so if `A_src`, `B_src`, and `V`
have leading coefficients `a_0`, `b_0`, and `v_0`, then

```text
lc(T)=a_0^2/v_0,
lc(d)=-b_0^2/(4v_0),
lc(s)=a_0b_0/(2v_0).
```

Consequently `r=-1` identically, and already `f(-1)=-64`.  The resultant
`(20)` is retained as an independent certificate that remains valid even
after forgetting this leading-coefficient sidecar.

## 6. The sole two-scale boundary `(a,b)=(0,1)`

Here

```text
v(T)=-3,                         v(d)=-1,
v(s)=-2.                                               (21)
```

For the eliminant calculation introduce two leading ratios:

```text
r=lc(Td)/lc(s^2),                zeta=lc(T^2)/lc(s^3). (22)
```

After removing the nonzero common monomials and scalar denominators, the
least coefficients of `Phi_Q`, `Psi_Q`, and `H_26` are

```text
P =
 143zeta^2-13728r^3-20592zeta r+288288r^2
 +48048zeta-480480r+96096,

G =
 1287zeta^2+5148zeta r^2-96096r^3-72072zeta r
 +480480r^2+60060zeta-288288r+13728,

J =
 143[(9-r)zeta+36r^2-168r+84].                       (23)
```

The first two vanish by the same flux argument as above.  Moreover
`v(K)=1` and `v(T)=-3`, so `(9)` requires `v(H_26)=4`; in particular its
least coefficient `J` must vanish.

The same source calculation as above actually forces `r=-1`.  Then

```text
J(-1,zeta)=286(5zeta+144).
```

Thus `J=0` forces `zeta=-144/5`, at which

```text
P=-24490752/25,                  G=-50189568/25.
```

This already empties the boundary.  The following two-variable
elimination is a stronger independent check which deliberately forgets
the forced ratio.  Exact elimination gives

```text
<P,G,J>=<1> in Q[r,zeta].                             (24)
```

For a transparent resultant certificate, `r=9` makes `J` nonzero.
Otherwise `J=0` gives

```text
zeta=12(3r^2-14r+7)/(r-9).                           (25)
```

Substitution in `P,G` leaves, up to nonzero factors, the two quintics

```text
p=2r^5+3r^4-488r^3+2842r^2-6930r+4011,
g_5=13r^5-182r^4+1974r^3-8776r^2+13173r-5130,

Res_r(p,g_5)=-2^36 3^2 31^9 37.                     (26)
```

Thus the exceptional lane is empty.  Sections 4--6 exhaust all
nonnegative `a,b`.

## 7. Application to both sextic carriers

For THM-2817's power carrier,

```text
V=(16/9)x^5(x^3-1)^3,
M=x(x^3-1)(x^3-1/2).                                 (27)
```

At every cube root of one, `(27)` has

```text
(v(V),v(M))=(3,1).                                   (28)
```

For its centered Chebyshev carrier,

```text
V=(256/9)(y^2-1/4)^4(y^2-1)^3,
M=(T_3(y)/4)(y^2-1/4)(y^2-1).                        (29)
```

At `y=1` and `y=-1`, `(29)` has the same pair `(28)`.
The local lemma therefore excludes both carriers from the normalized
degree-twenty-six nonsplit polynomial exact-square-prefix chart.

Unlike THM-2822, this conclusion permits arbitrary polynomial
`A_src,B_src,C_0,E_0`; centeredness is not assumed.

## 8. Scope

This is a local chart-entry obstruction at one reduced degree.  It does
not invalidate THM-2817's abstract rational response maps, exclude their
entry at another Faber degree or outside the polynomial exact-square-prefix
chart, classify all degree-twenty-six response carriers, or prove `JC(2)`
or `DC(2)`.

The mechanism is reusable: any future response carrier with one finite
local pair `(ord V,ord M)=(3,1)` fails the same degree-twenty-six chart.
Other multiplicity pairs require a new Newton polygon.

## 9. Exact companion

The primary companion:

1. reconstructs all six rows in `(4)` by recurrence and independent
   generalized multinomial expansion;
2. verifies the translation coefficient in `(3)` and identity `(7)`;
3. computes every weighted initial form in Sections 4--6;
4. verifies `(20)`, the unit ideal `(24)`, and resultant `(26)` exactly;
5. checks the unique monomial `(15)` and both regular `H_26` bounds; and
6. verifies the exact multiplicities in `(27)--(29)`.

It contains no Python `assert` node.  Run

```text
python 04-computation/jc_degree26_triple_pole_faber_valuation_obstruction_thm2823.py
python -O 04-computation/jc_degree26_triple_pole_faber_valuation_obstruction_thm2823.py
```

The hostile audit uses the logarithmic-derivative Faber recurrence rather
than either primary extraction path, reconstructs every weighted initial
form, verifies the exhaustive valuation partition and both resultants,
checks `<P,G,J>=<1>`, and proves the forced `r=-1` shortcut directly.  Run

```text
python 04-computation/jc_degree26_triple_pole_faber_valuation_independent_audit_thm2823.py
python -O 04-computation/jc_degree26_triple_pole_faber_valuation_independent_audit_thm2823.py
```

Normal, optimized, and stored transcripts agree exactly for both
companions.
