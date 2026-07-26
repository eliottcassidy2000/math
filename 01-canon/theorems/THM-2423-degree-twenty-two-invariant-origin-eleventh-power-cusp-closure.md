---
id: THM-2423
title: "Degree-twenty-two invariant-origin eleventh-power cusp closure"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  In the open first-flux chart of the genuine nonsplit polynomial
  exact-square-prefix degree-twenty-two branch, the invariant coefficient
  origin B=C=D=E=W=0 is empty. The first two fluxes make
  v=u/y^2 a root of one fixed quintic and reconstruct the nonzero constant
  zeta=Z/y^3. Hence y=h^2 and T=t_0 h^3. The third flux is a nonzero
  multiple of h^11, so the Keller one-form leaves only an eleventh-power
  monomial cusp. The full polynomial Faber sidecar has a unique nonzero
  h^11 pole there and rules it out. This removes one weighted stratum of
  the A!=0 chart; it does not close degree twenty two, JC(2), or DC(2).
source: klein-2026-07-26-degree-twenty-two-origin-cusp
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
related:
  - THM-2247-nonsplit-terminal-quartic-degree-fourteen-closure
  - THM-2262-degree-eighteen-trigonal-spectral-discriminant-reduction
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2406-degree-eighteen-H4-weighted-pole-deep-wall-collapse
script: 04-computation/jc2_degree22_invariant_origin_cusp_thm2423.py
output: 05-knowledge/results/jc2_degree22_invariant_origin_cusp_thm2423.out
script_sha256: 16b1223a60e8d4e7e418b11c9dd7f6e34a777cc9c58ea3c31d8b67c22e808ff7
output_sha256: 418cb934de555014d6c5f94356bc23a60bded32e66cd154527bcec6ced100d77
hash_basis: working-tree bytes (LF)
---

# THM-2423 -- the degree-twenty-two invariant origin is empty

**PROVED CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

THM-2411 empties the divisor on which the first Faber flux loses its
`Z=T^2` coefficient. This theorem begins the complementary chart. Its exact
conclusion is

```text
genuine degree-22 trajectory,
mathcal A!=0,
B=C=D=E=W=0
    => contradiction.                                           (1)
```

The mechanism is not the degree-eighteen `H_4` elimination itself. What
transfers from THM-2297 and THM-2406 is the operation:

```text
constant-field flux quotient
  -> exact weighted cusp
  -> Keller rational primitive
  -> whole-polynomial sidecar pole.                              (2)
```

The fixed coefficient and the whole-polynomial sidecar are both
load-bearing. Neither the projected spectral curve nor its weighted scalar
quotient is a sufficient carrier.

## 1. Inherited normalized flux chart

Use the target-translated degree-twenty-two coordinates of THM-2411:

```text
y=11s,                  u=dT,                 Z=T^2,

wt(y,u,Z,B,C,D,E,W)=(1,2,3,2,3,4,5,6).             (3)
```

The first two flux equations are

```text
N_1=1331 mathcal A Z+4 mathcal K=0,
N_2=0,                                                     (4)
```

with the exact polynomials in THM-2411, equations (12)--(16). Work on

```text
mathcal A!=0,             B=C=D=E=W=0.                     (5)
```

If `y=0`, then `mathcal A=-1089u!=0`, while (4) gives
`uZ=0`. Thus `Z=0`, contradicting `T!=0`. Hence

```text
y!=0.                                                        (6)
```

Put

```text
v=u/y^2,                 zeta=Z/y^3.                         (7)
```

Then

```text
mathcal A=9y^2(7-121v),
```

so the open chart also gives `121v-7!=0`. Dividing (4) by their nonzero
weighted powers of `y` gives

```text
f_1
 =-1229844v^2+483153v zeta+33880v-27951zeta-84=0,            (8)

f_2
 =-396829664v^3+49193760v^2-54113136v zeta-406560v
   +5314683zeta^2+745360zeta+224=0.                          (9)
```

Their exact resultant is

```text
Res_zeta(f_1,f_2)=-595244496 L_5(v),                         (10)
```

where

```text
L_5
 =155624547606v^5+3215383215v^4-1700698560v^3
   +58124770v^2-855470v+2583.                                (11)
```

Therefore `v` is algebraic over the algebraically closed constant field
`C`. Since `v in C(x)`, it follows that

```text
v in C.                                                       (12)
```

Equation (8) reconstructs

```text
zeta
 =28(121v-3)(363v-1)/[3993(121v-7)].                         (13)
```

There is no excluded root hidden here:

```text
L_5(7/121)=-44800,
L_5(3/121)=-6144,
L_5(1/363)=51200/81.                                         (14)
```

Thus `zeta in C*`.

## 2. The weighted curve is one eleventh-power cusp

Choose `t_0 in C*` with `t_0^2=zeta`. From

```text
T^2=zeta y^3
```

and coprimality of `2` and `3`, comparison of the divisors in `C(x)`
gives one `h in C(x)^*` such that, after absorbing a constant,

```text
y=h^2,                    T=t_0h^3,
u=vh^4,                   d=(v/t_0)h.                       (15)
```

This is a lossless parametrization, not a root choice in a larger
function field.

Let `F=R_Q/q` be the third Faber flux. A direct degree-twenty-two Laurent
calculation gives

```text
F=3y^7 H(v,zeta)/(14992384T),                              (16)
```

where

```text
H(v,zeta)
 =-43923zeta^2-1449459zeta v^2+139755zeta v-770zeta
   +1229844v^3-33880v^2+84v.                               (17)
```

On (13),

```text
H(v,zeta)
 =-56(121v-3)(363v-1)
    (5314683v^3-307461v^2+22869v-203)
   /[363(121v-7)^2].                                      (18)
```

The gcd of `L_5` and the numerator in (18) is one. Hence (14) and (18)
show that the coefficient is nonzero at every possible `v`. Equations
(15)--(16) become

```text
F=f_0h^11,                     f_0 in C*.                    (19)
```

## 3. The Keller one-form leaves only a monomial cusp

Write `A_src` for the polynomial in the exact-square-prefix
factorization, to distinguish it from the flux coefficient `mathcal A`.
The inherited Keller one-form is

```text
A_src(2T F'+F T')=2kappa T,       kappa in C*.               (20)
```

Substituting (15) and (19), then cancelling the nonzero `T`, gives

```text
A_src(h^11)' in C*.                                      (21)
```

The rational-primitive lemma of THM-2214 now has only two cases.

If `A_src` is constant, then `h^11` is a nonconstant affine polynomial.
That is impossible because its finite zero is simple.

Otherwise, after translating `X=x-xi`,

```text
A_src=a_0X^m,              m>=2,

h^11=c_0+c_1X^(1-m),       c_1!=0.                        (22)
```

If `c_0!=0`, the nonzero finite roots of
`c_0X^(m-1)+c_1` are simple, again impossible for an eleventh power.
Consequently, for some integer `ell>=1`,

```text
c_0=0,               m=11ell+1,
h=h_0X^-ell.                                               (23)
```

Thus the flux quotient leaves exactly one eleventh-power monomial cusp.

## 4. The full polynomial sidecar kills the cusp

At the invariant origin, the normalized mate is the single Faber seed
`E_22`. Write

```text
P=H_0^2+mathcal L
```

for the polynomial approximate root and its polynomial remainder. The full
degree-twenty-two sidecar of THM-2411 is now

```text
E_22-mathcal R_6(P,H_0)
 =33T/2048 [
      14mathcal L^4-28mathcal L^3s-14mathcal L^2Td
      +42mathcal L^2s^2-mathcal L T^2
      +56mathcal L Tds-56mathcal Ls^3
      +14T^2d^2+6T^2s-140Tds^2+70s^4
   ] in C[x,z].                                            (24)
```

The part independent of `mathcal L` is

```text
33T y^4 J(v,zeta)/2048,                                  (25)
```

where

```text
J(v,zeta)
 =14v^2+(6/11)zeta-(140/121)v+70/14641.
```

On (13),

```text
J(v,zeta)
 =14(1771561v^3-73205v^2+4235v-23)
    /[14641(121v-7)].                                     (26)
```

The gcd of `L_5` and the numerator in (26) is one. Thus (25) has a
nonzero coefficient and, on (23), exact valuation

```text
-11ell.                                                    (27)
```

Because `mathcal L` is polynomial, the terms of (24) containing,
respectively, one, two, three, or four copies of `mathcal L` have valuation
at least

```text
-9ell,             -7ell,             -5ell,             -3ell.
                                                                    (28)
```

Hence the pole (27) is unique and cannot cancel. This contradicts the
polynomiality in (24). The monomial cusp is empty.

If `y` had been constant instead, (12)--(13) would make `u,Z,T`, and then
`q`, constant. The genuine deck fixes the constant field but sends
`q` to `-q`, contradicting `q!=0`. This covers both cases and proves (1).

## 5. Scope, connection ledger, and next obstruction

The theorem removes exactly the invariant weighted origin from the
complementary `mathcal A!=0` degree-twenty-two chart. It does not remove a
one-sparse axis, classify the remaining weighted coefficient cone, or treat
split/even and integral order-raising branches. It proves neither `JC(2)`
nor `DC(2)`.

The lawful cross-thread transfer is methodological:

```text
source:
  continuation profiles and affine-origin sidecars;

map:
  replace a scalar quotient by its fixed constant-field coefficient
  and whole-polynomial response;

preserved predicate:
  constancy under the deck and absence of finite poles;

destroyed information in the scalar quotient:
  the polynomial approximate-root remainder;

restoring sidecar:
  equation (24);

hostile boundary:
  the flux equations really do leave the cusp (23), so the sidecar is
  necessary rather than decorative.
```

No knot invariant, tournament orientation, additive graph, or Prony count is
identified with a Keller trajectory. Their reusable lesson is only that a
quotient must retain the coordinate consumed by the next operation.

The next degree-twenty-two target is consequently precise: classify the
nonzero weighted parameter strata in `mathcal A!=0`, starting with the five
one-sparse axes, and retain the third flux and full-mate sidecar before
projecting to any spectral square class.

## 6. Exact verification

Run

```bash
python3 04-computation/jc2_degree22_invariant_origin_cusp_thm2423.py
python3 -O 04-computation/jc2_degree22_invariant_origin_cusp_thm2423.py
```

The companion independently reconstructs the three degree-twenty-two
Laurent observables by recurrence and finite multinomial expansion. It
verifies (8)--(18), the nonzero wall and zero-factor controls (14), both gcd
certificates, and the complete sidecar (24)--(26). All truth-bearing checks
use explicit exceptions and remain active under optimized Python.

The computation verifies the exact algebra. The constant-field,
coprime-divisor, rational-primitive, deck, and pole-valuation arguments are
the proof above rather than delegated computer conclusions.
