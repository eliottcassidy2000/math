---
id: THM-4294
title: "Lambda-zero central Keller extinction"
status: >
  PROVED RELATIVE TO THM-4103/4230/4272/4292 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On W=Lambda=0, U!=0, the pullback of the good-target
  invariant differential has exact vertical order nine at the generic point
  of the genus-seven central component, so that component is Keller-constant.
  THM-4292 makes every positive-genus exceptional tail constant, while every
  remaining component is rational. Proper-flat degree conservation therefore
  contradicts every positive generic response degree. The W-and-Lambda-zero
  slice is excluded; exact-M=12 seam entry and JC(2) remain open.
source: root/jc2-signal-20260830
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
  - THM-4292-lambda-zero-repeated-face-keller-extinction
related:
  - THM-4289-a23-blowdown-observer-kahler-dualizing-quotient
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
  - THM-4291-lambda-zero-genus-five-tail-degree-forty-two-equivariant-hostile
  - THM-4293-lambda-zero-normalized-response-degree-stratification
primary_script: 04-computation/jc23_lambda_zero_central_keller_extinction_thm4294.py
primary_output: 05-knowledge/results/jc23_lambda_zero_central_keller_extinction_thm4294.out
primary_script_sha256: 0b85b5df9ea6297df84844bc69b921299ede7806a04bfaeca3b3a2bc4e72634a
primary_output_sha256: 1e22c4bd9c697ff31dec7a408585d92780b3a678343c226ae31682691190993f
independent_audit_script: 04-computation/jc23_lambda_zero_central_keller_extinction_independent_audit_thm4294.py
independent_audit_output: 05-knowledge/results/jc23_lambda_zero_central_keller_extinction_independent_audit_thm4294.out
independent_audit_script_sha256: 2d98db2d8c1691f2b70a3e397d9972356d42b9a7a29c6022e31de24f531dd25c
independent_audit_output_sha256: 555a79f5e737d2e66daa5460fb55daa4219aa98792b62d7b0eb12dc87a9b5f17
sidecar_script: 04-computation/jc23_lambda_zero_repeated_discriminant_ladder_thm4294.py
sidecar_output: 05-knowledge/results/jc23_lambda_zero_repeated_discriminant_ladder_thm4294.out
sidecar_script_sha256: 4a31565fae680a18c238ea34ca52fc3b942d093e76d741d9005f91838b5c43bb
sidecar_output_sha256: d4db9fb3609e16660633bc4b28c3655922939e2088e2fed7e93c6d79b6139e16
hash_basis: raw LF bytes
audit: >
  PASS. The primary SymPy path reconstructs the scaled exponent ledger,
  literal wall component, restriction of G_P, generic-unit tests, and the
  visible r=3 hostile. A dependency-free sparse-polynomial path independently
  reconstructs the special-fibre and G_P identities and audits the remaining
  finite ledgers with geometric inputs explicitly typed as imports. Normal,
  optimized, and fixed-hash-seed streams byte-match the frozen outputs. A
  non-load-bearing SymPy sidecar reconstructs the repeated-discriminant ladder
  by discriminant and moving-critical-value paths, then checks the imported
  residue/genus ledgers; its three streams also byte-match.
---

# THM-4294 -- Lambda-zero central Keller extinction

**PROVED RELATIVE TO THM-4103/4230/4272/4292 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. THE `W=Lambda=0`, `U!=0` SLICE IS EXCLUDED.
EXACT-`M=12` SEAM ENTRY AND `JC(2)` REMAIN OPEN.**

## 1. Statement

Work over `C` in the exact-weight-twelve reduced `(2,3)` seam and assume

```text
W=0,                   Lambda=U+Z=0,                   U!=0.       (1)
```

Let `sigma^12=Q=q_target^(-1)` and let `E_sigma` be THM-4230's good elliptic
target after the weight-twelve base change.

> **Theorem.** No positive-degree Keller response exists under `(1)`.
> Equivalently, the `W=Lambda=0`, `U!=0` codimension-two slice inside the
> exact-`M=12` reduced seam is empty of nonautomorphic Keller candidates.

The inheritance pass is:

- closest proved mechanism: THM-4291's positive vertical order on one
  genus-five tail;
- canonical hostile: THM-4293's `r=3` degree-`36/28` eigenmaps, which pass
  every degree and Eisenstein-norm test;
- corrected near miss: stratifying the normalized response degree is valid
  but is not the final obstruction;
- least-used decisive sidecar: the divisorial order of the **good** invariant
  differential at the generic point of the central component.

The live concept board was `{residue scaling, central component, repeated
tails, degree conservation, deck eigenline}`. The new connection is

```text
residue scaling -> central vertical order -> central constancy
                -> all-component extinction -> degree contradiction.       (2)
```

Its source is THM-4103's residue identity, its target is the actual pullback
line bundle on a regular graph model, its preserved predicate is generic map
degree, and the information it deliberately discards is the detailed pole
packet. THM-4292 supplies the needed tail sidecar. The cheapest decisive test
is whether `G_P` is a unit at the central generic point.

## 2. Exact scaling of the good invariant differential

THM-4230 uses

```text
Q=sigma^12,          s=sigma^(-1)S,          p=sigma^(-2)P,
F_Q=sigma^(-2)G.                                               (3)
```

At fixed `sigma`, `P=sigma^2p`; hence

```text
(F_Q)_p=G_P,                  ds=sigma^(-1)dS.          (4)
```

THM-4103's exact Keller residue identity is

```text
phi^*(dA/(2C_target))=Q ds/(F_Q)_p.                    (5)
```

For the integral target coordinates

```text
A=sigma^(-4)X,             C_target=sigma^(-6)Y,
eta_0=dX/(2Y),                                               (6)
```

one has `dA/(2C_target)=sigma^2 eta_0`. Substitution into `(3)--(6)` gives
the exact identity

```text
phi^*eta_0=sigma^9 dS/G_P.                               (7)
```

No continuity from another coefficient gate and no choice of response packet
enters `(7)`. After a ramified base change `sigma=unit*pi^e`, its vertical
factor has order `9e`.

## 3. Exact order nine on the central component

The special source is

```text
G_0=R C,
R=S^2-P,
C=1-U P^6-W S^2P^5-Z S^4P^4.                            (8)
```

On `(1)`, `Z=-U` and

```text
C=1-U P^6+U S^4P^4,
C_P=-2U P^3(3P^2-2S^4).                                 (9)
```

Differentiating `(8)` and reducing at the generic point `xi_C` gives

```text
G_P=-C+R C_P = R C_P              in k(C).             (10)
```

Every factor in `(10)` is nonzero in `k(C)`:

- `P` cannot vanish identically because `C|_(P=0)=1`;
- `R` cannot vanish identically because `C|_(P=S^2)=1` on `Lambda=0`;
- if `3P^2-2S^4` held identically, reducing `C` by that relation would give
  `1+(4U/27)S^12`, not the zero polynomial.

Equivalently, the exact certificate checks

```text
gcd_P(C,C_P)=gcd_P(C,S^2-P)=1,
Res_P(C,3P^2-2S^4)=(4US^12+27)^2.                      (11)
```

Thus `G_P` is a unit in the DVR at `xi_C`. The implicit-function chart makes
`dS` a nonzero relative differential generator there. Equation `(7)` now
proves

```text
ord_C(phi^*eta_0)=9                                      (12)
```

exactly.

## 4. Why order nine forces the central map to be constant

The good target has equation of the form

```text
Y^2=X^3+1-a_8 sigma^8X-a_12 sigma^12,                  (13)
```

so its special fibre is the smooth elliptic curve

```text
E_0:Y^2=X^3+1                                           (14)
```

and `eta_0` is a nowhere-zero relative invariant differential. Properness of
the target extends the generic Keller map over the DVR at `xi_C`; no raw
codimension-two graph descent is needed for this generic-point statement.
Reducing `(12)` gives

```text
phi_C^*eta_(0,special)=0.                               (15)
```

A nonconstant morphism between smooth characteristic-zero curves induces a
separable function-field extension and has nonzero differential. Therefore
`(15)` proves that `phi_C:C->E_0` is constant. Normalizing or resolving the
graph does not change the valuation at `xi_C`, so the strict-transform map
is constant as well.

This is a genuine good-model calculation. The vanishing factor in `(7)` is
not the original singular-target differential discarded in THM-4289; it is
the pullback of the nowhere-zero differential on the good elliptic scheme.

## 5. Every special-fibre component is constant

THM-4272/4292 give the complete component inventory under `(1)`:

1. the raw toric special fibre is `R union C`, both of multiplicity one;
2. the family is smooth away from the unique `A_23` contact;
3. every new component from normalization, regularization, or semistable
   reduction is centered above that contact;
4. fixed toric-subdivision and ordinary graph-blowup components are rational;
5. THM-4292 proves every noncentral divisorial component above the contact,
   including all positive-genus repeated Newton faces, Keller-constant.

The rational component `R` and every rational exceptional component map
constantly to an elliptic curve. Section 4 handles `C`; THM-4292 handles all
remaining positive-genus components. Hence **every** component of a proper
regular graph model has constant special map.

This inventory is the needed sidecar. Central order nine alone would not
exclude degree migration to a positive-genus exceptional tail; THM-4292 is
what removes that loophole.

## 6. Proper-flat degree contradiction

Let `L=O_E(O)` be the relative degree-one origin bundle on the good elliptic
scheme, and let `f:X->E` be the resolved morphism after a finite base change.
Write the special fibre as

```text
X_0=sum_i m_i X_i.                                      (16)
```

For the actual pullback line bundle `M=f^*L`, proper-flat intersection gives

```text
deg(M_generic)=sum_i m_i deg(M|X_i).                    (17)
```

Every restriction `f|X_i` is constant by Section 5, so every term on the
right of `(17)` is zero. Thus the generic map has degree zero. But the Keller
response is nonconstant and therefore has positive degree. This
contradiction proves the theorem.

There is no vertical-twist loophole: `(17)` uses the actual line bundle
`f^*L`, not an arbitrarily chosen extension. Fibre multiplicities are retained
explicitly. Resolving codimension-two base points preserves generic degree
and introduces only components already covered in Section 5.

## 7. The `r=3` hostile survives arithmetic but not the differential

THM-4293's unique noncritical arithmetic survivor has degree

```text
36=4N(3),                     28=4N(3+omega).           (18)
```

The exact deck eigenline writes its central candidate as `a u`. In
THM-4272's normalized wall coordinates

```text
x=alpha P,                y=beta SP,                b=1/S,
u=(-x^2,y^2).                                                (19)
```

For the target uniformizer `t_E=-X/Y` at the origin,

```text
u^*t_E=x^2/y^2=(alpha^2/beta^2)b^2.                    (20)
```

Thus every nonzero multiple `a u` has local index two and pullback-
differential order one at the fixed contact. Multiplication on `E_0` is
etale in characteristic zero. This is the positive control: the norm-`9/7`
maps are genuine nonconstant eigenmaps, but they cannot equal a Keller
specialization whose good invariant differential reduces to zero on `C`.

Consequently THM-4293's noncritical degree/genus stratification remains a
valid corrected diagnostic, while its `r=3` and repeated-discriminant OPEN
statuses are superseded on this slice by this theorem. The repeated-locus
square classes and pole packets are no longer load-bearing for its exclusion.

## 8. Niche signal: the repeated-discriminant ladder

An orthogonal exact probe completed the repeated locus even though Section 6
makes it unnecessary for exclusion. On

```text
c_1=...=c_5=0,                  c_6^2+4U=0,             (21)
```

shift the double root by `u=z^6(-1/c+v)`, where `c=c_6`. The local
Weierstrass discriminant has exactly the successive possible orders

```text
m in {13,14,15,16,18};          m=17 is absent.         (22)
```

The first splitter coefficients are proportional to

```text
alpha_11,
8c+3upsilon_5,
eta,
405c-5152,
0,
Q(541Q+322),                                             (23)
```

respectively. The maximal stratum is THM-4292's deepest witness

```text
c=5152/405,                    Delta=4672/135.           (24)
```

Its residue square class is `Q(541Q+322)`, nonsquare in `C(Q)`. Before the
central-order obstruction is applied, normalization and the Eisenstein
eigenline sieve leave only

```text
m=14: d=28,                  m=18: d=16 or 12.          (25)
```

The `m=14` survivor has the same norm-seven degree `28` as the noncritical
`r=3` finite response; exactly eight pole units move between the collision
and quartic packets. This is useful structural signal, not a dependency of
the slice exclusion. At `m=18` the quadratic residue field already supports
nonzero target points, so a bare rationality argument would not determine its
origin label. Section 6 bypasses that label completely.

The sidecar reconstructs `(22)--(24)` from the literal response polynomial,
then verifies `(25)` and canonical-divisor saturation against the inherited
residue/genus formulas. It records that distinction in its local scope.

## 9. Scope and reproduction

This theorem closes exactly

```text
exact-M=12 reduced seam,        W=Lambda=0,        U!=0.       (26)
```

It does not prove entry into the exact-`M=12` seam, cross `U=0`, `Z=0`, or
`D=0`, or prove `JC(2)` or `DC(2)`. The same exponent-nine formula suggests a
wider central-extinction test on other coefficient gates, but their complete
component inventories must be audited separately before using `(17)`.

Run

```bash
python3 -B 04-computation/jc23_lambda_zero_central_keller_extinction_thm4294.py
python3 -B -O 04-computation/jc23_lambda_zero_central_keller_extinction_thm4294.py
python3 -B 04-computation/jc23_lambda_zero_central_keller_extinction_independent_audit_thm4294.py
python3 -B -O 04-computation/jc23_lambda_zero_central_keller_extinction_independent_audit_thm4294.py
python3 -B 04-computation/jc23_lambda_zero_repeated_discriminant_ladder_thm4294.py
```

The first path uses SymPy; the second is a dependency-free sparse-polynomial
audit. **QED.**
