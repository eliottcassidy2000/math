---
id: THM-4293
title: "Lambda-zero normalized response-degree stratification"
status: >
  PROVED RELATIVE TO THM-4103/4120/4230/4290/4292 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. On W=Lambda=0, U!=0, the double top-edge root lowers
  the two collision indices from 11 to 11-r. The normalized response degrees
  are 42-2r or 34-2r, not 42 or 34, and the normalized genus is 18-r. Uniform
  tail extinction and the exact deck eigenlattice exclude every noncritical
  stratum r=1,...,6 except r=3. THM-4294 subsequently excludes both r=3 and
  the repeated-discriminant locus on this `W=Lambda=0` slice by
  central-component extinction. The general coefficient wall, exact-M=12
  entry, and JC(2) are not proved.
source: root/jc2-signal-20260830
depends_on:
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion
  - THM-4292-lambda-zero-repeated-face-keller-extinction
related:
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
  - THM-4291-lambda-zero-genus-five-tail-degree-forty-two-equivariant-hostile
  - THM-4294-lambda-zero-central-keller-extinction
primary_script: 04-computation/jc23_lambda_zero_normalized_response_degree_stratification_thm4293.py
primary_output: 05-knowledge/results/jc23_lambda_zero_normalized_response_degree_stratification_thm4293.out
primary_script_sha256: 810ca8458fe56ae413f89cd76f900eecf58a7470d982a4e9c25354d944da7708
primary_output_sha256: 0e0fcdaa0120674747a4729a4fa4018bc4b17285870d895a58d5318629c35c3e
independent_audit_script: 04-computation/jc23_lambda_zero_normalized_response_degree_stratification_independent_audit_thm4293.py
independent_audit_output: 05-knowledge/results/jc23_lambda_zero_normalized_response_degree_stratification_independent_audit_thm4293.out
independent_audit_script_sha256: 3b6138762e2dc32834c3eeb021657dd362b26f66fa72e615d36fbfca450a6106
independent_audit_output_sha256: 990c1217dc13e28a88d29ce056670103a2893aa9cd62e158f26480e2b8211bfa
hash_basis: raw LF bytes
audit: >
  PASS. The primary exact path reconstructs the literal infinity chart,
  normalized branch orders, residue indices, global packet, genus and
  differential-zero ledgers, and Eisenstein-norm sieve. A standard-library
  clean-room path independently reconstructs the chart and all six strata.
  Normal, optimized, and fixed-hash-seed streams byte-match frozen outputs.
---

# THM-4293 -- Lambda-zero normalized response-degree stratification

**PROVED RELATIVE TO THM-4103/4120/4230/4290/4292 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. FIVE OF THE SIX NONCRITICAL `Lambda=0` COLLISION
STRATA FAIL THE DEGREE SIEVE; `r=3` SURVIVES IT. THM-4294 LATER CLOSES THE
SURVIVOR AND REPEATED LOCUS ON THIS SLICE. THE GENERAL `Lambda=0` WALL AND
`JC(2)` REMAIN OPEN.**

## 1. Statement

Work over `C` in the exact-weight-twelve reduced `(2,3)` seam and assume

```text
W=0,              Lambda=U+Z=0,              U!=0.                 (1)
```

Put

```text
c_1=alpha_11+beta_11,          c_2=upsilon_5+xi_10,
c_3=eta+zeta_3,                c_4=Delta+Theta,
c_5=Phi,                       c_6=7168/135-(7/6)Delta.            (2)
```

Define a **noncritical collision order** `r` as follows:

- for `1<=r<=5`, `c_1=...=c_(r-1)=0` and `c_r!=0`;
- for `r=6`, `c_1=...=c_5=0` and `c_6^2+4U!=0`.

> **Theorem.** The normalized source response has genus `18-r`, and the
> Keller map has exactly the two possible degrees
>
> ```text
> d_full=42-2r,                 d_finite=34-2r.                    (3)
> ```
>
> Every exceptional component in its exact `sigma`-degeneration is
> Keller-constant by THM-4292, so all of `(3)` lands on the genus-seven
> central component `C`. Deck equivariance forces its degree to be `4N(a)`
> for an Eisenstein integer `a`. This excludes every `r` except `r=3`.

The complete table is

| `r` | collision indices | genus | full / finite | conclusion |
|---:|:---:|---:|:---:|:---|
| `1` | `10,10` | `17` | `40 / 32` | excluded: `10,8` are not Eisenstein norms |
| `2` | `9,9` | `16` | `38 / 30` | excluded modulo four |
| `3` | `8,8` | `15` | `36 / 28` | survives: `9,7` are Eisenstein norms |
| `4` | `7,7` | `14` | `34 / 26` | excluded modulo four |
| `5` | `6,6` | `13` | `32 / 24` | excluded: `8,6` are not Eisenstein norms |
| `6` | `5,5` | `12` | `30 / 22` | excluded modulo four |

Thus the inherited `34/42` interface does **not** specialize through
`Lambda=0`. The closest proved mechanism is THM-4103's residue-degree
response. The canonical hostile is THM-4291's abstract degree-`42` tail map:
it exists, but the literal wall specialization carrying that tail has total
degree only `30` or `22`. The corrected near miss is counting a double edge
root twice with the simple-root index. The least-used sidecar is the
normalized derivative order `ord_z(Fbar_w)=r`.

## 2. The literal top infinity chart

Let `Q=q_target^(-1)` and use THM-4230's generic response curve

```text
F_Q=(s^2-p)(1-QH(p,sp))-Q s^2/2.                                (4)
```

At top infinity put

```text
z=1/s,               w=s^2/p,               u=w-1,
Fbar=z^14 w^7 F_Q(z^-1,z^-2w^-1).                              (5)
```

The top edge is

```text
Q(1-w)(U+Ww+Zw^2).                                             (6)
```

On `(1)` this becomes

```text
QU(1-w)^2(1+w).                                                (7)
```

Thus `w=1` is double, while `w=-1` remains simple. Direct expansion of the
literal polynomial gives

```text
Fbar(z,0)=-(Q/2)z^12,
Fbar(0,u)=QUu^2(u+2),

[u]Fbar=-Q(c_1z+c_2z^2+c_3z^3+c_4z^4+c_5z^5+c_6z^6
             +(8/3)z^8-3z^10+...).                            (8)
```

The six coefficients in `(8)` are exactly `(2)`. No continuity from
`Lambda!=0` is used.

## 3. Normalization of the double edge root

For `1<=r<=5`, Newton-Hensel normalization gives two smooth
`K=C(q_target)`-rational formal places (not extra global components) with orders

```text
ord_z(u_1)=r,                 ord_z(u_2)=12-r,
ord_z(u_1-u_2)=r.                                                (9)
```

For `r=6`, both have order six:

```text
u_i=lambda_i z^6+...,
2U lambda^2-c_6 lambda-1/2=0.                                  (10)
```

The discriminant condition in the definition of `r=6` makes the roots
distinct. Because the constant field is `C`, both branches are `K`-rational.
Their intersection multiplicity, and hence the local delta invariant, is
exactly `r`.

The toric closure has arithmetic genus eighteen and no other new singularity,
so

```text
g(normalization)=18-r.                                          (11)
```

## 4. The corrected residue indices

THM-4103 gives

```text
phi^*(dA/(2C_target))=Q ds/(F_Q)_p.                             (12)
```

Differentiating `(5)` on the curve and using `ds=-z^-2dz` yields the exact
identity

```text
phi^*(dA/(2C_target))
 =Q z^10 w^5 dz/Fbar_w.                                        (13)
```

On either normalized collision branch,

```text
ord_z(Fbar_w)=r,
ord_z(phi^*eta)=10-r.                                          (14)
```

If `tau=A/C_target` is the target parameter at `O`, then
`eta=-d tau+...`. Hence each collision branch has local map index

```text
e=11-r.                                                        (15)
```

The branch is `K`-rational, and THM-4120 proves `E_q(K)={O}`; both branches
therefore lie over `O`.

The other boundary contributions do not change:

```text
w=-1:                         one rational point, e=11;
AB boundary:                  one rational point, e=1;
quartic carrier:              four conjugate points, e=2 each. (16)
```

Over `C(q_target)`, before any Kummer extension, the quartic remains an
irreducible separable degree-four point: its leading coefficient is
`Z=-U!=0`, and it is the generic minimal polynomial of its degree-four
polynomial carrier (including on the composition wall). Conjugacy puts all
four points simultaneously over `O` or away from it, so the finite alternative
still subtracts exactly eight. A later field extension may split the closed
point but cannot change that total degree. Summing `(15)--(16)` proves `(3)`.

There is an independent completeness check. The listed differential zeros
have total degree

```text
2(10-r)+10+4=34-2r=2(18-r)-2.                                (17)
```

This saturates the canonical divisor of the normalized genus-`18-r` curve;
the listed packet exhausts all ramification, affine and boundary.

## 5. All response degree lands on `C`

After the Kummer base change and graph resolution, the good elliptic target
extends as an abelian scheme. Write the special fibre as

```text
X_0=C+sum_i e_i Gamma_i.                                      (18)
```

For a degree-one target bundle `L`, proper-flat intersection gives

```text
deg(phi_K^*L)=deg(phi_C^*L_0)
              +sum_i e_i deg(phi_(Gamma_i)^*L_0).             (19)
```

Every rational component maps constantly in characteristic zero. THM-4292
proves the same for every exceptional positive-genus component, including
all repeated Newton faces. The central component `C` has multiplicity one,
so `(19)` puts the entire degree `(3)` on `C`. Base-change multiplicities,
dual-graph genus, gluing, and vertical twists cannot absorb map degree.

## 6. The deck eigenline and the five exclusions

The exact deck action is

```text
(S,P)->(xi_12 S,xi_12^2P),
(X_E,Y_E)->(xi_12^4X_E,-Y_E)=[-omega].                          (20)
```

Graph specialization preserves equivariance. THM-4290's intrinsic wall
corollary, equivalently THM-4249's integral projector, identifies the required
eigenline as

```text
ker(T+omega)=O*u,                   O=Z[omega],
deg(a u)=4N(a).                                                (21)
```

For even `r`, both degrees in `(3)` are `2 mod 4`, proving exclusion. For
`r=1`, division by four asks for norms `10` or `8`; for `r=5`, it asks for
norms `8` or `6`. In `Z[omega]`, every rational prime `p=2 mod 3` has even
valuation in a norm. The prime `2` occurs to odd order in `10,8,6`, so none
is a norm. These strata are excluded as well.

For `r=3`, the quotients are

```text
36/4=9=N(3),                 28/4=7=N(3+omega).              (22)
```

The deck lattice therefore supplies no contradiction. This is the unique
noncritical survivor.

## 7. Repeated-discriminant boundary and THM-4291 repair

The theorem makes no exclusion when

```text
c_1=...=c_5=0,                  c_6^2+4U=0.                  (23)
```

There the local Weierstrass discriminant has order strictly greater than
twelve. Its successive square classes and normalized pole divisor require a
separate audit; one may not insert `r=6` into `(3)`. THM-4292 still proves
that every exceptional tail is Keller-constant, but the generic normalized
response degree remains to be stratified.

For THM-4291's displayed specialization,

```text
c_1=...=c_5=0,                  c_6=7168/135.                (24)
```

Its smoothness condition

```text
18225U+12845056!=0                                           (25)
```

is exactly `c_6^2+4U!=0`. Hence it lies in `r=6`: the normalized genus is
twelve and the response degrees are `30` or `22`. Its genus-seven central
component plus the displayed genus-five tail exhaust that genus. The
abstract equivariant degree-`42` tail map is therefore degree-mismatched to
the literal Keller family even before its differential order eight is used.

The first failed implication in the old `34/42` transport was precise: the
edge-index rule used for a simple root requires `Fbar_w!=0` at leading order.
At `Lambda=0`, the double root makes `ord_z(Fbar_w)=r`, subtracting `r` from
each of two local indices.

## 8. Scope and reproduction

This theorem's degree sieve closes five noncritical strata of the
`W=Lambda=0`, `U!=0` wall. It does not itself close `r=3` or the
repeated-discriminant locus `(23)`; THM-4294 later closes both by a different
central-differential mechanism on this slice. The general `Lambda=0` wall,
the walls `U=0`, `Z=0`, and `D=0`, exact-`M=12` seam entry, `JC(2)`, and
`DC(2)` remain open.

Run

```bash
python3 -B 04-computation/jc23_lambda_zero_normalized_response_degree_stratification_thm4293.py
python3 -B -O 04-computation/jc23_lambda_zero_normalized_response_degree_stratification_thm4293.py
PYTHONHASHSEED=0 python3 -B 04-computation/jc23_lambda_zero_normalized_response_degree_stratification_independent_audit_thm4293.py
```

The first path is an exact symbolic reconstruction. The second is a
standard-library clean-room audit. **QED.**
