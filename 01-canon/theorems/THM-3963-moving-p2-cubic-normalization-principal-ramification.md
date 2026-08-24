---
id: THM-3963
title: "Moving P-squared cubic normalization has a principal ramification arm"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For every nonzero
  c(t), the nonnormal cubic
  T^3-3PT-c(t)P^2 has exact finite normalization B=A[w], where generically
  w=c(t)v. The two charts D(w+3) and D(w) are smooth even over multiple zeros
  of c, so B is the actual regular normalization. The prime E2=V(w+2) is a
  reduced principal ramification divisor with coordinate ring k[t,1/c]. Any
  same-field Keller affine-plane open must delete E2, making w+2 a forbidden
  nonconstant unit. This closes the moving scalar P^2 conductor-debt family,
  not arbitrary P^2q2(P,t), the other repeated-hidden-factor debt, or JC(2).
source: jc-cohn3709 + jc-extra-debt-local / post-THM-3961 moving conductor-debt lane, 2026-08-24
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS (jc_zero_debt_lift and
  jc_degree6_one_place, 2026-08-24). Both audits reconstructed the cubic
  field, integral moving coordinate, determinantal relations, exact chart
  presentations and their primitive-linear injectivity, smoothness at
  arbitrary multiple zeros of c, and the reduced principal E2 ramification
  prime. A reproducibility audit replaced a chart determinant check that
  divided by c with the cancellation-free identity
  det+3T(w+2)=cP-wT, so the exact companion also covers c=0 fibres without
  localization. Normal and optimized 36-gate runs byte-match the frozen
  output; hashes and documentation checks pass.
depends_on:
  - THM-3961-arbitrary-q-hidden-repetition-normality-and-conductor-debt
related:
  - THM-3960-natural-one-parameter-cubic-normal-monogenic-closure
  - THM-3862-russell-finite-completion-nonmonogenic-branch-contract
script: 04-computation/jc2_moving_p2_normalization_principal_ramification_thm3963.py
output: 05-knowledge/results/jc2_moving_p2_normalization_principal_ramification_thm3963.out
script_sha256: b19847721eed96e5d86696a935584186b8add9fbbcf9ba27796c1acd5085b54c
output_sha256: c840f7802c74369bf48f04a2e2f59a5031f57c09ee590536d28630865c6fa990
semantic_sha256: a7338a3a38fa1fa510734def07cd2f74a4cf6337836114e17345d8e9a8ca4e30
hash_basis: raw LF bytes
---

# THM-3963 -- the branch that escapes to infinity supplies the missing normalization chart

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of
characteristic zero. Let

```text
c=c(t) in k[t],              c!=0,
R=k[P,t],
F=T^3-3PT-cP^2,
A=R[T]/(F),                  X_0=Spec A.                 (1)
```

THM-3961 detects `(1)` as `P^2` conductor debt: the monogenic surface is
singular along `(P,T)`. The present theorem computes its actual finite
normalization, including every zero of the moving coefficient, and excludes a
same-function-field Keller plane atlas there.

## 1. The cubic field and the integral moving coordinate

Put `K=k(t)` and introduce a transcendental parameter `v` over `K`. Define

```text
P=v^2(3+cv),                 T=v(3+cv).                  (2)
```

Then `F=0` and `P=vT`. Since `c` is nonzero in `K`, the polynomial map

```text
v -> P=v^2(3+cv)                                           (3)
```

has degree three. Hence `[K(v):K(P)]=3`. Moreover `v=P/T`, so
`K(P,T)=K(v)`. The monic cubic `F` is therefore the minimal polynomial of
`T` over `K(P)`. In particular

```text
A is a domain,                 Frac(A)=K(v).              (4)
```

The coordinate `v` need not be integral across a zero of `c`: the branch
`v=-3/c` runs to infinity there. The correct finite coordinate is

```text
w=cv=cP/T.                                                (5)
```

It satisfies the monic equations

```text
w^2+3w-cT=0,
w^3+3w^2-c^2P=0.                                         (6)
```

and is therefore integral over `A`. Put

```text
B=A[w] subset K(v),                 X=Spec B.             (7)
```

Then `B` is finite and birational over `A`. The parameter identities give

```text
wT=cP,
w(w+3)=cT,
P(w+3)=T^2,
c^2P=w^2(w+3).                                           (8)
```

The cancellation used for the third identity is valid in the domain `(4)`;
alternatively all four identities follow directly from `(2),(5)`.

## 2. Two charts prove that B is the regular normalization

The opens `D(w)` and `D(w+3)` cover `X`, because their defining elements
differ by the unit `3`.

On `D(w+3)`, equations `(8)` eliminate `P` and give

```text
B_(w+3)
 =k[t,T,w,(w+3)^(-1)]/[cT-w(w+3)],
P=T^2/(w+3).                                             (9)
```

The displayed hypersurface is a domain: its equation is primitive and linear
in `T`, and its fraction field is `K(v)`. Its decisive partial derivatives are

```text
partial_T=c,                    partial_w=-(2w+3).        (10)
```

A singular point would have `c=0,w=-3/2`, but the equation at those values is
`9/4=0`, impossible. Thus this chart is smooth, independently of the
multiplicities of the roots of `c`.

On `D(w)`, eliminate `T` instead:

```text
B_w
 =k[t,P,w,w^(-1)]/[c^2P-w^2(w+3)],
T=cP/w.                                                  (11)
```

Again the equation is primitive and linear in `P`, and the chart has fraction
field `K(v)`. Its decisive partials are

```text
partial_P=c^2,                  partial_w=-3w(w+2).       (12)
```

Because `w` is a unit, a singular point would require `c=0,w=-2`; the equation
would then read `-4=0`. This chart is also smooth over arbitrary multiple
zeros of `c`.

Equations `(9)` and `(11)` are not merely birational models. Each maps
surjectively to the indicated localization of `B`; both source rings are
domains with the same fraction field, so the maps are injective as well.
Consequently they are the exact chart presentations. They cover `X`, hence

```text
B is regular and normal.                                  (13)
```

Since `B` is finite and birational over `A`, `(13)` identifies it with the
integral closure of `A` in `K(v)`.

### 2.1 What happens at a zero of c

The two charts retain both vertical limits that the single `v` coordinate
cannot see. On the first chart, setting `c=0` forces `w=0` and leaves

```text
P=T^2/3.                                                  (14)
```

On the second, setting `c=0` forces `w=-3`, leaves `P` free, and gives
`T=0`. Thus the moving point `v=-3/c` is not discarded; it becomes the
`w=-3` vertical component. This explains geometrically why simply declaring
the normalization to be `k[t,v]` would fail when `c` is nonconstant.

Over the target zero section `P=T=0`, equation `(6)` gives

```text
w^2(w+3)=0.                                               (15)
```

The addresses `w=0` and `w=-3` are the ramified graph point and its etale
companion, respectively. The second-chart different below has value `9` at
`w=-3`.

## 3. A principal ramification prime closes every coefficient c(t)

Regard `(11)` as a finite algebra over the target ring `R=k[P,t]`. Its monic
equation in `w` is `(6)`, so its relative differential module on `D(w)` is

```text
Omega_(B_w/R)=B_w dw/[3w(w+2)dw].                        (16)
```

The factor `w` is a unit on this chart. The other factor defines

```text
E_2=V(w+2).                                               (17)
```

It is a prime divisor, because

```text
B/(w+2)
 =k[t,P]/(c^2P-4)
 =k[t,1/c].                                               (18)
```

In particular `E_2` lies wholly in `D(w)`, even over a coefficient with many
or multiple zeros. Equation `(16)` makes its generic point ramified, and

```text
d[3w(w+2)]/dw at w=-2 = -6                              (19)
```

shows that it is reduced in the ramification divisor. Most importantly,
`E_2` is the complete zero divisor of the global element `w+2`.

Suppose the same cubic function field admitted a planar Keller map with
target coordinates `(P,t)`. Since `X` is the finite normalization,
normalization-form Zariski Main would identify its source with an open

```text
U isomorphic to A2 -> X                                  (20)
```

on which `X -> A2_(P,t)` is etale. The ramification prime `E_2` must lie in
`X minus U`. Hence `(w+2)|_U` has no zero and is a unit. But `w+2` is
nonconstant in `K(v)`, whereas

```text
Gamma(A2,O)^*=k*.                                        (21)
```

This contradiction excludes every `c(t)!=0`, including constants and
polynomials with arbitrarily multiple zeros.

## 4. Scope and handoff

For constant nonzero `c`, the coordinate `v=w/c` is global and `(17)` is the
familiar ramification line `v=-2/c`. The theorem shows exactly how that line
survives when `c` moves: the correct global equation is `w+2=0`, while the
second normalization chart adds the branch that escaped through `v=-3/c`.

The case `c=0` is reducible and is already closed componentwise by THM-3960.
The present result closes only

```text
q(P,t)=c(t)P^2.                                          (22)
```

It does not yet close arbitrary `P^2q_2(P,t)`, repeated nonzero adjusted
hidden factors from THM-3961, a nonmonogenic cubic field outside that grammar,
or JC(2). **QED.**

## Reproduction

```bash
python3 04-computation/jc2_moving_p2_normalization_principal_ramification_thm3963.py
python3 -O 04-computation/jc2_moving_p2_normalization_principal_ramification_thm3963.py
```

The frozen transcript is
`05-knowledge/results/jc2_moving_p2_normalization_principal_ramification_thm3963.out`.
