---
id: THM-2725
title: "Split even-Faber nonzero-first-flux unified closure"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  The entire
  chosen-sheet split polynomial exact-square-prefix degree-22 even-Faber
  chart with nonzero first flux is physically empty.  If y is nonzero, the
  arbitrary-b prime-23 geometry has five q-pole points while THM-2723 permits
  at most two on the source.  If y is identically zero, THM-2724's degree-23
  resultant forces the response into the constant field and contradicts the
  third flux.  This includes every B_0, including zero, but not lambda=0,
  odd Faber seeds, the broader split branch, or JC(2).
source: root/split-even-nonzero-flux-unified-closure-2026-07-28
audit: thm2704-hostile-audit-2026-07-28-unified-even-closure (independent two-chart proof and pole-pullback audit)
depends_on:
  - THM-2713-split-prime23-component-divisor-budget-and-perfect-power-normal-form
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2724-split-even-faber-y-zero-resultant-closure
related:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2718-split-prime23-five-pole-rational-primitive-closure
  - THM-2726-a21-transverse-integral-split-response-three-pole-closure
---

# THM-2725 -- the nonzero-flux split even-Faber chart is empty

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**

The two apparent holes in the prime-23 chart were artifacts of one scale
choice.  Keeping the scale independent of `B_0` covers `B_0=0` whenever
`y!=0`; direct elimination covers the complementary identity `y=0`.  The
third Faber observable supplies the continuation obstruction in both cases.

## 1. Exact statement

Let a polynomial Keller pair over `C` lie on a chosen split sheet of the
polynomial exact-square-prefix terminal chart.  After a target shear, suppose
its reduced mate is

```text
Q=E_22+B_0 E_14+C_0 E_10+D_0 E_6+E_0 E_2,            (1)
```

so all eleven odd Faber coefficients vanish.  Assume its first two exact
fluxes are

```text
Phi_Q=lambda in C*,                 Psi_Q=W_0 in C.   (2)
```

Then no such physical Keller trajectory exists.

There is no hypothesis on `B_0`, and the centered scale coordinate `y` may
be zero or nonzero.  The conclusion is confined to `(1)--(2)`.

## 2. Nonzero `y`: five target poles exceed source capacity

First suppose that `y` is a nonzero element of the source function field.
Choose any `rho in C*`; it need not satisfy `rho^2=B_0`.  Put

```text
b=B_0/rho^2,       c=C_0/rho^3,       d=D_0/rho^4,
e=E_0/rho^5,       w=W_0/rho^6,

t=rho/y,           v=u/y^2,           zeta=Z/y^3.     (3)
```

THM-2713's arbitrary-`b` addendum applies, including when `b=0`.  The
response lies on the geometrically integral projective curve

```text
C_(b,c,d,e,w,eta):
F2_b=0,                  zeta F1_b^4=eta t^23,
eta=7496192^4 lambda^4/rho^23!=0.                     (4)
```

On its normalization `X`, the chosen-sheet physical coefficient satisfies

```text
div(q_X)=5N-3O,                                      (5)
```

where `O=O_1+...+O_5` consists of five distinct normalization points.
Thus `q_X` has five distinct pole points.

The physical coefficient functions give a rational map

```text
gamma_0:P1_x ---> C_(b,c,d,e,w,eta).                  (6)
```

It is nonconstant.  Indeed, if its affine generic value were constant, then
`t,v,zeta` would be constant.  Equation `(3)` would make `y,u,Z` constant,
and the exact chosen-sheet reconstruction

```text
q_source=-7496192 lambda t^5/(rho^5 f1_b)             (7)
```

would make `q_source`, hence `T=q_source^2`, constant as well.  Therefore
`d_ctr=u/T`, `s_ctr=y/11`, and every input of `R_Q` would be constant.  This
contradicts the split Hamiltonian identity

```text
R_Q'=kappa/U!=0.                                      (8)
```

Geometric integrality lets `(6)` factor through the normalization.  Since
the response and its normalization are projective, the nonconstant map
extends to a finite surjective morphism

```text
gamma:P1_x -> X.                                      (9)
```

Every fibre `gamma^(-1)(O_i)` is nonempty and the five fibres are disjoint.
Pullback of `(5)` therefore gives at least five distinct pole points of

```text
gamma^*(q_X)=q_source=A_src/U.                        (10)
```

THM-2723 proves that `(10)` has at most two pole points on `P1_x`.  This
contradiction closes every `y!=0` specialization, uniformly in `B_0`.

## 3. Zero `y`: constant-field resultant trap

It remains to suppose

```text
y=0 identically.                                      (11)
```

This locus is not a point at `t=infinity` in `(3)` and is not inferred by
closure of the prime-23 family.  THM-2724 treats it directly.  The first two
flux equations at `(11)`, with `Z=q^4`, have an exact eliminant

```text
P_23(q)=0,
[q^23]P_23=96059601,
[q^0]P_23=-992137445376 lambda^3.                     (12)
```

The fixed leading coefficient makes `P_23` nonzero under every specialization
of `B_0,C_0,D_0,E_0,W_0`, while `lambda!=0` excludes `q=0`.  Hence
`q in C(x)` is algebraic over `C` and lies in `C*`.  The second flux retains
a fixed nonzero cubic leading coefficient in `u`, so `u in C` as well.
Consequently all inputs of `R_Q` are constant, contradicting `(8)`.

The alternatives `y!=0` and `y=0` exhaust the source function field.  This
proves the theorem.

## 4. Scope and exact residual

The closed chart is exactly

```text
split polynomial exact-square prefix;
reduced degree 22;
even-Faber bank (1);
nonzero chosen-sheet first flux lambda.               (13)
```

It includes arbitrary `B_0` and both scale cases.  It uses no generic genus
claim and does not require the abstract response curve itself to have positive
genus.

The theorem does not treat `lambda=0`, any of the eleven odd Faber seeds, a
nonsplit or nonpolynomial exact prefix, another reduced degree, or a source
chart not represented by `(13)`.  The broader split branch, `JC(2)`, and
`DC(2)` remain open.

No new finite computation is used in this synthesis.  THM-2713 supplies the
uniform arbitrary-`b` divisor geometry, THM-2723 supplies the all-degree
source capacity, and THM-2724 supplies the complementary constant-field
boundary.

An independent hostile audit rechecked the exhaustive `y=0`/`y!=0` field
dichotomy, the arbitrary nonzero scale and `B_0=0` specialization, the exact
chosen-sheet `q` identification, the constant-map contradiction, and the
normalization/surjectivity pullback of all five pole fibres.  It found no
additional hypothesis and certified the residual in Section 4.
