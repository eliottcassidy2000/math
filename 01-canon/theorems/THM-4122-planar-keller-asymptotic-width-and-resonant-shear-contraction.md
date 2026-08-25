---
id: THM-4122
title: "Planar Keller intrinsic asymptotic-width bound"
status: >
  PROVED FROM CITED NGUYEN VAN CHAU, JELONEK--LASON, AND LANG + THM-3544 +
  VERIFIED-EXACT arithmetic controls. For every irreducible component of the
  nonproperness curve of a hypothetical planar Keller counterexample, the
  normalization is A1 and the intrinsic target pole pair is
  (rho*d,rho*e), where (d,e) is the reduced polynomial-degree pair. In every
  affine-linear source chart rho*max(d,e) is at most both D-1 and the chart's
  narrow source-direction width. Hence gcd(deg P,deg Q)>=2 and every chart
  has width at least max(4,d,e). The rho<=2 conclusion for (72,108) is only
  CONDITIONAL on an actual chart with narrow width six. JC(2) remains OPEN.
source: planar-jacobian-squeeze / 2026-08-25
audit: >
  PASS. An independent agent rederived the normalization-factorization
  argument, separated intrinsic pole multiplier rho from the cover-inflated
  parameter m printed in Nguyen's theorem, checked the D-1 and fixed-chart
  width bounds against the primary statements, and rejected an unconditional
  width-six reading of THM-3586. Normal and optimized executions of the exact
  arithmetic/hostile certificate byte-match the frozen output.
depends_on:
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
related:
  - THM-3586-nodal-cylinder-cap38-width-period-and-second-conductor-keller-gates
  - THM-4124-planar-keller-integral-degree-ratio-all-vertex-shear
external:
  - "Nguyen Van Chau, Non-proper value set and the Jacobian condition, Ann. Polon. Math. 84 (2004), Theorem 1, arXiv:math/0305088."
  - "Z. Jelonek and M. Lason, Quantitative properties of the non-properness set of a polynomial map, Manuscripta Math. 156 (2018), Theorem 3.2 and Corollary 3.5, arXiv:1411.5011v2."
  - "J. Lang, Newton polygons of Jacobian pairs, J. Pure Appl. Algebra 72 (1991), 39--51, DOI:10.1016/0022-4049(91)90128-O."
script: 04-computation/jc2_asymptotic_width_shear_controls_thm4122_4124.py
output: 05-knowledge/results/jc2_asymptotic_width_shear_controls_thm4122_4124.out
script_sha256: 7fd2f9a2e32b199489506c22bf0f34c849ab5b0f4750c20e78d40038bd491f5a
output_sha256: d1bf27ab87eb78b79667854c3c57181942e248d83145a729e3460e641c8e2ca4
hash_basis: raw LF bytes
---

# THM-4122 -- intrinsic asymptotic width of a planar Keller map

**PROVED FROM CITED RESULTS + THM-3544 + VERIFIED-EXACT controls; JC(2)
OPEN.** Let

```text
F=(P,Q):C^2 -> C^2,                 Jac(P,Q) in C*,       (1)
m=deg P>1,   n=deg Q>1,   G=gcd(m,n),
(d,e)=(m/G,n/G),                    E=max(d,e),
D=max(m,n)=GE.                                           (2)
```

Assume that `F` is not proper, and let `C` be any irreducible component of
its nonproperness set `S_F`. In a fixed affine-linear source chart
`X=(x_1,x_2)`, put

```text
w_X(F)=min_j max(deg_(x_j) P,deg_(x_j) Q).               (3)
```

Then the normalization of `C` is `A1`. If `t` is its coordinate, the two
target coordinates have their unique poles at `t=infinity`, of orders

```text
(poleord_infinity P_C,poleord_infinity Q_C)
                              =(rho_C d,rho_C e)          (4)
```

for an integer `rho_C>=1` depending on the component. Moreover,

```text
rho_C E <= min(D-1,w_X(F)).                              (5)
```

In particular every hypothetical planar Keller counterexample satisfies

```text
G>=2,                       w_X(F)>=max(4,E)              (6)
```

in every affine-linear source chart.

## 1. The multiplier in Nguyen's parametrization is not intrinsic

Make a generic source-linear change, independently of the fixed chart `X`, so
that both `P,Q` are monic in the chosen second variable. Nguyen's Theorem 1
then gives a **surjective** polynomial parametrization

```text
gamma:A1 -> C                                               (7)
```

whose coordinate degrees are `(M d,M e)` for a positive integer `M`.
The extension of `gamma` to `P1` lifts to the projective normalization
`tilde(Cbar)`. The lift is nonconstant, hence finite and surjective. By Luroth
(equivalently Riemann--Hurwitz), `tilde(Cbar)` has genus zero and is `P1`.
Every finite source point maps into affine `C`, so the inverse image of the
normalization boundary lies in the singleton `{infinity}`. Surjectivity forces
every boundary point to have a preimage; there is exactly one normalization
place at infinity. Thus the affine normalization is `P1-{infinity}=A1`. Let

```text
nu:A1 -> C                                                  (8)
```

be that normalization. Write Nguyen's parametrization as
`gamma=nu o h_N`. Its lift `h_N:A1->A1` is polynomial because finite source
points stay in affine `C`. If the intrinsic pole orders of the coordinate
functions are `(alpha,beta)`, then

```text
(M d,M e)=(deg h_N)(alpha,beta).                          (9)
```

Thus `alpha/beta=d/e`. Coprimality gives `(4)` for a positive integer
`rho_C`, and `(9)` gives

```text
M=(deg h_N)rho_C.                                         (10)
```

The raw integer in an arbitrary parametrization may therefore be inflated by
a polynomial cover. Only `rho_C` is an invariant of the normalized target
curve. This distinction is essential in applying any degree bound.

## 2. Two independent bounds on the intrinsic degree

Jelonek--Lason Theorem 3.2 covers `S_F` by polynomial curves of degree at
most `D-1`. Choose `p` in `C` but on no other component. A nonconstant
irreducible parametric curve through `p` lies in `S_F`; its closure lies in
one irreducible component, and `p` forces that component to be `C`. Equal
dimension makes its closure all of `C`. Its parametrization lifts uniquely as
`nu o h_J` with polynomial `h_J:A1->A1`, and its degree is

```text
(deg h_J)rho_C E.                                         (11)
```

Consequently `rho_C E<=D-1`.

Their Corollary 3.5 applied in the fixed chart `(3)` instead supplies such a
parametrization of degree at most `w_X(F)`. The same factorization gives
`rho_C E<=w_X(F)`, proving `(5)`. This is a target-curve statement obtained
from escaping source directions; it is not a claim that a source fibre has
one place at infinity.

If `G=1`, `(5)` would read `E<=E-1`; hence `G>=2`. Also `(5)` and
`rho_C>=1` give `w_X(F)>=E`. Finally, THM-3544's all-direction fibre gate
gives

```text
max(deg_(x_j)P,deg_(x_j)Q)>=4                            (12)
```

for both `j`, and `(6)` follows.

## 3. The exact conditional squeeze at `(72,108)`

For `(m,n)=(72,108)` one has

```text
G=36,                    (d,e)=(2,3),          E=3.      (13)
```

If an independently established affine-linear source chart has

```text
w_X(F)=6,                                                    (14)
```

as would happen in an attained narrow width-`(4,6)` chart, then `(5)` gives

```text
rho_C in {1,2},              pole pair (2,3) or (4,6).    (15)
```

Lang's Newton similarity adds an exact directional sidecar. For each source
variable in `X`,

```text
(deg_(x_j)P,deg_(x_j)Q)=(2r_j,3r_j),
w_X(F)=3r_min,                       r_min=min_j r_j.      (16)
```

THM-3544 forces `r_j>=2`, so no nonautomorphic reduced `2:3` chart has
`w_X<=5`. Combining `(5)` and `(16)` gives

```text
rho_C<=min(35,r_min).                                     (17)
```

Thus `(14)` is exactly `r_min=2` and yields `(15)`. THM-3586 says that
`(4,6)` is the first unexcluded width in its cited reduced cell; it does
**not** prove that every surviving `(72,108)` pair attains `(14)`. Without a
small-width chart, the global `D-1` term in `(17)` remains `rho_C<=35`.

## 4. Hostiles and exact controls

- Intrinsic pole pair `(rho*d,rho*e)` followed by a cover of degree `h`
  produces `(h*rho*d,h*rho*e)`. Thus Nguyen's displayed `M` cannot replace
  `rho_C` in `(5)`.
- A hypothetical chart with `w_X=30` at `(72,108)` permits
  `rho_C=1,...,10`; the set `{1,2}` is not an unconditional degree-pair
  consequence.
- Different components may have different `rho_C`. Neither cited theorem
  synchronizes them.

The companion certificate checks all divisibility, width, inflation, and
hostile rows in exact integer arithmetic. It does not certify the cited
algebraic-geometric theorems themselves.

## 5. Boundary

The theorem supplies necessary conditions for a nonproper planar Keller map.
It does not prove that `S_F` is nonempty for an arbitrary displayed formal
scaffold, produce a width-six chart, bound the number of components, identify
their lower Puiseux jets, or prove `JC(2)`.
