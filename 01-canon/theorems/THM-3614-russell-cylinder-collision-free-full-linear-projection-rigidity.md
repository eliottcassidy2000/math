---
id: THM-3614
title: "Russell-cylinder collision-free full linear-projection rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every polynomial cylinder graph, no rank-two linear projection of all
  four Russell coordinates (B,C,Y,S) has nonzero constant ordinary
  Jacobian; no collision hypothesis is needed.  THM-3610 already closes
  planes meeting span(B,C).  A transverse plane is injective on x=0, so
  Gwozdziewicz makes any Keller projection an automorphism.  Ordinary degree
  profiles then violate Jung--van der Kulk divisibility, except at
  h=-xq+n, where the source Jacobian has the parameter-independent
  coefficient [q]J=8.  No nonlinear projection, implicit source plane, or
  JC(2) counterexample is claimed.
source: root / full-cylinder degree-divisibility follow-up, 2026-08-21
audit: >
  PASS -- independent hostile audit reconstructed the transverse line inverse,
  every ordinary-degree branch and quadratic cancellation, and the exceptional
  coefficient [q]J=8; normal, optimized, and stored 14,535-gate transcripts
  are byte-identical.
depends_on:
  - THM-3610-russell-cylinder-full-linear-projection-collision-rigidity
related:
  - THM-3605-russell-cylinder-graph-slice-puncture-no-filling
  - THM-3607-russell-cylinder-mixed-projection-degree-seven-gate
  - THM-3611-russell-cylinder-arm-separating-nonlinear-first-coordinate-rigidity
  - THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity
external:
  - "Gwozdziewicz, Injectivity on one line, arXiv:alg-geom/9305008, Theorem 1.1."
  - "Jung--van der Kulk: degree divisibility for polynomial automorphisms of the affine plane."
script: 04-computation/jc2_russell_cylinder_collision_free_full_projection_thm3614.py
output: 05-knowledge/results/jc2_russell_cylinder_collision_free_full_projection_thm3614.out
script_sha256: 268d4c23d8db7715dee23d1d89732eee13df05ad1b077e9d1622802d0fc733ba
output_sha256: 391abbc78d5dd1063ccd39cc3c5a8033d4ddafbd3ba9460e6c6991b787d995bf
hash_basis: raw LF bytes
---

# THM-3614 -- Russell-cylinder collision-free full linear-projection rigidity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This result removes the only collision hypothesis from THM-3610.  It still
concerns polynomial source graphs and linear target projections; nonlinear
target functions and implicit source planes remain separate problems.

All rings and derivatives are over `C`, and `deg` means ordinary total
degree in `(x,q)`.

## 0. Statement

Retain the THM-3561 and Russell functions

```text
D=1+x^2q,
c=xD(D+2),
b=(D-1)(D+2)^2,
e=q(D+3),                                               (1)

B=b+ch,
C=c,
Y=ce+(2b+4)h+ch^2,
S=((b+2)(e+3h^2)+ch(3e+h^2))/8                        (2)
```

for an arbitrary polynomial `h in C[x,q]`.  For every rank-two matrix
`Lambda in Mat_(2x4)(C)`, the theorem claims

```text
(F_Lambda,G_Lambda)^T=Lambda(B,C,Y,S)^T,
Jac_(x,q)(F_Lambda,G_Lambda) notin C*.                  (3)
```

No equality among graph values and no collision condition appears in `(3)`.

## 1. Only the transverse row-space branch remains

Put

```text
V_0=span_C(B,C),                   U=span_C(Y,S),
W=rowspan_C(Lambda).                                      (4)
```

THM-3610 proves, without using any collision, that `(3)` holds whenever
`W intersect V_0` is nonzero.  Indeed, its three branches are the fixed-`C`
Hamiltonian gate, the `(B,C)` factor `c`, and the formal-arm transport for
`B+rho C` paired with `Y+tau C` or `S+sigma Y+tau C`.

It remains to take

```text
W intersect V_0=0.                                      (5)
```

Projection `W -> U` is then an isomorphism.  After a target `GL_2(C)`
change of basis, the outputs are

```text
F=Y+alpha B+beta C,
G=S+gamma B+delta C.                                    (6)
```

On the source line `x=0`, write `phi(q)=h(0,q)`.  Since

```text
B=C=0,                 Y=4phi(q),
S=q+3phi(q)^2/4,                                      (7)
```

the restriction of `(F,G)` is

```text
q |-> (4phi(q),q+3phi(q)^2/4).                         (8)
```

It is injective: the first coordinate recovers `phi`, and the second then
recovers `q`.  Thus, if `(6)` had nonzero constant Jacobian,
Gwozdziewicz's cited one-line theorem would make `(F,G)` a polynomial
automorphism of `A2`.

For a plane polynomial automorphism, Jung--van der Kulk implies that the
ordinary degrees of its two nonconstant components are equal or the smaller
divides the larger.  Sections 2--3 contradict this necessary condition, or
directly contradict constancy in the sole equal-degree boundary.

## 2. Ordinary degree profiles

Write

```text
c_7=x^5q^2,                        k=xq.                 (9)
```

These are the top homogeneous form of `c` and the degree-two quotient that
organizes every cancellation.

### 2.1 Graph degree at least three

Let `d=deg h>=3`, with top form `h_d`.  In `(2),(6)`, the unique top forms
are

```text
F_(2d+7)=c_7 h_d^2,
G_(3d+7)=c_7 h_d^3/8.                                  (10)
```

The `B,C` corrections have strictly smaller degree.  Hence

```text
(deg F,deg G)=(2d+7,3d+7).                             (11)
```

The difference is `d`, strictly between zero and `2d+7`; therefore the
smaller degree does not divide the larger.

### 2.2 Graph degree at most one

If `d<=1`, the base terms `ce` and `be/8` dominate:

```text
F_11=c_7 k^2,
G_13=c_7 k^3/8,
(deg F,deg G)=(11,13).                                 (12)
```

Again the smaller degree does not divide the larger.

### 2.3 Quadratic graph

Let `d=2` and write `h=h_2+h_1+n`, where `h_i` is homogeneous of degree
`i` and `n in C`.  Direct grouping of the degree-eleven and degree-thirteen
layers gives

```text
F_11=c_7(k+h_2)^2,
G_13=c_7(k+h_2)^3/8.                                  (13)
```

If `h_2!=-k`, the profile is again `(11,13)`.  If `h_2=-k` but `h_1!=0`,
the next nonzero layers are

```text
F_9=c_7 h_1^2,
G_10=c_7 h_1^3/8,
(deg F,deg G)=(9,10),                                  (14)
```

which also violates degree divisibility.

Only one shape is not decided by degrees:

```text
h=-xq+n.                                               (15)
```

## 3. The final cancellation has an unavoidable `q` coefficient

For `(15)`, restrict the first derivatives of `(6)` to `x=0`.  Exact
expansion of `(1),(2)` gives

```text
F(0,q)=4n,
G(0,q)=q+3n^2/4,

F_x(0,q)=8q+3n^2+3alpha n+3beta,
F_q(0,q)=0,
G_q(0,q)=1.                                           (16)
```

Consequently

```text
Jac(F,G)(0,q)=8q+3n^2+3alpha n+3beta,
[q] Jac(F,G)=8.                                       (17)
```

The Jacobian is not constant, independently of `gamma,delta` and of all
other coefficients.  This closes `(15)`, completes the transverse branch,
and proves `(3)`.

## 4. Preservation/loss ledger and boundaries

| item | content |
|---|---|
| source | every polynomial graph `w=h(x,q)` in the stabilized source |
| target | every rank-two linear projection of `(B,C,Y,S)` |
| retained | polynomiality; no collision or symmetry is assumed |
| first sidecar | injectivity of the transverse projection on `x=0` |
| second sidecar | Jung--van der Kulk degree divisibility |
| final exceptional test | the parameter-independent coefficient `[q]J=8` |
| still lost | nonlinear target coordinates, `Y`-dependent shears, and implicit non-graph source planes |

The theorem is strictly stronger than THM-3610, but the older result remains
the useful formal-row-space decomposition.  The even-fold non-graph planes
of THM-3612 are not covered: they are not polynomial graphs over `(x,q)`.
THM-3611 covers a broad nonlinear first-coordinate family, but not arbitrary
nonlinear target pairs.  No planar Keller counterexample and no result about
all polynomial maps `A2->A2` is claimed.

## 5. Exact companion contract

The companion must verify without truth-bearing `assert` statements:

- the Russell polynomial identities and line restriction `(7),(8)`;
- the top homogeneous forms `(10),(12)--(14)` over explicit all-degree and
  generic quadratic coefficient universes;
- the degree-divisibility failures for a stated integer grid, with the
  symbolic inequalities supplied by the proof;
- the exact exceptional formulas `(16),(17)` with symbolic
  `alpha,beta,gamma,delta,n`;
- hostile cancellations `h_2=-xq`, `h_1=0`, and nonzero `h_1`; and
- finite row-space controls confirming that `(6)` is exactly the branch
  transverse to `span(B,C)`.

The companion checks the displayed algebra and finite support census.  The
all-degree inference uses the cited one-line and plane-automorphism theorems;
no finite degree box is substituted for those results.
