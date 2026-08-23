---
id: THM-3723
title: "Automorphic Cohn complete E-plus/E-minus two-right nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For arbitrary nonzero polynomial parameters
  u,v, neither alternating two-left decoration of M0 E_+(v)E_-(u) is a
  Jacobian matrix.  If either right parameter is nonconstant, leading-form
  directional-curl gates prevent the first row closure.  In the remaining
  constant cell, a row closes precisely at uv=-1; the right matrix then has
  R^3=-I, and both closed rows are scalar forms of the Broughton polynomial
  L+L^2S.  THM-3716 prevents the second closure.  This closes one complete
  nonidentity two-right order, not E_-(u)E_+(v), longer words, or JC(2).
source: root + jc-sparse-direct-search / 2026-08-22
audit: >
  PASS.  An independent audit rederived all three nonconstant leading-form
  regimes, the positive-degree left-parameter valuation descent, both
  constant curls and their classification, the projective C3 relation, the
  common Broughton potential, and both complementary debt signs.  The exact companion
  checks the arbitrary-polynomial exposed curls, both constant
  classifications, R^3=-I and R^2-R+I=0, the Broughton source coordinates
  and both remaining debts, 110 nonconstant and 200 constant hostile cells,
  the resonant homogeneous kernel through degree six, and finite Hamiltonian
  cokernels through degree 12.  Normal and optimized runs byte-match.
depends_on:
  - THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
related:
  - THM-3709-cohn-alternating-two-by-two-elementary-decoration-nonentry
script: 04-computation/jc2_automorphic_cohn_c3_constant_resonance_thm3723.py
output: 05-knowledge/results/jc2_automorphic_cohn_c3_constant_resonance_thm3723.out
script_sha256: d968ecc5b3daf25eec1b501fff65feaa94e81bd04f9b14b2dfe020ee83d5c5ff
output_sha256: 4cf5eaa1cd9600eb35a1172d43dd0ac65f50594b46bf87f0d6c8dd80a38a0715
semantic_sha256: 3d3b50d3a9a571026aaf5819c258d3879367b3d29fd991fda52db4070262902b
hash_basis: raw LF bytes
---

# THM-3723 -- one complete two-right word closes only at C3

**PROVED + VERIFIED-EXACT.**  THM-3721 closes both single-right-factor cells
around the automorphic Cohn core

```text
M_0=[4T^2,2XT-1; 1+2XT,X^2].                           (1)
```

The first genuinely reduced successor has two alternating right factors.
It does contain closed rows, but only on a sharply defined modular resonance.

Let `k` be a characteristic-zero field and take nonzero
`u,v in k[X,T]`.  Put

```text
R=E_+(v)E_-(u)=[1+uv,v;u,1],
N=M_0R,                                                det N=1.    (2)
```

For arbitrary `f,g,h in k[X,T]`, neither

```text
E_+(f)E_-(h)N,                    E_-(g)E_+(h)N        (2a)
```

is a Jacobian matrix.  Write the rows of `N` as `alpha,beta`.  If either
`u` or `v` is nonconstant, neither `beta+h alpha` nor `alpha+h beta` is
closed.  If `u,v` are constants, the lower exposed row is closed if and only
if

```text
uv=-1,                         h=-v.                   (3)
```

The upper exposed row `alpha+h beta` is closed if and only if

```text
uv=-1,                         h=u.                    (4)
```

In either case no second polynomial left shear closes the complementary row.
Consequently neither alternating two-left decoration of `N` is a Jacobian
matrix.

## 1. Every nonconstant right parameter dies at the first curl

All degrees in this section are total degrees, and subscripts denote leading
homogeneous forms.

First suppose `q=deg u>=1` and `p=deg v>=1`.  The leading rows of `N` are

```text
alpha_(p+q+2)=(4T^2u_qv_p,0),
beta_(p+q+2) =(2XTu_qv_p,0).                           (5a)
```

If `deg h=m>=1`, the leading first component of either exposed row is the
corresponding entry in `(5a)` times `h_m`; its `T` derivative is nonzero,
because the product is nonzero and divisible by `T`.  Every possible
second-component curl is at least `q` degrees lower.  If `h` is constant,
the two leading first components are

```text
2Tu_qv_p(X+2hT),                  2Tu_qv_p(2T+hX).     (5b)
```

Their `T` derivatives are again nonzero.  Thus neither row closes.

Next suppose `u` is nonconstant and `v=c in k*`.  The leading first
components of `alpha,beta` are

```text
2T(X+2cT)u_q,                     X(X+2cT)u_q,         (5c)
```

while both second components have degree two.  A positive-degree `h`
multiplies one of `(5c)` and makes it uniquely highest.  For constant `h`,
the leading component is one of

```text
u_q(X+2cT)(X+2hT),       u_q(X+2cT)(2T+hX).            (5d)
```

No nonzero polynomial divisible by `X+2cT` lies in `k[X]`; hence every
displayed `T` derivative is nonzero.

Finally suppose `u=c in k*` and `v` is nonconstant.  The leading rows are

```text
4T^2v_p(c,1),                         2XTv_p(c,1).     (5e)
```

The leading curl of a scalar multiple `F(c,1)` is

```text
(c partial_T-partial_X)F.                              (5f)
```

The kernel of `(5f)` is `k[T+cX]`.  For positive-degree `h`, the relevant
`F` is `4T^2v_ph_m` or `2XTv_ph_m`; for constant `h`, it is one of the two
expressions in `(5b)` with `u_qv_p` replaced by `v_p`.  Every such `F` is
nonzero and divisible by `T`, whereas no nonzero member of `k[T+cX]` is
divisible by `T`.  Again neither row closes.

It remains only to justify that when `u,v in k*`, a closing `h` must itself
be constant.  Put

```text
w=(uX+2(1+uv)T, X+2vT),
D=w_1 partial_T-w_2 partial_X.                          (5g)
```

The quadratic leading rows are `2Tw` and `Xw`.  If `h_m` has positive
degree, closure of the lower or upper exposed row would respectively make

```text
F=2Th_m                 or                 F=Xh_m       (5h)
```

satisfy

```text
DF+(1+2uv)F=0.                                         (5i)
```

For the first choice, take the least `T`-valuation of `F`.  Since
`w_1(X,0)=uX`, the unique term one valuation lower is a nonzero multiple of
`uX`; contradiction.  For the second, use the least `X`-valuation and
`w_2(0,T)=2vT`.  Thus `h` is constant in every surviving constant cell.

## 2. Exact classification of the two constant exposed rows

Direct multiplication and differentiation give

```text
curl(beta+h alpha)
 =2{u(h+v)X+[4huv+3h-v]T},                            (5)

curl(alpha+h beta)
 =2{u(hv+1)X+[-hv+4uv+3]T}.                           (6)
```

Because `u,v` are nonzero, the `X` coefficient in `(5)` first forces
`h=-v`.  Its `T` coefficient then becomes `-4v(uv+1)`, proving `(3)`.
Likewise, the `X` coefficient in `(6)` forces `hv=-1`; the `T` coefficient
then becomes `4(uv+1)`, proving `(4)`.  This is a complete symbolic
classification, not a bounded search.

The zero-parameter boundaries reduce to the single-factor cells in THM-3721:
for example `u=0` gives the upper-shear lower survivor `h=v/3`, while `v=0`
gives the lower-shear lower survivor `h=0`.  They are deliberately excluded
from the nondegenerate statement above.

## 3. The resonance is the modular C3 factor

At `u=-1/v`, the right matrix is

```text
R_v=[0,v;-1/v,1].                                      (7)
```

It obeys

```text
R_v^2-R_v+I=0,                    R_v^3=-I.             (8)
```

Thus its image in `PSL_2(k)` has order three.  The first positive
two-right closure is not an accidental coefficient cancellation: it is the
`C_3` free-factor relation in the modular-group grammar.  The parameter `v`
is only a diagonal gauge on that projective orbit.

## 4. Both closed rows are the same Broughton line

Introduce source coordinates

```text
L=2vT-X,
M=X+vT,
S=-M/(3v).                                             (9)
```

Their Jacobian is

```text
J_(X,T)(L,S)=1.                                        (10)
```

Set

```text
Q=L+L^2S
 =-(2vT-X)[2v^2T^2+vTX-X^2-3v]/(3v).                 (11)
```

At the resonance, direct differentiation yields

```text
dQ=beta-v alpha,
alpha-(1/v)beta=-(1/v)dQ.                              (12)
```

So `(3)` and `(4)` expose the same gradient line.  The two apparently
different row closures are one Broughton component seen with different
target scalings.

The complementary curls, expressed in `(L,S)`, are

```text
curl_(X,T)(alpha)=6S,
curl_(X,T)(beta) =6vS.                                 (13)
```

After the lower exposure, a final left shear would therefore require

```text
J_(L,S)(f,L+L^2S)=6S.                                  (14)
```

After the upper exposure, its scaled form would require the same equation
with a nonzero scalar multiple of `S`.  THM-3716 proves that no such
polynomial exists: the target begins an isolated monomial chain whose forced
coefficients never terminate.  This proves the second-stage nonentry.

## 5. Counterexample-search meaning and boundary

Every elementary decoration of `M_0` remains outside `E_2`.  Hence a
Jacobian matrix in this cell would be a planar Jacobian counterexample by
THM-3652.  The result instead gives a precise three-level anatomy:

```text
nonconstant parameters:  no first closed row,
generic constants:       no first closed row,
projective C3 resonance: one closed gradient line,
Hamiltonian cokernel:    no complementary closed row.              (15)
```

This is the first place in the automorphic-Cohn search where the repo's
recurring ternary grammar appears as an exact group relation.  It also shows
why the resonance is not itself the missing counterexample: the `C_3` move
changes the elementary word but leaves the closed component inside the
Broughton orbit.

The theorem does not classify the opposite word `E_-(u)E_+(v)`, three right
factors, or general automorphic images.  Those are the live successors.
Reproduce the exact checks with

```bash
python3 -B 04-computation/jc2_automorphic_cohn_c3_constant_resonance_thm3723.py
python3 -B -O 04-computation/jc2_automorphic_cohn_c3_constant_resonance_thm3723.py
```

The 110-cell nonconstant grid and 200-cell rational constant grid are
independent finite guards; the leading-form and exact equations above prove
the stated classification uniformly.  **QED.**
