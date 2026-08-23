---
id: THM-3725
title: "Automorphic Cohn opposite two-right hyperbolic resonance nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For arbitrary
  nonzero polynomial parameters
  u,v, neither alternating two-left decoration of M0 E_-(u)E_+(v) is a
  Jacobian matrix.  If either right parameter is nonconstant, leading-form
  curl kernels prevent the first row closure.  In the constant cell a row
  closes precisely at uv=2; the right matrix then has trace four and Pell
  law R^2-4R+I=0, while both closed rows are scalar forms of the Broughton
  polynomial -(L+L^2S).  THM-3716 prevents the second closure.  This closes
  the opposite nonidentity two-right order, not longer words or JC(2).
source: jc-sparse-direct-search / 2026-08-22
audit: >
  PASS.  An independent audit rederived the exact row decomposition, all
  three nonconstant leading-kernel regimes, both positive-degree
  left-parameter valuation descents, the constant resonance, trace-four
  Pell law, shared Broughton potential, and both debt signs.  The exact companion
  checks the arbitrary-polynomial exposed curls, both constant
  classifications, the trace-four/Pell matrix law and diagonal gauge, the
  Broughton source coordinates and both remaining debts, 110 nonconstant
  and 200 constant hostile cells, resonant homogeneous kernels through
  degree six, and finite Hamiltonian cokernels through degree 12.  Normal
  and optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
  - THM-3716-monomial-broughton-hamiltonian-obstruction-family
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
related:
  - THM-3723-automorphic-cohn-c3-two-right-resonance
script: 04-computation/jc2_automorphic_cohn_opposite_two_right_hyperbolic_thm3725.py
output: 05-knowledge/results/jc2_automorphic_cohn_opposite_two_right_hyperbolic_thm3725.out
script_sha256: fd9688113ad5328614c503687d7da3c43658e6f800e8c07efcd7c7b69d1a65ef
output_sha256: aa3221da084925a8ff80bd5044c124f463e617f1e5324162c5195d6efab2a367
semantic_sha256: d45eb546df8a1a3a2e5195375ec7ae32e75868d3f9f31828a38a5add41b2050b
hash_basis: raw LF bytes
---

# THM-3725 -- the opposite two-right word has one hyperbolic resonance

**PROVED + VERIFIED-EXACT.**  Work over a characteristic-zero field `k`.
Recall the automorphic Cohn core

```text
M_0=[4T^2,2XT-1; 1+2XT,X^2],                 det M_0=1.       (1)
```

THM-3723 closes the right word `E_+(v)E_-(u)`, whose exceptional constant
cell is projective `C_3`.  The opposite word has a different and initially
less visible resonance.  For nonzero `u,v in k[X,T]`, put

```text
R=E_-(u)E_+(v)=[1,v;u,1+uv],
N=M_0R,                                           det N=1.    (2)
```

For arbitrary `f,g,h in k[X,T]`, neither

```text
E_+(f)E_-(h)N,                    E_-(g)E_+(h)N          (3)
```

is a Jacobian matrix.  If either `u` or `v` is nonconstant, neither exposed
row closes.  If both are constants, writing the rows of `N` as `alpha,beta`,
the lower exposed row closes exactly when

```text
uv=2,                         h=v,                         (4)
```

and the upper exposed row closes exactly when

```text
uv=2,                         h=1/v.                       (5)
```

These two closures are the same gradient line.  No final polynomial left
shear closes the complementary row.

## 1. Structural row form

Put

```text
a=4T^2,       b=2XT-1,       c=1+2XT,       d=X^2,
A=a+bu,                                      B=c+du.       (6)
```

Multiplication in `(2)` gives the useful exact decomposition

```text
alpha=(A,b+vA),                         beta=(B,d+vB).       (7)
```

This makes the two asymmetrical degree kernels below visible.  All degrees
are total degrees, and a subscript denotes a leading homogeneous form.

## 2. Nonconstant right data cannot expose a gradient

First let `q=deg u>=1` and `p=deg v>=1`.  The degree-`p+q+2` parts of the
second components of `alpha,beta` are

```text
2XT u_qv_p,                              X^2u_qv_p.         (8)
```

If `deg h=m>=1`, the uniquely highest second component of the lower or upper
exposure is respectively

```text
2XT u_qv_ph_m,                           X^2u_qv_ph_m.      (9)
```

Its negative `X` derivative is the leading curl; the `T` derivative of the
first component is at least `p` degrees lower.  Neither derivative in `(9)`
vanishes: its source is nonzero and divisible by `X`, whereas a polynomial
with zero `X` derivative lies in `k[T]`.  If `h` is constant, the two
highest second components instead are

```text
Xu_qv_p(X+2hT),                          Xu_qv_p(2T+hX),    (10)
```

and the same argument applies.

Next let `u` be nonconstant and `v=c_0 in k*`.  At top degree the two rows
are scalar one-forms

```text
2XT u_q(1,c_0),                          X^2u_q(1,c_0).     (11)
```

The leading curl of `F(1,c_0)` is

```text
D_(c_0)F=(partial_T-c_0 partial_X)F.                    (12)
```

Its polynomial kernel is `k[X+c_0T]`.  For positive-degree `h`, the
relevant `F` is `2XT u_qh_m` or `X^2u_qh_m`; for constant `h`, it is

```text
Xu_q(X+2hT),                              Xu_q(2T+hX).     (13)
```

Every displayed nonzero form is divisible by `X`, while no nonzero member
of `k[X+c_0T]` is divisible by `X`.  Thus no leading curl vanishes.

Finally let `u=c_0 in k*` and `v` be nonconstant.  The quadratic parts

```text
A_2=2T(2T+c_0X),                         B_2=X(2T+c_0X)   (14)
```

show that the uniquely highest second component of an exposure is, for
positive-degree `h`,

```text
2T(2T+c_0X)v_ph_m,                       X(2T+c_0X)v_ph_m,
```

and, for constant `h`,

```text
(2T+c_0X)v_p(X+2hT),        (2T+c_0X)v_p(2T+hX).       (15)
```

Closure would make the relevant expression lie in `k[T]`.  But every one
is divisible by the linear prime `2T+c_0X`, which is not associated to `T`.
This is impossible.  Equations `(8)--(15)` exhaust all cells in which at
least one of the two nonzero right parameters is nonconstant.

## 3. Positive-degree left data do not survive a constant right cell

Suppose now `u,v in k*`.  The quadratic leading rows of `N` are

```text
alpha_2=2T w,                    beta_2=Xw,
w=(2T+uX, 2vT+(1+uv)X).                               (16)
```

If the leading form `h_m` has positive degree, closure of the lower or upper
exposed row would respectively make

```text
F=2Th_m                       or                       F=Xh_m
```

satisfy

```text
DF+(1-uv)F=0,
D=(2T+uX)partial_T-[2vT+(1+uv)X]partial_X.             (17)
```

For the first choice, take the least `T`-valuation of `F`.  The term
`uX partial_T F` is the unique term one valuation lower, since `u!=0`; it
cannot cancel.  For the second choice, take the least `X`-valuation.  The
term `-2vT partial_X F` is uniquely one valuation lower.  Hence any closing
`h` in a constant right cell must itself be constant.

## 4. Exact constant classification

For constant `u,v,h`, direct differentiation gives

```text
curl(beta+h alpha)
 =2{u(h-v)X+[-huv+3h-v]T},                            (18)

curl(alpha+h beta)
 =2{u(1-hv)X+[3-hv-uv]T}.                             (19)
```

Since `u,v` are nonzero, the `X` coefficient in `(18)` forces `h=v`; the
`T` coefficient then becomes `2v(2-uv)`.  Similarly `(19)` forces `h=1/v`
and then `uv=2`.  This proves `(4)--(5)` uniformly.

If `u=0` or `v=0`, the word shortens to a single right shear; THM-3721
closes both such boundaries for arbitrary polynomial surviving parameter.

## 5. The resonance is hyperbolic/Pell, not C3

At `u=2/v`, the right matrix is

```text
R_v=[1,v;2/v,3].                                      (20)
```

It satisfies

```text
tr(R_v)=4,                       R_v^2-4R_v+I=0.       (21)
```

Moreover, with `D_v=diag(v,1)`,

```text
R_v=D_v[1,1;2,3]D_v^(-1).                             (22)
```

Thus over the reals its eigenvalues are `2+sqrt(3)` and `2-sqrt(3)`; the
integer base in `(22)` generates the corresponding Pell recurrence.  This
is a hyperbolic trace-four orbit, sharply distinct from the projective
order-three resonance of THM-3723.

## 6. The resonant gradient is again Broughton

Introduce affine source coordinates

```text
L=X+2vT,                       S=-(2X+vT)/(3v).        (23)
```

They have `J_(X,T)(L,S)=1`.  Put

```text
Q=-(L+L^2S).                                          (24)
```

At the resonance, exact differentiation gives

```text
dQ=beta+v alpha,
alpha+(1/v)beta=(1/v)dQ.                              (25)
```

Thus `(4)` and `(5)` expose the same gradient line.  In the same source
coordinates the complementary curls are

```text
curl(alpha)=-6S,                         curl(beta)=6vS. (26)
```

After the lower exposure, closure of `alpha+f dQ` would require

```text
J_(L,S)(f,Q)=-6S.                                     (27)
```

After the upper exposure, closure of `beta+(g/v)dQ` would require

```text
J_(L,S)(g,Q)=6v^2S.                                   (28)
```

Since `Q=-(L+L^2S)`, either equation would give a polynomial solution of
`J(F,L+L^2S)=lambda S` with `lambda!=0`.  Scaling `F` would then solve the
normalized `6S` equation excluded by THM-3716.  Therefore the second row
never closes.

## 7. Counterexample-search meaning and scope

The automorphic Cohn core is outside `E_2(k[X,T])`, and elementary
decorations preserve its nonentry.  By THM-3652, a Jacobian matrix found in
this cell would therefore yield a planar Jacobian counterexample.  The
complete anatomy is instead

```text
nonconstant right data: no first closed row,
generic constants:      no first closed row,
trace-four resonance:   one Broughton gradient line,
Hamiltonian cokernel:   no complementary closed row.              (29)
```

Together THM-3723 and the present theorem close both alternating
nonidentity two-right orders around `M_0`.  Identity-factor boundaries are
already THM-3721.  Three or more right factors, other non-elementary cores,
general automorphic images, and `JC(2)` remain open.

Reproduce the exact checks with

```bash
python3 -B 04-computation/jc2_automorphic_cohn_opposite_two_right_hyperbolic_thm3725.py
python3 -B -O 04-computation/jc2_automorphic_cohn_opposite_two_right_hyperbolic_thm3725.py
```

The 110-cell nonconstant grid and 200-cell rational constant grid are
independent finite guards; Sections 2--6 prove the all-degree statement.
**QED.**
