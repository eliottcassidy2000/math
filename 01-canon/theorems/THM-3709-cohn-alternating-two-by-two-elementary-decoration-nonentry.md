---
id: THM-3709
title: "Cohn alternating two-by-two elementary-decoration nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Let C be Cohn's non-elementary determinant-one matrix.  Decorate it by
  either alternating pair of right elementary factors whose parameters are
  both nonzero, and by either alternating pair of left elementary factors
  with arbitrary polynomial parameters.  The resulting matrix is never a
  Jacobian matrix.  Positive-degree right data die by an ordinary or
  directional leading-curl obstruction.  The remaining constant resonance
  is a linear form of the Broughton polynomial X+X^2 T, whose Hamiltonian
  cokernel kills the second correction.  This closes the complete
  nonidentity alternating two-by-two Cohn repair cell, not an identity
  shortening, longer reduced words, general non-elementary classes, or JC(2).
source: root / 2026-08-22
audit: >
  PASS.  An independent audit rederived both right products, both exposed
  left rows, the arbitrary positive-degree leading forms, all four mixed
  nonzero-constant cases, and the all-constant resonance.  It checked the
  constant and positive-degree left-parameter split, the ordinary- and
  directional-derivative kernels, the Broughton normal form and homogeneous
  Hamiltonian-kernel descent, and the E_2/Wright counterexample consequence
  without a normality assumption.
depends_on:
  - THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
related:
  - THM-3655-cohn-alternating-two-source-factor-row-obstruction
script: 04-computation/jc2_cohn_alternating_two_by_two_decoration_thm3709.py
output: 05-knowledge/results/jc2_cohn_alternating_two_by_two_decoration_thm3709.out
script_sha256: 4ddf624d34b0503f33f0e29e30a60ce694f658ed6bccbcbb76deb68642fee14f
output_sha256: 5bde3ec902181f77f9bf2d30a70e348b50917c37d3e19bded155ab44ec0caf52
semantic_sha256: 3f6ba7550d4acdffcae88d4aaaa412f43266d121b76a90a77d765aee0d7c7ab8
hash_basis: raw LF bytes
---

# THM-3709 -- no alternating two-by-two decoration repairs the Cohn core

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This attacks a positive counterexample target left open by THM-3653 and
THM-3655.  Elementary decoration preserves Cohn's nontrivial
`SL_2/E_2` class.  Consequently, if a decorated matrix were the derivative
of a polynomial map, its determinant would be one and Wright's criterion
would certify that the map is not an automorphism: it would be a planar
Jacobian counterexample.  The theorem rules out the entire nondegenerate
alternating two-left/two-right cell with no identity right factor.

Work over a characteristic-zero field `k`.  Put

```text
A=1+xy,       B=x^2,       G=-y^2,       D=1-xy,

C=[A B;G D],                                  det C=1,       (1)
E_+(h)=[1 h;0 1],              E_-(h)=[1 0;h 1].       (2)
```

Let `u,w in k[x,y]` be nonzero.  First suppose both are nonconstant, and write

```text
p=deg u>=1,                q=deg w>=1,
u_p=leading homogeneous form of u,
w_q=leading homogeneous form of w.                    (3)
```

For arbitrary `f,g in k[x,y]`, neither matrix

```text
E_+(f)E_-(g) C R,             E_-(g)E_+(f) C R       (4)
```

is a Jacobian matrix when

```text
R=E_+(w)E_-(u)       or       R=E_-(u)E_+(w).         (5)
```

Thus all four alternating two-left/two-right orders are impossible whenever
both right factors are nonidentity, with no degree bound.

## 1. The row-combination gate

Write the rows of `N=CR` as `rho_1,rho_2`, regarded as polynomial one-forms.
For a row `rho=(p,q)`, use

```text
curl(rho)=p_y-q_x.                                    (6)
```

The two left words in `(4)` expose one row before the other:

```text
bottom row of E_+(f)E_-(g)N = rho_2+g rho_1,
top row of E_-(g)E_+(f)N    = rho_1+f rho_2.           (7)
```

It therefore suffices to prove

```text
curl(rho_2+h rho_1)!=0,
curl(rho_1+h rho_2)!=0                 for every h in k[x,y].       (8)
```

Notice that `(8)` is stronger than a bounded coefficient search: the left
parameters may have arbitrary degree and arbitrary lower terms.

## 2. Right order `E_+(w)E_-(u)`: a forced y-derivative

Set

```text
s=B+Aw,                 r=D+Gw.                        (9)
```

Direct multiplication gives

```text
rho_1=(A+us,s),              rho_2=(G+ur,r).          (10)
```

Their leading homogeneous pieces are

```text
s_(q+2)=xyw_q,             r_(q+2)=-y^2w_q,
(rho_1)_(p+q+2)=(u_pxyw_q,0),
(rho_2)_(p+q+2)=(-u_py^2w_q,0),                       (11)
```

where the displayed pairs retain only the top component; both second
components have degree `q+2`, lower by `p`.

Let `h!=0` have degree `m` and leading form `h_m`.  If `m>=1`, the top
homogeneous parts of the two curls in `(8)` are respectively

```text
[curl(rho_2+h rho_1)]_(m+p+q+1)
       =partial_y(u_p x y w_q h_m),

[curl(rho_1+h rho_2)]_(m+p+q+1)
       =partial_y(-u_p y^2 w_q h_m).                  (12)
```

The `q_x` half of the curl has degree at most `m+q+1`, lower by `p`, so it
cannot cancel either expression.  Each product inside `partial_y` is nonzero
and divisible by `y`.  A polynomial with zero `y`-derivative lies in `k[x]`;
a nonzero such polynomial cannot be divisible by `y`.  Both expressions in
`(12)` are therefore nonzero.

If `h=c in k`, including `c=0`, the degree-`p+q+1` curl pieces are

```text
partial_y[u_p y w_q(-y+cx)],
partial_y[u_p y w_q( x-cy)],                           (13)
```

in the same order as `(12)`.  The polynomials inside the derivatives are
again nonzero and divisible by `y`, so `(13)` cannot vanish.  This proves
`(8)` for the first right word.

## 3. Right order `E_-(u)E_+(w)`: a forced x-derivative

Now put

```text
t=A+Bu,                 z=G+Du.                        (14)
```

Then

```text
rho_1=(t,B+wt),              rho_2=(z,D+wz),          (15)
```

and the leading pieces are

```text
t_(p+2)=x^2u_p,             z_(p+2)=-xyu_p,
(rho_1)_(p+q+2)=(0,x^2u_pw_q),
(rho_2)_(p+q+2)=(0,-xyu_pw_q).                        (16)
```

Here the `-q_x` half dominates.  For `m>=1`, the leading curl pieces are

```text
[curl(rho_2+h rho_1)]_(m+p+q+1)
       =-partial_x(x^2u_pw_qh_m),

[curl(rho_1+h rho_2)]_(m+p+q+1)
       =-partial_x(-xyu_pw_qh_m).                     (17)
```

The products are nonzero and divisible by `x`, hence cannot have zero
`x`-derivative.  The `p_y` half has degree at most `m+p+1`, lower by `q`.
For `h=c in k`, including zero, the degree-`p+q+1` pieces are

```text
-partial_x[xu_pw_q(-y+cx)],
-partial_x[xu_pw_q( x-cy)].                           (18)
```

They are nonzero for the same reason.  Equations `(17)--(18)` prove `(8)`
for the second right word.

Combining `(7)--(8)` proves that every matrix in `(4)--(5)` has a nonclosed
row when both right parameters are nonconstant.

## 4. A mixed nonzero-constant right word also fails

It remains to allow exactly one of `u,w` to be constant.  The constant is
nonzero by hypothesis.  Only the leading vector field changes.

First take the right order `E_+(w)E_-(u)`.

- If `w=c in k*` and `deg u=p>=1`, the leading first components of the two
  rows are

  ```text
  u_p x(x+cy),                -u_p y(x+cy).            (19)
  ```

  Their second components have degree two.  For positive-degree `h`, the
  leading curl is the `y`-derivative of `(19)` times `h_m`; for `h=d in k`,
  it is the `y`-derivative of

  ```text
  u_p(x+cy)(-y+dx),           u_p(x+cy)(x-dy).         (20)
  ```

  None can vanish: a polynomial with zero `y`-derivative lies in `k[x]`,
  whereas a nonzero polynomial in `k[x]` cannot be divisible by `x+cy`.

- If `u=c in k*` and `deg w=q>=1`, each leading row is `(cF,F)`.  The two
  positive-degree choices of `F` are

  ```text
  h_mxyw_q,                   -h_my^2w_q,              (21)
  ```

  and the constant-`h=d` choices are `yw_q(-y+dx)` and `yw_q(x-dy)`.
  The leading curl is

  ```text
  (c partial_y-partial_x)F.                            (22)
  ```

  Its kernel is `k[y+cx]`.  No nonzero member of that kernel is divisible by
  `x` or by `y`, while every displayed `F` is.  Hence `(22)` is nonzero.

For the right order `E_-(u)E_+(w)`, exchange `x` with `y` and the two curl
halves.  If `w=c in k*`, each leading row is `(F,cF)` and the operator is
`partial_y-c partial_x`, whose kernel is `k[x+cy]`; the displayed factors
carry `x`.  If `u=c in k*`, the leading second components carry `y+cx`, and
their `x`-derivatives cannot vanish.  This proves `(8)` in all four mixed
cases.

## 5. The all-constant resonance is a Broughton obstruction

It remains to take `u,w in k*`.  Both constant right words have the form

```text
R=[a b;c d] in SL_2(k),                b=w!=0, c=u!=0.    (23)
```

Put

```text
ell_1=ay+cx,                  ell_2=by+dx.                (24)
```

The quadratic leading rows of `N=CR` are

```text
x(ell_1,ell_2),              -y(ell_1,ell_2).             (25)
```

First expose `rho_2+h rho_1`.  If `deg h=m>=1`, write
`h_m=x^m p(t)`, where `t=y/x` and `deg p<=m`.  Vanishing of its top curl
would give

```text
[bt^2+(a+d)t+c]p'
    =[(m+1)bt+(m+2)d-a]p.                                (26)
```

If `r=deg p`, comparison of the top powers of `t` in `(26)` gives
`r=m+1`, impossible because `b!=0`.  Thus a closed exposed row forces
`h=g in k`.  Its curl vanishes exactly when

```text
g(a-2d)=c,                         gb=d-2a.              (27)
```

These equations and `ad-bc=1` imply

```text
2(a-d)^2=1.                                                (28)
```

Away from `(28)`, even the exposed row cannot close.  At the resonance put

```text
s=a-d,       X=gx-y,       Y=(a-2d)x-by.                 (29)
```

Then `X,Y` are independent because `J(X,Y)=3s!=0`, and direct integration
of the exposed row gives the potential

```text
Q=2sX+(1/3)X^2Y.                                         (30)
```

Thus the only survivor is, after a linear rescaling, the Broughton
polynomial `X+X^2T`.  The complementary row `rho_1` has

```text
curl(rho_1)=Y.                                           (31)
```

If the remaining left parameter `f` closed `rho_1+f dQ`, exterior
differentiation would therefore force

```text
J(f,Q)=Y.                                                (32)
```

Equation `(32)` has no polynomial solution.  Indeed, the cubic part is
`Q_3=X^2Y/3`.  For a homogeneous polynomial `H`, expansion in monomials
shows

```text
J(H,X^2Y)=0     iff     H is a scalar power of X^2Y.     (33)
```

More explicitly, `J(X^iY^j,X^2Y)=(i-2j)X^(i+1)Y^j`, up to the fixed
nonzero coordinate determinant.  If `f_m` is the top homogeneous part of a
hypothetical solution of `(32)`, degree comparison gives
`J(f_m,Q_3)=0`.  Hence `m=3j` and `f_m=lambda Q_3^j`.  Subtracting
`lambda Q^j` leaves `(32)` unchanged and strictly lowers the degree.
Iterating reduces `f` to a constant, whose bracket with `Q` is zero,
contradicting the nonzero right side `Y`.

For the other exposed row, write `h_m=y^m p(x/y)`.  Since `c!=0`, the same
top-coefficient comparison again forces `deg p=m+1`.  A constant `h=f`
can close only when

```text
fc=a-2d,                         f(d-2a)=b,              (34)
```

which gives the same resonance `(28)`.  Now take

```text
X=x-fy,                 Y=cx+(2a-d)y.                  (35)
```

The exposed potential is again `(30)`, while the complementary row has curl
`-Y`.  Its final correction would require `J(g,Q)=-Y`, excluded by the same
descent.  This proves `(8)` for every all-constant nonzero right word.

## 6. Counterexample-search meaning and boundary

Cohn's matrix is determinant one and non-elementary.  Multiplication by
elementary matrices preserves both facts.  If all rows of one of the
decorations above were closed, integration over characteristic zero would
produce a Keller map whose derivative remains outside `E_2`; by THM-3652 it
would be a nonautomorphism and hence a counterexample.  The proof identifies
the precise failed implication: two alternating right corrections cannot
cancel the coordinate factor in the top curl, no matter how high either
side's degree is.

The theorem does **not** treat a zero right parameter (an identity shortening),
more than two elementary factors on either side, another non-elementary
`SL_2/E_2` core, or `JC(2)`.  Its positive successor is therefore sharply
typed: settle the shortened boundary, or make the reduced word longer.

Reproduce the exact leading-form audit with

```bash
python3 -B 04-computation/jc2_cohn_alternating_two_by_two_decoration_thm3709.py
python3 -B -O 04-computation/jc2_cohn_alternating_two_by_two_decoration_thm3709.py
```

The symbolic companion reconstructs all four products before extracting the
leading pieces, checks `(11)--(22)` on hostile grids of independent right and
left degrees, and verifies `(23)--(35)` together with the Hamiltonian descent
on a degree grid.  **QED.**
