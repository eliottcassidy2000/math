---
id: THM-3709
title: "Cohn alternating two-by-two elementary-decoration nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Let C be Cohn's non-elementary determinant-one matrix.  Decorate it by
  either alternating pair of right elementary factors whose two parameters
  are nonconstant polynomials, and by either alternating pair of left
  elementary factors with arbitrary polynomial parameters.  The resulting matrix is
  never a Jacobian matrix.  The highest curl is a nonzero x- or y-derivative
  of a polynomial carrying a forced coordinate factor.  This closes the
  complete positive-degree alternating two-by-two Cohn repair cell, not a
  constant right boundary, longer reduced words, general non-elementary
  classes, or JC(2).
source: root / 2026-08-22
audit: >
  PASS.  An independent audit rederived both right products, both exposed
  left rows, and the arbitrary-degree leading forms.  It checked the constant
  and positive-degree left-parameter split, proved that the opposite curl
  component is strictly lower degree, verified the forced-coordinate
  derivative cannot vanish, and checked the E_2/Wright counterexample
  consequence without a normality assumption.
depends_on:
  - THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
related:
  - THM-3655-cohn-alternating-two-source-factor-row-obstruction
script: 04-computation/jc2_cohn_alternating_two_by_two_decoration_thm3709.py
output: 05-knowledge/results/jc2_cohn_alternating_two_by_two_decoration_thm3709.out
script_sha256: b282ce1922759f6eeb4121e5f0e92ee58f9e7f0d567e10516b00f822db208010
output_sha256: 2ea501226f9822fa67882dbbfe61c0771319502695ce92acb3fee6c9b33bc292
semantic_sha256: 3487770566278b4fc17a8f95c14a80198255f747d14c0e475a9d2545a1ddda1e
hash_basis: raw LF bytes
---

# THM-3709 -- no alternating two-by-two decoration repairs the Cohn core

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This attacks a positive counterexample target left open by THM-3653 and
THM-3655.  Elementary decoration preserves Cohn's nontrivial
`SL_2/E_2` class.  Consequently, if a decorated matrix were the derivative
of a polynomial map, its determinant would be one and Wright's criterion
would certify that the map is not an automorphism: it would be a planar
Jacobian counterexample.  The theorem rules out the entire positive-degree
alternating two-left/two-right cell of that construction grammar.

Work over a characteristic-zero field `k`.  Put

```text
A=1+xy,       B=x^2,       G=-y^2,       D=1-xy,

C=[A B;G D],                                  det C=1,       (1)
E_+(h)=[1 h;0 1],              E_-(h)=[1 0;h 1].       (2)
```

Let `u,w in k[x,y]` be nonconstant, and write

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
both right parameters are nonconstant, with no degree bound.

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
row, so none is a Jacobian matrix.

## 4. Counterexample-search meaning and boundary

Cohn's matrix is determinant one and non-elementary.  Multiplication by
elementary matrices preserves both facts.  If all rows of one of the
decorations above were closed, integration over characteristic zero would
produce a Keller map whose derivative remains outside `E_2`; by THM-3652 it
would be a nonautomorphism and hence a counterexample.  The proof identifies
the precise failed implication: two alternating right corrections cannot
cancel the coordinate factor in the top curl, no matter how high either
side's degree is.

The theorem does **not** treat a constant right parameter, more than two
elementary factors on either side, another non-elementary `SL_2/E_2` core,
or `JC(2)`.  Its positive successor is therefore sharply typed: shorten a
constant boundary by separate identities, or make the reduced word longer.

Reproduce the exact leading-form audit with

```bash
python3 -B 04-computation/jc2_cohn_alternating_two_by_two_decoration_thm3709.py
python3 -B -O 04-computation/jc2_cohn_alternating_two_by_two_decoration_thm3709.py
```

The symbolic companion reconstructs all four products before extracting the
leading pieces, and checks `(11)--(18)` on a hostile grid of independent
right and left degrees.  **QED.**
