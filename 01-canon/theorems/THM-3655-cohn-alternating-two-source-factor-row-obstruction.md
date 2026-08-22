---
id: THM-3655
title: "Cohn alternating two-source-factor row obstruction"
status: >
  PROVED + SYMBOLIC-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over every
  characteristic-zero field, neither row of Cohn's determinant-one matrix
  becomes a closed polynomial one-form after either alternating pair of
  right elementary factors.  Consequently no additional single left
  elementary factor can make the matrix a Jacobian matrix.  No statement is
  made about longer words, other non-elementary classes, Keller pairs, or
  JC(2).
source: kps-s191 / THM-3653 non-elementary-core continuation, 2026-08-21
audit: >
  PASS -- an independent hostile audit rederived all four row formulas and
  supplied characteristic-zero leading-coefficient proofs in both factor
  orders, including every zero-polynomial and degree-zero edge case.  It
  verified that one left elementary factor leaves a proved nonclosed row
  untouched and found no algebraic-closure assumption or scope leak.
depends_on: []
related:
  - THM-3653-cohn-factorial-repair-and-weighted-rectangle-holonomy
  - THM-3652-wright-elementary-jacobian-criterion-reduced-word-reproof
script: 04-computation/jc2_cohn_alternating_two_source_factor_obstruction_thm3655.py
output: 05-knowledge/results/jc2_cohn_alternating_two_source_factor_obstruction_thm3655.out
script_sha256: 46950a9492b2991a19c2740d776fdf2d9e0ee306ae322901c816d7f8e084e008
output_sha256: 50144414d166e710d98f5cc36d7fb7b4cd7e76e863b8357945c80cd9649f0517
semantic_sha256: d1f05c8841156820f7ab6573e16315995dfbd197dc2e02e403f27179dc7298e7
hash_basis: raw LF bytes
---

# THM-3655 -- Cohn alternating two-source-factor row obstruction

**PROVED + SYMBOLIC-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This closes the
first two-factor source-side escape left open by THM-3653.  Right elementary
factors are treated only as polynomial matrix factors; no claim is made that
they are themselves derivatives of source automorphisms.

Let `k` be a characteristic-zero field and put

```text
A=1+xy,       B=x^2,       G=-y^2,       D=1-xy,
C=[A B;G D],                                      det C=1,       (1)
E_+(h)=[1 h;0 1],             E_-(h)=[1 0;h 1].       (2)
```

A row `(P,Q)` is closed exactly when `P_y=Q_x`.

## 1. A logarithmic-derivative gate

We repeatedly use the elementary fact

```text
0!=H in k[z], H'=qH with q in k[z]
       implies q=0 and H is constant.                 (3)
```

Indeed `deg H'<deg H`, whereas a nonzero `qH` has degree at least `deg H`.
In particular, a nonzero polynomial divisible by `z` cannot have zero
derivative or satisfy `(3)`.  Characteristic zero is load-bearing here.

## 2. The order `E_+(w)E_-(u)`

Direct multiplication gives

```text
C E_+(w)E_-(u)=
 [A+u s   s]
 [G+u r   r],
s=B+Aw=x^2+(1+xy)w,
r=D+Gw=1-xy-y^2w.                                    (4)
```

### Top row

Closure would require

```text
(A+us)_y=s_x.                                         (5)
```

First suppose `w!=0`.  Let `m=deg_y w` and let `b(x)` be its leading
coefficient.  Then `deg_y s=m+1` and the leading coefficient of `s` is the
nonzero polynomial `xb(x)`, which is divisible by `x`.  If `u!=0`, write
`n=deg_y u` and let `a(x)` be its leading coefficient.

- If `n>=2`, then `(us)_y` has degree `n+m>m+1`, larger than every other
  term in `(5)`.
- If `n=1`, the coefficient of `y^(m+1)` in `(5)` is

  ```text
  (xb)'=(m+2)a(xb),                                   (6)
  ```

  contradicting `(3)`.
- If `n=0`, or if `u=0`, the same coefficient gives `(xb)'=0`, also
  impossible by `(3)`.

If `w=0`, equation `(5)` reduces to

```text
x^2 u_y=x,                     equivalently x u_y=1,  (7)
```

which has no polynomial solution.  Thus the top row in `(4)` is never
closed.

### Bottom row

Closure of the bottom row in `(4)` is

```text
(ur)_y-r_x=2y.                                         (8)
```

As a polynomial in `x`, `r` has degree `M>=1`; its leading coefficient
`R_M(y)` is nonzero and divisible by `y`.  Explicitly it is `-y`,
`-y(1+yb)`, or `-y^2b` according as `deg_x w` is at most zero, equal to one,
or at least two.

If `u!=0`, let `a(y)` be its leading coefficient in `x`.  The highest
`x`-coefficient in `(8)` forces

```text
(a R_M)'=0.                                            (9)
```

But `aR_M` is nonzero and divisible by `y`, contradicting `(3)`.  If `u=0`,
then modulo `y^2`, equation `(8)` reads

```text
y=2y mod y^2,                                         (10)
```

again impossible.  Hence neither row in `(4)` is closed.

## 3. The order `E_-(u)E_+(w)`

Now

```text
C E_-(u)E_+(w)=
 [t   B+w t]
 [z   D+w z],
t=A+Bu=1+xy+x^2u,
z=G+Du=-y^2+(1-xy)u.                                  (11)
```

### Top row

Closure would give

```text
t_y=2x+(wt)_x.                                        (12)
```

If `w!=0`, view the equation in `y`.  The leading coefficient `T(x)` of
`t` is always nonzero and divisible by `x`; in the only potentially delicate
linear case it is `x(1+xa_1)`, which cannot vanish for polynomial `a_1`.
The highest `y`-coefficient in `(12)` therefore gives

```text
(bT)'=0,                                               (13)
```

where `b(x)` is the leading `y`-coefficient of `w`.  This contradicts `(3)`.
If `w=0`, equation `(12)` is once more `xu_y=1`.  The top row is never
closed.

### Bottom row

Closure is

```text
z_y+y-(wz)_x=0.                                       (14)
```

If `u=0`, equation `(14)` becomes `y^2w_x=y`, impossible in `k[x,y]`.
Suppose `u!=0`, put `n=deg_xu`, and let `a(y)` be its leading coefficient.
Then `N=deg_xz=n+1` and the leading coefficient

```text
Z(y)=-y a(y)                                          (15)
```

is nonzero and divisible by `y`.  Let `m=deg_xw` when `w!=0`.

- If `m>=2`, `(wz)_x` has degree `m+N-1>N`, impossible in `(14)`.
- If `m=1`, the coefficient of `x^N` gives

  ```text
  Z'=(N+1)bZ,                                         (16)
  ```

  contradicting `(3)`.
- If `m=0`, or if `w=0`, the same coefficient gives `Z'=0`, also
  contradicting `(3)`.

Thus neither row in `(11)` is closed.

## 4. One additional left factor cannot help

Left multiplication by `E_+(f)` changes only the top row and leaves the
bottom row untouched.  Left multiplication by `E_-(f)` changes only the
bottom row and leaves the top row untouched.  Sections 2--3 prove that both
rows of `CR` are nonclosed for either alternating word `R`.  Therefore

```text
L C R is not a Jacobian matrix                           (17)
```

for `L in {I,E_+(f),E_-(f)}` and

```text
R in {E_+(w)E_-(u), E_-(u)E_+(w)}.                    (18)
```

## 5. Verification and boundary

Reproduce the symbolic identities and hostile controls with

```bash
python3 04-computation/jc2_cohn_alternating_two_source_factor_obstruction_thm3655.py
python3 -O 04-computation/jc2_cohn_alternating_two_source_factor_obstruction_thm3655.py
```

The assertion-free companion checks both exact matrix products, determinant
one, all four curl reductions, the boundary jets used above, and 144 exact
hostile pairs in each order.  Normal and optimized streams match the stored
transcript.  The finite bank guards signs and order; the universal proof is
the degree argument in Sections 1--4.

Same-sign adjacent right factors merge and return to THM-3653's one-factor
regime.  This theorem does not classify two factors on the left, words with
two factors on both sides, longer alternating words, other non-elementary
`SL_2/E_2` classes, row-integrable matrices in general, Keller pairs,
collisions, nonproperness, counterexamples, or JC(2).  **QED.**
