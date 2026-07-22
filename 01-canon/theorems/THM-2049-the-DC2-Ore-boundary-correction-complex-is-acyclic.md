---
id: THM-2049
title: "The DC(2) Ore boundary correction complex is acyclic"
status: >
  PROVED. The beta filtration beta(sum a_k ell^k)=min(v_x(a_k)-2k)
  is multiplicative on the HYP-8803 Ore algebra, its associated-graded
  commutator is the scaled canonical bracket {ell,q}_0=2, and the linearized
  simultaneous (S,T) correction map is surjective in every relevant grade.
  Hence the grade-six anomaly is not a cohomology obstruction. An exact
  recursive lift exists in the beta-completion; polynomial termination remains
  open, and the currently chosen first eight corrections advance grades 6
  through 13 without terminating.
source: codex-2026-07-21-DC2-LRC14-termination
related:
  - THM-2044
  - THM-2046
  - HYP-8802
  - HYP-8803
  - HYP-8841
script: 04-computation/dc2_ore_descent_codex_20260721.py
output: 05-knowledge/results/dc2_ore_descent_codex_20260721.out
---

# THM-2049 -- acyclicity of the boundary correction complex

In the Ore algebra

```text
O=Q[x,q][ell;delta],
delta=3x^2 partial_x+(2-6xq) partial_q,              (1)
```

write every element uniquely in coefficient-left normal form and define

```text
beta(sum_k a_k(x,q) ell^k)=min_k(v_x(a_k)-2k).       (2)
```

Here `v_x` is the ordinary `x`-adic valuation of a polynomial in `Q[x,q]`.

## 1. The filtration and its bracket

The derivation in (1) does not lower `v_x`, since

```text
delta=2 partial_q+x(3x partial_x-6q partial_q).       (3)
```

Ore multiplication is

```text
(a ell^i)(b ell^j)
 =sum_(r=0)^i binom(i,r) a delta^r(b) ell^(i-r+j).  (4)
```

The term indexed by `r` has beta at least

```text
beta(a ell^i)+beta(b ell^j)+2r.                     (5)
```

Consequently `beta(AB)>=beta(A)+beta(B)`.  The `r=0` terms cancel in a
commutator, so

```text
beta([A,B])>=beta(A)+beta(B)+2.                     (6)
```

The degree-two associated-graded bracket is obtained from the `x`-initial
part `delta_0=2 partial_q`:

```text
{A,B}_0=2(partial_ell A partial_q B
           -partial_q A partial_ell B),             (7)
```

with `x` central.  Thus `{ell,q}_0=2`.

## 2. Initial canonical pair

Let `T_W,S_W` be the Weyl-ordered Ore lifts from HYP-8803, and put

```text
u=x^2 ell.                                           (8)
```

Their lowest beta symbols are

```text
tau=in_-1(T_W)=(2/3)x ell(4-u),
sigma=in_-2(S_W)=ell(1/2-(5/18)u+(1/27)u^2).        (9)
```

Differentiating with respect to `ell` gives

```text
partial_ell tau=(4/3)x(2-u),
partial_ell sigma=(2u^2-10u+9)/18.                 (10)
```

## 3. Surjectivity in every grade

Suppose a residual has homogeneous beta-`g` symbol

```text
x^g H(q,u),                 g>=1.                   (11)
```

Use corrections

```text
C_S=x^(g-1) Atilde(q,u),
C_T=x^g Btilde(q,u).                                 (12)
```

Writing `A=partial_q Atilde`, `B=partial_q Btilde`, the linearized grade-`g`
change is

```text
d_g(C_S,C_T)/x^g
 =(8/3)(u-2)A+(2u^2-10u+9)B/9.                    (13)
```

The two coefficient polynomials are coprime, because

```text
(2u^2-10u+9)|_(u=2)=-3.                             (14)
```

More explicitly, for a desired right side `K(q,u)`, choose

```text
B=-3K(q,2),
A=(3/8)[K-(2u^2-10u+9)B/9]/(u-2).                 (15)
```

The numerator in (15) vanishes at `u=2`, so `A` is polynomial.  Since
`partial_q:Q[q,u]->Q[q,u]` is surjective in characteristic zero, polynomial
antiderivatives `Atilde,Btilde` exist.  Therefore every `d_g` is surjective.

## 4. The affine syzygy gauge

The solution chosen in (15) is one section of a larger affine family.  Put

```text
D(u)=2u^2-10u+9.                                      (16)
```

If `(A_0,B_0)` is any solution of the grade equation

```text
(8/3)(u-2)A+D(u)B/9=H,                                (17)
```

then every solution, and no others, is

```text
A=A_0+D(u)C,
B=B_0-24(u-2)C,             C in Q[q,u].              (18)
```

The displayed change is in the kernel by direct substitution.  Conversely,
after multiplying the homogeneous equation by nine, coprimality of `D(u)`
and `u-2` forces `D(u)` to divide `A-A_0`, giving (18).  Thus the exact
descent script's nonzero grade-fourteen residual after correcting grades six
through thirteen concerns only its chosen section; it is not invariant under
the syzygy gauge `C`.
Moreover, after integrating `A=partial_q Atilde` and
`B=partial_q Btilde`, the two independent additions

```text
Atilde |-> Atilde+a(u),       Btilde |-> Btilde+b(u)          (19)
```

do not change the current linearized grade but can change later exact grades.
A termination search should retain `C_g`, `a_g`, and `b_g` at every grade and
optimize support or future residual size.

## 5. Consequences

HYP-8803's Weyl residual begins in grade six with

```text
2x^6(u-1)(u-2).                                    (20)
```

Equation (13) kills it with the particularly simple fixed-`T` correction

```text
C_S=-(3/4)q x^5(u-1).                              (21)
```

The exact residual then begins in grade seven.  Iterating (15), while now
allowing both `S` and `T` to move, advances the exact residual through

```text
6 -> 7 -> 8 -> 9 -> 10 -> 11 -> 12 -> 13 -> 14.    (22)
```

The script verifies every cancellation symbolically.  In the beta-completion,
the same recursion constructs a formal pair satisfying `[S,T]=1`: nonlinear
products of the new corrections lie in strictly later grades and are handled
at subsequent steps.

This proves a negative and a positive result.  Negative: grade six is not the
DC(2) obstruction, and no finite associated-graded class can be.  Positive: a
formal simultaneous quantum descent exists for the `(S,T)` relation.  What is
not proved is termination in the polynomial Weyl algebra, nor compatibility
with the remaining `D` relations.  Those are now the only honest gates for
this route.
