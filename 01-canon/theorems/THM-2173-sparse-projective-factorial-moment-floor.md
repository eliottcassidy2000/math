---
id: THM-2173
title: "Sparse projective factorial-moment floor"
status: >
  PROVED. In every prescribed t-dimensional polynomial support envelope,
  there is a nonzero member whose first t-1 factorial moments vanish.
  For envelopes on exactly two or three monomials, every such projective
  witness has full support. If the exponents are positive, the two-charge
  lift W+Z h(ZW) has its first 2t-1 Gaussian moments zero and a later even
  moment nonzero. Thus fixed charge span two has detection depth at least
  2t on every t-slot radial envelope; for t=2,3 the witness has exactly
  t+1 monomials. This is a lower-bound theorem, not the Strong Factorial
  Conjecture's matching upper cutoff.
source: codex-2026-07-24-sparse-projective-moment-floor
depends_on:
  - THM-1510
related:
  - THM-1790
  - THM-2022
  - HYP-8765
external: >
  Eric Edo and Arno van den Essen, The Strong Factorial Conjecture,
  arXiv:1304.3956v2, Definition 2.7 and Conjecture 2.8.
---

# THM-2173 -- sparse projective factorial-moment floor

Let

```text
L:C[s]->C,                 L(s^n)=n!.                  (1)
```

THM-1790 used the full degree envelope `C[s]_(<=d)`. The mechanism only
uses dimension, so it survives on arbitrary sparse or lacunary envelopes.

## 1. The support-envelope theorem

Let `V` be any `t`-dimensional complex linear subspace of `C[s]`, with
`t>=2`. Then there is a nonzero `H in V` such that

```text
L(H^j)=0,                  1<=j<=t-1.                 (2)
```

In particular, for any distinct nonnegative exponents

```text
E={e_1,...,e_t},
V_E=span_C{s^(e_1),...,s^(e_t)},                      (3)
```

some nonzero polynomial in the prescribed support envelope `V_E` satisfies
(2). The theorem does not yet say that all `t` allowed coefficients are
nonzero.

### Proof

Choose coordinates `x_1,...,x_t` on `V` and put

```text
F_j(H)=L(H^j),              1<=j<=t-1.                (4)
```

Each `F_j` is a homogeneous polynomial of degree `j`. The ideal

```text
I=(F_1,...,F_(t-1)) subset C[x_1,...,x_t]             (5)
```

has height at most `t-1` by Krull's height theorem. If its affine zero set
were only the origin, then

```text
rad(I)=(x_1,...,x_t),                                 (6)
```

whose height is `t`, a contradiction. Therefore (4) has a nonzero common
zero. QED.

This is a dimension obstruction, not a factorial peculiarity: the
factorial functional enters only when interpreting the witness and invoking
eventual nonvanishing.

## 2. Exact support for two and three slots

For `t=2`, take `E={a,b}`, `a<b`. The equation `L(H)=0` has the explicit
projective solution

```text
H=b!*s^a-a!*s^b.                                     (7)
```

Both coefficients are nonzero. Moreover `H` is a nonzero real polynomial,
so

```text
L(H^2)=integral_0^infinity H(s)^2 exp(-s) ds>0.       (8)
```

Thus the two-slot lower bound is exact: the first factorial moment can
vanish, but the second does not for (7).

For `t=3`, let `E={a,b,c}` with `a<b<c`. Section 1 gives a nonzero common
zero of

```text
L(H)=L(H^2)=0                                         (9)
```

in `P(V_E)`. Such a zero cannot lie on a coordinate boundary. A one-term
polynomial has nonzero first moment. If a two-term polynomial

```text
H=A*s^u+B*s^v,                u<v,                    (10)
```

has `L(H)=0`, then

```text
B/A=-u!/v! in R_<0.                                  (11)
```

Hence `H=A R` for a nonzero real polynomial `R`, and

```text
L(H^2)=A^2 integral_0^infinity R(s)^2 exp(-s) ds!=0.  (12)
```

This contradicts (9). Therefore every witness supplied by (9) uses all
three prescribed monomials.

We do not assert that its third factorial moment is nonzero. That assertion
is already the `n=1`, three-term upper side of the Strong Factorial
Conjecture.

## 3. Exact two-charge lift

Assume now that every exponent in `E` is positive. Write

```text
H=s h(s)
```

and define, with `s=ZW`,

```text
P=W+Z h(s).                                           (13)
```

The two summands have charges `-1` and `+1`. Every odd moment vanishes by
charge. In an even moment only the equal-choice term is balanced, giving

```text
E[P^(2j)]=binom(2j,j)L(H^j).                          (14)
```

Combining (2) and (14),

```text
E[P^m]=0,                  1<=m<=2t-1.                (15)
```

THM-1510's eventual moment property applied to the nonzero `H` gives some
`j>=t` with `L(H^j)!=0`; then (14) gives a later nonzero even Gaussian
moment. Thus `P` is a genuine two-sided near-nullcone witness, not a
one-sided polynomial.

Consequently, on every prescribed `t`-slot radial envelope

```text
W+Z*span{s^(e_1-1),...,s^(e_t-1)},                   (16)
```

the Gaussian moment detection depth is at least

```text
2t.                                                    (17)
```

For `t=2` and `t=3`, Section 2 makes the lifted witness have exactly `t+1`
monomials. At `t=2`, (8) and (14) show that moment `4` is nonzero, so the
depth is exactly `4` for the explicit projective witness.

## 4. Strong Factorial boundary

For a polynomial with exactly `t` nonzero monomials, the `n=1` case of the
Strong Factorial Conjecture asks that at least one of

```text
L(H),L(H^2),...,L(H^t)                               (18)
```

be nonzero. The present theorem shows:

1. on every `t`-dimensional support envelope, `t-1` homogeneous moment
   equations necessarily have a nonzero common projective zero;
2. for `t=2,3`, that zero can be forced to have exactly `t` terms;
3. for larger `t`, excluding all coordinate-boundary solutions is precisely
   where a support-sensitive upper theorem is needed.

Thus the proposed HYP-8765 cutoff `(k-1)R` is already sharp from below on
the exact two- and three-slot radial families: with `k=t+1` monomials and
primitive return `R=2`, the witnesses in (13) vanish through moment

```text
(k-1)R-1=2t-1.                                        (19)
```

Neither (18) nor the general HYP-8765 upper cutoff is proved here.

## 5. Transfer boundary

The projective argument is over `C` and finds a point in a linear support
envelope. It does not automatically find a point in the coefficient torus
or over `Q`:

```text
x*y=0
```

has only coordinate-boundary zeros in `P^1`, while

```text
x^2+y^2=0
```

has complex projective zeros but no rational one. Section 2 supplies the
extra positivity needed to exclude the boundary only for two and three
factorial slots. Any transfer to an exact-support Jacobian locus, an
integer LRC carrier, a CRT residue class, or a sign-constrained tournament
therefore needs an additional torus, rational-point, or sign sidecar.

QED.
