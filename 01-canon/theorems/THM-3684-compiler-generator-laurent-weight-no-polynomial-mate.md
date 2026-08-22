---
id: THM-3684
title: "Compiler-generator Laurent-weight no-polynomial-mate theorem"
status: >
  PROVED + VERIFIED-EXACT.  For the THM-3561 compiler, none of the three
  source polynomials b, c, e has a nonzero-constant Jacobian mate even in the
  full source ring C[x,q].  More generally, no nonconstant polynomial in one
  of b,c,e has such a mate.  Consequently a Russell-cylinder construction
  cannot use a one-generator surface output and repair only the other output
  by arbitrary stable-coordinate mixing.  This is an all-degree no-go for
  that lane, not a counterexample or a general mixed-pair obstruction.
source: root / affine-mixed Hamiltonian-slice reframe, 2026-08-22
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
related:
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
script: 04-computation/jc2_compiler_generator_laurent_weight_no_mate_thm3684.py
output: 05-knowledge/results/jc2_compiler_generator_laurent_weight_no_mate_thm3684.out
script_sha256: 8d81dd6def597d36c493101dffe535f7a14e4eb1751a9eb0acba2d37203e7eb4
output_sha256: f0ab98356d890ed38f2a99e72674100d5015fc0ce316c5e83c8d9ba01480041e
semantic_sha256: 285475b809b60618be6c799a35ff0c1dd28caeb5ed24273ef6052d04480b4250
hash_basis: raw LF bytes for files; displayed generator/weight/operator ledger for semantic hash
---

# THM-3684 -- compiler generators have no polynomial Jacobian mate

**PROVED + VERIFIED-EXACT.**  A stable variable enlarges the target ring, so
the one-generator slice obstructions inside the Danielewski surface do not by
themselves exclude a mate after pullback.  The Laurent-weight calculation
below closes that gap in the strongest possible ambient ring: the proposed
mate may be any polynomial in the source coordinates.

All rings are over `C`.

## 1. Compiler and Laurent chart

Put

```text
D=1+x^2q,
b=(D-1)(D+2)^2,
c=xD(D+2),
e=q(D+3)=x^(-2)(D-1)(D+3).                            (1)
```

Localizing at `x` gives the injective chart

```text
C[x,q] -> C[x,x^(-1),D],             q=x^(-2)(D-1).  (2)
```

Every source polynomial therefore has a unique finite Laurent expansion

```text
G=sum_m x^m K_m(D),                    K_m in C[D].    (3)
```

The coefficient polynomials in `(3)` are not arbitrary, but only their
polynomiality will be needed.  Since the coordinate change `(x,q)->(x,D)`
has determinant `x^2`, one has

```text
Jac_(x,q)(F,G)=x^2 Jac_(x,D)(F,G).                    (4)
```

## 2. The pure-weight mate equation

Let

```text
F=x^a S(D),                 G_m=x^m K(D).              (5)
```

Direct differentiation in `(4)` gives

```text
Jac_(x,q)(F,G_m)
 =x^(a+m+1)[a S K'-m S'K].                            (6)
```

The powers of `x` in the Laurent ring are independent over `C[D]`.  Hence a
nonzero constant in `Jac(F,G)` can only come from

```text
m=-a-1,                                                (7)
```

and its coefficient must solve the polynomial ODE

```text
a S K'+(a+1)S'K=kappa,              kappa in C*.       (8)
```

Suppose `deg S=d>=1` and `K` is nonzero of degree `n`.  Unless its displayed
coefficient vanishes, the left side of `(8)` has degree `d+n-1` and leading
coefficient

```text
[a n+(a+1)d] lc(S)lc(K).                              (9)
```

This one-variable leading term is the entire obstruction; no bounded support
assumption has entered.

## 3. The three compiler generators

For `b`, use

```text
a=0,       S=(D-1)(D+2)^2,       d=3.                 (10)
```

The only possible Laurent coefficient has `m=-1`, and `(8)` is

```text
S'K=kappa.                                             (11)
```

Its degree is `n+2` and its leading coefficient is `3lc(K)`, so `(11)` has
no polynomial solution.

For `c`, use

```text
a=1,       S=D(D+2),             d=2.                 (12)
```

The only possible coefficient has `m=-2`, and `(8)` is

```text
S K'+2S'K=kappa.                                       (13)
```

Its degree is `n+1` and its leading coefficient is `(n+4)lc(K)`, again
excluding a polynomial solution.

For `e`, use

```text
a=-2,      S=(D-1)(D+3),         d=2.                 (14)
```

Now `m=1`, and `(8)` becomes

```text
-2S K'-S'K=kappa.                                     (15)
```

Its degree is `n+1` with leading coefficient `-(2n+2)lc(K)`.  It also has no
polynomial solution.  Therefore

```text
Jac(b,G), Jac(c,G), Jac(e,G) notin C*
for every G in C[x,q].                                (16)
```

This conclusion is stronger than nonexistence inside the target ring.

## 4. Polynomial functions of one generator

Let `z` be any one of `b,c,e`, and suppose `h in C[T]` is nonconstant.  If

```text
Jac(h(z),G)=h'(z)Jac(z,G)=kappa in C*,                (17)
```

then both factors on the left are units of `C[x,q]`.  Since `z` is
nonconstant and the characteristic is zero, this forces `h'` to be a nonzero
constant.  Thus `h` is affine, and `(17)` contradicts `(16)`.  Hence no
nonconstant polynomial in a single compiler generator has a polynomial
Jacobian mate.

## 5. Consequence for the counterexample search

Let a Russell-cylinder target pair be pulled back through any polynomial
source graph, affine or critical.  Its second output is still a polynomial in
`C[x,q]`.  If the first output is `h(b)`, `h(c)`, or `h(e)`, equation `(17)`
rules out a nonzero constant source Jacobian even when the second output has
arbitrary degree and arbitrary stable-coordinate dependence.  The same holds
after swapping the outputs.

Thus the live mixed-pair lane has a sharper entry condition:

```text
each output must mix at least two genuinely independent target directions;
stable mixing in only the opposite output cannot rescue a generator slice. (18)
```

The theorem does **not** exclude a nonlinear surface coordinate involving
several generators, a coordinate already mixing a surface generator with the
stable variable, the open Darboux system, or an arbitrary planar polynomial
map.  It constructs no Keller pair and no counterexample; `JC(2)` remains
**OPEN**.

## 6. Exact verification

The companion checks the compiler relation, Laurent conversion, the master
identity `(6)` on `360` pure-weight controls, the exact degree and leading
coefficient in all three ODEs for degrees `0` through `12`, and direct
source/chart agreement on a rectangular monomial bank.  It contains no
Python assertion statements.

Reproduce with

```bash
python3 -B 04-computation/jc2_compiler_generator_laurent_weight_no_mate_thm3684.py
python3 -O -B 04-computation/jc2_compiler_generator_laurent_weight_no_mate_thm3684.py
```

Both transcripts must agree byte-for-byte with the stored output.  **QED.**
