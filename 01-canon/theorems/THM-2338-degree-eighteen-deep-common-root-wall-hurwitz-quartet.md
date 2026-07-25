---
id: THM-2338
title: "Degree-eighteen deep common-root wall and the Hurwitz quartet"
status: >
  PROVED + VERIFIED-EXACT. On the simultaneous wall 126D=25B^2 and
  21W=-20BC, the structured Mordell polynomial 4P^3+49Q^2 has
  square-class degree at most four at exactly three weighted-projective
  parameter orbits. Every one has square-class degree exactly four;
  hence this wall contains three H_4 S_4^2 orbits and no H_2 S_5^2
  orbit. These are square-class candidates only. No Keller trajectory,
  scalar branch, JC(2), or DC(2) closure is asserted.
source: codex-2026-07-25-deep-common-root-wall
depends_on:
  - THM-2332-degree-eighteen-genus-zero-square-class-and-dessin-trap
  - THM-2335-degree-eighteen-cyclic-square-class-stratum-empty
related:
  - THM-2297-degree-eighteen-target-translation-normal-form
  - THM-2314-degree-eighteen-bd-ratio-bank-closure
script: 04-computation/jc2_degree18_deep_common_root_wall_thm2338.py
output: 05-knowledge/results/jc2_degree18_deep_common_root_wall_thm2338.out
script_sha256: 2c183329f2338277555051147466b49865cdee8549ef55d45b8c86faf5c40c95
output_sha256: 4006502b1f23d739be4dfdc6c74c93061ba282ee6ab5cd54e611ae3a6b21932f
hash_basis: working-tree bytes (LF)
---

# THM-2338 -- three Hurwitz quartets on the deep common-root wall

**PROVED + VERIFIED-EXACT.**

Retain THM-2332's structured polynomials

```text
P
 =245y^4+1890By^2-24300B^2+122472D,

Q
 =539y^6+11340By^4+183708Cy^3
  +(72900B^2-367416D)y^2
  +(2361960BC+2480058W)y,

F=4P^3+49Q^2.                                      (1)
```

Consider the simultaneous common-root wall

```text
126D=25B^2,                21W=-20BC.               (2)
```

On its nonzero weighted-projective parameter space, `F` has
square-class degree at most four at exactly three orbits. All three have
degree exactly four. Thus (2) contains no cyclic or mixed low-branch
point and exactly three four-simple square-class candidates.

## 1. The wall exposes a residual sextic

Conditions (2) make both depressed-cubic covariants vanish deeply at
`y=0`:

```text
P=35y^2(54B+7y^2),

Q=7y^3(1620By+26244C+77y^3).                       (3)
```

Consequently,

```text
F=9261y^6 G(y),                                    (4)
```

where

```text
G
 =7889y^6+211680By^4+1047816Cy^3
  +1814400B^2y^2+22044960BCy
  +2916000B^3+178564176C^2.                        (5)
```

The constant in (4) is nonzero and `y^6` is a square. Over `C`, the
square class of `F` is therefore exactly the square class of `G`.

The discriminant of the sextic factors completely:

```text
Disc_y(G)
 =-2^29 3^30 5^9 7^15 23

  *(361B^3+30618C^2)^3

  *(62500B^6+7654500B^3C^2+199644669C^4).          (6)
```

This factorization is the wall's complete collision ledger.

## 2. Weighted-projective reduction

If `B=0`, projectivity forces `C!=0`. Both nonconstant factors in (6)
are then nonzero:

```text
30618C^2 !=0,               199644669C^4 !=0.       (7)
```

So `G` is squarefree and its square-class degree is six.

Suppose `B!=0`. The weighted action

```text
(B,C,D,W) ->
(lambda^2B,lambda^3C,lambda^4D,lambda^5W)          (8)
```

preserves the square class after the corresponding rescaling of `y`.
Normalize `B=1` and put

```text
t=C^2/B^3=C^2.                                     (9)
```

The two factors in (6) become

```text
A(t)=361+30618t,

K(t)=62500+7654500t+199644669t^2.                  (10)
```

Their three roots are

```text
t_0=-361/30618,

t_+=-250/13041 + (500/117369)sqrt(3),

t_-=-250/13041 - (500/117369)sqrt(3).              (11)
```

They are distinct and nonzero. Indeed,

```text
K(t_0)=383/108 !=0,

Disc(K)=8680203000000 !=0.                         (12)
```

Outside (11), equation (6) makes `G` squarefree, so its square-class
degree is six.

The two choices of `C` over a fixed `t` give the same
weighted-projective orbit: the residual scaling `lambda=-1` fixes
`B,D` and negates `C,W`. Hence (11) represents exactly three orbits,
not six.

## 3. Every exceptional sextic has one double root

It remains to distinguish square-class degree four from degree two.
The exact Euclidean algorithm for `G` and `dG/dy` does this without
factoring the residual quartic.

At `t=t_0`,

```text
gcd(G,G')=y+(486/19)C.                              (13)
```

For `t=t_+` and `t=t_-`, respectively, put

```text
k_+=-81/5-(27/5)sqrt(3),

k_-=-81/5+(27/5)sqrt(3).                           (14)
```

Then

```text
gcd(G,G')=y-k_+C                 at t=t_+,

gcd(G,G')=y-k_-C                 at t=t_-.          (15)
```

The quadratic-branch root scale also has the radical-free
characterization

```text
2500k+3168963t+101250=0,

25k^2+810k+4374=0.                                (16)
```

Every gcd in (13)--(15) has degree one. Therefore each exceptional
sextic has exactly one double root and four other simple roots. If `r`
is the double root, then

```text
G=(y-r)^2 H_4(y),                 H_4 squarefree,

F=9261 [y^3(y-r)]^2 H_4(y),      deg(H_4)=4.        (17)
```

Thus all three exceptional orbits lie in the `H_4 S_4^2` stratum, and
none lies in `H_2 S_5^2`.

All coordinates are nonzero at these orbits: (11) gives `B,C!=0`, and
(2) then gives `D,W!=0`. They are genuine four-support points, not
another disguised two-sparse face.

## 4. Frontier effect

THM-2335 showed that the abstract cyclic `(3,3)` signature is never
realized anywhere in the coefficient cone. This theorem shows a
different phenomenon: the four-transposition square class is realized,
but the deepest common-root wall compresses its entire realization to
the three exact ratios (11).

The residual degree-eighteen search now has a concrete three-point
four-support bank on (2). Each point should be tested next against:

```text
local normality and actual field discriminant,
the Faber equations,
the Keller one-form and first flux,
the nonsplit deck condition,
and the retained rational trajectory.              (18)
```

Away from (2), both `H_2 S_5^2` and `H_4 S_4^2` remain open. Equation
(17) is a square-class certificate, not a Keller trajectory or even a
claim that every even discriminant zero is genuine cover
ramification. Other Newton edges, `JC(2)`, and `DC(2)` remain open.

## 5. Exact reproduction

Run

```bash
python3 04-computation/jc2_degree18_deep_common_root_wall_thm2338.py
python3 -O 04-computation/jc2_degree18_deep_common_root_wall_thm2338.py
```

Both transcripts are byte-identical to the stored output. The companion
verifies (3)--(6), both projective charts, all three exact ratios and
their separation, each number-field gcd in (13)--(15), the
radical-free relations (16), and generic squarefree hostile controls.
No executable check uses Python `assert`.

Independent audit is pending. QED.
