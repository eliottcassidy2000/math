---
id: THM-3900
title: "F-zero generic y-polynomial root-color response classification"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit of the assembled candidate.  Over k(x)[y], the
  complete THM-3881 f=0 square chart has only the two T-coordinates T=0 and
  T=-2K/(3a^2).  After the proved THM-3895 y-degree cutoff, quadratic T is
  classified by the intrinsic colors G=+/-L^2 at its roots.  Same-color
  roots force the nonzero hostile point.  Split colors have one sharp
  constant-parameter solution, but require f^2=169/512 and are excluded by
  the actual nonconstant x-family.  Linear T has odd leading degree, and a
  duplicated fifteen-cell rational-root exhaustion closes constant T.
  No elliptic-surface or Mordell--Weil claim is used.
source: tournament-jc-breakthrough / post-THM-3895 root-response scout, 2026-08-23
audit: >
  SELF-AUDITED PROVISIONAL CANDIDATE, preceded by an independent hostile
  derivation of the mathematical packet.  The focused exact companion
  verifies the normalization, both known solutions, the same-color scalar
  equation, all split Hermite responses, the complete split response ideal,
  the sharp split hostile, the actual-family mismatch, the linear odd-degree
  obstruction, and all fifteen constant-T rational-root coefficient gcds in
  47 active gates.  Normal and optimized runs byte-match the frozen output.
  A final independent audit of this assembled file must recheck the generic
  scope of the THM-3895 cutoff, exhaustion of root multiplicities and colors,
  the Hermite interpolation quotient, and the rational-root theorem before
  status promotion.
depends_on:
  - THM-3895-f-zero-quartic-covariant-and-high-y-degree-emptiness
related:
  - THM-3885-cusp-residual-f-zero-arm-dichotomy-and-quadratic-closure
  - THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality
  - THM-3897-f-zero-residual-all-degree-global-emptiness
script: 04-computation/jc2_f_zero_generic_y_root_color_classification_thm3900.py
output: 05-knowledge/results/jc2_f_zero_generic_y_root_color_classification_thm3900.out
script_sha256: 366059dda300d6ec657592bd0708399bda365b0762f8777b1ec9fcaa216fcd95
output_sha256: 1778b9b63eff819ec289279a930f6f3889a53826e110c0efe459391e6ffd8525
semantic_sha256: 9e63b16c78b79eb2033b94d78951be2662d27a73c66bf43131f1c64e75d39638
hash_basis: raw LF bytes
---

# THM-3900 -- the generic f-zero chart has two T-coordinates

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit of the assembled candidate.**  Work over an
algebraically closed field `k` of characteristic zero.  Put

```text
C=k(x),                     a=x+1,
L=9x+4,                     F=15x^2+15x+4,
K=y^2-F.                                                    (1)
```

Consider the complete generic `f=0` square equation

```text
G^2=L^4-6aL^2T^2-8KT^3-3a^2T^4,          T,G in C[y].     (2)
```

The provisional classification is

```text
T=0,                         G=+/-L^2,                    (3)
```

or

```text
T=T_*=-2K/(3a^2),
G=+/-G_*,                    G_*=4K^2/(3a^3)-L^2.         (4)
```

Both alternatives are exact solutions.  There are no others.  The result is
over `k(x)[y]`; `(4)` is the canonical hostile showing that generic
classification and polynomial descent in `x` are different problems.

## 1. Normalized chart and the proved degree cutoff

Pass temporarily to the algebraic closure of `C`.  Choose nonzero `s,e`
with

```text
s^2=a,                       e^2=s^3L,                    (5)
```

and set

```text
Y=y/e,                       u=sT/L,
g=G/L^2,                     z=K/(s^3L)=Y^2-f,
f=F/(s^3L).                                                 (6)
```

Then `(2)` becomes exactly

```text
g^2=1-6u^2-8zu^3-3u^4.                                   (7)
```

The proof of THM-3895 works over `overline{k(x)}[y]` and gives

```text
deg_Y(u)<=2.                                               (8)
```

The scaling `(5),(6)` is invertible, so it remains to classify the three
degrees in `(8)`.  Notice also that

```text
f^2=F^2/(a^3L^2),                                         (9)
```

a nonconstant rational function of `x`.

## 2. Quadratic roots carry two intrinsic colors

Factor `(7)` as

```text
(g-1)(g+1)=-u^2(6+8zu+3u^2).                             (10)
```

The two left factors are coprime in `overline C[Y]`.  At every root of `u`,
the last factor on the right is `6`, so a root of multiplicity `r` in `u`
occurs with exact multiplicity `2r` in exactly one of `g-1` and `g+1`.
This assigns every root one of two colors, including a repeated quadratic
root.  Since the right side of `(7)` has degree at most eight,

```text
deg_Y(g)<=4.                                              (11)
```

### 2.1. Same color gives the hostile point

Suppose all roots of a quadratic `u` have the same color `epsilon`, where
`epsilon^2=1`.  Then degree `(11)` and the exact multiplicities in `(10)`
give a scalar `c` with

```text
g=epsilon+c u^2.                                         (12)
```

Substitution in `(7)` and division by the nonzero `u^2` gives

```text
(c^2+3)u^2+8zu+(2epsilon*c+6)=0.                         (13)
```

The first two terms in `(13)` are divisible by the nonconstant polynomial
`u`.  Hence the scalar last term vanishes:

```text
c=-3epsilon.                                              (14)
```

Equation `(13)` then reduces to

```text
u=-2z/3,                   g=epsilon(1-3u^2).             (15)
```

Undoing `(6)` gives exactly `(4)`, including the two signs of `G`.

### 2.2. Split colors have one sharp second-response family

Now suppose the two roots are distinct and receive opposite colors.  Write

```text
X=Y-t,
u=A(X^2-q^2),             A*q!=0.                        (16)
```

After changing the sign of `g` if necessary, the color and first-response
conditions are

```text
g(q)=1,   g'(q)=0,         g(-q)=-1,   g'(-q)=0.          (17)
```

The unique cubic Hermite interpolant is

```text
g_0=3X/(2q)-X^3/(2q^3).                                  (18)
```

Every polynomial of degree at most four satisfying `(17)` differs from
`g_0` by a scalar multiple of the fourfold root polynomial.  Absorbing the
nonzero factor `A^2` into the scalar, write

```text
g=g_0+c u^2.                                              (19)
```

Substitute `(16),(19)` and

```text
z=(X+t)^2-f                                               (20)
```

in `(7)`, and multiply by `4q^6`.  After removing nonzero common factors,
the coefficients in degrees seven, five, and three are

```text
16Aq^3t-c,
48Aq^3t-5c,
48Aq^3t-7c.                                               (21)
```

The first two force

```text
t=0,                         c=0.                         (22)
```

The degree-eight response is then `3A+8`, so

```text
A=-8/3.                                                   (23)
```

The remaining even response ideal has the exact reduced lexicographic basis

```text
[3f+13q^2, 512q^4-9].                                    (24)
```

If `u=AY^2+BY+C`, equations `(22)-(24)` give

```text
A=-8/3,                    B=0,
C^2=1/8,                   f^2=169/512.                  (25)
```

These conditions are sharp.  Whenever

```text
512q^4=9,                  f=-13q^2/3,                   (26)
```

the explicit pair

```text
u=-8(Y^2-q^2)/3,
g=3Y/(2q)-Y^3/(2q^3)                                     (27)
```

satisfies `(7)`.  Thus split colors are not formally contradictory; their
second-response payment is a genuine hostile family.

They do not occur in the actual `x`-family.  Combining `(9),(25)` would
require

```text
512F^2=169a^3L^2.                                        (28)
```

The left side of `(28)` has degree four in `x`, while the right side has
degree five.  More explicitly, their difference has leading coefficient
`-13689`.  Hence `(28)` is impossible, closing the split-color branch.

## 3. The linear channel has odd degree

Let

```text
u=A Y+B,                    A!=0.                         (29)
```

In `(7)`, the term `-8zu^3` has the unique leading form

```text
-8A^3Y^5.                                                 (30)
```

All other terms have degree at most four.  A nonzero polynomial square has
even degree, so the linear channel is empty.

## 4. Constant T: a self-contained fifteen-cell exhaustion

It remains to take `T=t in C`.  If `t=0`, equation `(2)` gives `(3)`.  If
`t!=0`, the right side of `(2)`, as a polynomial in `y`, is

```text
-8t^3y^2+C(t),
C(t)=L^4-6aL^2t^2+8Ft^3-3a^2t^4.                        (31)
```

There is no linear `y` term.  A square root must be linear with nonzero
leading coefficient; its constant term is therefore zero.  Thus a survivor
would require

```text
C(t)=0.                                                   (32)
```

The polynomial `C(t)` is primitive in `k[x][t]`.  Its constant and leading
coefficients are `L^4` and `-3a^2`, and `gcd(a,L)=1`.  The UFD rational-root
theorem therefore leaves exactly

```text
t=c L^i/a^j,
c in k^*,                    0<=i<=4, 0<=j<=2.            (33)
```

There are fifteen exponent pairs.  The focused THM-3900 companion duplicates
the complete calculation locally: for each pair, substitute `(33)`, clear
the denominator, collect the coefficients in `x`, and take their gcd in
`Q[c]`.  Every gcd is the unit `1`:

```text
(0,0):1 (0,1):1 (0,2):1
(1,0):1 (1,1):1 (1,2):1
(2,0):1 (2,1):1 (2,2):1
(3,0):1 (3,1):1 (3,2):1
(4,0):1 (4,1):1 (4,2):1.                                (34)
```

A unit coefficient ideal has no common zero after extension to `k`, so none
of `(33)` solves `(32)`.  This closes every nonzero constant.

## 5. Exhaustion, dependencies, and scope

The alternatives quadratic, linear, nonzero constant, and zero constant
exhaust `(8)`.  Sections 2--4 therefore prove the provisional classification
`(3),(4)` once the THM-3895 cutoff is imported.

The proof depends only on:

1. THM-3895's proved and independently audited generic-y cutoff;
2. the exact residual identity `(2)` inherited there; and
3. elementary UFD, root-multiplicity, Hermite-interpolation, and finite exact
   polynomial calculations reproduced in the focused companion.

THM-3888 is related provenance for `(4)` and the fifteen candidates, but no
elliptic surface, Mordell--Weil rank, bad-fibre, integrality, or other
provisional claim from THM-3888 is a dependency.

Over `k[x,y]`, the hostile `(4)` has an `a=x+1` denominator, so the generic
classification is consistent with the global `f=0` closure.  This theorem
does not address `f!=0`, construct a Keller atlas, or settle JC(2).

Reproduce with

```bash
python3 04-computation/jc2_f_zero_generic_y_root_color_classification_thm3900.py
python3 -O 04-computation/jc2_f_zero_generic_y_root_color_classification_thm3900.py
```

Both streams must byte-match
`05-knowledge/results/jc2_f_zero_generic_y_root_color_classification_thm3900.out`.
