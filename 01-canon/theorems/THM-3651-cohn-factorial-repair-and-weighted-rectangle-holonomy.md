---
id: THM-3651
title: "Cohn factorial repair obstruction and weighted rectangle holonomy"
status: >
  PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED; CITED for Cohn
  non-elementarity; CONDITIONAL only for the motivating Weihrauch corollary.
  The Cohn core has no polynomial first repair, and no elementary decoration
  with at most one factor on each side is a Jacobian matrix.  Size-n
  rectangles in 2 by n are odd-rigid/even-balanced, but a coefficient cycle
  closes only when its multiplier product is (-1)^n.  The actual Cohn
  factorial holonomy never closes.
source: boxeph / parity recovery and Cohn-core hostile audits, 2026-08-21
depends_on: []
related:
  - THM-3650-wright-elementary-jacobian-criterion-reduced-word-reproof
  - THM-3613-three-by-four-size-seven-ray-parity-gate
script: 04-computation/jc2_cohn_parity_cycle_repair.py
output: 05-knowledge/results/jc2_cohn_parity_cycle_repair.out
script_sha256: 1a9ff3b9cb8190e434307400f658ff8bf476d1629924af91b61b8bfcc81091ee
output_sha256: 34f2ff8c8e9d78c1d48d37e1931d912c4713b1a05f5db5daa37e2b23b5a2a8f1
semantic_sha256: 5aa8a03f0a3530fd6c03161b19488f302b1c9ea348c5acd0f0d0f4a5f57510aa
external:
  - https://www.numdam.org/item/PMIHES_1966__30__5_0/
---

# THM-3651 -- Cohn repair tails and weighted odd/even holonomy

**PROVED + FINITE-EXACT + INDEPENDENTLY HOSTILE-AUDITED; CITED for the
non-elementarity statement.**

Work over a characteristic-zero field `k` and put

```text
C=[1+xy   x^2]
  [-y^2  1-xy].                                        (1)
```

Cohn proved that `(1)` lies in `SL_2(k[x,y])` but not in
`E_2(k[x,y])`.  The theorem below concerns its row integrability and a finite
odd/even transport mechanism; Cohn's literature result is not reproved.

## 1. The first repair is an infinite factorial tail

The determinant and row curls are

```text
det C=1,                 curls(C)=(-x,-y).              (2)
```

Left multiplication by `E_-(Z)` could close the lower row only if

```text
L(Z)=y,
L=(1+xy)partial_y-x^2 partial_x-x.                     (3)
```

Give `x^i y^j` weight `j-i`.  The operator `L` lowers weight by one, so a
solution of `(3)` must have weight two:

```text
Z=sum_(i>=0) a_i x^i y^(i+2).                          (4)
```

Coefficient comparison gives

```text
2a_0=1,                 (i+2)a_i+a_(i-1)=0,
a_i=(-1)^i/(i+2)!.                                      (5)
```

The degree-`D` truncation satisfies exactly

```text
L(Z_D)=y+(-1)^D x^(D+1)y^(D+2)/(D+2)!,                (6)
```

so no polynomial solution exists.  The unique regular formal solution is

```text
Z_hat=(exp(-xy)-1+xy)/x^2.                             (7)
```

The symmetric upper repair is

```text
R_hat=-(exp(xy)-1-xy)/y^2,                             (8)
```

and is likewise nonpolynomial.

## 2. One elementary factor on each side is still impossible

For arbitrary `W,Z in k[x,y]`, write the lower row of
`E_-(Z) C E_+(W)` as `(P,Q)`, where

```text
A=1+xy,       P=-y^2+AZ,
Q=1-xy+x^2Z+WP.
```

Then

```text
P_y-Q_x=-y+L(Z)-partial_x(WP).                         (9)
```

If `deg W=d>=1` and `deg Z=N>=1`, the unique top term in the closure equation
is

```text
-partial_x(xy W_d Z_N),                               (10)
```

of degree `d+N+1`.  It is nonzero.  Thus one of `W,Z` is constant.

For `Z=z`, closure becomes

```text
(y^2-z-zxy)W_x-zyW=xz+y.                              (11)
```

If `z=0`, this requires `W_x=1/y`.  If `z!=0`, `x`-degree at least two has
top coefficient `-zy(m+1)w_m`; degree one forces `w_1=-1/(2y)`; degree zero
cannot provide `xz`.  All are impossible.

For `W=c!=0`, the leading homogeneous equation with
`Z_N=x^N p(t)`, `t=y/x`, is

```text
t(2+ct)p'=(N+1)(1+ct)p.                               (12)
```

For even `N` the nonzero solution has a half-integral exponent.  For odd
`N` it is proportional to `[t(2+ct)]^((N+1)/2)`, of `t`-degree `N+1>N`.
The case `c=0` is `(3)--(7)`, and a constant `Z` fails directly.  Therefore
the lower row never closes.  Source/target exchange proves the symmetric
statement for `E_+(R) C E_-(U)`.

The same-sign layouts die from their untouched rows:

```text
E_-(Z) C E_-(U):       U_y=1/x,
E_+(R) C E_+(W):       W_x=1/y.                       (13)
```

Hence no decoration with at most one elementary factor on each side of `C`
is a Jacobian matrix.

## 3. Rectangle classification

Let a nonempty rectangle `A x B` in `g x n` have cardinality `n`.  Its height
`|A|` divides `g` and `n`, hence divides `gcd(g,n)`.  In particular, a
size-`n` rectangle in `2 x n` is

```text
n odd:   one full row;
n even:  one full row or 2 x S with |S|=n/2.           (14)
```

Suppose further that two rectangular size-`n` tooth blocks and `n`
rectangular size-two stripe blocks both partition `2 x n`, with every
tooth/stripe intersection a singleton.  For odd `n`, the teeth are the two
rows and the stripes the `n` columns, up to relabeling.  For `n=2m`, there is
one additional type:

```text
teeth=2 x S, 2 x S^c;
stripes=rowwise perfect matchings S -> S^c.            (15)
```

Conversely every pair of such matchings realizes `(15)`.  With fixed target
coordinates and unlabeled tooth blocks, the additional type has

```text
binomial(2m,m)/2 * (m!)^2=(2m)!/2                     (16)
```

realizations.  Its invariant sidecar is the cycle type of the relative
matching permutation.

This finite theorem does not prove the motivating Weihrauch implication.
That implication additionally requires a raw rectangular answer set of
exactly `n` cells and a committed decoder that preserves tooth/stripe
incidence bijectively.

## 4. Weighted cycle holonomy

For an abstract unweighted cyclic transport, the equations

```text
b_i+b_(i-1)=0                                           (17)
```

have matrix `I+S_n`, and

```text
det(I+S_n)=2          for odd n,
rank(I+S_n)=n-1       for even n.                       (18)
```

The even kernel is the alternating vector.  This is the same support-level
odd/even split as `(14)`, with overall sign the one-bit gauge.

Multiplier holonomy cannot be discarded.  The weighted equations

```text
alpha_i c_i+c_(i-1)=0                                  (19)
```

have determinant

```text
product_i alpha_i-(-1)^n.                              (20)
```

Thus even cycles need multiplier product `+1`, while odd cycles need `-1`.
Factorial rescaling of `(5)` flattens only the interior path: a wrap seam
retains Cohn weights `2,3,...,n+1`, whose product is `(n+1)!`.  Hence every
actual cyclic Cohn closure has determinant

```text
(n+1)!-(-1)^n !=0.                                     (21)
```

The smallest even rational kernel uses reciprocal weights `2,1/2`.  Parity
is therefore a support filter, never a substitute for edge-gain holonomy.

## 5. Reproduction and boundary

Run

```bash
python3 04-computation/jc2_cohn_parity_cycle_repair.py
python3 -O 04-computation/jc2_cohn_parity_cycle_repair.py
```

The exact companion verifies `(2)--(8)`, the generic identity `(9)`, cycle
ranks and rectangle counts through the displayed finite ranges, every actual
Cohn weighted cycle, and the reciprocal positive kernel.  Normal and
optimized transcripts reproduce the stored output byte for byte.

The theorem supplies no row-integrable non-elementary matrix, no polynomial
map, no collision, and no nonproperness.  Its positive construction target is
an elementary--Cohn--elementary decoration with at least two corrections on
one side, an exact reciprocal-weight transport cycle, both curls zero, and
the non-elementary class retained.  **QED.**
