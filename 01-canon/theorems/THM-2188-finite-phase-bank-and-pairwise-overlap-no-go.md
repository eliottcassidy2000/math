---
id: THM-2188
title: "Finite phase banks and pairwise overlaps do not determine LRC safety"
status: >
  PROVED. For every finite denominator bound B and septimal window depth h,
  there are two labelled primitive thirteen-speed rows with identical
  residues modulo every Q<=B and modulo 7^h, identical complete singleton
  and pair danger-overlap matrices, and a common saturated rank-eleven
  relation carrier of full rank modulo every base, but one has zero safe
  measure and the other has positive safe measure. Their full relation
  lattices differ by one primitive quotient slope, congruent in every
  prescribed finite bank. They even have the same endpoint phase word and
  unscaled numerator. The separating coordinate is the Archimedean lift
  `1/W`. Thus no finite bank/window, weighted pair graph, fixed high-rank
  relation carrier, or unscaled endpoint label determines the LRC target.
source: codex-2026-07-24-LRC-phase-bank-hostile-family
depends_on:
  - THM-541
  - THM-1166
  - THM-2174
related:
  - THM-2161
  - THM-2167
  - THM-2184
  - THM-2187
---

# THM-2188 -- finite phase and pairwise data remain blind

At radius `1/14`, write

```text
D_v={t in R/Z:||vt||<1/14},
G_E=(R/Z) \ union_(v in E) D_v.                       (1)
```

The theorem gives one uniform hostile family for finite residue banks,
radix windows, pairwise Tournament Analysis, and scale-free endpoint data.

## 1. The two labelled rows

Put

```text
C={1,2,3,4,5,6,8,9,10,11,12,13},
A=C union {7}={1,...,13},                             (2)

L_0=lcm{14c:c in C}=720720.                           (3)
```

Fix arbitrary integers `B,h>=1`, and define

```text
L_*=lcm(L_0,1,2,...,B,7^h),
W=7+L_*,
A_(B,h)=C union {W}.                                  (4)
```

The star label is `7` in `A` and `W` in `A_(B,h)`; all labels in `C`
retain their positions. Both rows are positive, distinct, and primitive.

THM-541 gives

```text
mu:=|G_C|=4319/51480.                                 (5)
```

The complete AP row `A` has

```text
|G_A|=0.                                              (6)
```

We now compute the second measure exactly.

## 2. The same endpoint word, a different scale

All endpoints of `G_C` have denominator dividing `L_0`. Because

```text
W=7 mod L_0,                                          (7)
```

the star speeds `7` and `W` induce the same complete endpoint phase word.
In THM-2174's notation,

```text
C_C(W)=C_C(7).                                        (8)
```

Apply its exact phase-scale formula first at speed `7`. Equations (5)--(6)
give

```text
0=6mu/7-C_C(7)/7,
C_C(7)=6mu=4319/8580.                                (9)
```

Apply the same formula at `W` and use (8):

```text
|G_(A_(B,h))|
 =6mu/7-6mu/W
 =(6mu/7)(1-7/W)
 =4319(W-7)/(60060W)>0.                              (10)
```

Thus the two rows have opposite zero-safety outcomes even though the
complete endpoint labels and the **unscaled** numerator in (9) agree.
The separating coordinate is exactly the scale divisor `1/W`.

## 3. Every prescribed finite bank and window agrees

By construction, `L_*` is divisible by every `Q<=B` and by `7^h`. Hence

```text
W=7 mod Q                 for every Q<=B,
W=7 mod 7^h.                                            (11)
```

The complete labelled residue word of the two rows therefore agrees at
every denominator in the prescribed bank and through the first `h`
septimal digits. This includes every danger clutter computed from those
residue words, not merely their zero masks.

The quantifiers are important: one family defeats each **fixed** pair
`(B,h)`. It does not assert that all unbounded residues agree.

## 4. The full pair-overlap matrix agrees

Pairs internal to `C` are unchanged, and every individual danger band has
measure `1/7`. For `c in C`, one has

```text
gcd(c,7)=gcd(c,W)=1,             W=7 mod 14.          (12)
```

Here `gcd(c,W)=gcd(c,7)` because `L_*` is divisible by every `c in C`.
THM-1166's exact two-comb formula, together with the symmetry
`c <-> 14-c`, gives

```text
|D_7 intersection D_c|
 =|D_W intersection D_c|
 =1/49.                                               (13)
```

Consequently all thirteen singleton masses and the entire labelled
`13 x 13` pair-overlap matrix are identical for `A` and `A_(B,h)`.
No weighted pair graph or tournament extracted from those values can
distinguish the target.

This conclusion is deliberately pairwise. The full Boolean atom law and
higher intersections need not agree; THM-2184 identifies that joint
continuation profile as the faithful finite object on a fixed affine ray.

## 5. A common rank-eleven carrier and one quotient slope

Let `x_*` denote the coefficient at the star label and define

```text
pi:Z^13 -> Z^2,
pi(x)=(x_*,sum_(c in C)c x_c).                        (14)
```

This map is onto: the star coordinate maps to `(1,0)` and the core
coordinate at speed `1` maps to `(0,1)`. Therefore

```text
L_C=ker(pi)                                           (15)
```

is a saturated rank-eleven lattice. Every element of `L_C` annihilates
both labelled rows, and saturation makes its reduction have rank eleven
modulo every integer base. Thus the hostile pair shares a universal-radix
relation carrier far larger than one chosen relation plane.

The full relation lattice of a row `C union {w}` has the exact quotient
description

```text
Lambda(C union {w})
 =pi^(-1)(Z(1,-w)).                                   (16)
```

Hence the two rank-twelve lattices differ only by the primitive transverse
slopes

```text
(1,-7),                    (1,-W)                    (17)
```

in `Z^13/L_C`. Equation (11) makes those slopes identical modulo every
prescribed `Q`, but their Archimedean lifts differ. The endpoint formula
reads that difference as `1/7` versus `1/W`.

The complete relation lattices do not literally agree. For example, the AP
row has additional small star relations; (16) says exactly what is shared.

## 6. Consequence and carrier boundary

For every finite `(B,h)`, the two rows agree on:

```text
- all labelled residue words modulo Q<=B;
- the first h base-seven digits;
- all singleton and pair danger masses;
- one fixed saturated rank-eleven universal-radix relation carrier;
- the deletion core C, its mass and components;
- the complete core-endpoint phase word;
- the unscaled signed endpoint numerator.             (18)
```

Yet (6) and (10) give different LRC targets. Hence none of these data, even
in combination, is a target-determining finite quotient.

The result sharpens rather than contradicts THM-2187: saturation makes any
finite radix bank free, but the final phase sidecar must retain higher-order
continuation and the scale-sensitive current. A pairwise tournament is
intrinsically too low-arity. QED.
