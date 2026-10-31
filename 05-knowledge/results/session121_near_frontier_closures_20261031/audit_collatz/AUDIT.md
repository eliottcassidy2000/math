# Independent audit: one E14 native-price obligation closes

## Verdict and inherited scope

**PASS, conditional on the pinned source's E10 certificate and native-return identities.**
Pinned source commit: `0578f16726e4e535189fb6074cc2ddb024c95d82`.
Read `precision_lifts_ramsey_deficits_20261030/COLLATZ.md`, C5–C6, and
`NATIVE14.md` directly from that commit. This audit does not replay the E10
certificate and does not establish universal Collatz or all-layer control.

The stronger universal argument below proves, without requiring the five-form
grammar, that every non-one primitive13 excursion has fee at least
`r=1062011676/1220703125 > 4/5`. In particular it proves

\[
 g_{13}(-(4,4))\ge 4/5.
\]

It also excludes all four padded forms for the other obligation. The only
remaining possible strict failure of `g13(0,0)>=1` is a **bare primitive11 core**
which is also primitive13 and has the required identity endpoint modulo `5^14`,
with one of the following necessary score/clock pairs:

| Core eighth-score j | Maximum core clock A |
|---|---:|
| 23 | 48,828,467 |
| 24 | 97,656,592 |

These are two score classes, not two words or existence claims. Thus E14 itself
remains open.

## 1. Integral score and exact small-clock reduction

Write `eta=183105426/48828125`, `D12=8*5^11=390625000`, and
`D14=8*5^13=9765625000`. A non-one primitive11 core c has

\[
 H_c=\frac{5A_c}{8}-S_{11}(c)=\frac{j_c}{8},\qquad j_c\in\mathbb Z.
\]

Its inherited native fee at slope lambda12 is at least `eta-1`. Since
`lambda12=5/8-1/D12`, this says exactly

\[
 \frac{j_c}{8}\ge\eta-1+\frac{A_c}{D12}.
\]

Positive A gives `j_c>=22`. For `j_c=22`, the displayed inequality gives

\[
 A_c\le D12(22/8-\eta+1)=342.
\]

For scores23 and24 the same exact calculation gives the caps in the table.
No rounding of eta is used.

## 2. Exact charge accounting, not a surrogate price

Let `T11=19531250`, `T12=97656250`. For an all-one primitive e block,

\[
 S_e(1^{T_e})=\sum_{k=1}^e\frac{T_e}{2\,5^{k-1}},\qquad
 \frac{5T_e}{8}-S_e(1^{T_e})=\frac14.
\]

For a bare primitive11 core c which is also primitive13, no proper prefix
returns at levels11,12,13. Its endpoint contributes exactly two additional
units to S13 relative to S11. Therefore

\[
 8H_{13}(c)=j_c-16.
\]

For a T11 decoration, cut the unique primitive12 constituent into the core and
the all-one primitive11 block. Their S11 charges add because the entering maps
are identity modulo5^11. There is exactly one extra level12 return and one
extra level13 return. The all-one block adds H=1/4, so `8H13=j_c-14`.

For a T12 decoration, the undecorated core is the unique non-one primitive12
constituent and is itself primitive11. Its S12 is S11+1. The all-one primitive12
block adds H=1/4, and the full primitive13 word adds its single level13 charge.
Again `8H13=j_c-14`. Prefix and suffix placements are separate words; these
charge identities are valid for either only under the native conditions of
the inherited grammar. They do not declare arbitrary concatenations native.

Consequently the exact fee formulas are

\[
 \delta_{13}(c)=\frac{j_c-16}{8}-\frac{A_c}{D14},
\]

\[
 \delta_{13}(1^{T_e}c)=\delta_{13}(c1^{T_e})
 =\frac{j_c-14}{8}-\frac{A_c+T_e}{D14}
\]

for each valid decorated candidate. In the second formula equality means the
fee determined by the native cut identities, not an assertion that both
placements exist with the same core.

## 3. Strict price cutoffs

The inherited five-form bounds are independently recalculated by the checker
from `(target-eta+3-padding_fee)/(lambda14-lambda12)`, with
`lambda14-lambda12=3/5^13`. The strict rational bounds and integer maxima are

| Form / target | Strict upper bound on A_c | Integer maximum |
|---|---:|---:|
| bare / 1 | 305176850/3 | 101725616 |
| T11 / 1 | 2442475/3 | 814158 |
| T12 / 1 | 12208100/3 | 4069366 |
| bare / 4/5 | 61036225/3 | 20345408 |

Inserting the integer maxima into the exact fee formulas shows:

- At the low target, `j_c>=23` implies fee at least
  `7/8-20345408/D14 > 4/5`.
- At either T11 high-target placement, `j_c>=23` implies fee at least
  `9/8-(814158+T11)/D14 > 1`.
- At either T12 high-target placement, `j_c>=23` implies fee at least
  `9/8-(4069366+T12)/D14 > 1`.
- At the bare high target, `j_c>=25` implies fee at least
  `9/8-101725616/D14 > 1`.

Thus the first three rows of exclusions reduce to `j_c=22`, hence `A_c<=342`.
The bare high target also loses score22 by the same obstruction below, leaving
exactly the necessary score possibilities23 and24.

## 4. Literal multiplier obstruction

Put `q=5^14=6103515625`. A positive valuation word of length L and clock A has
multiplier `3^L/2^A mod q`. A non-one positive word satisfies `1<=L<=A-1`.
The endpoint required by tangent zero is `(1,0)`, and that required by
`-(4,4)=(1,1)` is `(1+5^13,5^13)` modulo q.

An all-one block has affine map `(a,a-1)`, where `a=(3/2)^T mod q`.
For either identity-target placement, the core must have multiplier `a^-1`.
Indeed its full required affine map is `(a^-1,a^-1-1)`, checked by composing
on both sides; no translation constraint is silently relaxed to imply
existence. The four required multipliers are

| Form | Required core multiplier modulo q |
|---|---:|
| bare, high target | 1 |
| bare, low target | 1220703126 |
| T11, either placement | 2246093751 |
| T12, either placement | 5126953126 |

The independent checker directly evaluates powers modulo q and tests

\[
 3^L\equiv m\,2^A\pmod q
\]

for every one of the 58,311 pairs `1<=L<A<=342` against all four multipliers.
There are **zero matches**. This excludes the core already at the multiplier
level; no translation, realizability, or native-word assumption is needed
for the exclusion. There is no reliance on discrete logarithms, source
floor-sum counts, or producer imports.

## 5. Stronger universal route: independently audited after the candidate proof

The producer's `collatz/QUANTIZED_NATIVE14.md` gives a stronger route. This
was independently rechecked, rather than inferred from the four target-specific
modular scans above. A primitive11 core must have identity multiplier modulo
`5^11`. The added direct scan tests `3^L=2^A mod5^11` for every
`1<=L<A<=342`, again 58,311 pairs, and finds zero matches. Consequently
**every** non-one primitive11 core has `H>=23/8`, independently of its later
endpoint or any price cutoff.

Put `gamma=eta-1` and `D=D12`. The inherited inequality implies
`A/D<=H-gamma`. Hence

\[
 \lambda_{14}A-S_{11}
 =H-\frac{A}{25D}
 \ge\frac{24}{25}H+\frac{\gamma}{25}
 \ge 2+r,\qquad r=\frac{1062011676}{1220703125}.
\]

For an arbitrary non-one primitive13 excursion, cut at all its actual12
returns. Let B be the number of non-one primitive12 blocks, and a12 the
number of all-one primitive12 blocks. Cut only those B non-one blocks at
actual11 returns, obtaining C non-one cores and a11 all-one primitive11
blocks. Since a concatenation of all-one words is all-one, `C>=B>=1`.
All clocks are retained. The total extra level12 charge over those B inner
cuts is exactly B, and the full word has exactly one extra level13 charge.
Therefore

\[
 \delta_{13}\ge(2+r)C-B-1+\frac{31}{125}a_{11}+\frac6{25}a_{12}.
\]

The universal lower bound follows from the **identity**

\[
 (2+r)C-B-1-r=(1+r)(C-1)+(C-B)\ge0.
\]

This is not a finite check of possible C,B. It proves the inequality for all
counts. Thus the uniform fee is at least r, and strictly above4/5.
For a strict fee below1, C>=2 is impossible because `2(1+r)-1>1`.
So C=B=1. Both positive decoration counts are impossible because
`r+31/125>1` and `r+6/25>1`. This recovers the bare-core reduction without
using the source five-form grammar. For a bare core with j>=25, the same
linear bound gives

`delta13 >= (24/25)*(25/8)+gamma/25-2 > 1`,

also without the old coarse clock cutoff. This leaves exactly the same
score23 and24 clock bands established above. All separating inequalities and
the rational r are independently checked in `check_independent.py`.

## 6. Remaining precise target and reproducibility

The unresolved identity-target core must satisfy all the original native,
affine-translation and positive-word conditions, not just the table's clock
and score conditions. Integrality also gives the harmless extra filters
`A_c=5*j_c mod8`: score23 has `A_c=3 mod8`, and score24 has `A_c=0 mod8`.
These necessary conditions do not prove existence or nonexistence.

Run:

```sh
python check_independent.py > receipt.json
python -O check_independent.py > receipt_optimized.json
cmp receipt.json receipt_optimized.json
```

Both executions pass and their receipts are byte-identical. Checks use explicit
exceptions, not removable Python assertions. The first proof route assumes the pinned five-form grammar and native11 fee.
The stronger route assumes only the native11 fee and native-cut identities.
The checker independently checks the numerical bounds, exact all-one prices,
endpoint inverses, and every small-clock multiplier pair needed for both routes.
The native charge accounting is established in the written proof, not by a
claim to enumerate all words.
