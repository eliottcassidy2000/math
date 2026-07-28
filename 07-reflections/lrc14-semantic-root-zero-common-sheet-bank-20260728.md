# The natural 81-sheet clutch has universal units, not a universal gain

> **STATUS: FINITE-EXACT CANDIDATE; AWAITING INDEPENDENT HOSTILE AUDIT.**
> On canonical rail 8, identify the present and delayed clock labels, insert
> the full `E3` mask and delayed `Q_(3,{1,2})` word, and scan the `9 x 9`
> target sheets common to the two root-zero endpoint cylinders.  Every one of
> the 81 natural one-sided source/target pairs is a constant-line private
> unit at both roots.  The normalized target/source gain is nevertheless `7`
> on only 22 sheets and takes five values overall.  The first hostile is
> `(s,t)=(0,8)`, with `12 -> 3` rather than `12 -> 6`.  A type-correct
> unipotent shear exists only after adjoining a pair-derived reference table;
> no physical reference packet, endpoint current, or LRC(14) conclusion is
> constructed.

## 1. Inheritance and question

`THM-2744-relative-present-unit-repair-and-root-zero-overlap-clutch.md`
constructs the open physical rechart

```text
(carry 12, right root 12) -> (carry 6, left root 1),
```

under translation by `7/13^6`.  Before semantic/target cuts, canonical rail
8 has equal raw source and target vectors and private units.  The exact
single natural sheet `(s,t)=(0,3)` at commit `ad375e51e` remains a unit on
both sides but has normalized profiles `5 -> 9`; its gain is

```text
9/5 = 7 = 6/12 mod 13.                                  (1)
```

`THM-2749-fully-marked-root-zero-clutch-and-target-character-profile.md`
instead intersects the marked sheet with its translated pullback.  That
two-sided common section has an exact raw intertwiner and normalized mirror
gain `-1`.  The natural one-sided object and the two-sided fibre product must
therefore remain distinct.

The proposed THM-2754 addendum asks the next exact question: after identifying
the present and delayed clock indices, is `(1)` a genuine clutch law on all
common target sheets, a sheet character, or an accident of one sheet?

## 2. Exact bank and typing

The source and target semantic points have common label banks

```text
s in (0,1,2,3,8,9,10,11,12),
t in (3,4,5,6,7,8,9,10,11).                            (2)
```

For each of the 81 pairs, the companion independently:

1. restricts the rail-8 overlap by both relative-present complements;
2. inserts the full unshifted `E3` mask in the source and forward target
   coordinates separately;
3. inserts the four shifted safe factors defining the natural `U_(s,t)`;
4. rebuilds delayed sector zero with both lower blockers dangerous, namely
   `Q_(3,{1,2})`, and uses the same clock `ell` in the present and delayed
   factors;
5. evaluates carry `12` at root `12` and carry `6` at root `1`, then divides
   by content `26`, root-normalizes, and reduces modulo `Phi_7`.

Delayed prefix values are cached only within a fixed clock.  Source and target
interval banks and weights remain separate.  Thus the scan neither freezes a
source-coordinate mask on the target nor reuses a coefficient across charts.

## 3. Universal survival and exact failure of gain seven

All 81 source and all 81 target vectors have shape

```text
(0,a,a,a,a,a,a),
```

have content `26`, and are private units.  Their scalar profile-pair census is

```text
(5,9)   x10,      (12,3) x19,      (9,8)  x2,
(4,12)  x 9,      (7,12) x 2,      (8,6)  x5,
(11,6)  x 2,      (7,5)  x20,      (2,1)  x12.         (3)
```

The target/source gains are

```text
gain 3: 9 sheets,     gain 4: 5 sheets,
gain 7:22 sheets,     gain10:41 sheets,
gain11: 4 sheets.                                      (4)
```

The lexicographically first failure of `(1)` is

```text
(s,t)=(0,8):       source=12, target=3,
                   target/source=10 != 7.              (5)
```

Thus gain `7` is only the predecessor ratio on a 22-sheet subbank.  It is not
a uniform multiplicative clutch scalar on the natural 81-sheet carrier.
Because `|F_13^*|=12`, every homomorphism `C_13 -> F_13^*` is trivial; none of
the five nontrivial gains can be an additive-sheet character either.

The strongest survivor is consequently sharper and differently typed:
**every common natural sheet carries a unit on both sides, but phase/gain
transport is sheet dependent.**

## 4. The smallest type-correct additive model

Let `p_-` and `p_+` be the two scalar profiles and note that the carry
difference is

```text
6-12=7 mod 13.                                         (6)
```

On `F_13^2`, put

```text
N(r,p)=(0,r),       U(a)=I+aN,       N^2=0.             (7)
```

Then `U(a)U(b)=U(a+b)` and `U(13)=I`, so `(7)` is the minimal nontrivial
characteristic-13 representation of the additive carry group.  For every
sheet there is a unique reference

```text
r_(s,t)=(p_+-p_-)/7                                   (8)
```

such that

```text
U(7)(r_(s,t),p_-)=(r_(s,t),p_+).                       (9)
```

The reference census is

```text
r=3:11,     r=8:29,     r=9:25,     r=10:2,     r=11:14. (10)
```

Hence the observed sheet `(0,3)` does admit

```text
U(7)(8,5)=(8,9),                                      (11)
```

but `8` is not a global reference.  On the 22 gain-seven sheets only,
`p_-/12=p_+/6`, so the reference can also be read from either endpoint by
carry normalization.  On the other 59 sheets this endpointwise
normalization disagrees, and `(8)` genuinely uses the ordered endpoint pair.

This is the exact type boundary.  The table `(8)` is uniquely determined
once the oriented carry difference is fixed; it is not an arbitrary numeric
basis choice.  But it is still a **derived coefficient table**, not an
independently supplied physical packet or common-endpoint amplitude.  Using
it as a reference extension requires a new sidecar.

Flattened as a `9 x 9` target-label tensor, the source profile, target
profile, gain, and reference tables each have rank exactly three over
`F_13`.  The reference table is constant on the three `t`-blocks

```text
(3,4,5,6,7),       (8,9),       (10,11),               (12)
```

and its three corresponding column vectors are independent.  Therefore
three separable labelled references are algebraically necessary and
sufficient to realize `(8)`.  This rank-three statement is a coefficient
factorization, not yet three physical arms or a `C3` action.

The rank lower bounds have explicit `3 x 3` witnesses, using columns
`t=(3,8,10)`.  Modulo `13`, their determinants are

```text
source rows s=(0,1,2):       9,
target rows s=(0,1,2):       4,
gain rows s=(0,1,2):         5,
reference rows s=(0,2,3):    6.                         (13)
```

Column constancy on `(12)` gives rank at most three, and `(13)` gives rank at
least three.  No numerical rank inference is used.

This gives a literal but sharply scoped `2/3` interface: additive carry first
requires the two-dimensional Jordan extension `(7)`, while the resulting
global sheet reference has separable tensor rank three.  The three target
blocks have sizes `5,2,2` and are not cyclically permuted, so this is not yet
the `C2*C3` or modular-group grammar; a genuine `C3` interpretation would
need an operator rotating three physical reference packets, not merely rank
three.

## 5. Minimal hostile to the scalar interpretation

Multiplication by `7` on one `F_13` line has order twelve, not thirteen:

```text
7^12=1,       7^13=7 !=1.                              (14)
```

Two additive carry-seven moves compose to carry `14=1`, whereas two scalar
gain-seven maps compose to gain `49=10`.  In contrast,

```text
U(7)^2=U(14)=U(1).                                     (15)
```

Thus even on the 22 sheets satisfying `(1)`, the scalar ratio is at most a
single chart-transition coordinate.  It cannot represent the additive
odometer class by itself.  The minimal legal object is the unipotent
reference/profile doublet `(7)--(9)`.

## 6. Consequence and next test

The natural sheet problem no longer has a survival debt: all 81 coefficients
are units.  It has a reference/phase debt.  The cheapest decisive continuation
is to decompose the pair-derived rank-three table into the three blocks in
`(12)` and construct those three separable references directly on the same
diagonal-clock carrier.  MISTAKE-313 and the repaired/refuted THM-2751 forbid
obtaining them by subtracting THM-2749's fixed-present-clock common section:
that comparison omits or folds incompatible clock factors, and the lawful
full-unclocked target wing augments to zero.  A direct positive construction
would produce a genuine marked unipotent clutch.  A negative result would
identify the exact missing common-endpoint reference rather than another
support or unit obstruction.

## 7. Exact reproduction

Companion:

```text
04-computation/lrc14_semantic_root_zero_common_sheet_bank_probe_20260728.py
```

Run:

```bash
python 04-computation/lrc14_semantic_root_zero_common_sheet_bank_probe_20260728.py
python -O 04-computation/lrc14_semantic_root_zero_common_sheet_bank_probe_20260728.py
```

The two transcripts must byte-match the stored output.  The computation uses
exact integer interval arithmetic, explicit exceptions, and no truth-bearing
Python `assert`.

LF-normalized SHA-256:

```text
script  5894f964809164298e4bb5b28183f177d06df0d4773e690ccad8f03c9bf1792d
output  59d34deab60d2e797d296e64c170a1a024d315728db698e8eac3a5b2fd2b5b9e
```
