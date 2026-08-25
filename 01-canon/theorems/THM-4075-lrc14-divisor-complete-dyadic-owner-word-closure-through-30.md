---
id: THM-4075
title: "LRC(14) divisor-complete dyadic owner-word closure through 30"
status: >
  PROVED RELATIVE TO THM-2061/2066 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. Every eleven-element positive core C with max(C)<=30 that
  divisibility-covers 2,...,14 closes the dyadic seam
  2C union {x,y} for every two distinct positive odd tails. The new exact
  layers max(C)=25,...,30 contain 1,764,356 cores, including all 389
  nonprimitive cores, and every one has an empty owner-word obstruction at
  some clock 15<=N<=38. Together with THM-2066 this closes 1,824,236 cores.
  Hence a THM-4041/4070 residual has max(H)>=31. In the full THM-3818
  physical producer this numerical floor is dominated by the inherited
  crossing-height bound; this is a finite dyadic-seam closure, not LRC(14).
source: codex-frontier-synthesis-creative-20260825c / recovered owner-word lane
audit: >
  PASS. A primary owner-word bitset implementation and a no-import literal
  two-lift implementation independently enumerate each exact-maximum layer,
  all odd residue pairs with repetition, and clocks 15,...,40. They agree on
  every candidate count, primitive/nonprimitive split, first-certificate
  histogram, four clock-38 boundary cores, and semantic digest. Both scripts
  have zero assert nodes and zero float literals; normal and optimized
  outputs are byte-identical.
depends_on:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-4070-lrc14-d2-mod14-two-bank-affine-ray-firewall
related:
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-4052-lrc14-affine-component-width-escape-cones
script: 04-computation/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075.py
output: 05-knowledge/results/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075.out
independent_audit_script: 04-computation/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075_independent_audit.out
script_sha256: d62779e46d8bf64d394ffb61daf598820e2fd4d2125e0da2b16726dd3649b911
output_sha256: 3a5b4b0400a376cc4a61559fda16fe40e1a6313ab7a637f190c5d4b0f57c53ef
independent_audit_script_sha256: b01b61037c6478cf5d3dd01e72243c50cdcf4010e97fcee44887b7fe65204ce0
independent_audit_output_sha256: 3fbf8385f86fe211ab3cf8a70fbaf0a678ec2b0737931ee172da99f5bed7d0bd
hash_basis: raw LF bytes
---

# THM-4075 -- LRC(14) divisor-complete dyadic owner-word closure through 30

**PROVED RELATIVE TO THM-2061/2066 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** The small-denominator residual recovered by THM-4070 is exactly
THM-2061's older divisor-complete dyadic seam. Reconnecting that condition
to THM-2066's owner-word clocks closes six new maximum layers, including the
nonprimitive cores that were outside THM-2066's normalized finite census.

## 1. Inheritance and exact type match

For an eleven-element set `C` of distinct positive integers and distinct
positive odd integers `x,y`, put

```text
S(C;x,y)=2C union {x,y}.                                (1)
```

THM-4041's typed `d=2,c_2=9` row is literally `(1)` with

```text
C=H,                 (x,y)=(alpha,beta).               (2)
```

The eleven divided pack speeds are distinct because the original even pack
speeds are distinct; the two exceptions are distinct and odd. A row is not
`1/14`-lonely exactly when its attained maximum clearance is strictly below
`1/14`, which is THM-2061/2066's strict-counterexample convention. Thus
there is no open/closed or distinctness mismatch in `(2)`. Primitivity of
`H` is not inherited and is not assumed below.

THM-4070 proves that a nonlonely row `(2)` must satisfy

```text
for every d=2,...,14, some c in C is divisible by d.    (3)
```

Condition `(3)` is precisely THM-2061 Section 3. THM-4070 gives a shorter
explicit bank proof and new `q=14` consequences; it does not introduce a
new residual condition. The recovered downstream mechanism is THM-2066's
labelled owner word, not another scalar divisor count.

## 2. The all-height owner-word obstruction

Write

```text
|a|_N=min_{k in Z}|a-kN|,
```

and label every class modulo `N` by its canonical representative
`0<=r<N`. For a clock `N`, define the closed pack-safe labelled packet

```text
A_N(C)={0<=r<N: 14|cr|_N>=N for every c in C}.          (4)
```

For an odd residue `z mod 2N`, call `z` eligible when

```text
7|zr|_N<N                    for every r in A_N(C),     (5)
```

and on the labelled packet put

```text
omega_z(r)=nint(zr/N) mod 2.                            (6)
```

Here `nint` is the unique nearest integer; `(5)` puts `zr/N` within `1/7`
of an integer, so there is no tie. The canonical label is load-bearing:
replacing `r` by `r+N` flips every odd tail's owner bit, although it preserves
complementarity after the common gauge flip.

Let `E_N(C)` be the eligible odd residue classes and define

```text
R_N(C)={(u,v) in E_N(C)^2:
          omega_v(r)=1-omega_u(r) for every r in A_N(C)}.
```

Repeated residue classes are retained in this finite relation because two
distinct odd tails may coincide modulo `2N`.

THM-2066 proves that two odd tails spoil both lifts above every point of
`A_N(C)` exactly when both are eligible and their words in `(6)` are
pointwise complementary. Consequently every strict counterexample must have
a complementary eligible pair at every clock. In particular,

```text
R_N(C)=empty for some N  ==>  S(C;x,y) is 1/14-lonely
for every two odd tails x,y.                            (7)
```

This criterion does not require `gcd(C)=1`. It is an all-height rational
bank sidecar retaining tail residues and sheet ownership, strictly stronger
than `(3)`.

## 3. Complete exact extension through maximum 30

Every eleven-element `C subset {1,...,30}` satisfying `(3)` closes by `(7)`.
For the already proved range `max(C)<=24`, THM-2066 enumerates `59,880`
primitive divisor-complete cores. Divisor completeness actually makes
primitivity automatic there: `13 in C`, and a common divisor greater than
one would therefore be `13`, impossible for eleven distinct numbers at
height at most `24`.

For each new exact maximum `M`, both companions enumerate

```text
C=C_0 union {M},       C_0 in binom({1,...,M-1},10),  (8)
```

filter by all thirteen conditions `(3)`, and test clocks `15,...,40` in
increasing order. The complete counts are

```text
M    searched       divisor-complete   primitive   nonprimitive   certified
25    1,961,256          22,953           22,953          0          22,953
26    3,268,760         112,241          112,207         34         112,241
27    5,311,735         149,354          149,354          0         149,354
28    8,436,285         457,165          457,062        103         457,165
29   13,123,110         252,939          252,939          0         252,939
30   20,030,010         769,704          769,452        252         769,704.   (9)
```

Thus all `1,764,356` new cores close, including all `389` nonprimitive
cores. With THM-2066's earlier range, this is `1,824,236` closed cores. Every
new core has a certificate by `N=38`; clocks `39,40` are unused.

The endpoint of the ordered clock bank is genuine. Exactly four cores first
certify at `N=38` after surviving every clock `15,...,37`:

```text
{1,3,11,23,24,25,26,27,28,29,30},
{1,3,22,23,24,25,26,27,28,29,30},
{2,3,11,23,24,25,26,27,28,29,30},
{2,3,22,23,24,25,26,27,28,29,30}.                    (10)
```

This is necessity only relative to the ordered bank `15,...,38`, not a
claim that clock `38` is universally minimal among all possible clocks.

## 4. Independent literal-lift audit

The primary companion computes `(4)--(6)` as bitsets and tests word
complementarity. The independent companion never forms a nearest-integer
word. For each odd `z mod 2N`, packet residue `r`, and label `j in {0,1}` it
computes the literal strict-danger bit

```text
L_(z,j)(r)=1[14|z(r+jN)|_(2N)<2N].                    (11)
```

It then exhausts unordered odd residue pairs with repetition and declares a
clock certificate exactly when no pair's two masks cover `A_N(C)` on both
labels. Formula `(11)` evaluates the actual times

```text
(r+jN)/(2N),                                           (12)
```

so it independently audits eligibility, parity ownership, strictness, and
residue repetition. Both paths obtain the same first-certificate histogram
in every row of `(9)`, the same four cores `(10)`, and the same semantic
digest.

As a small explicit residual killed beyond the divisor sieve, take

```text
C_*={1,2,3,4,8,9,10,11,12,13,14}.                    (13)
```

It satisfies `(3)`. At `N=23`,

```text
A_23(C_*)={4,9,10,13,14,19}.                          (14)
```

The only eligible odd class modulo `46` is `23`, whose owner word in the
displayed order is `(0,1,0,1,0,1)`. There is no complementary word. Hence
the twelve explicit times `(r+23j)/46`, with `r` in `(14)` and `j=0,1`,
close `2C_* union {x,y}` for every odd pair. This shows concretely why
divisor completeness is not the terminal obstruction.

## 5. Consequence and honest physical scale

Combining `(2)`, THM-4070, and Section 3 gives the typed necessary condition

```text
nonlonely THM-4041/4070 residual  ==>  max(H)>=31.      (15)
```

This is a genuine owner-word strengthening of the divisor pins, but it is
not a competitive numerical height gain in the full THM-3818 producer.
There, let `b` be any divided body-pack owner and `c` either divided pair
owner. The primitive two-speed crossing relation between actual speeds
`2b,2c` has coefficient height

```text
max(b,c)/gcd(b,c).                                    (16)
```

THM-3818's unresolved `W=V_dec` branch forbids every crossing row of height
at most

```text
Q=91^6=567869252041.                                  (17)
```

Thus `(16)>Q` for every such cross-component pair, and in particular
`max(H)>Q`. The useful content of this theorem is therefore the reusable
clock obstruction `(7)` and the complete abstract seam closure, not the
small physical floor `(15)`. Rows above height `30`, the all-height
intersection problem, the other component shapes, and LRC(14) remain open.

## 6. Replay

```text
python3 -B 04-computation/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075.py
python3 -B -O 04-computation/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075.py
python3 -B 04-computation/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075_independent_audit.py
python3 -B -O 04-computation/lrc14_divisor_complete_dyadic_owner_word_closure_through30_thm4075_independent_audit.py
```

All four runs reproduce the frozen raw-LF outputs. **QED.**
