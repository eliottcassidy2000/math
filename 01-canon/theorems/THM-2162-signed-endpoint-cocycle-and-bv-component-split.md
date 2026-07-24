---
id: THM-2162
title: "Signed endpoint cocycle, Euler/BV component split, and exact one-swap stability"
status: >
  PROVED (general endpoint identity and analytic tail) + FINITE-EXACT (232
  far replacements, 12 small replacements, two independent exact interval
  paths for the equality row). For a finite union of positive-length
  intervals, adding one danger comb has discrepancy equal to a signed
  endpoint cocycle of a continuous periodic primitive. At danger density p,
  only positive-length components pay the sharp oscillation p(1-p)/W;
  isolated weak-safe points cancel identically. Applied at p=1/7, the
  drop-6 OPEN-Q-108 extremizer is a strict local minimizer under every proper
  one-element replacement. The exact neighboring floor is
  3859/420420=7/858+1/980, uniquely at replacement 10->20. This corrects the
  false incoming claim that the next local value was 426/35035. The theorem
  is local around the extremizer and does not prove global OPEN-Q-108 or
  LRC(14).
source: opus-2026-07-24-puzzle-atlas
depends_on: []
related:
  - THM-541-lrc14-ap-window-single-hole-collar
  - THM-816-quartic-order-three-s3-uniform-looseness
  - THM-2047-phase-height-toric-arrangement-for-lrc
  - THM-2048-fiber-quantization-tax-for-lrc-discrepancy
script: 04-computation/lrc_q108_signed_endpoint_ramanujan_scout_opus_20260724.py
output: 05-knowledge/results/lrc_q108_signed_endpoint_one_swap_thm2162_opus_20260724.out
script_sha256: 25f7e8c41d59a3f2c86b78c25cf40a93862dad1416819009e32f2008c4b764e4
output_sha256: 970cc92d1a401eb67ba1eb5e507fe300e5245f1f47e9c12f476487b3f44578a7
---

# THM-2162 -- the signed endpoint cocycle and exact local stability

Put

```text
theta=1/14,
D_W={t in R/Z: ||Wt||<theta},
S_W=(R/Z)\D_W,
G(C)=intersection_(v in C) S_v.                       (1)
```

The first theorem separates two component notions which must not be
interchanged:

```text
Euler components: positive intervals and isolated weak-safe points;
BV components:    positive-length intervals only.     (2)
```

THM-2047 needs the first list to detect a weak lonely witness. Measure
discrepancy needs the second.

## 1. Exact endpoint cocycle at arbitrary density

More generally let `D` be a centered circular interval of measure
`0<p<1`, put

```text
h=1_D-p,
```

and let `H` be a continuous one-periodic primitive of `h`. It may be
normalized so that

```text
max H-min H=p(1-p).                                  (3)
```

Indeed `H` rises with slope `1-p` across an interval of length `p`, and
falls with slope `-p` across its complement.

Let `G` be a finite union of `K` positive-length circular intervals and
finitely many isolated points. Choose interval representatives `[l_i,r_i]`,
splitting at zero if necessary. For every positive integer `W`, define

```text
epsilon_W(G)=|G intersection {t:Wt mod 1 in D}|-p|G|.
```

Then exactly

```text
epsilon_W(G)
 =1/W sum_(i=1)^K [H(Wr_i)-H(Wl_i)],                 (4)

|epsilon_W(G)|<=K p(1-p)/W.                          (5)
```

**Proof.** On one interval, substitute `x=Wt` and use `H'=h` away from
the two harmless threshold points:

```text
integral_[l,r] h(Wt)dt=[H(Wr)-H(Wl)]/W.
```

Sum over the positive intervals. Isolated points have measure zero. If an
implementation represents one as `[x,x]`, its two endpoint terms cancel
literally:

```text
H(Wx)-H(Wx)=0.                                       (6)
```

Finally each difference in (4) is at most the oscillation (3) in absolute
value. This proves (4)--(5). The constant is sharp for one arbitrary
interval because its endpoint phases may approach a maximum and a minimum
of `H`. QED.

At `p=1/7`, take `H(0)=0`. On `[0,1]`,

```text
H(x)=
  (6/7)x,                              0<=x<=1/14,
  1/14-x/7,                            1/14<=x<=13/14,
  -3/49+(6/7)(x-13/14),                13/14<=x<=1.
                                                               (7)
```

Its range is `[-3/49,3/49]`. Therefore, if `G(C)` has measure `mu`
and `K` positive-length components, (4) gives the exact signed formula

```text
|G(C union {W})|
 =6mu/7-epsilon_W(G(C)),                              (8)
```

and the corrected one-far bound

```text
|G(C union {W})|>=6mu/7-6K/(49W).                    (9)
```

The older estimate `6mu/7-N/(7W)` pays a larger primitive variation and,
when `N` is obtained by `len` of a closed-interval list, may also charge
isolated points. Both losses are removed in (9).

## 2. The binding core has six topological points but only eight BV arcs

Let

```text
A={1,2,3,4,5,7,8,9,10,11,12,13}.                    (10)
```

This is the drop-6 extremizer from THM-541, with

```text
|G(A)|=7/858.                                        (11)
```

For `v in A`, put `C_v=A\{v}`. The empirically binding eleven-core is

```text
C_10={1,2,3,4,5,7,8,9,11,12,13}.
```

Exact interval arithmetic gives

```text
|G(C_10)|=313/9702.                                  (12)
```

Its closed interval list has fourteen entries, but six are the isolated
unit fractions

```text
1/14,3/14,5/14,9/14,11/14,13/14.                    (13)
```

Thus

```text
Euler component count=14,
positive-length BV component count=8.                (14)
```

Formula (9) charges `48/(49W)`, not `14/(7W)`. Relative to the target
`7/858`, this closes the core analytically for every `W>=51`; the old
count would only close at `W>=103`.

This is not a cosmetic convention. The points (13) are indispensable for
the weak-safe Euler predicate but invisible to every measure integral.

## 3. Exact one-swap stability of the sharp extremizer

**Theorem.** If `S` is a proper one-element replacement of `A`, namely

```text
S=(A\{v}) union {w},             v in A, w notin A,  (15)
```

then

```text
|G(S)|>=3859/420420
       =7/858+1/980.                                  (16)
```

Equality holds only for

```text
v=10, w=20,
S_*={1,2,3,4,5,7,8,9,11,12,13,20}.                  (17)
```

Consequently `A` is a strict local minimizer in the one-swap graph, and
its exact local stability gap is `1/980`.

### Analytic tail

Write

```text
q_*=3859/420420.
```

For each `v`, exact interval arithmetic supplies the core mass `mu_v`,
the number `K_v` of positive-length components, and the number `I_v` of
isolated points. Formula (9) proves `|G(C_v union {W})|>q_*` whenever

```text
W>A_v:=[6K_v/49]/[6mu_v/7-q_*].                       (18)
```

All denominators in (18) are positive. The exact ledger is:

| removed `v` | `mu_v` | `K_v` | `I_v` | `A_v` |
|---:|---:|---:|---:|---:|
| 1 | `239/3003` | 6 | 4 | `308880/24821` |
| 2 | `461/12012` | 6 | 6 | `23760/767` |
| 3 | `1159/18018` | 8 | 4 | `411840/19321` |
| 4 | `389/12012` | 8 | 6 | `411840/7811` |
| 5 | `379/8580` | 12 | 4 | `617760/12059` |
| 7 | `4739/51480` | 10 | 6 | `257400/14657` |
| 8 | `2243/42042` | 8 | 6 | `2882880/107567` |
| 9 | `11191/168168` | 10 | 4 | `900900/35213` |
| 10 | `313/9702` | 8 | 6 | `2882880/54367` |
| 11 | `1951/25480` | 12 | 4 | `540540/20767` |
| 12 | `1546/35035` | 6 | 6 | `2162160/84299` |
| 13 | `907/17640` | 14 | 4 | `315315/6418` |

Thus only the integer windows

```text
14<=W<=floor(A_v)                                    (19)
```

remain. They contain 232 rows in total.

### Finite exact window and the small replacement

The companion computes all 232 rows with `Fraction` arithmetic. For every
row it:

1. intersects the eleven core safe-comb interval lists;
2. evaluates the signed endpoint cocycle (4), retaining all active owners;
3. intersects with the `W`-safe comb directly; and
4. checks exact agreement of the two values.

Every row satisfies (16), and equality occurs only at `(v,W)=(10,20)`.

Since `w` in (15) is a positive integer not in `A`, either `w=6` or
`w>=14`. The twelve `w=6` rows are the other AP-drop rows. An exact replay
finds their minimum

```text
426/35035,                  at v=12,                  (20)
```

which is strictly larger than `q_*`. This completes the proof of (16).

As an independent path for the equality classifier, the companion clips
and merges every open danger interval for (17), then subtracts its rational
union length from one. The seven merged danger components again give

```text
|G(S_*)|=3859/420420.                                 (21)
```

Normal and optimized Python executions match the stored transcript byte for
byte.

## 4. Correction and scope

The exact inequalities are

```text
7/858 < 3859/420420 < 426/35035,

3859/420420-7/858=1/980,
426/35035-3859/420420=179/60060.                      (22)
```

Therefore the incoming claim that the two smallest global twelve-core
measures were `7/858` and `426/35035` was false: its bounded scan did not
contain the speed `20`. What survives is stronger and correctly typed:
`7/858` is an isolated strict minimum on the entire one-swap neighborhood,
with the new exact local gap (22).

This theorem does **not** classify cores at swap distance two or more and
does not prove that (16) is the global second value. It does not prove
OPEN-Q-108 or LRC(14).

The supplied puzzle lenses clarify why this carrier works:

- the connected-sum “common intermediate” move says to retain the whole core
  union before charging a new operation; (8) does exactly that;
- Ramanujan/Fourier expansion becomes absolutely summable only after the
  primitive `H` contributes the extra frequency denominator, but the endpoint
  cocycle is already the exact finite version;
- partial-cube/Euler topology must retain isolated witnesses, while BV
  discrepancy must quotient them out;
- a permanent, Asano contraction, or unsigned norm would forget the signs in
  (4), which are the cancellation being measured.

The reusable lesson is not “discard isolated points.” It is:

```text
retain them for existence topology; cancel them for measure variation.
```

QED.
