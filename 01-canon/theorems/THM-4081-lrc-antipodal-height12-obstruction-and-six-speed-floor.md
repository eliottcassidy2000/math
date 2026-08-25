---
id: THM-4081
title: "LRC antipodal height-twelve obstruction and six-speed floor"
status: >
  PROVED + PROVED RELATIVE TO CITED LOWER-CASE LRC(7) + VERIFIED-EXACT +
  INDEPENDENTLY VERIFIED-EXACT. The doubled-circle parity reduction is exact.
  D*={1,3,4,5,7,8,9,11,12} is an inclusion-minimal antipodal obstruction,
  and a subset S of {1,...,12} obstructs iff D* is contained in S. Odd
  dilations preserve obstruction and minimality; even dilation does not.
  Relative only to the cited six-speed LRC theorem, the global minimum
  obstruction cardinality kappa satisfies 7<=kappa<=9. The cases kappa=7,8
  remain OPEN.
source: codex-frontier-synthesis-creative-20260825e / antipodal-cover anchor
audit: >
  PASS. The primary Fraction-exact transformed-circle audit constructs 36
  closed odd-layer intervals, checks component counts 1,3,5,5,5,5 and all
  five strict even-cover maxima, verifies 99 strong-deletion survivor gates
  plus nine strict failures, classifies all 4,096 subsets of {1,...,12}, and
  retains the six isolated endpoint solutions of an eight-speed hostile.
  The independent literal-theta audit constructs every equality wall and
  open cell: 185+185 tests for D*, 221+221 tests for the full height-twelve
  universe, and 301+301 tests for the endpoint hostile; it checks 6,660,
  10,608, and 9,632 literal phase inequalities respectively, plus 198 direct
  deletion-survivor inequalities. Normal and optimized outputs byte-match;
  both scripts have zero assert nodes and zero floating literals.
depends_on:
  - LRCUpTo13 # used only through its classical lower-case LRC(7) consequence
related:
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists
  - THM-526-lrc-arc-width-lemma-large-stranger-discharge
  - MISTAKE-126
  - MISTAKE-273
  - MISTAKE-274
  - MISTAKE-382
  - MISTAKE-464
script: 04-computation/lrc_antipodal_height12_obstruction_thm4081.py
output: 05-knowledge/results/lrc_antipodal_height12_obstruction_thm4081.out
script_sha256: 92edbb90e237f079a0694266b6ffe12d3ddab5dcbe5325e21f799dea19147d7f
output_sha256: 8dbb911cfeb358833d22fed6c2e079a3ad917b28caf1325c7d99418f02967cde
independent_audit_script: 04-computation/lrc_antipodal_height12_obstruction_thm4081_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_antipodal_height12_obstruction_thm4081_independent_audit.out
independent_audit_script_sha256: 9f750af317ca89e4a645923488b41a9fc43d48d4746079ee3ea91746044bcc9b
independent_audit_output_sha256: cda0903b58ebfc96aa4f18ae195e96c4202c2fbecff647b9f2bdcd078dd8a8ef
hash_basis: raw LF bytes
---

# THM-4081 -- the first-height antipodal obstruction and the six-speed floor

Put `delta=1/14`. For a finite set `D` of positive integer speeds, call
`theta in R/Z` an **antipodal safe point** when

```text
||d theta||>=delta and ||d(theta+1/2)||>=delta   for every d in D. (1)
```

Call `D` an **antipodal obstruction** when no such `theta` exists, and let

```text
kappa=min{|D|: D is an antipodal obstruction}.                       (2)
```

The explicit obstruction in this theorem is

```text
D*={1,3,4,5,7,8,9,11,12}
  ={1<=d<=12: d is not congruent to 2 modulo 4}.                       (3)
```

## 1. Exact parity reduction

Set `t=2 theta`. The elementary identity

```text
min(||x||,||x+1/2||)=(1/2)||2x||                              (4)
```

gives the exact equivalence

```text
d odd:      (1) for d  iff  ||d t||>=1/7,
d=2e even:  (1) for d  iff  ||e t||>=1/14.                    (5)
```

No endpoint is discarded: every inequality in `(1)` and `(5)` is weak, so
all safe sets below are closed. In particular, if odd `o` and its double
`2o` are both present, the `2o` condition is a strict-threshold shadow of the
`o` condition: both use the multiplier `o`, but `1/14<1/7`.

This explains the residue anatomy in `(3)`. Inside `{1,...,12}`, the speeds
`2,6,10` shadow `1,3,5`; the non-shadow core consists of every odd speed and
the multiples of four `4,8,12`.

## 2. Symbolic interval obstruction

For odd `o`, its safe set in the `t` circle is exactly

```text
S_o = union over 0<=k<o of
      [(7k+1)/(7o),(7k+6)/(7o)].                              (6)
```

Intersecting `(6)` successively for `o=1,3,5,7,9,11` gives the following
closed unions. This display is also an endpoint certificate; singleton
intersections are retained throughout.

```text
after 1:
  [1/7,6/7]

after 3:
  [1/7,2/7] U [8/21,13/21] U [5/7,6/7]

after 5:
  [1/7,6/35] U [8/35,2/7] U [3/7,4/7]
  U [5/7,27/35] U [29/35,6/7]

after 7:
  [8/49,6/35] U [8/35,13/49] U [22/49,27/49]
  U [36/49,27/35] U [29/35,41/49]

after 9:
  [8/49,6/35] U [5/21,13/49] U [29/63,34/63]
  U [36/49,16/21] U [29/35,41/49]

after 11:
  I_1=[8/49,13/77],       I_2=[5/21,20/77],
  I_3=[36/77,41/77],      I_4=[57/77,16/21],
  I_5=[64/77,41/49].                                         (7)
```

The even speeds `12,8,4,8,12`, respectively, strictly kill the five closed
intervals in `(7)`. Their transformed multipliers are `6,4,2,4,6`; the exact
maximum danger norms on the whole intervals are

| interval | killing speed | maximum transformed norm | strict comparison |
|---|---:|---:|---:|
| `I_1` | `12` | `1/49` | `<1/14` |
| `I_2` | `8` | `1/21` | `<1/14` |
| `I_3` | `4` | `5/77` | `<1/14` |
| `I_4` | `8` | `1/21` | `<1/14` |
| `I_5` | `12` | `1/49` | `<1/14` |

These are maxima over the entire closed intervals, not endpoint samples.
They follow by splitting `||mt||` at its integral and half-integral break
points; the primary audit performs exactly that symbolic rational operation.
Thus every point surviving all six odd speeds fails one of `4,8,12`, proving

```text
D* is an antipodal obstruction.                               (8)
```

## 3. Strong deletion witnesses and the complete height-twelve classification

The obstruction is strongly inclusion-minimal. For every `d in D*`, the row
below gives `t_d` that is safe not merely for `D*\{d}` but for the larger set
`{1,...,12}\{d}`. The omitted speed strictly fails. The `theta` column is the
literal original-circle witness `theta_d=t_d/2`.

| omitted `d` | `t_d` | `theta_d` | omitted transformed norm | a tight survivor |
|---:|---:|---:|---:|---:|
| `1` | `1/14` | `1/28` | `1/14<1/7` | `2` at `1/14` |
| `3` | `15/49` | `15/98` | `4/49<1/7` | `7` at `1/7` |
| `4` | `36/77` | `18/77` | `5/77<1/14` | `11` at `1/7` |
| `5` | `8/21` | `4/21` | `2/21<1/7` | `3` at `1/7` |
| `7` | `1/7` | `1/14` | `0<1/7` | `1` at `1/7` |
| `8` | `5/21` | `5/42` | `1/21<1/14` | `9` at `1/7` |
| `9` | `8/35` | `4/35` | `2/35<1/7` | `5` at `1/7` |
| `11` | `29/63` | `29/126` | `4/63<1/7` | `9` at `1/7` |
| `12` | `8/49` | `4/49` | `1/49<1/14` | `7` at `1/7` |

Direct substitution in `(5)` checks all `9*11=99` survivor inequalities.
The table deliberately records a tight survivor in every row: replacing
`>=` by `>` would invalidate the certificate. The independent audit instead
checks the `9*11*2=198` literal inequalities in `(1)`, so it does not inherit
the doubled-coordinate reduction.

Consequently, for **every** subset `S of {1,...,12}`,

```text
S is an antipodal obstruction  iff  D* is a subset of S.      (9)
```

The forward direction is immediate from a strong deletion witness: if `S`
misses `d in D*`, then `theta_d` is safe for `S`. The reverse direction
follows from `(8)` by monotonicity. Therefore the complete list consists of

```text
D* union E,       E any subset of {2,6,10}.                  (10)
```

There are exactly eight such sets, with cardinality histogram

```text
size 9: 1,     size 10: 3,     size 11: 3,     size 12: 1.   (11)
```

In particular `D*` is the unique inclusion-minimal obstruction of height at
most twelve, the unique nine-speed obstruction in that universe, and no
obstruction has maximum speed at most eleven. Hence the least possible
obstruction height is exactly twelve, attained minimally only by `D*`.

## 4. Odd-dilation symmetry and the even-dilation boundary

Let `q` be odd. Parity is preserved, and the multiplier and threshold in `(5)`
satisfy

```text
m(qd)=q m(d),                  rho(qd)=rho(d),                (12)
```

where `rho` is `1/7` on odd speeds and `1/14` on even speeds. Thus `t` is safe
for `qD` exactly when `qt` is safe for `D`. Multiplication by `q` is onto on
`R/Z`, so

```text
qD obstructs iff D obstructs.                                (13)
```

Moreover `t_d/q` lifts every deletion witness. Hence `qD*` is
inclusion-minimal for every positive odd `q`.

Even dilation is genuinely different. For `2D*`, take

```text
t=1/13,                    theta=1/26.                        (14)
```

Every speed `2d` is even, its multiplier is `d`, and for `1<=d<=12`

```text
||d/13||>=1/13>1/14.                                           (15)
```

Thus `2D*` has an antipodal safe point. Obstruction is not invariant under
even dilation; this is the precise parity boundary.

## 5. Cited six-speed floor and the global cardinality window

This paragraph alone uses an external citation. For an arbitrary finite speed
set `D`, deduplicate the transformed multipliers

```text
A(D)={d: d in D odd} union {d/2: d in D even}.                (16)
```

If `|D|<=6`, then `|A(D)|<=6`. The cited classical lower-case `LRC(7)` theorem
(six nonzero speeds at separation `1/7`, routed in the repository through the
`LRCUpTo13` citation node and THM-526's Barajas--Serra reference) supplies
`t in R/Z` such that

```text
||a t||>=1/(|A(D)|+1)>=1/7       for every a in A(D).         (17)
```

Equation `(5)` then makes every odd speed safe at its required `1/7` threshold
and every even speed safe at the weaker `1/14` threshold. Hence no set of at
most six speeds is an antipodal obstruction. Combining this cited floor with
`D*` gives the global bound

```text
boxed: 7<=kappa<=9.                                             (18)
```

This theorem does **not** decide whether `kappa` is `7`, `8`, or `9`. In
particular, the exact height-twelve classification must not be globalized into
a claim that every eight-speed set has an antipodal safe point.

## 6. Endpoint-only eight-speed hostile

The warning is substantive. The eight-speed set

```text
H_8={1,3,5,8,11,13,23,36}                                    (19)
```

is not an obstruction, but its transformed safe set is exactly

```text
{1/7,2/7,3/7,4/7,5/7,6/7}.                                  (20)
```

All six components in `(20)` are isolated closed points; there is no positive
length safe interval. On the original `theta` circle their two lifts are

```text
{k/14: k=1,2,3,4,5,6,8,9,10,11,12,13}.                      (21)
```

The primary closed-interval engine obtains `(20)` exactly. Independently, the
literal-theta arrangement finds twelve safe equality walls and zero safe open
cells. Any audit that samples only cell interiors or midpoints would therefore
misclassify `H_8` as an eight-speed obstruction.

## 7. Two independent exact audits

The primary script uses `(5)` and closed `Fraction` interval algebra. It
reconstructs every stage of `(7)`, maximizes the five even danger functions on
whole intervals, verifies the strong witnesses, and independently enumerates
all `2^12=4096` subsets to recover `(9)`--`(11)`. It also retains the six
singleton components of `(20)` and checks exact odd- and even-dilation
certificates.

The second script does not use `t=2theta`, transformed multipliers, or the
interval formulas `(5)`--`(7)`. For each literal inequality it constructs the
walls

```text
theta=(k +/- 1/14)/d-h/2 mod 1,       h in {0,1},             (22)
```

tests every wall, and tests one exact midpoint in every complementary open
cell. Between consecutive walls all literal inequality signs are constant,
so this is exhaustive. It then classifies a subset by OR-ing its speed-failure
masks at all wall/cell tests. The independent path again obtains exactly the
eight masks in `(10)`, the nine strong deletion certificates, and the twelve
endpoint-only solutions in `(21)`.

Both outputs are frozen beside their scripts. No floating arithmetic, solver,
randomness, or positive-width assumption is used. No height-120 or global
eight-speed claim is part of this theorem.
