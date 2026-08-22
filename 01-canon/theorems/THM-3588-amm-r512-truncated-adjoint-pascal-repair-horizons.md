---
id: THM-3588
title: "AMM R=512 truncated-adjoint Pascal repair horizons"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the five
  failing golden-floor Rule-A
  prefixes at R=512 and offsets D0=0,...,4, positive integer
  truncated-adjoint Pascal multipliers give exact Farkas certificates forcing
  every surviving alternative to differ from Rule A by rows
  (35,38,39,42,44), respectively.  These improve THM-3577's capacity bounds
  (49,52,53,56,59).  This is a fixed-epoch prefix obstruction, not a
  feasible alternative, a uniform extractor, or a value of the AMM constant.
source: kps-s188 / exact top-light-cone attack, 2026-08-21
audit: >
  ACCEPT.  The difference-inequality direction, fatal-top sign,
  truncated-adjoint recursion, cut indexing, exact D0=4 path sum, boundary
  signs, multiplier ledger, replays, and hashes were independently rebuilt.
  The positive preceding cuts remain controls, not feasibility claims.
depends_on:
  - THM-3577-amm-r512-offset-transition-and-causal-horizon
related:
  - THM-3374-amm-r512-rule-a-causal-repair-horizon
  - THM-3371-r8-integer-entry-section-and-hostile-face-separation
script: 04-computation/amm12592_r512_truncated_adjoint_pascal_horizons_thm3588.py
output: 05-knowledge/results/amm12592_r512_truncated_adjoint_pascal_horizons_thm3588.out
script_sha256: 692ac6411f9a99c2734fe85995915b144a9752528816bae6df73cd0943c53985
output_sha256: a71c4b1a0503b6ee79f5bde7c0d0d4448cb34c0e7f0a245e26bd11d8b5c5bdd8
hash_basis: LF-normalized bytes
---

# THM-3588 -- AMM R=512 truncated-adjoint Pascal repair horizons

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The crude
causal-capacity triangle inequality
can be replaced by an exact positive dual functional on the fatal top cell's
backward Pascal light cone.  At the five failing offsets this moves the last
possible Rule-A agreement row fourteen or fifteen rows earlier.

## 1. Statement

Use the exact `R=512` golden-floor degree profiles

```text
d_i(D0)=floor(log_5(phi^(2(512+i))))+D0.                (1)
```

THM-3577 proves that Rule A first dies at rows

```text
j=(107,110,113,116,121)                                (2)
```

for offsets `D0=0,1,2,3,4`.  Its coefficient-capacity argument showed only
that any surviving alternative must differ from Rule A by rows
`(49,52,53,56,59)`.

The stronger exact conclusion is

| `D0` | death row `j` | first negative dual cut `s` | every survivor differs by row | old bound | gain |
|---:|---:|---:|---:|---:|---:|
| 0 | 107 | 36 | `<=35` | `<=49` | 14 |
| 1 | 110 | 39 | `<=38` | `<=52` | 14 |
| 2 | 113 | 40 | `<=39` | `<=53` | 14 |
| 3 | 116 | 43 | `<=42` | `<=56` | 14 |
| 4 | 121 | 45 | `<=44` | `<=59` | 15 |

Thus no repair confined to the final `72,72,74,74,76` causal rows,
respectively, can even survive the Rule-A death cell.

## 2. Halved causal recurrence

Write `y_(i,t)` for the halved junk coefficient and `u_(i,t)` for the
halved correction chosen at row `i`, cell `t`.  If

```text
epsilon_i=d_i-d_(i-1) in {0,1},
K_i(x)=(1+x)^(1+epsilon_i),                            (3)
```

then the recurrence before truncation is

```text
y_i=K_i y_(i-1)+f_i-u_i.                              (4)
```

The independent Lucas-box coordinates satisfy

```text
L_(i,t)=-binom(d_i-1,t),
U_(i,t)= binom(d_i-1,t-1),
L_(i,t)<=u_(i,t)<=U_(i,t).                            (5)
```

The top cell has no state variable.  Survival requires its new load to lie
in `[0,1]`.

Let `y^A,u^A` be the exact Rule-A prefix and let another admissible prefix
agree with it through row `s-1`.  Put

```text
z_i=y_i-y_i^A.                                        (6)
```

For every state cell at rows `s<=i<j`, equations `(4)--(5)` give the valid
one-sided inequality

```text
z_(i,t)-(K_i z_(i-1))_t
   =u^A_(i,t)-u_(i,t)
   <=u^A_(i,t)-L_(i,t).                               (7)
```

At `i=s` the previous difference is zero.  At the fatal top cell, the lower
survival bound gives

```text
-(K_j z_(j-1))_(d_j) <= load^A_(j,d_j).               (8)
```

The right side of `(8)` is the negative fatal Rule-A load.

## 3. The truncated adjoint multiplier

Start at the death cell and propagate backwards.  If

```text
K_i(x)=sum_q k_(i,q)x^q,                               (9)
```

define nonnegative integer weights on actual state cells by

```text
lambda_(j-1,t)=k_(j,d_j-t),                           (10)

lambda_(i-1,u)=sum_t lambda_(i,t) k_(i,t-u),          (11)
```

where terms outside `0<=u<d_(i-1)` are discarded.  This truncation is the
essential improvement over the scalar capacity bound: Pascal ancestry that
leaves the actual state width cannot help repair the death cell.

Multiply `(7)` by `lambda_(i,t)`, sum over `s<=i<j`, and add `(8)`.  Equations
`(10)--(11)` cancel every coefficient of every free `z_(i,t)` exactly.
All multipliers are nonnegative integers, so one obtains the Farkas
inequality

```text
0 <= B_s,                                              (12)

B_s=load^A_(j,d_j)
    +sum_(i=s)^(j-1) sum_t
       lambda_(i,t)(u^A_(i,t)-L_(i,t)).                (13)
```

There is no floating-point optimization in `(12)--(13)`.  It is a positive
integer combination of literal causal inequalities whose left side is zero.

## 4. Five exact sign walls

For each offset, the value immediately before the sign wall is positive and
the next value is negative:

```text
D0=0: B_35=
 41086809435714651554305675711968187478323158540573926649141693358061153498971915247,
       B_36=
-76874861671160752918838113835252414677773118569312018596735213281104563872967120639.

D0=1: B_38=
 132432154406083499404833844855611470852311846908930910441711037452812262877844835930,
       B_39=
-67299131366770904358812237887063009272989246608815495819036140457947571682105501209.

D0=2: B_39=
 3338985983972956972751398645882828535963003571887262875884142255942013894214739377351,
       B_40=
-9349426258946131715675734429741281552564767664183924828959904668783436277847131542484.

D0=3: B_42=
 12585302820800871834254225868388244300005120240474964547720968195252363811275385411024,
       B_43=
-6065349518170092926836849962397980226126874773437073169325318976224259558021484964106.

D0=4: B_44=
 7877782680146233652648282350697035268590370134970386735210385525966333032831210948230500,
       B_45=
-552094518178517855711445625652899715444284968095857944791483565103784726945679420728711.
                                                               (14)
```

The negative members contradict `(12)`.  For example, an alternative that
agreed through row `44` at `D0=4` would have to satisfy `0<=B_45<0`.
Therefore it must already differ by row `44`.  The other four rows follow
identically.

The preceding positive values are sharp controls for this particular
adjoint certificate family.  They do **not** construct a feasible prefix at
the earlier cut and do not show that the new bounds are globally optimal.

## 5. Mechanism and scope

THM-3577 bounded the possible repair at the fatal coordinate by adding the
absolute capacities of an entire backward triangle.  The new functional is
the dual operation: it transports the fatal top evaluation backwards by the
adjoint Pascal kernels, clips it at the real row widths, and charges only the
available lower-box slack encountered by that transported functional.  Its
exact cancellation is why it can beat the triangle inequality without
solving the full entry polytope.

This theorem proves no:

- surviving alternative prefix at any offset;
- infeasibility when changes are allowed before the displayed row;
- first-feed-free `F1--F3` entry certificate;
- extractor uniform in `R`; or
- bound, lower bound, or value for the uniform AMM constant `C*`.

## 6. Reproduction

Run

```bash
python3 04-computation/amm12592_r512_truncated_adjoint_pascal_horizons_thm3588.py
python3 -O 04-computation/amm12592_r512_truncated_adjoint_pascal_horizons_thm3588.py
```

The standard-library companion reconstructs the exact golden floors, feeds,
Rule-A rows and five fatal loads.  At every cut it constructs the positive
integer multipliers, sums their literal right sides, and verifies by a sparse
integer coefficient ledger that the weighted left side is identically zero.
It checks that each offset has exactly one sign wall, pins both boundary
integers, and replays identically with optimization enabled.

**QED.**
