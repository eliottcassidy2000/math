---
id: THM-4136
title: "Fixed-body universal odd-tail LRC(14) completion"
status: >
  PROVED ELEMENTARY UNIVERSAL FIXED-BODY COMPLETION + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. Let U=(1,4,6,8,10,12,14,15,16,18,22). For every
  two distinct positive odd integers a,b, the thirteen distinct speeds
  2U union {a,b} have a common phase of clearance at least 1/14. This
  strictly generalizes THM-4132 from tails {t,9t} to the whole odd-tail
  parity class of this fixed body. It proves no arbitrary-body entry or
  LRC(14).
source: codex-frontier-synthesis-creative-20260825u
depends_on:
  - THM-4129-universal-two-speed-completion-of-the-eleven-speed-lrc14-body
related:
  - THM-3910-lrc14-auxiliary-center-erosion-and-t-sheet-variance-response
  - THM-4132-fixed-body-exceptional-scale-two-lrc14-completion
script: 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136.py
output: 05-knowledge/results/lrc14_fixed_body_universal_odd_tail_completion_thm4136.out
independent_audit_script: 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_fixed_body_universal_odd_tail_completion_thm4136_independent_audit.out
script_sha256: 50c40de4756362d05ffe7f6064a0990e062a08b8b317708febd554e780de8c23
output_sha256: d8b01ff355561f69a4b0bcfa3d25fe065eb9b253f1367ef8d5457e31a620c179
independent_audit_script_sha256: 6d1f308fc009e67ad9c9e5ba8c544db56cc09375a778bfe7b454efdf7afbc034
independent_audit_output_sha256: 7a2905f0e33c027db0c61845ace430d9ae80a38a9e5e49d79d86534a350ee418
semantic_sha256: 4e62dbd9b4cd87a7044687d48668449b7129ac11adb29512803ae7d01be61dbf
independent_semantic_sha256: 92b303d4737b146a0ff31de4c13fed920e70531d486495e31793128e5ca4482e
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone exact-Fraction implementation rebuilds the body's 244
  walls and 28 positive components, computes every low quotient component,
  stress-tests the q-tooth bound on 1,053 primitive ratios through q=101,
  verifies all 68 residual ratios against three clocks, and directly solves
  eight full rows including nonprimitive and q=18,009 controls. Normal,
  optimized, and hash-seeded replays byte-match the frozen output.
independent_audit: >
  ACCEPT. A clean-room standard-library referee imports no primary code. It
  joins quotient cells only across strictly active walls, independently
  recovers the low topology and all 68 residual rows, checks 1,047 primitive
  q>=9 ratios through q=101 plus 70 unbalanced and near-diagonal hostiles
  through q=1,001, and verifies the literal thirteen-speed clock clearances.
  Normal, optimized, and two hash-seeded replays byte-match its frozen output.
---

# THM-4136 -- fixed-body universal odd-tail completion

**PROVED ELEMENTARY UNIVERSAL FIXED-BODY COMPLETION + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**

THM-4132 used the two physical lifts of one body-safe quotient phase only for
the exceptional ratio `(1,9)`. The same lift carrier has more structure than
that special row suggests: after removing the common odd scale, every
two-sheet obstruction is a cross-intersection of two danger combs. Its
components shrink with the larger primitive tail. A three-clock bank closes
the finitely many ratios below that shrinkage threshold.

## 1. Statement

Put

```text
U=(1,4,6,8,10,12,14,15,16,18,22),       delta=1/14.       (1)
```

> **Theorem.** For every two distinct positive odd integers `a,b`, the
> thirteen distinct positive speeds
>
> ```text
> 2U union {a,b}                                             (2)
> ```
>
> have a phase `x in R/Z` such that
>
> ```text
> min_(v in 2U union {a,b}) ||vx|| >= 1/14.                 (3)
> ```

Parity separates the two tails from the eleven even body speeds. Interchange
the tails if necessary and write

```text
a=pt,       b=qt,       t=gcd(a,b),       0<p<q,            (4)
```

where `p,q,t` are odd and `gcd(p,q)=1`.

## 2. The body arc and its two physical sheets

THM-4129 proves that the closed interval

```text
J=[33/70,27/56] subset G_U,             |J|=3/280,          (5)
```

is `delta`-safe for all of `U`. Runner `15` uniquely owns its left equality
and runner `4` uniquely owns its right equality. Thus for every `y in J`,
both half-lifts

```text
x_0=y/2,                       x_1=(y+1)/2                 (6)
```

keep all of `2U` safe.

Put `w=ty` in `R/Z`, and let `z` be either solution of `2z=w`. Because
`p,q,t` are odd, the tail packets at the two phases in `(6)` are, up to
interchanging the sheets,

```text
(pz,qz),                   (p(z+1/2),q(z+1/2)).            (7)
```

The possible interchange records the parity of the integer removed when
`ty` is reduced modulo one; it loses no physical sheet.

For odd `r`, define the strict danger comb

```text
D_r={z in R/Z: ||rz||<delta}.                              (8)
```

Let `C_(p,q)` be the set of quotient phases `w` for which both packets in
`(7)` are bad. Equality is safe, so `C_(p,q)` is open.

## 3. Cross-comb component bound

For every odd `r`,

```text
D_r intersect (D_r-1/2)=empty,                            (9)
```

because the two values of `rz` differ by `1/2` while two strict danger radii
sum to only `1/7`. Expanding the two-sheet bad predicate therefore leaves
only the cross terms

```text
B_(p,q)=D_p intersect (D_q-1/2),
B_(p,q)+1/2=D_q intersect (D_p-1/2).                      (10)
```

Doubling identifies these antipodal copies. Hence

```text
C_(p,q)=2B_(p,q).                                         (11)
```

Each component of `B_(p,q)` lies in one shifted `q`-tooth, whose length in
the `z` circle is `1/(7q)`. After doubling, the tooth centres form the full
`1/q` grid and the image teeth have length `2/(7q)<1/q`; distinct teeth do
not concatenate. Consequently every component of `C_(p,q)` has length at
most

```text
beta(p,q) <= 2/(7q).                                      (12)
```

For the primitive ratios with `q<9`, direct strict-wall arithmetic gives

| `(p,q)` | open components of `C_(p,q)` | `beta(p,q)` |
|---:|:---|---:|
| `(1,3)` | empty | `0` |
| `(1,5)` | empty | `0` |
| `(3,5)` | `(13/35,8/21)`, `(13/21,22/35)` | `1/105` |
| `(1,7)` | `(6/49,1/7)`, `(6/7,43/49)` | `1/49` |
| `(3,7)` | `(2/7,15/49)`, `(34/49,5/7)` | `1/49` |
| `(5,7)` | `(20/49,3/7)`, `(4/7,29/49)` | `1/49` |

Combining this table with `(12)` yields the uniform bound

```text
beta(p,q) <= 2/63,                                        (13)
```

sharp at `(p,q)=(1,9)`. Also `0 notin C_(p,q)`: at `z=0` the first packet is
bad but the antipodal odd packet is exactly at distance `1/2`.

## 4. Compact-to-open scale closure

Suppose first that `t>=3` and, contrary to `(3)`, both lifts in `(6)` are
tail-bad for every `y in J`. Then the circular image `tJ` is a compact
connected subset of the open set `C_(p,q)`. If `t|J|>=1`, this image is the
whole circle, contradicting `0 notin C_(p,q)`. If `t|J|<1` but the image
wraps, it again contains zero. Otherwise it is a compact interval inside one
open component, which forces the strict inequality

```text
t|J| < beta(p,q).                                         (14)
```

But already at the first admissible odd scale,

```text
3|J|-2/63=9/280-2/63=1/2520>0.                           (15)
```

Thus `(14)` is impossible for every odd `t>=3`.

It remains to take `t=1`. If `q>=27`, then `(12)` gives

```text
|J|-2/(7q) >= 3/280-2/189=1/7560>0,                      (16)
```

and the same compact-to-open argument closes the row.

## 5. The 68 residual ratios

The only cases not yet covered have `t=1` and `q<=25`. There are exactly
`68` coprime odd pairs `p<q` in that range. The following three literal
physical clocks close all of them.

| condition | rows | phase `x` | body clearance | minimum full-row clearance |
|:---|---:|---:|---:|---:|
| `13 notin {p,q}` | `56` | `89/1176` | `9/98` | `89/1176` |
| one tail is `13`, the other is `3,5,7,9,11,15,17,19,21,23` | `10` | `181/4704` | `15/196` | `15/196` |
| `{p,q}={1,13}` or `{13,25}` | `2` | `431/4480` | `17/224` | `17/224` |

Every displayed full-row minimum is strictly larger than `1/14`. Exact
residue evaluation of the thirteen speeds in every one of the `68` rows
proves the residual statement. Together with Sections 2--4 this proves the
theorem. **QED.**

## 6. Sharp hostile and information ledger

A fixed physical sheet cannot replace the two-sheet carrier. For the sharp
ratio `(p,q)=(1,9)`, take

```text
w=1/9,                         z=w/2=1/18.                (17)
```

The two tail-gap packets are

```text
sheet 0: (1/18,1/2),          sheet 1: (4/9,0).          (18)
```

Tail `1` kills the first sheet and tail `9` kills the second. Thus both the
closed body arc and the sheet choice are load-bearing; the quotient width
`2/63` is attained on this same ratio.

The connection contract is

```text
source:       THM-4129's closed U-safe interval J
target:       every fixed-body row 2U union {a,b}, a,b odd
map:          y -> {y/2,(y+1)/2}, then w=gcd(a,b)y
preserved:    all eleven body clearances and one common physical phase
destroyed:    the parity label of the physical half-lift after reducing w
sidecar:      the cross-comb component C_(p,q), not a single fixed sheet
decisive test: component width versus gcd(a,b)|J|, then three low clocks. (19)
```

## 7. Scope and exact replay

This theorem strictly generalizes THM-4132: every row `2U union {t,9t}` with
odd `t` is the primitive ratio `(1,9)`, but no special ratio is now required.
It closes the whole two-odd-tail parity class over the fixed body `2U`. It
does **not** close an arbitrary eleven-speed body, physical entry into the
`11+2` branch, mixed/even tail parities, another collision type, or LRC(14).

The primary audit reconstructs all `244` body walls and `28` positive
components, the low quotient topology, `1,053` primitive quotient controls
through `q=101`, the complete clock bank, and eight direct full-row controls
including `(a,b)=(2001,18009)`. Replay with

```text
python3 -B 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136.py
python3 -B -O 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136.py
```

The three streams byte-match the primary frozen output. The clean-room
referee independently joins only across active strict walls, scans the same
primitive range plus `70` hostile ratios through `q=1001`, and directly
checks all thirteen clearances in every residual row. Replay with

```text
python3 -B 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136_independent_audit.py
python3 -B -O 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136_independent_audit.py
PYTHONHASHSEED=8675309 python3 -B 04-computation/lrc14_fixed_body_universal_odd_tail_completion_thm4136_independent_audit.py
```

These streams byte-match the independent frozen output. **QED.**
