---
id: THM-4142
title: "Common safe-arc/clock pool universal odd-tail LRC(14) completion"
status: >
  PROVED ELEMENTARY COMMON-CERTIFICATE FAMILY + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. There is an explicit maximal 26-speed alphabet P for
  the THM-4136 arc and three-clock certificate. For every eleven-element
  subset B of P and every two distinct positive odd integers a,b, the thirteen
  speeds 2B union {a,b} are LRC(14)-safe. This closes 7,726,160 bodies at once
  and strictly generalizes THM-4136. Maximality is relative to this fixed arc
  and clock bank; arbitrary bodies, physical entry, mixed tail parity, and
  LRC(14) remain open. The namespace was renumbered from the colliding
  THM-4140 reservation by MISTAKE-508 before this promotion.
source: codex-frontier-synthesis-creative-20260825y
depends_on:
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
related:
  - THM-4129-universal-two-speed-completion-of-the-eleven-speed-lrc14-body
  - THM-4132-fixed-body-exceptional-scale-two-lrc14-completion
script: 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142.py
output: 05-knowledge/results/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142.out
independent_audit_script: 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142_independent_audit.out
script_sha256: 27e57b9a117cbf21c31da1e1c500d0dcf1b892b17f95bdd3c711c14a5f126e59
output_sha256: 31d805faa67c818a1b60b3ee705e4338934423ec871d5f2782136d3513563b11
independent_audit_script_sha256: 6cdfe76483895c6bd515894b7ac34c9bf28b67222d73c1c2180c316fed8585d5
independent_audit_output_sha256: 8c2dff6db1c74236ced2d9dff23e48cf46f383fde8a4c23b0be30c86165cfbd7
semantic_sha256: fe13698f821b5328bf7fc56267c0745280a7607261770f1b99424627f6686a3e
independent_semantic_sha256: dfddeafede97f47d6b8c38a4676cd6bc24c16559833690f51400addae8a5801b
hash_basis: raw LF bytes
primary_audit: >
  PASS. A pinned THM-4136 kernel and standalone exact-Fraction layer compute
  the complete speed universe, all arc and clock margins, 1,053 primitive
  quotient hostiles through q=101, all 68 residual ratios, and six complete
  thirteen-speed rows including new bodies, nonprimitive tails, and a tail
  of height 18,009. Normal, optimized, and hash-seeded outputs byte-match.
independent_audit: >
  ACCEPT. A clean-room standard-library referee imports neither primary code
  nor THM-4136. It reconstructs the pool by real-interval/integer separation,
  rebuilds strict two-sheet wall cells, checks 2,350 primitive ratios through
  q=151 and four far hostiles through q=1,001, verifies every residual clock,
  and directly solves four disjoint full rows. Normal, optimized, and
  hash-seeded outputs byte-match.
---

# THM-4142 -- common safe-arc/clock pool odd-tail completion

**PROVED ELEMENTARY COMMON-CERTIFICATE FAMILY + VERIFIED-EXACT +
INDEPENDENTLY AUDITED.**

THM-4136 closes every odd-tail pair over one fixed eleven-speed body. Its
infinite argument does not use the body's individual labels: it uses one
closed safe arc and three residual clocks. Intersecting those four
certificates before choosing the body exposes a much larger family.

## 1. Statement

Put

```text
delta=1/14,
J=[33/70,27/56],                         |J|=3/280,

P={1,3,4,6,8,9,10,11,12,14,15,16,18,22,24,28,30,
   32,35,37,41,43,45,49,60,64}.                         (1)
```

> **Theorem.** For every eleven-element subset `B subset P` and every two
> distinct positive odd integers `a,b`, the thirteen distinct speeds
>
> ```text
> 2B union {a,b}                                          (2)
> ```
>
> have a phase `x in R/Z` satisfying
>
> ```text
> min_(v in 2B union {a,b}) ||vx|| >= 1/14.               (3)
> ```

There are exactly

```text
binom(26,11)=7,726,160.                                  (4)
```

The THM-4136 body
`(1,4,6,8,10,12,14,15,16,18,22)` is one member.

## 2. The complete common certificate alphabet

For a positive integer `s`, define

```text
A(s):  min_(y in J)||sy|| >= delta.                       (5)
```

The exact minimum is zero if the real interval `sJ` contains an integer and
otherwise is the smaller endpoint distance. Since

```text
94|J|=141/140>1,                                         (6)
```

every `s>=94` fails `(5)`. Exact rational evaluation for `1<=s<=93` gives
the complete `J`-safe alphabet

```text
P_J={1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,
     26,28,30,32,33,35,37,39,41,43,45,47,49,60,62,64,66}.
                                                               (7)
```

Retain THM-4136's clocks

```text
xi_0=89/1176,       xi_1=181/4704,       xi_2=431/4480.   (8)
```

Intersect `(7)` with the three strict conditions

```text
||2s xi_j||>delta,                         j=0,1,2.        (9)
```

The result is exactly `P` in `(1)`. The ten `J`-safe exclusions are

```text
5,7,13,20,26,33,39,47,62,66.                            (10)
```

Each has a literal failure in `(9)`; speed `26` fails all three. Thus `P` is
maximal among positive integer speeds sharing this particular arc and all
three clocks. This is certificate maximality, not a claim that a body using
one of `(10)` is unsafe.

For later use, the common body minima and their unique owners are

| clock | `min_(s in P)||2s xi_j||` | owner |
|:---|---:|---:|
| `xi_0` | `4/49` | `60` |
| `xi_1` | `11/147` | `64` |
| `xi_2` | `17/224` | `10` |

All are strictly above `1/14`. Equations `(5)` and `(9)` therefore hold for
every subset `B`, without enumerating the `7,726,160` choices.

## 3. The body arc and its two physical sheets

Fix `B subset P`, `|B|=11`. Equation `(5)` gives

```text
J subset {y: min_(s in B)||sy||>=delta}.                 (11)
```

Hence both half-lifts

```text
x_0=y/2,                    x_1=(y+1)/2                  (12)
```

keep every speed in `2B` safe.

Interchange the odd tails if necessary and write

```text
a=pt,        b=qt,        0<p<q,                         (13)
```

where `t=gcd(a,b)` and `p,q,t` are odd with `gcd(p,q)=1`. Put `w=ty` and
let `z=w/2` denote one square-root sheet. The two tail packets at `(12)` are,
up to sheet interchange,

```text
(pz,qz),              (p(z+1/2),q(z+1/2)).               (14)
```

The interchange loses only the parity of the integer removed when reducing
`ty`; retaining both physical phases makes it harmless.

## 4. The body-independent cross-comb bound

For odd `r`, put

```text
D_r={z in R/Z: ||rz||<delta}.                             (15)
```

The strict sets `D_r` and `D_r-1/2` are disjoint because their values differ
by `1/2`, while two danger radii total only `1/7`. Thus simultaneous failure
of both sheets in `(14)` consists only of the two cross terms

```text
D_p intersect (D_q-1/2),
D_q intersect (D_p-1/2).                                 (16)
```

After doubling to the `w`-circle they are antipodal copies. Every component
lies in one doubled `q`-tooth and therefore has length at most

```text
2/(7q),                                      q>=9.        (17)
```

The six primitive ratios with `q<9` have exact maximum component widths

```text
(1,3):0, (1,5):0, (3,5):1/105,
(1,7):1/49, (3,7):1/49, (5,7):1/49.                      (18)
```

Consequently every two-sheet bad component has width at most `2/63`, sharp
at `(p,q)=(1,9)`, and zero is outside the open bad quotient.

If `t>=3`, compact containment of the connected image `tJ` inside one open
bad component would force `t|J|<2/63`, but

```text
3|J|-2/63=1/2520>0.                                      (19)
```

If `t=1` and `q>=27`, `(17)` instead contradicts

```text
|J|-2/(7*27)=1/7560>0.                                   (20)
```

Wrapping contains zero and is impossible as well. Thus only `t=1` and
`q<=25` remain.

## 5. The 68 residual ratios close uniformly over every body

There are exactly `68` coprime odd pairs `p<q<=25`. Use `(8)` according to
the same three categories as THM-4136. Because `(9)` holds for every
`s in P`, the body can now be chosen *after* the clock.

| condition | rows | minimum tail gap | minimum `2P` body gap | full gap |
|:---|---:|---:|---:|---:|
| `13 notin {p,q}` | `56` | `89/1176` | `4/49` | `89/1176` |
| one tail `13`, excluding `(1,13),(13,25)` | `10` | `541/4704` | `11/147` | `11/147` |
| `(1,13)` or `(13,25)` | `2` | `431/4480` | `17/224` | `17/224` |

Every full gap is strictly greater than `1/14`. Sections 3--5 prove `(3)`,
and `(4)` gives the body count. **QED.**

## 6. Boundary, losses, and hostiles

Closed endpoint safety is load-bearing: speeds `4`, `15`, and `60` attain
exactly `1/14` somewhere on `J`. Strictifying `(5)` would discard valid
bodies. Conversely, allowing weak clock safety would be invalid because the
residual tail proof needs a strict full-row margin.

The sharp fixed-sheet hostile survives unchanged. At primitive ratio `(1,9)`
and quotient phase `w=1/9`, the two packets are

```text
sheet 0: (1/18,1/2),             sheet 1: (4/9,0).        (21)
```

Tail `1` kills the first sheet and tail `9` the second. Neither a single
sheet nor a naked interval-length comparison proves the theorem.

The connection contract is

```text
source:       positive speeds sharing J and all three residual clocks
target:       7,726,160 eleven-bodies with arbitrary distinct odd tails
map:          choose B subset P, then y -> {y/2,(y+1)/2}
preserved:    eleven body clearances, three low clocks, one physical phase
destroyed:    larger body-specific arcs and alternative clock certificates
sidecar:      the strict two-sheet cross-comb component, with sheet parity
decisive test: component width versus gcd(a,b)|J|, then the 68 clocks. (22)
```

The operation is intersection of compatible certificates before
scalarization. Intersecting only the three clock-safe sets would lose the
compact arc; intersecting only body-safe labels would lose the finite
residual closure.

## 7. Scope and exact replay

This theorem promotes a fixed-body result to a large common-certificate
family. It does not say that `P` contains every safe body label, that bodies
outside `P` fail, or that these are all odd-tail-completable bodies. It proves
no physical entry into this parity class, no mixed/even-tail theorem, and no
general LRC(14).

The primary audit pins THM-4136's source hash, reconstructs `(7)--(10)`,
checks `1,053` primitive quotient ratios through `q=101`, all `68` residual
ratios, and six complete full rows. Replay with

```text
python3 -B 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142.py
python3 -B -O 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142.py
```

The independent referee imports neither implementation. It reconstructs the
pool by real-interval/integer separation, rebuilds strict quotient cells,
checks `2,350` ratios through `q=151`, four far hostiles through `q=1,001`,
all residual clocks, and four disjoint direct rows. Replay with

```text
python3 -B 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142_independent_audit.py
python3 -B -O 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142_independent_audit.py
PYTHONHASHSEED=8675309 python3 -B 04-computation/lrc14_common_safe_arc_clock_pool_odd_tail_completion_thm4142_independent_audit.py
```

Each replay family byte-matches its frozen output. **QED.**
