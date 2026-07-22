---
id: THM-2078
title: "Bounded terminal guard-containment closure through height 24"
status: >
  PROVED by exact integer-bitset census. No nontrivial THM-2073 dyadic tower
  can have terminal maximum at most 24. THM-2076 originally restricts
  terminal size to 6 through 10 and THM-2080 subsequently removes rank 6;
  the terminal maximizer interval bounds the preceding odd
  guard by (13-s)h<2(s+1)max(Q). Among all 4,484,931 cores in those ranks,
  exactly 30,594 are hereditarily primitive and divisor-complete through 14.
  Every allowed core/guard pair already fails containment at a rational phase
  of denominator 8192; zero pairs survive to continuous atomization. This
  closes a low-terminal/unbounded-top slice, not LRC(14).
source: codex-2026-07-21-LRC-bounded-terminal-closure
depends_on:
  - THM-2073
  - THM-2076
  - THM-2077
related:
  - THM-2066
  - THM-2068
  - THM-2072
script: 04-computation/lrc14_terminal_guard_containment_referee_codex_20260721.py
output: 05-knowledge/results/lrc14_terminal_guard_containment_referee_codex_20260721.out
script_sha256: 4e0d36e38aa2fabf74bf7ba012f24c764e95cc79dd44e4a538d12814b318de53
output_sha256: b542a94aebec338ca4b91253b205429f986b881b93f276ac0036237e1bd888d2
hash_basis: normalized repository blobs (LF)
---

# THM-2078 -- bounded terminal guard-containment closure

Put `delta=1/14`. Retain a nontrivial hypothetical tower from THM-2073:

```text
C=Q_0,
Q_i=2Q_(i+1) union {h_i},    0<=i<r,    r>=1.           (1)
```

Every quotient is primitive and divisor-complete through `14`, the terminal
`Q_r` is hereditarily primitive, and the last guard satisfies

```text
G_(Q_r) subset {t:||h_(r-1)t||<1/7}.                    (2)
```

Then

```text
max(Q_r)>=25.                                           (3)
```

## 1. The finite terminal box

Let

```text
s=|Q_r|,    B=max(Q_r),    h=h_(r-1).                   (4)
```

THM-2076 gives `6<=s<=10`; THM-2080 subsequently removes `s=6`, so the live
range is `7<=s<=10`. Let `mu=M(Q_r)` and
`rho=(mu-delta)/B`. Settled LRC for the `s` terminal speeds gives

```text
mu>=1/(s+1),
rho>=(13-s)/(14(s+1)B).                                 (5)
```

The maximizer interval of radius `rho` lies in `G_(Q_r)`. Condition (2) puts
that closed interval strictly inside one guard tooth of radius `1/(7h)`.
Therefore

```text
h<1/(7rho)<=2(s+1)B/(13-s),
(13-s)h<2(s+1)B.                                       (6)
```

If `B<=24`, the complete search box is

```text
Q subset {1,...,24},    6<=|Q|<=10,
Q hereditarily primitive and divisor-complete through 14,
h positive odd,    (13-|Q|)h<2(|Q|+1)max(Q).            (7)
```

The largest possible guard occurs at `s=10,B=24` and is at most `175`.

## 2. Exact census

The original exhaustive run used the larger pre-THM-2080 range, so its total
core counts and arithmetic survivors are

| `s` | total `C(24,s)` | hereditary divisor-complete |
|---:|---:|---:|
| 6 | 134,596 | 6 |
| 7 | 346,104 | 131 |
| 8 | 735,471 | 1,198 |
| 9 | 1,307,504 | 6,375 |
| 10 | 1,961,256 | 22,884 |
| **sum** | **4,484,931** | **30,594** |

For every surviving core and every odd `h` allowed by (6), the exact rational
grid

```text
T_8192={j/8192:0<=j<8192}                               (8)
```

contains a witness `t` such that

```text
||qt||>=1/14 for every q in Q,
||ht||>=1/7.                                            (9)
```

Thus `t in G_Q` but the guard is not eligible, contradicting (2). The stored
census records

```text
sampled necessary-condition survivors: 0,
exact continuous-containment survivors: 0.              (10)
```

This is a rigorous finite proof, not a mesh approximation. A true continuous
containment would hold at every rational grid point and hence would have to
pass the bitset filter. Since no pair passes, resolution cannot hide a
continuous survivor. The companion also implements complete rational-boundary
atomization, including endpoints, for any prefilter survivor; that stage is
dormant because the exact necessary grid closes every pair. This proves (3).
QED.

## 3. Relation to the bounded owner atlas

THM-2066/2068 close original quotient cores `C` with `max(C)<=24`. The present
theorem is different: only the terminal quotient is bounded, while the
original core contains scaled terminal speeds and guards and may be much
larger. The last guard alone contradicts terminal safe-set containment before
the original odd tails or earlier address bits are enumerated. Thus (3)
closes a low-terminal/unbounded-top rectangle in the THM-2073 atlas.

THM-2077's terminal interval is the lawful reason the guard search is finite.
Searching an arbitrary cutoff such as `h<=175` without (6) would be evidence
only; equations (5)--(7) make the audit exhaustive.

## 4. Scope and next target

Every nontrivial tower now satisfies `max(Q_r)>=25`, but no upper bound is
known. THM-2072 rules out a fixed clock bank uniform in the core. A scalable
completion should instead prove an adaptive terminal conductor:

```text
for every hereditary divisor-complete Q of size 7,...,10,
and every odd h satisfying (6), find a rational t with (9).               (11)
```

THM-2075 says that the component and endpoint word to be hit is entirely
terminal, while THM-2077 supplies a quarter escape and relative-height box.
Those sidecars are better candidates for a scale-dependent clock than merely
increasing the fixed denominator `8192`.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that bounding the original seam core is the
only meaningful finite slice. Safe-child homeomorphism makes the terminal
maximum a separate coordinate, and the last guard tests it before any
large-scale reconstruction. The rational-grid quotient preserves the sole
implication needed for exclusion—continuous containment implies sampled
containment. Its loss of widths and owner labels is harmless only because no
sampled pair survives.

Tournament Analysis is not a proof carrier. Orienting guards by the number of
sampled safe residues they miss gives a search heuristic but cannot certify
set containment. The exact objects are core obligations, guard eligibility
bitsets, and their rational witness phases. QED.
