---
id: THM-3594
title: "Berggren positive two-cube slope atlas through denominator 201"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
  AUDIT.  Among all 2,072 primitive parity-correct reduced slopes through
  denominator 201, an exact screen by every prime at most 499 leaves 30
  rows.  Complete generalized-Pell class enumeration and unit-orbit analysis
  modulo 2n leaves exactly 15 parity-admissible slopes.  Every admissible
  class is moved into an invariant positive cone and gives an infinite ray
  of distinct positive two-cube points on the Berggren U-spine.  Eleven rays
  are new beyond THM-3580.  This is bounded at n<=201; no parametrization,
  density, or asymptotic is claimed.
source: kps-s188 / Berggren arithmetic wildcard, 2026-08-21
audit: >
  Author proof audit plus optimization-safe exact replay.  The finite screen,
  LMM completeness mechanism, all orbit addresses, fifteen compiled seeds,
  and positive-cone recurrences are pinned.  Independent hostile audit is
  still required before the stronger audit label is used.
depends_on:
  - THM-3547-positive-two-cube-slope-atlas-through-101
  - THM-3580-berggren-positive-cube-slope-atlas-completion-through-101
related:
  - THM-3370-berggren-two-cube-biquadratic-norm-collision
  - THM-3375-berggren-positive-two-cube-pell-ray
script: 04-computation/berggren_positive_cube_slope_atlas_201_thm3594.py
output: 05-knowledge/results/berggren_positive_cube_slope_atlas_201_thm3594.out
script_sha256: d3466965aa2ccaca81c3f4f4369dbca8f54f03c8472c7d5094459e38ea8977bd
output_sha256: dd86e7548773cf59113969bb7f0f7d99d345f98adbdf1de8b722533d15a02856
hash_basis: LF-normalized bytes
---

# THM-3594 -- Berggren positive two-cube slope atlas through 201

**PROVED + FINITE-EXACT + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE
AUDIT.**  This doubles the denominator range of THM-3580 while retaining its
complete global norm-class and parity-orbit sidecars.  It is not another
local screen presented as a classification.

## 1. Universe and compiler

For coprime integers

```text
3<=n<=201,        n odd,        m even,        n/2<m<n,       (1)
```

put

```text
d=n^2 u^2+2,
e=uW,
a=mn^2u^3+(2m+n)u,
q=m^2n^2u^4+2m(m+n)u^2+1.                              (2)
```

There are exactly `2,072` slopes `(m,n)` in `(1)`.  As in THM-3547, the
load-bearing equation is

```text
3W^2=n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2).                (3)
```

For every survivor below, neither `m` nor `n` is divisible by three.  Write

```text
K=(4m^2-n^2)/3,
C=4(2m^2+2mn-n^2)/3,
h=nu.                                                   (4)
```

Then `(3)` is exactly

```text
W^2-Kh^2=C.                                             (5)
```

The original odd compiler address is the iff condition

```text
W odd,                 h == n (mod 2n).                (6)
```

Indeed `(6)` is equivalent to `h=nu` with `u` odd.  Whenever `(5)--(6)`
hold and `0<W<n h`, formulas `(2)` give

```text
x=(d+e)/2 > y=(d-e)/2 >0,
r=(a-1)/2,
x^3+y^3=a^2+2=(2r+1)^2+2.                              (7)
```

The companion verifies both polynomial identities behind `(7)`:

```text
a^2+2=dq,                    4q-d^2=3e^2.              (8)
```

## 2. Exact local reduction

For every prime `p<=499`, the verifier exhausts

```text
3W^2 == n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2) (mod p).     (9)
```

The 95-prime screen excludes 2,042 rows and leaves exactly

```text
(14,23), (26,29), (26,47), (38,47), (38,53),
(50,53), (50,71), (62,95), (74,95), (74,101),
(98,101), (86,125), (98,125), (86,149), (110,149),
(122,149), (86,167), (98,167), (110,167), (146,167),
(98,173), (110,173), (122,173), (170,173), (98,191),
(110,191), (170,191), (110,197), (146,197), (182,197). (10)
```

Nothing is inferred merely from surviving `(9)`.  The screen is only the
cheap first half of the proof.

## 3. Complete generalized-Pell and parity-orbit decision

For every row of `(10)`, the exact Lagrange--Mollin--Matthews `P,Q,a,A,B,G`
algorithm inherited and independently audited in THM-3580 enumerates one
representative of every solution class of `(5)`.  Let

```text
epsilon=P+Q sqrt(K),                  P^2-KQ^2=1        (11)
```

be the fundamental positive unit.  Every integer solution comes from a
listed representative under powers of `epsilon`, signs, and conjugation.
The action

```text
(W,h) |-> (PW+KQh,QW+Ph) (mod 2n)                     (12)
```

is finite, so every class is followed through its exact full orbit and
tested against `(6)`.

The complete decision ledger is:

| slope | LMM classes | orbit period(s) | admissible classes |
|---|---:|---|---:|
| `(14,23)` | 3 | `11` | 2 |
| `(26,29)` | 6 | `10` | 2 |
| `(26,47)` | 9 | `46` | 6 |
| `(38,47)` | 0 | -- | 0 |
| `(38,53)` | 9 | `2` | 0 |
| `(50,53)` | 0 | -- | 0 |
| `(50,71)` | 0 | -- | 0 |
| `(62,95)` | 1 | `60` | 0 |
| `(74,95)` | 1 | `10` | 0 |
| `(74,101)` | 6 | `34` | 0 |
| `(98,101)` | 6 | `34` | 4 |
| `(86,125)` | 3 | `50` | 2 |
| `(98,125)` | 0 | -- | 0 |
| `(86,149)` | 6 | `50` | 2 |
| `(110,149)` | 6 | `25` | 2 |
| `(122,149)` | 5 | `150` | 4 |
| `(86,167)` | 6 | `83` | 4 |
| `(98,167)` | 2 | `83` | 2 |
| `(110,167)` | 12 | `166` | 8 |
| `(146,167)` | 6 | `83` | 4 |
| `(98,173)` | 0 | -- | 0 |
| `(110,173)` | 0 | -- | 0 |
| `(122,173)` | 0 | -- | 0 |
| `(170,173)` | 6 | `58` | 0 |
| `(98,191)` | 0 | -- | 0 |
| `(110,191)` | 0 | -- | 0 |
| `(170,191)` | 6 | `95` | 4 |
| `(110,197)` | 4 | `198` | 4 |
| `(146,197)` | 3 | `66` | 0 |
| `(182,197)` | 6 | `11` | 2 |

Hence the exact admissible set is

```text
(14,23), (26,29), (26,47), (98,101),
(86,125), (86,149), (110,149), (122,149),
(86,167), (98,167), (110,167), (146,167),
(170,191), (110,197), (182,197).                       (13)
```

The first four are THM-3580's rows.  The remaining eleven are new.

## 4. Every admissible row gives an infinite positive ray

Choose one admissible class and a phase satisfying `(6)`.  Independent sign
and conjugation let us take `W,h>0`.  Since

```text
sqrt(K)<n                                                   (14)
```

is equivalent to `m<n`, positive unit iteration eventually gives `W<n h`.
The verifier constructs such a point explicitly in every row of `(13)`.

The positivity is then permanent, not a one-point accident.  Put

```text
L=n h-W,                       H=nW-Kh.                 (15)
```

At the constructed seed, `L,H>0`.  If `P_T+Q_T sqrt(K)` is the positive
unit power for one complete orbit period modulo `2n`, then

```text
L'=P_T L+Q_T H,
H'=P_T H+KQ_T L.                                      (16)
```

Thus the cone is invariant, while `(6)` is restored after every period.
The `h` coordinate grows strictly.  Equations `(7)--(8)` therefore compile
infinitely many distinct positive ordered cube pairs for every slope in
`(13)`.

The exact selected seeds range from 18-bit to 7,286-bit cube coordinates.
Rather than print several thousand decimal digits, the output pins each
integer tuple by its bit lengths and SHA-256 digest.  Their combined exact
seed ledger has digest

```text
8ff45843ce97b3ea6bdeeeab691bf1c06b53f7b58ea82abb1fa28edc14952a40. (17)
```

## 5. Scope and next target

The conclusion is exactly bounded by `(1)`.  It proves no arithmetic
parametrization or density of admissible slopes and no asymptotic for all
two-cube intersections.  The jump from four rows through `101` to fifteen
through `201` is evidence against extrapolating sparsity from the first
atlas, but it is not itself an infinitude theorem for slopes.

The sharp unbounded question is now whether the parity-admissible LMM classes
in `(13)` lie on finitely many Vieta/Pell recurrences in the **slope**
variables `(m,n)`.  That is a different recurrence from the already proved
unit action within each fixed slope.

## 6. Reproduction

Run

```text
python -B 04-computation/berggren_positive_cube_slope_atlas_201_thm3594.py
python -B -O 04-computation/berggren_positive_cube_slope_atlas_201_thm3594.py
```

The normal and optimized transcripts agree line-for-line with the maintained
output.  The verifier hash-pins both audited parent programs, all 2,072
candidate slopes, the 95-prime screen, every class and orbit period, all
fifteen positive seeds and their next period iterate, and the cone
recurrences `(16)`.  Its semantic digest is

```text
2389ac8015a1f1394f938fe5ebf87933175dcdc7b4480901bdeba2c3a6b57de0.
```

**End of proof.**
