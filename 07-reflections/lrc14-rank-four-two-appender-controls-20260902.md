# LRC(14) rank-four two-appender controls, 2026-09-02

**Status: FINITE-EXACT on two named THM-4333 residual pairs; analytic
consequence relative to THM-4170/4150. LRC(14) remains open.**

## Inheritance and exact target

THM-4333 proves a strict rank-three retained-mass floor `L_3/D>2/27` for
nine-pool bodies on every pair in the frozen 181,194-pair remainder. Its
canonical exact-fallback hostile is `(50,70)`; its weakest degree certificate
is `(509,640)`. The theorem explicitly leaves higher failure ranks open.

For a residual pair `Q={q,r}`, retain the THM-4333 pool `P`. For a pair-safe
wall cell `C`, let `F(C)` be its pool-failure set and let `w_C` be its integer
width on the exact grid `D`. For every `K in binom(P,8)`, define

```text
L_4(K;Q)=sum_(C: |F(C)|<=4, F(C) intersect K=empty) w_C.       (1)
```

This is retained truncation mass, so

```text
0 <= L_4(K;Q)/D <= mu(G_(K union Q)).                          (2)
```

Two successive THM-4170 discrepancy losses dictate the strict target

```text
(6/7)(7/81)=2/27,        (6/7)(2/27)=4/63.                    (3)
```

Thus the precise finite question is

```text
for every K in binom(P,8),     81 L_4(K;Q)-7D > 0.            (4)
```

The source is the pair-safe wall arrangement; the target is every eight-body
on each named residual pair; the map discards all failure ranks at least five.
It preserves exact rank-at-most-four mass and overlap incidence, but destroys
the omitted safe mass, cyclic address, owner, and endpoint data.

## Exact result

The literal-midpoint wall program flat-enumerates all
`binom(30,8)=5,852,925` bodies without pruning. A second exact optimizer uses
direct weighted coverage of the retained wall cells, not the rank-four
inclusion-exclusion tensors, and proves the same optimum and least-mask tie.

| residual pair | exact grid `D` | minimum `L_4` | `81L_4-7D` | least body |
|:--|--:|--:|--:|:--|
| `(50,70)` | `91,205,797,082,400` | `11,108,019,386,090` | `261,308,990,696,490` | `{10,80,95,120,143,145,168,170}` |
| `(509,640)` | `74,278,001,143,906,560` | `11,595,224,392,491,506` | `419,267,167,784,466,066` | `{8,80,88,95,145,168,170,286}` |

Both named pairs therefore satisfy `(4)` strictly. The normalized minima and
surpluses are

```text
(50,70):   L_4/D = 100981994419/829143609840,
           L_4/D-7/81 = 263948475451/7462292488560;

(509,640): L_4/D = 445970168941981/2856846197842560,
           L_4/D-7/81 = 1791740033266949/25711615780583040.    (5)
```

The root sum-of-largest-degrees screen is inconclusive on both pairs. Its
target ticks are respectively `-431,061,626,070,252` and
`-40,171,801,599,658,044`. Rank-four overlap data, rather than another
marginal-degree estimate, supplies the positive result.

## Retained truncation is not complete safe mass

On the two rank-four minimizers above, direct integration over **all**
pair-safe wall cells gives

| pair | retained `L_4` | complete safe mass on the same body | omitted safe mass |
|:--|--:|--:|--:|
| `(50,70)` | `11,108,019,386,090` | `14,349,106,369,428` | `3,241,086,983,338` |
| `(509,640)` | `11,595,224,392,491,506` | `14,171,259,234,541,482` | `2,576,034,842,049,976` |

These complete masses belong to the rank-four-minimizing bodies. They are
not asserted to be minima of complete safe mass over all bodies. No omitted
rank-at-least-five cell was used to prove `(4)`.

## Direct one-appender LRC(14) consequence

Write `m_Q` for the appropriate minimum retained ratio in `(5)`, and put

```text
C_Q=1891+q+r,                 gamma_Q=(6/7)m_Q-4/63,    (6)
```

where `1891` is the largest possible sum of eight labels of `P`. The exact
constants are

```text
Q=(50,70):   C_Q=2011,
              gamma_Q=118691847737/2902002634440,
              (6C_Q/49)/gamma_Q
                =714603342594960/118691847737<6021;

Q=(509,640): C_Q=3040,
              gamma_Q=703055796194263/9998961692448960,
              (6C_Q/49)/gamma_Q
                =3722062474903449600/703055796194263<5295.   (7)
```

Therefore, for every `K in binom(P,8)`, every integer `s>=6021` on the
first pair or `s>=5295` on the second, and `U=G_(K union Q)`, THM-4170's
discrepancy estimate gives

```text
mu(U intersect G_s)
 >= (6/7)mu(U)-6c(U)/(49s)
 >= (6/7)m_Q-6C_Q/(49s)
 >= 4/63.                                               (8)
```

Both rational thresholds in `(7)` are nonintegral, so the last inequality
is actually strict at the displayed integer cutoffs. At those cutoffs the
certified lower masses and their strict surpluses over `4/63` are

```text
(50,70):   41090163799831/647146587480120,
           surplus 497192957/215715529160040;
(509,640): 74714970194220413/1176544492478160960,
           surplus 41197729678199/3529633477434482880.  (9)
```

The cutoffs exceed every label of `K union Q`, so `K union Q union {s}`
has eleven distinct positive speeds. For every positive integer `d` and
every two distinct positive odd integers `a,b`, THM-4150 now proves safety
of

```text
2d(K union Q union {s}) union {a,b}.                   (10)
```

The doubled body consists of eleven distinct even speeds and cannot collide
with the two odd tails. Thus `(10)` has exactly thirteen relative speeds,
the correct LRC(14) count. This sharpens the cofinal newcomer bound on these
two named residual pairs only.

## Secondary two-appender wedge

The `7/81` target was generated by allowing two successive `6/7` losses.
Using an eight-body directly with two appenders and two odd tails would,
however, produce fourteen relative speeds, a structured LRC(15) row. A
same-count LRC(14) side consequence follows by heredity: take
`J in binom(P,7)`, extend it to an eight-body, and then drop the eighth
constraint. If `delta_Q=m_Q-7/81`, put

```text
C'_Q=1715+q+r,
epsilon_Q(s)=(6/7)delta_Q-6C'_Q/(49s).                 (11)
```

Then the ordered conditions

```text
s > C'_Q/(7 delta_Q),
t > max{s,(C'_Q+s)/(7 epsilon_Q(s))}                   (12)
```

give

```text
mu(G_(J union Q union {s,t}))>4/63.                   (13)
```

The first appender leaves mass at least `2/27+epsilon_Q(s)`, hence strictly
above `2/27`, and the second strict inequality gives `(13)`. The least
integer passing only the first condition is `7412` for `(50,70)` and `5872`
for `(509,640)`. Condition `(12)` makes all eleven core labels distinct;
doubling and adding two distinct odd tails again gives thirteen relative
speeds. This is an ordered wedge, not a cofinal rectangle, and it is
secondary to the simpler one-appender consequence `(10)`.

## Audit and next test

- Normal `-O2`, `-O3`, warnings-as-errors, and `-O0` undefined-behaviour
  sanitizer C++20 builds are byte-identical.
- The flat inclusion-exclusion census visits all `5,852,925` bodies per pair.
- Independent direct-cell coverage branch-and-bound proves the same rewards
  and least bodies, in `1,841` and `4,011` nodes.
- A no-import Python implementation uses the unfurled two-sided midpoint
  condition to reconstruct both wall/rank ledgers and directly reproduce the
  reported candidates' retained and complete masses. Normal, optimized, and
  hash-seeded outputs byte-match. This checks the wall semantics and candidates,
  while the independent global optimization check remains the C++ direct-cell
  branch-and-bound.
- Direct cell summation reproduces each minimizing retained mass and separately
  computes its complete safe mass.
- The rank-zero-through-three prefixes are exactly
  `22,084,503,304,648` and `21,053,923,224,812,998`, matching the audited
  THM-4333 packet.
- The frozen remainder hash remains
  `9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1`,
  and both named pairs occur in it.

Reproduce from the repository root:

```text
clang++ -std=c++20 -O2 -Wall -Wextra -pedantic \
  04-computation/lrc14_rank4_two_appender_pair_probe_20260902.cpp \
  -o /tmp/lrc14_rank4_two_appender_pair_probe
/tmp/lrc14_rank4_two_appender_pair_probe
python3 \
  04-computation/lrc14_rank4_two_appender_pair_probe_20260902_independent.py
```

The source/output SHA-256 hashes are respectively
`91f885d72bf989ad79964cee92e0b7da0b1e1a9ee16d6df4d60c3c7ffc8e2f3f`
and
`bd48c1d4a3d66146f7d3ded13c911240e664ee5bef39bd2b83aa96da44627894`.
The independent source/output hashes are
`4e69685eb7be35c4cdae81778a0c1d46217067893d4a425b4120caf03d085db1`
and
`fd342ff8ee0b4bb5c384209b2c779fae2fdd03d06d991dbe2bd591e73f7b3db7`.

The positive control suggests a complete rank-four scan of the frozen
remainder. The direct-coverage optimizer is especially promising: both
degree-hostile controls close in at most 4,011 search nodes. Nothing here
establishes `(4)` for the other 181,192 pairs, arbitrary outsiders, arbitrary
entry, or LRC(14).
