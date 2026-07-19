---
id: THM-1169
title: The first-safe owner window closes every three- and four-carrier r=6 row
status: PROVED — THM-1155's exact multi-comb density bound makes four carriers uniform immediately.  For three carriers, truncation at the second core speed is exactly THM-1154's d<=5 owner budget; density leaves x<=12,18,25,32,38,45 in the six possible core strata, and a complete endpoint-exact Fraction audit checks the resulting 473 two-comb candidates with zero covers and zero point-only residuals.  This theorem closes one through four thirteen-carriers; its explicit rho=5,6 remainder is subsequently closed by THM-1212
source: codex-2026-07-18-S75 multi-carrier owner attack
depends_on: [THM-1154, THM-1155]
related: [THM-1135, THM-1136, THM-1151, THM-1212]
script: 04-computation/lrc14_r6_three_four_carrier_owner_window_referee_codex_S75.py
output: 05-knowledge/results/lrc14_r6_three_four_carrier_owner_window_referee_codex_S75.out
---

# THM-1169 — the first-safe owner window closes three and four carriers

Let

```text
P={p1<p2<...<p7} subset {1,...,12},       M=p7,
13M<k1<...<k6.
```

Suppose exactly `rho` killers are divisible by `13`, and write their distinct
normalized values as

```text
O={x=o1<o2<...<o_rho},       carrier killers={13o:o in O}.       (1)
```

No Covering hypothesis is used below, and the `6-rho` noncarriers may have
arbitrary size.  Put

```text
D_o={u in R: ||ou||<1/14}.                              (2)
```

The variable `u` is the normalized carrier clock.  A safe boundary

```text
u=c/(14s),        s in O,                               (3)
```

feeds THM-1154's actual twelve-chart family

```text
t_a=a/13+u/13=a/13+c/(182s),       a=1,...,12.          (4)
```

The distinction between `u` and the final time `t_a` is essential.

## 1. Boundary readout from a carrier-safe window

The least carrier `x` is safe on its first closed safe window

```text
I_x=[1/(14x),13/(14x)],       |I_x|=6/(7x).             (5)
```

Suppose a closed subinterval `J subset I_x` is not covered by the open danger
sets of the later carriers.  Its carrier-safe subset is compact and nonempty.
Take its least point `u`.  If `u` is the left endpoint of `J`, that endpoint
is already of the form (3), owned by `x`.  Otherwise, points immediately to
the left are covered, so `u` is the boundary of some `D_s`, `s in O`.  Every
such boundary has the form (3) for a positive integer `c`.

Two upper endpoints turn this topological readout into THM-1154's arithmetic
conditions.

1. If `u<=13/(14x)`, then `cx<=13s`.  Since `M<x`, this gives the broad core
   cap `cM<13s`.
2. If additionally `u<=1/(14p2)`, then `cp2<=s`.  Both `p1` and `p2` are
   therefore un-crossed core owners, so

   ```text
   d=#{p in P:cp>s}<=5.                                 (6)
   ```

For `rho=4`, every possible `d<=7` already satisfies `d<2rho`.  For
`rho=3`, (6) is exactly `d<2rho`.  Once the later carrier combs fail to cover
the chosen window, THM-1154 absorbs the noncarriers: the obstruction ledger is

```text
rho=4: d+2(6-rho)<=7+4=11<12,
rho=3: d+2(6-rho)<=5+6=11<12.                           (7)
```

Thus at least one multiplier chart (4) survives.

## 2. Four carriers close uniformly by density

Assume `rho=4`.  If no point of `I_x` were safe for all four carriers, the
three later danger combs would cover `I_x`.  Apply THM-1155's formalized
multi-speed density bound with `lambda=1/14`, interval length `L=6/(7x)`,
and three speeds at least `x+1`:

```text
L(1-2*3/14) <= (1/7) sum_{i=2}^4 1/o_i
          <= 3/(7(x+1)).                                (8)
```

Hence a cover would force

```text
L<=3/(4(x+1)).                                          (9)
```

But the exact difference is

```text
6/(7x)-3/(4(x+1))=(3x+24)/(28x(x+1))>0.                (10)
```

This contradiction leaves a carrier-safe point in `I_x`.  Section 1 gives a
boundary chart with `cM<13s`; `d<=7<8` completes (7).

> **Four-carrier corollary.**  Every seven-core/six-killer row with exactly
> four `13`-divisible killers is `1/14`-lonely.

This is search-free and uniform in all killer magnitudes.

## 3. The three-carrier owner window

Now assume `rho=3` and abbreviate `p=p2`.  Use the owner-truncated window

```text
J_{p,x}=[1/(14x), min(13/(14x),1/(14p))].               (11)
```

It lies in `I_x`, and every boundary read from it satisfies (6).  A seven-set
inside `{1,...,12}` has

```text
2<=p<=7,       M>=p+5,       x>=M+1>=p+6.              (12)
```

If the two later carriers `y<z` covered (11), THM-1155 with two combs would
give

```text
5L<=1/y+1/z<=2/(x+1),       L=|J_{p,x}|.               (13)
```

When `x>=13p`, (11) is the full first-safe window and

```text
6/(7x)-2/(5(x+1))=(16x+30)/(35x(x+1))>0,               (14)
```

contradicting (13).  When `x<13p`, its length is

```text
L=(x-p)/(14xp).                                         (15)
```

Therefore density alone closes the row whenever

```text
5(x-p)(x+1)>28xp.                                       (16)
```

The six possible `p` values leave the following exact finite residual.

| `p=p2` | cores with this `p` | first `x` closed by (16) | residual `x` |
|---:|---:|---:|---:|
| 2 | 252 | 13 | 8..12 |
| 3 | 252 | 19 | 9..18 |
| 4 | 168 | 26 | 10..25 |
| 5 | 84  | 33 | 11..32 |
| 6 | 30  | 39 | 12..38 |
| 7 | 6   | 46 | 13..45 |

The core counts sum to `C(12,7)=792`.  Checking (16) over this table is a
finite integer calculation; the referee also verifies that it remains true
from the displayed first value through `13p`, after which (14) applies.

## 4. The finite two-comb lemma

> **Finite two-comb lemma.**  For every `(p,x)` in the residual table and all
> distinct integers `x<y<z`, the closed interval `J_{p,x}` is not contained
> in the open union `D_y union D_z`.

The reduction to the exact audit is itself finite and rigorous.  If a cover
existed, the strict inequality `y<z` sharpens (13) to

```text
y<2/(5L).                                               (17)
```

Thus `y` has an explicit finite integer cap.  Remove the open comb `D_y`
from `J_{p,x}`.  If the residual has a closed positive-length cell of maximum
length `ell`, then a single open `z`-tooth, of length exactly `1/(7z)`, can
contain that closed cell only if

```text
z<1/(7ell).                                             (18)
```

If only isolated residual points occurred, (18) would not bound `z`; this is
why the referee records that case separately instead of silently discarding
it.  There are exactly zero point-only rows.

The complete exact census is

| `p` | admissible `(p,x,y)` rows | bounded `(p,x,y,z)` candidates | covers |
|---:|---:|---:|---:|
| 2 | 17  | 0   | 0 |
| 3 | 80  | 0   | 0 |
| 4 | 192 | 3   | 0 |
| 5 | 355 | 34  | 0 |
| 6 | 576 | 115 | 0 |
| 7 | 854 | 321 | 0 |
| **total** | **2074** | **473** | **0** |

No first later comb covers the window by itself.  There are also zero
point-only complements after the first comb, so every remaining complement
contains a positive-length closed cell and (18) is exhaustive.

Every bounded pair candidate is checked twice.  The first checker splits at all rational
tooth boundaries and tests every wall plus one midpoint of every open cell.
The second constructs connected components of the union of open teeth,
merging strict overlaps but not touching intervals.  They agree on all 473
candidates.  Endpoint equality is therefore accepted exactly as required by
the LRC predicate.

The finite lemma shows that (11) always contains a carrier-safe point.
Sections 1 and (7) finish the owner chart.

> **Three-carrier corollary.**  Every seven-core/six-killer row with exactly
> three `13`-divisible killers is `1/14`-lonely.

Together with THM-1136 and THM-1154, all rows with one, two, three, or four
carriers are now uniformly closed.

## 5. Exact remainder for five and six carriers

The same first-safe window gives an honest narrowing, but not an immediate
closure, at the remaining cardinalities.  Write `y=o2`.

For `rho=5`, four later combs covering `I_x` would imply

```text
3L<=sum_{i=2}^5 1/o_i <4/y,
so y<4/(3L)=14x/9.                                     (19)
```

For `rho=6`, five later combs covering `I_x` would imply

```text
2L<=sum_{i=2}^6 1/o_i <5/y,
so y<5/(2L)=35x/12.                                    (20)
```

Consequently the owner-window method already closes the complementary cones

```text
rho=5: y>=14x/9,
rho=6: y>=35x/12.                                      (21)
```

The exact remaining carrier-compatibility problem is the close-first-pair
region opposite (21).  Neither (19) nor (20) claims that those residual cones
contain counterexamples.  They are the next interval-cell classification
target.  For `rho>=4`, THM-1154's owner budget is automatic once a compatible
carrier boundary is found.

THM-1212 subsequently closes both regions by iterating nested prefix windows;
thus (21) records the exact handoff from this theorem, not a current open
frontier.

## 6. Tournament Analysis and the underlying object

The six `p2` residual strata can be made into a telemetry tournament.  Orient
a pair toward the stratum with the larger exact `(p,x,y,z)` candidate count;
break the sole `0=0` tie by increasing `p`.  The resulting tournament is
transitive, with score multiset `{0,1,2,3,4,5}`, no directed cycles, six
singleton SCCs, and one Hamiltonian path.  Reversing the tie switch flips one
edge.

That tournament proves nothing by itself.  The proof-bearing vertices are the
labelled rational wall cells of `J_{p,x}`, with the `y`- and `z`-comb cover
obligations attached.  This object preserves endpoint openness, metric cell
length, the boundary owner, and the existence of the exact chart (4).  A
tournament on runners, carriers, carrier ratios, or residual strata destroys
those data.

Candidate vertex sets challenged here were runners, carriers, ratios, danger
teeth, wall events, rational cells, owner obligations, and cover obligations.
The labelled wall-cell complex is the smallest faithful choice used by the
proof.  It forgets killer order outside the carrier set and all noncarrier
magnitudes, precisely because THM-1154 absorbs those coordinates afterward.

## 7. Exact replay

The dependency-free referee uses only `fractions.Fraction` and integers.  It
checks all 792 cores' `p2` distribution, both density identities, all six
integer threshold rows, all 2,074 finite `y` rows, all 473 bounded two-comb
candidates, both open-cover implementations, the absence of one-comb covers and point-only
residuals, the five/six-carrier ratio cones, and the tournament fingerprint.
Normal and optimized executions are byte-identical.

```text
04-computation/lrc14_r6_three_four_carrier_owner_window_referee_codex_S75.py
SHA-256 d49e22bc7526d2355fd1d9f30a1e24bcab4c98b2cd372835ae02f01647020f3b

05-knowledge/results/lrc14_r6_three_four_carrier_owner_window_referee_codex_S75.out
SHA-256 37b829d87af131756cd8d855b154beb26043d38a1b148a9116802c589f808da3
```
