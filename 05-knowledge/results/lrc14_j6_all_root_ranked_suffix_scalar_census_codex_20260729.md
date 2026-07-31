# All-root ranked-suffix scalar census

**Status: PROVED + FINITE-EXACT + VERIFIED; THM-2899.  LRC(14) remains
open.**

## Result

THM-2896 gives every seven-body root

```text
E in C({1,...,14},7)
```

an ordered external six-cover hitting gate of least size `K(E)<=25`.  Assign
each hypothetical six-cover to its earliest gate label.  On the resulting
literal first-apex carrier, forbid the complete gate prefix and globally
seal the five largest allowed singleton coverages

```text
q_1>=...>=q_5.
```

The complete exact universe and scalar verdict are

```text
roots                                      3432
marked suffix branches                    41415
q_1+...+q_5<h closures                    26609
scalar-hard branches                      14806
whole roots closed by the scalar test         5
roots nonterminal under this scalar test   3427.
```

The five scalar-terminal bodies are

```text
(1,2,3,4,5,6,13),
(1,2,3,4,6,7,14),
(1,2,3,4,6,11,13),
(1,2,3,4,6,12,13),
(7,8,9,11,12,13,14).
```

They are disjoint from the four THM-2895 roots and the independently proved
THM-2898 root.  Those three theorems therefore close ten distinct roots and
leave `3422` in the official seven-body residual.

## Global tail seal

For a carrier of mass `h` with `r` interval components, THM-735(ii) gives

```text
c(w)<h/7+(99/70)r/(7w).
```

The verifier first scans every allowed `15<=w<=1600`, then extends to

```text
tail_first=max(
  1601,
  ceil((99/70)r/(7(q_5-h/7)))
).
```

It verifies `q_5>h/7` and

```text
h/7+(99/70)r/(7 tail_first)<=q_5.
```

Thus every omitted speed is strictly below the retained fifth rank.  The
order statistics are global, not finite-window proposals.  The largest
tail entry anywhere is `2923`; the largest speed retained in a top five is
`378`.

## Uniform five-slot parity entry

Every one of the `14806` scalar-hard branches satisfies

```text
q_1<3h/7.
```

There are zero failures in every THM-2896 stratum and every surviving apex
rank `1,...,13`.  The smallest strict margin is

```text
3h/7-q_1 = 369209/35315280
```

at

```text
E=(2,4,9,10,12,13,14), rank=1, apex=22,
h=62558/315315, q_1=c(16)=75245/1009008.
```

THM-2895 therefore forces two labels of every hypothetical five-cover into

```text
H_4={w:c(w)>=(h-q_1)/4}.
```

The exact discrepancy cutoff for `H_4` is at most `2782`.  Its nearest-rank
quantiles over all hard branches are

```text
percent       0   25   50   75   90   95    99   100
cutoff      215  427  508  612  729  816  1071  2782.
```

Hence every marked branch now has one common proof grammar:

```text
scalar closure, or finite 5 -> 3 -> 1 parity descent.
```

This is an entry theorem, not an all-root child closure.

## Adaptive ranked flags

On a scalar-hard branch put

```text
R_d=q_1+...+q_d,       1<=d<=4.
```

Subadditivity makes `R_d` a global `d`-label complement cap.  The condition

```text
R_d<(d+2)h/7
```

forces at least `d+1` labels of a hypothetical five-cover into the
corresponding finite high core.  Exact eligibility is

```text
d    eligible    failures
1       14806           0
2       14555         251
3       11699        3107
4        6180        8626.
```

Choosing the strongest eligible singleton-sum trigger branchwise gives

```text
strongest d       1       2       3       4
branches        251    2856    5519    6180.
```

These are sufficient scalar caps, not exact union caps.  In particular,
the `251` failures at `d=2` are the exact target for a globally sealed
pair-union cap; failure of `q_1+q_2<4h/7` does not imply failure of
`B_2<4h/7`.

## Root-level comparison and hostile controls

The competing root-level `p=6` parity condition

```text
q_1(E)<2|G_E|/7
```

holds on `3200/3432` roots and fails on `232`.  Its finite-core cutoff has
median `1531` but maximum `7775459`, whereas the uniform marked-suffix
route has maximum `2782`.  Proof cost, not only forced-label count, should
select the route.

The smallest positive scalar margin is

```text
5557/3988104120
```

at `E=(3,4,5,7,8,13,14)`, rank `4`, apex `17`.  The closest scalar failure
is

```text
-17/119819700
```

at `E=(3,5,6,7,10,12,13)`, rank `3`, apex `22`.  The even arithmetic
progression is the canonical hostile:

```text
E=(2,4,6,8,10,12,14), rank=1, apex=22,
margin=-27077/270270.
```

Hardness is not monotone in gate rank: `912` roots have a closed rank before
a later open rank.  The unique rank-13 hard branch is recorded in the
canonical transcript and hard ledger.

## Reproduction

The script hash-pins the THM-2896 gate source and transcript and the literal
residual engine.  It checks `165660` scalar/vector controls, `41415`
literal/direct carrier reconstructions, all locked counts and extrema, and
the complete digest

```text
9770a2ef0e7ae063a482e475f47047de08344cd1fc24880f742651d4c71d167d.
```

Ordinary and optimized full replays are byte-identical, with SHA-256
`e6261bb6c4dd13624ddccba32cac7a90f7973929c1c0dc4881743f4aec5d4206`.
The script uses explicit runtime guards rather than Python `assert`.

```text
04-computation/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.py
SHA-256 e0ff69252870f194549bba61289c1c5b15bef451e37a72836d71f9e71b1016e9

05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_census_codex_20260729.out
SHA-256 dbd3dc5a8c44a55957a6e1ce660ca0e89fcd70e6c0d06d5ba47dc3a22f40c680

05-knowledge/results/lrc14_j6_all_root_ranked_suffix_scalar_hard_ledger_codex_20260729.out
SHA-256 6be9a6c9218f3b42b2eea733c9050f5d35160664af0f19390337b3c5be57cb37
```

The hard ledger retains every body, stratum, gate size, branch rank, apex,
excluded prefix, carrier mass, component count, scalar margin, tail cutoff,
and globally sealed top five.  THM-2899 does not enumerate the resulting
H4 pair complexes, close the other `3422` roots, close the seven-body rung,
or prove LRC(14).
