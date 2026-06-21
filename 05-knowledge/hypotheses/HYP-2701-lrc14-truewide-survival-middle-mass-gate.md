---
id: HYP-2701
title: LRC14 true-wide cap is a survival middle-mass gate
status: OPEN; exact S64 scout stored
source: codex-2026-06-20-S64
tangent: T936
depends_on:
  - HYP-2695
  - HYP-2699
  - THM-556
  - THM-535
related:
  - HYP-2703
  - THM-548
  - HYP-2680
  - HYP-2679
  - HYP-2700
  - HYP-2698
  - HYP-2697
  - HYP-2696
  - HYP-2684
  - HYP-2675
---

# HYP-2701 - True-Wide Survival Middle-Mass Gate

## Claim

For a row `E`, let `N_E(x)` be the number of missed inner sectors in
`{1,...,6}` and write its distribution as `p_t=Pr(N_E=t)`.
THM-556 gives

```text
U4(E) = p0 + p5 + 5 p6.
```

Equivalently,

```text
U4(E) <= floor_k=(k-6)/7
```

is the exact survival middle-mass inequality

```text
p1 + p2 + p3 + p4 - 4 p6 >= (13-k)/7.
```

Thus HYP-2695's true-wide floor gate is not primarily a raw Bonferroni
statement.  It asks whether true-wide geometry forces enough middle missed
sector mass to pay for the fully missed tail `p6`.  The exact cap only changes
the right side from `(13-k)/7` to `1-cap_k`, explaining why `k=8` can spend
the THM-535 dividend while `k>=9` should not need it.

The useful proof obligation is therefore:

```text
true-wide, k>=9
  => p1+p2+p3+p4 - 4p6 >= (13-k)/7.

true-wide, k=8
  => p1+p2+p3+p4 - 4p6 >= 1-cap_8 = 3637/5880,
     with finite dividend rows allowed to miss 5/7.
```

## Exact Scout

The S64 scout
`04-computation/lrc14_truewide_survival_middle_mass_codex_s64.py` stores its
output at
`05-knowledge/results/lrc14_truewide_survival_middle_mass_codex_s64.out`.
It keeps all arithmetic as `Fraction`s, asserts the identity

```text
floor_k - U4(E)
  = p1+p2+p3+p4 - 4p6 - (13-k)/7
```

on every audited row, and compares the same survival currency with the exact
cap right side `1-cap_k`.

Cap/floor currency table:

```text
k=8:  floor required=5/7, cap required=3637/5880, dividend=563/5880
k=9:  floor required=4/7, cap required=2025/4004, dividend=263/4004
k=10: floor required=3/7, cap required=36/91,    dividend=3/91
k=11: floor required=2/7, cap required=25/91,    dividend=1/91
k=12: floor required=1/7, cap required=1/7,      dividend=0
```

Exact true-wide box scans:

```text
k=8,  B18: truewide=16359, floor failures=3, cap failures=0
k=9,  B18: truewide=27020, floor failures=0, cap failures=0
k=10, B17: truewide=13299, floor failures=0, cap failures=0
k=11, B17: truewide=12441, floor failures=0, cap failures=0
k=12, B17: truewide=9009,  floor failures=0, cap failures=0
```

The worst audited `k=8` floor failure is

```text
E=(0,3,6,9,12,14,15,18)
currency=6193/8820
floor slack=-107/8820
cap slack=295/3528
p6=1/126, tail45=8/105, far_count=2.
```

The `k>=9` audited floor leaders are all cap-safe and floor-safe:

```text
k=9:  E=(0,4,6,8,10,12,14,15,16), floor slack=29/980
k=10: E=(0,2,4,6,8,10,12,14,15,16), floor slack=29/588
k=11: E=(0,2,4,6,8,9,10,12,14,15,16), floor slack=17/126
k=12: E=(0,4,6,8,9,10,11,12,13,14,15,16), floor slack=85469/420420
```

The far-count ledger is the main structural signal:

```text
k=8,  B18: far_count=2 best=-107/8820;       far_count=3 best=5/504
k=9,  B18: far_count=2 best=29/980;          far_count=3 best=239719/5717712
k=10, B17: far_count=2 best=29/588;          far_count=3 best=100855/952952
k=11, B17: far_count=2 best=17/126;          far_count=3 best=5213/33320
k=12, B17: far_count=2 best=85469/420420;    far_count=3 best=155489/720720
```

Thus the risky true-wide rows are barely true-wide: exactly two speeds exceed
`14` in every tight audited family.  Rows with at least three far speeds have
a visibly larger survival-currency margin in the boxes scanned.

## Two-Far Boundary Addendum

The addendum scout
`04-computation/lrc14_twofar_survival_currency_boundary_codex_s64.py` stores
its output at
`05-knowledge/results/lrc14_twofar_survival_currency_boundary_codex_s64.out`.
It transports THM-548's boundary-value lens from `p0` to the HYP-2701 currency
`C=p1+p2+p3+p4-4p6=1-U4`.

For a bounded core `B` and two decorrelated far runners, the missed-count
process is an exact death chain: if `B` misses `t` inner sectors, each far
runner hits one of the seven sector colors uniformly and only hits to missed
sectors reduce `t`.  The exact transition is

```text
Pr(t -> s after r hits)
  = binom(t,s) * sum_{j=0}^{t-s} (-1)^j binom(t-s,j) ((7-s-j)/7)^r.
```

This gives a fully decorrelated boundary currency `C_boundary(B,2)`.  The
bounded-core scan over `B subset {0,...,14}` found it floor-safe with explicit
minimum margins:

```text
k=8:  min boundary margin=1879/10290,  core=(0,1,2,3,4,5)
k=9:  min boundary margin=569/3430,    core=(0,1,2,3,4,5,6)
k=10: min boundary margin=5717/36015,  core=(0,1,2,3,4,5,6,7)
k=11: min boundary margin=5317/24010,  core=(0,1,2,3,4,5,6,7,8)
k=12: min boundary margin=35543/123480, core=(0,1,2,3,4,5,6,7,8,9)
```

Actual two-far rows spend this positive boundary margin through a negative
deviation

```text
C(B union {u,v}) - C_boundary(B,2).
```

For the `k>=9` tight leaders, the boundary-margin/deviation ledger is:

```text
k=9  leader (0,4,6,8,10,12,14,15,16):
  boundary margin=18119/72030, deviation=-6395/28812, slack=29/980.

k=10 leader (0,2,4,6,8,10,12,14,15,16):
  boundary margin=5717/36015, deviation=-15763/144060, slack=29/588.

k=11 leader (0,2,4,6,8,9,10,12,14,15,16):
  boundary margin=28976/108045, deviation=-9599/72030, slack=17/126.

k=12 leader (0,4,6,8,9,10,11,12,13,14,15,16):
  boundary margin=5584429/13733720, deviation=-152315/749112,
  slack=85469/420420.
```

Thus the two-far lemma has a sharper form:

```text
C_boundary(B,2) - (13-k)/7
  >= -(C(B union {u,v}) - C_boundary(B,2)).
```

The left side is now an exact finite bounded-core margin.  The right side is a
signed two-far deviation concentrated in low relation-distance pairs: in every
audited layer the tight rows have relation distance `1`, `2`, or `3` for
`u-v` at coefficient height `<=4`, usually consecutive or near-consecutive far
pairs.  This is precisely where THM-548/HYP-2679 expects Freiman/scale finite
atlases and signed Abel bounds, not a generic absolute discrepancy estimate.

## Proof Route

1. **Two-far survival-currency lemma.**  Prove the exact inequality
   `p1+p2+p3+p4-4p6 >= (13-k)/7` for true-wide `k>=9` with exactly two far
   speeds.  The tested leaders are low-growth, low-height relation packets, so
   this is the place to combine HYP-2679's two-far curvature, HYP-2676/HYP-2677's
   finite Ruzsa/tournament packet atlas, and THM-558's transfer-tax ledger.
2. **Three-or-more-far margin lemma.**  Treat `far_count>=3` as the easier
   branch.  The audited margin suggests that additional far carriers increase
   middle survival mass faster than they increase the fully-missed tail.  This
   should connect to HYP-2680's signed multi-far Newton hierarchy and HYP-2684's
   BV/decorrelation error.
3. **Cap-dividend templates for k=8.**  The three floor failures in `B18` are
   all two-far and cap-safe.  Do not force them through the `k>=9` no-dividend
   floor lemma; classify them as finite dividend templates and charge the exact
   THM-535 margin.
4. **Lift near equalities before scalarizing.**  The survival currency forgets
   sector labels.  Rows near equality must be lifted back to state-word,
   transfer-tax, residual-profile, and relation-lattice coordinates before any
   final theorem is claimed.

## Tournament Analysis

The S64 quotient deliberately challenges the assumption that tournament
vertices should be runners, arcs, or even sector colors.  Vertices are exact
proof obligations / row-family leaders.  The pairwise observable is smaller
floor slack, and the switch/gauge orients toward the row spending more
cap-floor currency.

The audited tournament is transitive:

```text
score_hist={0:1,1:1,...,13:1}
directed_3cycles=0
```

Its Hamiltonian risk path starts with the three `k=8` floor failures, then the
`k=9` leader, then the `k=10` leader.  This says the proof has a single visible
risk spine, not several incomparable obstruction classes, at least in the
audited boxes.

The quotient preserves the cap/floor implication, because the scalar survival
currency is exactly equivalent to the `U4` comparison.  It destroys sector
labels, cyclic order, and relation-lattice phase, so it is only the gate
coordinate, not the full proof state.

## Relation To Incoming Work

HYP-2702/T937 landed during the same S64 run and is strong signal in the same
direction.  There the sparse residual-tail failures disappear only after
changing basis to miss-zeta generated-context words; here the true-wide
cap-floor gate becomes legible only after changing basis from `U4` to survival
middle mass.  Both routes say the next proof should choose the right coordinate
system before taking scalar bounds.

Later incoming HYP-2703 is a compatible warning from the `measS7` extremality
side.  Its 7-band Sturmian identity is exact, but per-band consec dominance and
greedy index compression both fail; the proof-bearing fact is an aggregate
signed seven-band balance.  That mirrors HYP-2701's two-far addendum: the
boundary margin is positive, but the actual row is boundary value plus a signed
deviation, and proving either term independently in an over-refined quotient is
not enough.

## Status

No LRC14 proof is claimed.  The sharp target is now the two-far
survival-currency lemma, followed by a separate `far_count>=3` margin lemma and
finite `k=8` dividend templates.
