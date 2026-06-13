---
id: THM-349
name: good-cut-interval-union-recurrence
status: PROVED
date: 2026-05-29
session: opus-2026-05-29-S15
depends_on:
  - THM-336
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/GoodCuts.lean
scripts:
  - 04-computation/goodcut_interval_union_s15.py
results:
  - 05-knowledge/results/goodcut_interval_union_s15.out
  - 05-knowledge/results/lean_goodcuts_interval_union_verify_s15.out
---

# THM-349: Good-Cut Buckets Are an Interval-Union Gas

## Statement

Let a base-path tiling on `n` vertices have `N=n-1` legal cuts
`{1,...,N}`. Each upward tile `(hi,lo)` contributes the interval

```text
{lo+1, ..., hi}
```

of length at least `2`. Thus the good-cut set is exactly a union of selected
intervals of the `N`-cut path, where the allowed intervals are all contiguous
intervals of length at least `2`.

For a connected run of `L` good cuts, let `c_L` be the number of subsets of
length-at-least-2 intervals inside that run whose union covers the whole run.
Then

```text
c_L = sum_{A subset [L]} (-1)^|A| 2^{sum_R binom(|R|,2)},
```

where `R` ranges over the contiguous runs of `[L] \ A`.

Let

```text
B_N(x) = sum_g #{tilings on N cuts with g good cuts} x^g.
```

With `B_{-1}(x)=B_0(x)=1`, the bucket polynomial satisfies

```text
B_N(x) = B_{N-1}(x) + sum_{L=2}^N c_L x^L B_{N-L-1}(x).
```

Equivalently, if a good-cut set decomposes into runs of lengths
`L_1,...,L_r`, the number of tilings with exactly that good-cut set is

```text
prod_i c_{L_i},
```

and it is zero if any run has length `1`.

## Proof

The Lean-level membership theorem now states:

```text
k in goodCuts(b) iff exists upward tile t with k in cutInterval(t).
```

So every tiling selects a subset of intervals of length at least `2`, and the
good-cut set is their union.

Fix a candidate good-cut set `S`. No selected interval may cross a gap of `S`,
otherwise the union would include a cut outside `S`. Therefore interval choices
factor over the connected runs of `S`.

For one run of length `L`, inclusion-exclusion over uncovered cut positions
gives the displayed formula for `c_L`: after declaring `A` uncovered, the
remaining available positions decompose into runs `R`; an interval of length at
least `2` may be chosen independently inside each run, and a run of length
`r` contains `binom(r,2)` such intervals.

The recurrence for `B_N` follows by looking at the leftmost cut. Either it is
absent, giving `B_{N-1}`, or it begins a first good run of length `L>=2`. That
run contributes `c_L x^L`; if it does not reach the end, the next cut is forced
absent and the remaining suffix contributes `B_{N-L-1}`.

## Computed Values

Connected run cover counts begin:

```text
c_2=1
c_3=5
c_4=50
c_5=903
c_6=30773
c_7=2032504
c_8=264271477
```

The recurrence recovers the known THM-336 small buckets:

```text
[x^2] B_{n-1} = n-2
[x^3] B_{n-1} = 5(n-3)
[x^4] B_{n-1} = 50(n-4) + binom(n-4,2)
```

and extends the full bucket table without enumerating all tilings. For example:

```text
n=8: {0:1, 2:6, 3:25, 4:206, 5:2739, 6:61671, 7:2032504}
n=9: {0:1, 2:7, 3:30, 4:260, 5:3672, 6:92695, 7:4067314, 8:264271477}
```

## Verification

`04-computation/goodcut_interval_union_s15.py` verifies three independent
calculations agree for `n=3..8`:

1. the recurrence above;
2. summing component-cover weights over every good-cut subset;
3. direct enumeration of all tilings.

It also checks that the polynomial totals equal `2^binom(N,2)` through `n=13`.
The Lean targets `TournamentH7.GoodCuts` and `TournamentH7.Verify` build
successfully with the new interval-union membership lemmas; see
`05-knowledge/results/lean_goodcuts_interval_union_verify_s15.out`.

## Engineering Consequence

The good-cut bucket distribution is now computable in `O(N^2)` polynomial
operations after precomputing the `c_L`. This makes `g`, bucket probabilities,
and interval-run profiles cheap features for `tournament_tda.py` or any future
tiling-to-quotient extractor, and provides a checksum for exact tiling-census
builders.
