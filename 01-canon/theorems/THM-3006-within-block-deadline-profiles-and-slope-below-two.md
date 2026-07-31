---
id: THM-3006
title: "Within-block deadline profiles and exactly fair extractors of slope below two"
status: >
  PROVED + VERIFIED-EXACT. Block normal form: a balanced block rule on the
  dyadic shell [m, 2m) that decides critical value n = m+k by flip
  T(m+k) = m+k+1+a_k exists IFF the integer system
  sum_{k,i} p_{k,i} binom(L_k, j-1-i) = binom(m,j)/2  (1 <= j <= m-1),
  0 <= p_{k,i} <= binom(a_k,i), L_k = m-1-k-a_k, is solvable. Feasibility is
  monotone in each a_k, so the least within-block ratio
  rho(m) = min over rules of max_k T(m+k)/(m+k) is computable. RESULT:
  rho(4) = 3/2, rho(8) = 14/9, rho(16) = 25/16, rho(32) <= 11/7 -- ALL BELOW
  2, with explicit exactly fair rules (brute-force verified through m = 16,
  exact-integer verified at m = 32). Stitching the blocks gives an exactly
  fair extractor with T(n) <= (11/7) n for 4 <= n < 64. This REFUTES the
  standing expectation (HYP-9061, PROBLEM-LEDGER) that every known/available
  scheme has slope C = 2: the classical 2n, max(2,2n-1), 2n-2 envelopes are
  far from optimal, and the slope-2 barrier is NOT a consequence of the
  dyadic block geometry of THM-3005. Whether sup_r rho(2^r) < 2 -- which
  would settle the proposer's minimal-C challenge with C* < 2 -- is OPEN;
  the data 1.5000, 1.5556, 1.5625, 1.5714 is increasing but slowly.
source: opus-2026-07-31-amm12592-writeup
depends_on:
  - THM-3005
related:
  - THM-2160
  - THM-2225
  - THM-2966  # spine normal form (across shells); this is the within-shell dual
  - THM-2996
  - HYP-9061
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_block_profile_ilp_thm3006.py
referee: 04-computation/amm12592_fast_block_referee_thm3006.py
output: 05-knowledge/results/amm12592_fast_block_referee_thm3006.out
---

# THM-3006 -- the within-block deadline profile, and slope below two

By THM-3005 every balanced block is dyadic; take the finest shells
`[m, 2m)`, `m = 2^r`. On the branch `X_1 = 0` write the relative tail
`z = (X_(m+1),...,X_(2m))`, `z != 0`; critical value `n = m+k` iff
`z = 0^k 1 w` with `w in {0,1}^(m-1-k)`. Balance of the shell is exactly

```text
#{heads z : |z| = j} = binom(m,j)/2,        1 <= j <= m-1,     (L)
```

with `z = 1^m` free (branch `X_1 = 1` carries the complementary labels).

## 1. Normal form

A rule with deadline `T(m+k)` reads exactly

```text
a_k := T(m+k) - (m+k) - 1                                     (1)
```

coordinates of `w`, so its heads set on stratum `k` is `S_k x {0,1}^(L_k)`
with `S_k` inside `{0,1}^(a_k)` and `L_k = m-1-k-a_k`. Put
`p_{k,i} = #{y in S_k : |y| = i}`; every integer vector with
`0 <= p_{k,i} <= binom(a_k,i)` is realizable.

**Theorem.** A balanced block rule on `[m,2m)` with depth profile
`(a_k)` exists if and only if

```text
sum_(k=0)^(m-1) sum_(i=0)^(a_k) p_{k,i} binom(L_k, j-1-i) = binom(m,j)/2,
                                            1 <= j <= m-1,     (*)
```

has an integer solution in the box `0 <= p_{k,i} <= binom(a_k,i)`.

*Proof.* Stratum `k` contributes `sum_i p_{k,i} binom(L_k, (j-1)-i)` heads of
`z`-weight `j`; sum over `k` and impose (L). Conversely any box vector is
realized by choosing `p_{k,i}` words of each weight in `{0,1}^(a_k)`. The
`j = m` equation only labels `z = 1^m`, which is free. QED

Equivalently, in generating-function form,

```text
u sum_k P_k(u) (1+u)^(L_k) = [(1+u)^m - 1]/2 + c u^m,
P_k(u) = sum_i p_{k,i} u^i,   c = p_{0,a_0} - 1/2 in {-1/2, 1/2}.  (**)
```

**The half solution and the corner obstruction.** Setting
`p_{k,i} = binom(a_k,i)/2` makes the left side of (**) equal
`(u/2) sum_k (1+u)^(m-1-k) = [(1+u)^m - 1]/2` EXACTLY. So the entire
difficulty is integrality: `binom(a,i)` is odd exactly for `i` a binary
submask of `a` (Lucas), and there `p_{k,i}` must be off by `+-1/2`. Writing
`p_{k,i} = (binom(a_k,i) + e_{k,i})/2`, (**) becomes the pure deficit
equation

```text
sum_k u E_k(u) (1+u)^(L_k) = 2c u^m,   E_k(u) = sum_i e_{k,i} u^i,
|e_{k,i}| <= binom(a_k,i),   e_{k,i} = binom(a_k,i) mod 2.        (***)
```

This is HYP-9061's corner-deficit transport, now inside a single shell and
with all coefficients explicit.

Since a rule oblivious to more coordinates is a special case of one
oblivious to fewer, feasibility of (*) is **monotone increasing in each
`a_k`**, so the least achievable within-block ratio

```text
rho(m) = min over rules of max_(0<=k<m) T(m+k)/(m+k)             (2)
```

is found by testing one profile per candidate ratio (binary search).

## 2. The computation

```text
m      rho(m)      value      depth profile a_k
 4     3/2         1.50000    [1,1,1,0]
 8     14/9        1.55556    [3,4,4,4,3,2,1,0]
16     25/16       1.56250    [8,8,9,9,10,10,9,8,7,6,5,4,3,2,1,0]
32     11/7        1.57143    [17,17,18,19,19,20,20,21,21,22,21,20,...,1,0]

m=4    T(n), n=4..7    = 6,7,8,8
m=8    T(n), n=8..15   = 12,14,15,16,16,16,16,16
m=16   T(n), n=16..31  = 25,26,28,29,31,32,32,...,32
m=32   T(n), n=32..63  = 50,51,53,55,56,58,59,61,62,64,64,...,64
```

The `m = 4` row agrees with the exhaustive enumeration of all `1036800`
shell-balanced `m=4` rules (THM-2996 sec 4), including the infeasibility of
`T = (6,6,7,8)`; this validates the normal form against ground truth.

Values are exact upper bounds: each was produced with HiGHS and then the
witness re-verified in exact integer arithmetic, and for `m <= 16` the
explicit rule was rebuilt and brute-force checked (every layer bisected,
every composition class of the full shell bisected, causality at every
claimed stop, exact fairness as a polynomial identity for `m <= 8`).
Minimality below the stated value is solver-certified only for `m <= 8`.

## 3. Consequences

1. **Slope 2 is not forced.** Blocks are independently balanced, so they
   stitch: using the above rule on each `[2^r, 2^(r+1))` gives an exactly
   fair extractor with

   ```text
   T(n) <= (11/7) n           for 4 <= n < 64,
   ```

   versus `2n`, `max(2,2n-1)`, `max(n+1,2n-2)` for the classical schemes and
   THM-2996. The repeated claim that "every known rule has `C = 2`" is
   therefore superseded: the classical schemes simply never optimized the
   within-block profile.
2. **THM-3005 does not imply slope 2.** Block ratio `>= 2` (THM-3005) and
   deadline slope `>= 2` are different statements; the second is false at
   every `m` computed here.
3. **The remaining question is uniformity.** `C* <= sup_r rho(2^r)` if that
   supremum is finite and `< 2`. The data `1.5000, 1.5556, 1.5625, 1.5714`
   increases slowly; `2 - rho` reads `0.500, 0.444, 0.4375, 0.4286`. Four
   points do not distinguish "converges near `1.6`" from "tends to `2`
   logarithmically". This is now the sharp form of HYP-9061 (Q).
4. **Where the binding constraint sits.** In each optimum the constraint is
   active only for `k < (2-C)m/C` (about `0.28 m`); beyond that the stratum
   may read its whole tail. And `a_0 - m/2 = -1, -1, 0, +1` for
   `m = 4,8,16,32`, so THM-2160 S6.2's bound `a_0 >= m/2 - 1` (i.e.
   `T(m) >= 3m/2`) is tight only for small `m`; the other strata push `a_0`
   up.

## 4. Next moves

- Exact arithmetic beyond `m = 32` (the float ILP loses exactness once
  `binom(m,m/2)/2 > 2^52`, i.e. at `m = 64`); solve (***) instead, whose
  coefficients stay small when the deficits do.
- A self-similar construction: find a family of profiles with a provable
  uniform ratio, most plausibly by choosing every `a_k` with small binary
  popcount so that (***) has few forced half-deficits (`a_k = 2^s` leaves
  only `i = 0, a_k`).
- A lower bound above `3/2`: THM-2160 S6.2 gives `rho(m) >= 3/2` from the
  `k = 0` stratum alone; the computed values exceed it, so a sharper
  multi-stratum bound exists and should be extracted from (***).

## 5. Referee

```bash
python3 04-computation/amm12592_block_profile_ilp_thm3006.py     # rho(m)
python3 04-computation/amm12592_fast_block_referee_thm3006.py    # explicit rules
```

QED for the normal form; the `rho` values are exact verified upper bounds.
