---
id: HYP-2052
status: PARTIALLY PROVED (negative direction rigorous; positive conjecture open)
source: opus-2026-06-01-S551
related:
  - THM-369
  - THM-360
  - HYP-2039
  - HYP-2040
---

<!-- ⚠ ID COLLISION (MISTAKE-053): two files share HYP-2052. THIS file
(opus-S551 sieve-no-finite-completeness) is the first claimant and KEEPS the
number; the duplicate is oracle-S552 `lrc-loneliness-spectral-gap` (→ HYP-2065 in
a future cleanup). This is the canonical HYP-2052. -->

# HYP-2052: the division-point sieve has no finite completeness threshold for LRC(n)

## Rigorous part (PROVED)

For LRC(n) (observer + `m=n-1` distinct integer speeds; lonely iff
`||v_i t|| >= 1/n` for all `i`), the **division-point sieve** certifies
loneliness when some `t=a/q` has every `(v_i a) mod q` in the safe band
`min(r,q-r)*n >= q` (THM-369).

> **Construction.** Fix any `Q`. Take one speed `= lcm(2,...,Q)` and the rest
> arbitrary distinct positive integers. For every `q <= Q`, that speed is `≡ 0
> (mod q)`, hence in the danger band (`r=0`) for **every** multiplier `a`. So the
> set has **no division-point witness of denominator `<= Q`** ("sieve-blind up
> to Q").

Such blind sets are nonetheless lonely. Verified exactly for `n=14`,
`{lcm(2..Q)} ∪ {1,...,12}` at `Q = 14,16,18,20,24`: blind up to `Q`, but the
explicit witness `2/27` (and beyond) is checked safe ⇒ **unconditionally lonely**.

**Therefore:** for every tested bound `Q` there is an explicitly verified-lonely
n=14 set with no witness of denominator `<= Q`. The minimal witness denominator
is **unbounded** over lonely sets; a **bounded-modulus division-point sieve
cannot prove LRC(n)** for `n >= 8`. This is the converse-limitation of THM-369
and the sieve analogue of the S1/HYP-2039 fact that the *measure* bound is
trivial: both cheap certificates are provably too weak, so the content of LRC@14
is exactly their residual.

## Positive conjecture (OPEN)

> **HYP-2052+.** The minimal witness denominator `q_min(V)` of a primitive
> n-set `V` is controlled by its divisibility loading: `q_min(V)` is within a
> small additive/multiplicative constant of the least `q` that is *not* killed by
> some speed (the next prime power above the largest `q` that divides a speed).

Evidence (`lrc_n14_*_s551.py`): hardest min-witness-modulus is `34` for speeds
`<= 60` and `35` for speeds `<= 300`, in both cases on **fully loaded** sets
(every `q <= 14` divides some speed); random/structured sets resolve at `q <= 25`.
If proved, `q_min(V)` is bounded **in terms of the speeds**, turning the
(astronomical) bounded-speed exact check into a bounded-modulus check for any
fixed speed cap — a genuine finite reduction.

## Why it matters

Redirects strategy: stop seeking a finite-`Q` sieve proof of LRC@14 (provably
impossible) and instead either (a) prove HYP-2052+ to bound witnesses by speeds,
or (b) attack the set-vs-measure residual on the round-tournament slice (S525,
HYP-2039).

**See:** `07-reflections/lrc-n14-sieve-has-no-finite-completeness-witness-denominator-unbounded-s551.md`;
`04-computation/lrc_n14_multiprime_sieve_s551.py`,
`lrc_n14_sieve_completeness_probe_s551.py`,
`lrc_n14_witness_tournament_s551.py` (+ `.out` in `05-knowledge/results/`);
THM-369 (sieve), THM-360 (divisor necessity), HYP-2039 (set-vs-measure), S525.
