---
id: HYP-9210
title: "The latest first lonely time. For a primitive set v of n positive speeds let tau(v) = min{t > 0 : min_i ||t v_i|| >= 1/(n+1)} (the first time the rider t v enters the central cube; 'slide until the obstruction is reached'). tau*(2) = 4/9 and tau*(3) = 7/16 are PROVED; conjecture: tau*(n) = sup_v tau(v) equals 11/25, 19/42, 19/42, 49/104 for n = 4..7, every extremizer contains 1 and a multiple of n+1 and has its lonely set in (0,1/2] a single interval just below 1/2, and 1/2 - tau*(n) = Theta(1/n)"
status: >
  OPEN, with PROVED parts:
  - n = 2 exactly: tau*(2) = 4/9, attained only at the camel (1,3);
    tau(1,3k) = 1/3 + 1/(9k), tau(a,b) = 1/(3a) for a < b < 2a,
    tau(a,b) <= 2/(3a) for a >= 2, b > 2a (THM-4550 when promoted;
    chessboard_weave_20261006.md Theorem 2.4).
  - n = 3 exactly: tau*(3) = 7/16, attained only at (1,3,12) (sets containing 1: lane proof
    split on the middle speed; sets with smallest speed a >= 2: tau <= 3/8 -- a = 2 by LRC for
    four runners (Betke-Wills 1972, Cusick 1974) plus time symmetry, a = 3 by interval covering
    and two computed triples, a >= 4 by a ratio argument giving tau <= 3/(4a); see the note).
  - Lower bound for every n: the progression with 2 removed,
    F_N = {1,3,4,...,N} (n = N-1), has lonely set in (0,1/2] equal to
    [tau_N, 1/2 - 1/(4N)] with tau_N = 1/2 - 1/(2N) + 1/N^2 (N odd),
    1/2 - 1/(2N) + 1/(2N(N-1)) (N even); so tau*(n) >= 1/2 - 1/(2(n+1)) + O(1/n^2)
    and tau*(n) -> 1/2.
  - Upper bound: tau(v) < 1/2 whenever v is lonely at all (time symmetry plus an
    all-odd neighbourhood of 1/2); so tau*(n) <= 1/2 for every n where LRC is known.
  - Swap lemma: if T is tight with lonely times in (1/N)Z and a speed s of T is
    replaced by a multiple of N, every lonely time t of the new set has ||s t|| < 1/N.
  FINITE-EXACT (exact rational arithmetic, all primitive n-subsets of {1..B}):
  n = 2 (B = 400): 4/9; n = 3 (B = 150): 7/16 at (1,3,12); n = 4 (B = 60): 11/25 at
  (1,3,4,5), (1,3,5,8); n = 5 (B = 48): 19/42 at (1,5,6,7,8); n = 6 (B = 36): 19/42 at
  (1,3,4,5,7,18); n = 7 (B = 26): 49/104 at (1,5,6,7,8,11,13). Records stable from
  B = 12, 5, 8, 18, 13 respectively up to the box limit. Top-ten sets all contain 1 and
  a multiple of n+1 (100 of 100).
  Not monotone in n (n = 5 and 6 share 19/42).
source: mac-mini-2026-10-06-chessboard session, LRC lane (owner's chessboard prompt; "move ... unlimited amounts until an obstruction is reached")
related:
  - 05-knowledge/results/chessboard_weave_20261006.md (sections 2 and 6)
  - 01-canon/theorems/THM-4550-chessboard-lonely-runner-dictionary-and-the-knight-weave-law.md
  - 01-canon/theorems/THM-386-lrc-lonely-central-box-grounding.md (the central-box predicate)
  - 01-canon/theorems/THM-764-covering-small-period-signed-pair-deck-and-q25-refutation.md (no uniform good period; the discrete-time counterpart)
  - 01-canon/theorems/THM-3043-lrc-tight-instances-are-not-only-APs-and-what-FC2-does-not-transfer.md (tight instances; denominators n+1)
scripts:
  - 04-computation/experiments/chessboard_weave_20261006_lrc_BC_tau_period.py
  - 04-computation/experiments/chessboard_weave_20261006_lrc_B_probe.py
  - 04-computation/experiments/chessboard_weave_20261006_lrc_B_proofs.py
  - 04-computation/experiments/chessboard_weave_20261006_tau2_check.py
---

# HYP-9210 -- the latest first lonely time

**Object.** The lonely runner asks whether the rider `t -> t v (mod 1)` ever enters the
central cube `[1/(n+1), n/(n+1)]^n`. Its **first lonely time** `tau(v)` asks *when*. It is
the lonely-runner counterpart of the Collatz first-descent (stopping) time, and the
comparison is sharp: `tau(v) < 1/2` for every lonely `v`, while Collatz stopping times are
unbounded.

**Why late arrivals look like "the progression with 2 removed".** The progression
`{1, ..., N-1}` is tight and is lonely only at the multiples of `1/N`. Adding the speed `N`
kills exactly those times; removing `2` lets the times within `1/(2N)` of `1/2` survive.
So the rider has to slide almost to the half-way point. The camel `(1,3)` is the case
`N = 3`. Every census extremizer for `n = 2, 3, 4, 7` is such a swap of a tight set; for
`n = 5, 6` the swaps reach only `79/180` and `22/49`, below the census maxima, so another
mechanism is also in play there.

**Open.** Exact `tau*(n)` for `n >= 4` (`n = 2, 3` are proved); whether `1/2 - tau*(n)` is of order exactly `1/n`; whether every extremizer
contains 1 and a multiple of `n + 1`.

**Hostile/limit.** Statements about `tau` say nothing about LRC itself: they presuppose a
lonely time. They are first-passage facts at bounded `n`.
