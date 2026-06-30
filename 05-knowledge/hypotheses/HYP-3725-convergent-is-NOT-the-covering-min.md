---
id: HYP-3725
title: REFUTED -- the convergent n/Phi_6(n) is NOT the covering-min for n>=7; the construction {1,..,n-2,n(n-1)} is NOT the extremal covering set. EXACT counterexamples (smaller-M primitive covering (n-1)-sets) at n=7,8,9: n=7 {1,2,5,6,7,8} M=2/13<7/43 (t=4/13, mod 2n-1=13); n=8 {1,4,5,6,7,11,16} M=2/15<8/57 (t=8/15, mod 2n-1=15); n=9 {1,3,4,5,7,11,18,32} M=4/33<9/73 (t=29/33, mod 33). So the GENERAL claim 'convergent = covering-min for all n>=7' (HYP-3701's n>=7 half; opus's '14/183 via 107-set scan') is FALSE. The mediant 2/(2n-1) (mod 2n-1) is the covering-min at n=7,8 (exact), and a sub-convergent value (4/33) wins at n=9. The whole recent arc (HYP-3703 tiling-optimality, 3704 three-routes, 3717 three-gap, 3722 observer-escape) attached the Kershner/Eisenstein/A2 (mod Phi_6) structure to the construction, which is NOT the minimal covering set. LRC NOT threatened (all candidates > 1/n; the covering-min being the mediant makes the floor margin TIGHTER, ~1/(n(2n-1)) instead of ~1/n^2). n>=10 very likely also refuted but my searches are too weak there to confirm; opus's 14/183 is suspect (likely a restricted-family min) but not DIRECTLY refuted at n=14
status: REFUTED (the convergent-as-covering-min claim) -- exact, reproducible counterexamples at n=7,8,9. The corrected covering-min trajectory at n>=10 and the n=14 value are OPEN. Court case CASE-convergent-not-covering-min filed; MISTAKE-087 logged.
source: mac-mini-2026-06-30-S47
related:
  - HYP-3701  # my own S42 claim: 'convergent for n>=7' -- the n>=7 half is now REFUTED (no transition at 7)
  - HYP-3722  # observer escape (convergent vs mediant) -- the premise (convergent=covering-min n>=7) is false
  - HYP-3717  # klein-S28 three-gap covering-min 14/183 -- restricted family; not the global min
  - HYP-3704  # the three routes -- aimed at the wrong (non-extremal) set
  - HYP-3703  # tiling-optimality -- the Kershner/Eisenstein frame is on the construction, not the covering-min
results:
  - 04-computation/covering_min_trajectory_macmini_20260630.py
  - 04-computation/covering_min_n14_attack_macmini_20260630.py
  - 04-computation/covering_min_hillclimb_macmini_20260630.py
  - 05-knowledge/results/covering_min_n14_attack_macmini_20260630.out
---

# HYP-3725 -- the convergent is NOT the covering-min (REFUTED at n=7,8,9)

The owner asked to think *past the projective plane*. Doing so computationally **overturns the premise of
the entire recent arc**: the construction `{1,..,n-2,n(n-1)}` (M = `n/Phi_6(n)`, the "convergent") is NOT
the covering-min for `n>=7`. The PG(2,6)-transition-at-`n=7` narrative is a mirage.

## The exact counterexamples (verified with `Fraction` arithmetic)
Each is a primitive covering `(n-1)`-set with `M` strictly below the construction's `n/Phi_6`:

| n | set | M | construction n/Phi_6 | witness t | mod |
|---|-----|---|----------------------|-----------|-----|
| 7 | {1,2,5,6,7,8} | **2/13** = 0.15385 | 7/43 = 0.16279 | 4/13 | 13 = 2n−1 |
| 8 | {1,4,5,6,7,11,16} | **2/15** = 0.13333 | 8/57 = 0.14035 | 8/15 | 15 = 2n−1 |
| 9 | {1,3,4,5,7,11,18,32} | **4/33** = 0.12121 | 9/73 = 0.12329 | 29/33 | 33 |

All cover `{2,..,n}`, are primitive, have exactly `n-1` elements, and `M < n/Phi_6`. `M` is exact (the
max-min is always attained at `t = k/d` with `d` a pairwise sum/difference -- a complete candidate set).

## What this kills
- **HYP-3701's `n>=7` half** ("construction = `n/Phi_6` is the covering-min for `n>=7`"): there is **no
  transition at `n=7`**. The mediant/sub-convergent keeps beating the construction past `n=6`.
- **opus's "14/183 covering-min via 107-set scan"**: that scan was a **restricted family** (near-construction
  variants); the spread-structured beaters above were not in it. `14/183` is a restricted-family min, not
  the global covering-min. (NOT directly refuted at `n=14` -- see "open" -- but strongly suspect.)
- The Kershner/Eisenstein/A2 (mod `Phi_6`) structure of HYP-3703/3704/3717/3722 is **real for the
  construction** but the construction is **not the minimal covering set**, so that structure is NOT the
  covering-min's structure. The "three routes," the "observer's convergent escape," the "tiling-optimality"
  were all aimed at the wrong (non-extremal) set.

## What this does NOT kill
The **LRC itself is untouched**. Every candidate -- mediant `2/(2n-1)`, `4/33`, convergent `n/Phi_6` -- is
`> 1/n`. The covering-min being *smaller* (the mediant) simply makes the floor margin **tighter**:
`2/(2n-1) - 1/n = 1/(n(2n-1)) ~ 1/(2n^2)`, versus the construction's `(n-1)/(n Phi_6) ~ 1/n^2`. The floor
`M >= 1/n` still holds with positive margin; only the *identity and value of the extremal set* change.

## The corrected open questions
1. **The true covering-min trajectory.** It is NOT a clean `2/(2n-1)`: `n=7,8` give the mediant `2/(2n-1)`
   (mod `2n-1`), but `n=9` gives `4/33` (mod 33), not `2/17`. Mixed moduli `{13,15,33}`. Needs a proper
   exhaustive search (hard: the `n=9` winner already has a speed of 32 ~ `3.5 n`, so low-speed exhaustion
   misses it) or a theory of the extremal family.
2. **The `n=14` covering-min.** Is `14/183` beaten? My random/greedy/hill-climb searches do NOT beat it at
   `n>=10`, but they demonstrably fail to reproduce even the `n=9` winner, so the failure is uninformative.
   The honest status: the GENERAL convergent claim is dead (n=7,8,9); the SPECIFIC `n=14` value is open.
3. **What governs the witness modulus** (`2n-1` at n=7,8; `33` at n=9)? `2n-1` is the signed-LRC modulus `C`;
   the covering-min may live there (or in a mix), NOT in the Eisenstein `Phi_6`.

## "Past the projective plane" -- the real lesson
Pushing past PG reveals the projective-plane story was a *post-hoc fit* to small `n<=6`, not the mechanism.
The covering-min is governed by a messier mod-`(2n-1)`-flavored structure (the signed-LRC modulus), not by
the clean Eisenstein/hexagonal/Kershner object. The construction's `n/Phi_6` is a genuine, beautiful covering
set -- the hexagonal AP with three-distance gaps is real -- but it is one nice covering, not the extremal
one. The team (myself included, HYP-3701) over-fit the elegant structure to the optimization.
