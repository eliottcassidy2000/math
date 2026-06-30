---
id: HYP-3615
title: SEARCH + SYNTHESIS of the small-measure / thin-coverage regime for the runners -- the regime is not a corner of the LRC but its HEART, anchored by a clean (verified n=3..8, proved general) fact: the consecutive extremal W={1,..,n-1} has lonely measure EXACTLY 0 at every n, its lonely "set" being the phi(n) touch-points {k/n : gcd(k,n)=1} = the UNITS mod n (in phi(n)/2 antipodal pairs); the union bound goes NEGATIVE ((2-n)/n for n>=3, =-0.857 at n=14: the n-1 danger combs each of measure 2/n can cover all of [0,1)), so measure-based arguments are structurally doomed at the extremal. The catalogue (klein-S16 inf R'=0; klein-S18 rho_0->0 = the existence question; HYP-3580 the cusp = m_R->0; HYP-3562 Lebesgue=0 but count=phi; HYP-3548 two razor-thin lines; THM-590 the finite cyclotomic min) all says the SAME thing: at the small-measure regime measure VANISHES and EXISTENCE/COUNTING carries it -- the units, the odd cycle, the discrete non-bipartiteness certificate (THM-590/HYP-3606). The small-measure search REVEALS why the proof had to abandon measure and move to the discrete certificate
status: SEARCH + SYNTHESIS (catalogues the prior small-measure thread) + VERIFIED anchor (consecutive extremal lonely measure=0, lonely set=phi(n) units, n=3..8; proved: t=k/n is lonely iff gcd(k,n)=1, else runner j=n/gcd sits at 0). Generalizes klein-S8's n=14 {1..13}->units fact to all n. Not a new proof of LRC.
source: mac-mini-2026-06-30-S38
related:
  - HYP-3597  # klein-S16: inf R'=0 over the infinite family (measure vanishes); existence != measure
  - HYP-3599  # klein-S18: rho_0 = the existence question, ->0 at the cusp
  - HYP-3580  # mac-mini S29: the cusp (m_R->0) = the binding regime
  - HYP-3562  # mac-mini S23: Lebesgue=0 but counting=phi at the extremal (the two measures)
  - HYP-3606  # mac-mini S37: the least-eigenvalue/non-bipartiteness certificate (the discrete content that survives)
  - HYP-3548  # the two razor-thin lines (gap 10% safe; measure floor the thin one)
results:
  - 04-computation/intransitivity_among_n_things_macmini_20260630.py  # union-bound failure (related)
---

# HYP-3615 -- the small-measure regime is the heart of the Lonely Runner

The owner asked to search the prior work on very small areas of coverage/measure for the runners. The
search (an Explore over ~25 files) plus a fresh computation say one thing: **the small-measure regime is
not a corner of the LRC -- it IS the LRC.**

## The anchor (verified n=3..8, proved): the extremal has measure-ZERO lonely set
For the consecutive extremal relative speeds `W = {1, 2, ..., n-1}` (gap `1/n`):
- **lonely measure = EXACTLY 0** at every `n` (fine-grid verified n=3..7; structurally exact).
- the lonely "set" is the `phi(n)` **touch-points** `{k/n : gcd(k,n)=1}` = the **units mod n**, in
  `phi(n)/2` antipodal pairs `{k/n, (n-k)/n}` (verified n=3..8: counts `2,2,4,2,6,4 = phi(n)`).
- **Proof:** at `t=k/n`, `||w_j t|| = ||jk/n||`; the values `{jk mod n}` over `j=1..n-1`. If `gcd(k,n)=1`
  they are all nonzero residues, so `min_j ||jk/n|| = 1/n` (boundary, lonely). If `gcd(k,n)=d>1` then
  `j = n/d` gives `jk ≡ 0`, a runner exactly at the origin (`||0||=0 < 1/n`, NOT lonely). So `t=k/n` is
  lonely iff `k` is a unit. This GENERALIZES klein-S8's `n=14` fact (`{1..13}` touch-points = `(Z/14)*` =
  `{1,3,5,9,11,13}`, `phi/2=3` antipodal pairs) to all `n`.

## Why measure is doomed at the extremal: the union bound goes NEGATIVE
Each danger comb `D_j = {t : ||w_j t|| < 1/n}` has measure `2/n`; there are `n-1` of them. The union
(lonely-complement) bound is `1 - 2(n-1)/n = (2-n)/n`, which is `< 0` for all `n >= 3` (`= -0.857` at
`n=14`, where `13` combs of measure `1/7` sum to `1.857`). So the `n-1` danger zones CAN cover all of
`[0,1)`, and at the extremal they DO (up to the measure-0 touch-points). Any measure-from-below argument is
structurally hopeless at the extremal -- the measure really is 0.

## The catalogue (all of it says the same thing)
| theme | statement | where |
|---|---|---|
| **inf measure = 0** | over the INFINITE covering family `inf R' = 0` (R'=0.34->0.25->...->0 adding high speeds) | klein-S16 HYP-3597 |
| **rho_0 IS existence** | top descent level `rho_0 = m_S/(m_{O_0}m_{S'})` contains `m_S`, ->0 at the cusp; `rho_0>0 <=> LRC(S)` | klein-S18 HYP-3599 |
| **the cusp** | binding worst-case = `m_R->0` = the 4 cusps of `X_0(14)` (d=1,2,7,14); apex-7 hardest | mac-mini HYP-3580 |
| **Lebesgue 0, count >0** | LRC4 `{1,2,3}`: Lebesgue floor 0 but lonely count `=phi(4)=2`; Borsuk-Ulam index detects it | mac-mini HYP-3562 |
| **the units** | extremal touch-points = `(Z/14)*` in `phi/2` antipodal pairs (the counting/BU witness) | klein-S8 HYP-3571 |
| **CV blows up** | `CV(N_R)^2` unbounded as `m_R->0` (sup 8.74), speed-7 amplified -- the WRONG coordinate | klein-S4 HYP-3554 |
| **razor-thin** | gap line 10% safe (`7/89`,`14/183`); measure-floor line the thin/open one | HYP-3548 |
| **finite min** | over the finite `Z_7`-cores `inf = 4cos^2(3pi/7)>0` (vs infinite-family `0`) | THM-590, HYP-3597 |
| **p_0 -> ceiling** | consec maximizes `p_0` (prob. of being lonely); stays closest to the lonely ceiling | HYP-2635 |

## The synthesis -- the regime IS the proof's pivot
At the small-measure regime the lonely MEASURE genuinely vanishes (extremal = exactly 0; cusp/over-covering
-> 0). So measure cannot be the proof. What survives, and what the whole recent program converged on, is the
DISCRETE / EXISTENCE / COUNTING side:
- **the units** `(Z/n)*` are the measure-0 lonely touch-points (the Borsuk-Ulam / `sigma`-odd witness);
- **the odd cycle** (the apex `C_p`) is present even where its measure is 0 (intransitivity, HYP-3602);
- **the non-bipartiteness certificate** `lambda_min(2I+A(C_p)) = 4sin^2(pi/2p) > 0` (THM-590, HYP-3606) is
  the discrete floor that does NOT see the vanishing measure;
- **the 2-adic descent** finitizes the (measure-0-infimum) infinite family into the finite `Z_7`-cores
  where a true positive minimum is attained.
So "very small areas of coverage/measure" is exactly the place the LRC lives: the extremal is measure-zero,
the union bound is negative, and the conjecture holds there by COUNTING the units / detecting the odd cycle,
never by bounding a measure. The search confirms the frame the project arrived at (klein-S16's existence !=
measure, HYP-3606's discrete certificate) is forced by the small-measure regime, and pinpoints the cleanest
witness: the `phi(n)` units, at the consecutive extremal, at every `n`.

## What it buys
A complete catalogue of the small-measure thread, anchored by a clean general fact (extremal lonely measure
= 0, lonely set = the `phi(n)` units) that makes the regime concrete and explains -- via the negative union
bound -- WHY measure must be abandoned for counting. It unifies klein-S16/S18, HYP-3562/3580, and the
certificate (HYP-3606) as one response to one regime: when the area is zero, you count the units / the odd
cycle.
