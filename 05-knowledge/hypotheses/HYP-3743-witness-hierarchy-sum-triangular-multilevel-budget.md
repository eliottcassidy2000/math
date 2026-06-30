---
id: HYP-3743
title: Summing the witness hierarchy + the constant-residue budget over ALL moduli. The witness hierarchy is MULTI-LEVEL (radius bands). RADIUS-1 LEVEL: a covering set missing small speed k exposes the unit-pair {k,D-k} mod every band modulus D in (n-2+k, 2n-2], i.e. E(k)=n-k exposed moduli; SUMMED over k, Sum_{k=1}^{n-2} E(k) = T_(n-1)-1 (a TRIANGULAR number -- the project's triangle, verified n=10..16: 44,65,90,119). E(k)=n-k DECREASES in k, ordering the M-jumps (missing 1 most-exposed -> M=2/17; missing n-2 least -> 2/25). CONSTANT-RESIDUE EXTENDS TO ALL MODULI: k mod D = k for every D>k (not just primes), so the n-2 dense-core speeds are universal pair-coverers mod EVERY band modulus (the (n-2)*|band| pair-covers done by n-2 speeds = the budget's efficiency). BUDGET/SUMMING: missing core speed k over-commits the n-1 speeds at radius-1 [re-cover q=k if k>=2] + [resonances n-1,n] + [E(k)=n-k exposed pairs] -- the forced killer 182 covers ~0 exposed (verified) -- UNLESS one huge CRT speed covers all E(k) band-1 pairs at once; but that speed is SCATTERED at HIGHER radius levels (S57: the missing-1 CRT escape's hole is at mod 85, radius 12, OUTSIDE the band). So summing the hierarchy ACROSS LEVELS (constant-residue => speed k is the universal coverer at EVERY level) forces {1,..,n-2} subset S => covering-min(n>=12)=n/Phi_6
status: FRAMEWORK + the radius-1 sum EXACT. E(k)=n-k and Sum=T_(n-1)-1 are exact (verified n=10..16); the killer covers ~0 exposed (verified); the constant-residue-all-moduli is elementary (k mod D=k). The MULTI-LEVEL closure (the CRT escape caught at higher radius levels for ALL k, not just k=1) is the mechanism + S57 evidence, NOT a complete proof -- the cross-level bookkeeping is open. Strongest for small k (radius-1 E(k) large); the tightest large-k cases lean on the higher levels.
source: mac-mini-2026-06-30-S58
related:
  - HYP-3741  # the constant-residue principle (this sums it across the hierarchy + extends to all moduli)
  - HYP-3740  # the LRC14 hard core = the lowness lemma (this is its budget/summing proof attempt)
  - HYP-3736  # klein-S39 the killer-or-transversal budget (band primes); this extends to all band moduli + multi-level
  - THM-523   # the q-witness (the radius-0 layer); the radius-1+ levels are the lowness hierarchy
results:
  - 04-computation/witness_hierarchy_sum_budget_macmini_20260630.py
  - 05-knowledge/results/witness_hierarchy_sum_budget_macmini_20260630.out
---

# HYP-3743 -- summing the witness hierarchy (triangular) + the multi-level constant-residue budget

Working the owner's two asks -- *sum the witness hierarchy across primes* and *extend the constant-residue
budget to all moduli* -- gives an exact triangular sum at the radius-1 level and a multi-level budget framework
for the lowness lemma.

## The witness hierarchy is multi-level
The radius-demand criterion (klein-S38) layers by radius: radius-0 (`D<=n-1`, the THM-523 resonances),
radius-1 band (`D in (n, 2n-2]`), radius-2 band, ... A covering set with `M<=n/Phi_6` must satisfy every level.
The **lowness** obstruction lives in the radius-`>=1` levels.

## Summing the radius-1 level -- a triangular number
At the radius-1 band, the dense core `{1,..,n-2}` covers the pair `{j, D-j} mod D` for `j=1,..,n-2` and every
`D in (n,2n-2]` (each speed `j` by its constant residue). A set **missing speed `k`** leaves `{k, D-k}`
uncovered exactly when `D-k > n-2`, i.e. `D in (n-2+k, 2n-2]` -- so

  **E(k) = #exposed band moduli = n - k.**

`E(k)` decreases in `k` (missing 1 exposes the most, missing `n-2` the least), which *orders the M-jumps*
(missing 1 -> `2/17`, ..., missing 12 -> `2/25`). Summed over the small speeds:

  **Sum_{k=1}^{n-2} E(k) = Sum_{j=2}^{n-1} j = T_(n-1) - 1** (a TRIANGULAR number).

Verified: `n=10,12,14,16 -> 44, 65, 90, 119`. The witness hierarchy's radius-1 weight is the `(n-1)`-th
triangular number minus 1 -- "everything is the triangle," now at the LRC binding.

## Extending the constant-residue budget to ALL moduli
The S57 constant-residue principle (`k mod p = k` for primes `p>k`) extends verbatim to **every** modulus:
`k mod D = k` for all `D > k`. So the `n-2` core speeds are **universal pair-coverers mod every band modulus**
(`(n-2)*|band|` pair-covers with `n-2` speeds). Replacing a core speed `k` loses all `E(k)=n-k` of its
exposed pair-covers; a scattered (large) speed restores only the `D | (w∓k)`.

## The budget / summing -- forcing {1,..,n-2}
Missing core speed `k` over-commits the `n-1` speeds at the radius-1 level: `[re-cover q=k if k>=2]` +
`[resonances n-1,n]` + `[E(k)=n-k exposed pairs]`. The forced killer `n(n-1)` covers `~0` exposed pairs
(verified: missing-1 -> 0 of 13; missing-12 -> 0 of 2). The only escape is **one huge CRT speed** covering all
`E(k)` band-1 pairs at once -- but that speed is **scattered at higher radius levels** (S57: the missing-1 CRT
escape's lonely hole sits at `mod 85`, radius 12, *outside* the band). By constant-residue, speed `k` is the
universal coverer at **every** level; **summing the hierarchy across levels leaves no escape**, forcing
`{1,..,n-2} \subseteq S`, hence (HYP-3740 Step 2) `covering-min(n>=12) = n/Phi_6`.

## Honest status
- **Exact:** `E(k)=n-k`, `Sum = T_(n-1)-1`, the killer-covers-~0, the constant-residue-all-moduli.
- **Framework:** the multi-level closure -- that the CRT escape is caught at a higher radius level for *every*
  `k` (not just the verified `k=1`) -- is the mechanism + S57 evidence, not a complete proof. The cross-level
  bookkeeping (how the radius-`r` exposure of a level-1-CRT speed is itself uncoverable) is the remaining gap.
- **Strongest** for small `k` (`E(k)` large at radius-1); the tightest large-`k` cases (`E` small) lean on the
  higher levels.

## What it buys
Turns "sum the witness hierarchy" into the exact triangular weight `T_(n-1)-1` at the radius-1 level, extends
the constant-residue budget from band *primes* (klein-S39) to *all* band moduli, and frames the lowness lemma
as a **multi-level summing** in which the constant-residue small speed is the universal coverer at every level
-- the precise remaining work being the cross-level closure. This is the budget half of the LRC14 hard core
(HYP-3740), now with an exact level-1 accounting.
