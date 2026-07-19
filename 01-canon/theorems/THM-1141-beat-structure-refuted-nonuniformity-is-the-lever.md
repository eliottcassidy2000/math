---
id: THM-1141
title: THE BEAT/ALIGNMENT ROUTE IS REFUTED FOR THE CLUSTERED MAJORITY — clustered phase configurations are FROZEN, not sweeping — AND THE REAL LEVER IS GAP NONUNIFORMITY. (I) MY OWN NAMED-NEXT FROM THM-1140 IS WRONG: I proposed that for k_{j+1}/k_j ≈ 1 the combs beat with period 1/(k_{j+1}−k_j), so alignment (and a long gap) must occur somewhere inside a core component. It does not. The phase pair (i,j) sweeps ℓ·(k_j − k_i) full cycles across a component of length ℓ, and for real clustered quadruples that product is **0.027 to 0.54** — well below 1. The configuration is FROZEN across the whole component; whichever phase arrangement holds, holds throughout, and alignment is never forced. (II) YET SPREAD PHASES ARE NOT FATAL. The worst case of THM-1140 (core [1,3,5,6,7,8,11,12], killers 371/374/377/379) has phases at distance 0.325, 0.457, 0.239, 0.240 from the nearest integer — genuinely spread, not aligned — and still gives 7·k₄·L = 2.358. Its surviving component has length 0.0008888 against an aligned prediction of 0.0022616 and a uniform-4-spread prediction of 0.0002827: it sits **3.1× above the spread prediction**, not at it. (III) SO THE NAIVE WORST CASE NEVER OCCURS, and that is the whole content. Uniform interleaving of r combs gives gap = mean gap = 3/(7Σk) ≈ 3/(28k) at r=4, which is **below** the 4/(28k) the four-comb theorem needs — the naive bound genuinely fails. But across 500 clustered quadruples the minimum L·k₄ is **0.3584**, i.e. **3.34× the mean gap**, with **zero** spread-like cases (L·k₄ ≤ 0.25) and 213 aligned-like. Combs with *distinct* k interleave irregularly, and irregular interleaving always has a maximum gap exceeding its mean. (IV) THE TARGET LEMMA IS THEREFORE A NONUNIFORMITY STATEMENT: max gap ≥ (4/3)·mean gap suffices for the four-comb theorem, and the measured ratio is 3.34 — a factor 2.5 of room. This is three-distance-theorem territory, not beat/resonance territory
status: (I) REFUTED — my THM-1140 named-next, killed by the sweep measurement ℓ·d ≤ 0.54 < 1. (II),(III) MEASURED on clustered samples. (IV) is the identified target and its measured constant, **not a proof**. Uniform r=5 remains OPEN, and nothing here discharges a universal quantifier
source: kind-pasteur-2026-07-18-S128 (cont.69; owner: work the beat structure for the clustered majority)
depends_on:
  - THM-1140    # whose named-next this refutes and replaces
related: [THM-1097, THM-1094, MISTAKE-163, MISTAKE-164]
script: 04-computation/beat_structure_kps_S128c69.py (+ .out)
---

# THM-1141 — the beat route refuted; nonuniformity is the lever

## (I) My own proposal was wrong

THM-1140 ended by proposing the beat structure: for k_{j+1}/k_j near 1 the two combs differ
in phase by frac((k_{j+1}−k_j)t), which drifts with period 1/(k_{j+1}−k_j), so somewhere
inside a core component the phases should align and leave a long gap.

**The phases do not drift far enough.** Across a component of length ℓ the pair (i,j) sweeps
ℓ·(k_j − k_i) full cycles:

| core | ℓ | killers | min sweep | max sweep | regime |
|---|---|---|---|---|---|
| [1,3,5,6,7,8,11,12] | 0.02706 | 371,374,377,379 | 0.054 | 0.216 | **FROZEN** |
| [2,3,4,5,6,7,8,9] | 0.06746 | 400,401,402,403 | 0.067 | 0.202 | **FROZEN** |
| [1,3,5,6,7,8,11,12] | 0.02706 | 400,405,410,415 | 0.135 | 0.406 | **FROZEN** |
| [2,3,4,5,6,7,8,9] | 0.06746 | 400,440,480,520 | 2.698 | 8.095 | SWEEPS |

Sweeping only begins once the killers are *spread* — which is the regime THM-1140's gap
recursion already covers. In the clustered regime the configuration is **frozen**, and
alignment is never forced. The route is dead as proposed.

## (II) But spread phases turn out not to be fatal

THM-1140's worst case, dissected: the surviving component sits at phases

> 0.325, 0.457, 0.239, 0.240 from the nearest integer

— genuinely spread, not aligned. And yet 7·k₄·L = 2.358. Its length 0.0008888 compares to

| prediction | value |
|---|---|
| aligned, 12/(14k₄) | 0.0022616 |
| uniform 4-spread, 3/(28k₄) | 0.0002827 |
| **actual** | **0.0008888** |
| needed, 1/(7k₄) | 0.0003769 |

The actual gap is **3.1× the uniform-spread prediction**. The naive worst case is not what
happens.

## (III) Why: irregular interleaving beats uniform

Uniform interleaving of r combs gives every gap equal to the mean, 3/(7Σk) ≈ 3/(28k) at
r=4 — which is genuinely below the 4/(28k) the theorem needs, so the naive bound really does
fail. But combs with *distinct* k do not interleave uniformly, and any nonuniformity puts
the **maximum** gap above the mean. Across 500 clustered quadruples:

> min L·k₄ = **0.3584** (core [1,4,5,6,7,9,10,11], killers 550/553/554/558)
> aligned-like (L·k₄ ≥ 0.6): 213 · middle: 287 · **spread-like (≤ 0.25): 0**

Not one case reached the uniform-spread value. The measured minimum is **3.34× the mean
gap**, and 2.5× what the theorem needs.

## (IV) The target lemma

The four-comb theorem needs max gap > 1/(7k₄) while the mean gap is 3/(28k₄), so it suffices
to prove

> **max gap ≥ (4/3) · mean gap**

for the union of four combs with distinct moduli over a core component. The measured ratio
is 3.34, so there is a factor of 2.5 of room. This is a **three-distance-theorem** style
statement about the gap distribution of several incommensurate combs — *not* a beat or
resonance statement, which is where I was looking and where the answer is not.

## Honest status

Nothing here proves the four-comb theorem, and per MISTAKE-163/164 the measured constants do
not discharge the quantifier. **Uniform r=5 remains open.** What has changed is the target:
one dead route eliminated with a measurement that says exactly why (frozen, not sweeping),
and the live route identified as a gap-distribution lemma with a 2.5× margin to work in.

## Named next
- Prove max gap ≥ (4/3)·mean gap for a union of r ≤ 4 combs with distinct moduli on an
  interval. The three-distance theorem gives the exact gap structure for ONE comb plus an
  interval; the multi-comb version is the needed generalisation, and only a factor 4/3 is
  required where 3.34 is observed.
- The 213 aligned-like versus 0 spread-like split suggests the frozen configurations are
  biased *toward* alignment rather than uniform. If that bias has an arithmetic cause
  (killers just above 13·max P sharing small residues), it may be provable directly and
  would be a cleaner route than the general gap lemma.
