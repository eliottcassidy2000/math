---
id: THM-1140
title: The three-comb discrepancy method stops structurally; exact interval transfer gives a 7/3 spread cone, while the old 7/6 bank remains reconnaissance
status: CORRECTED. Part I is proved: the discrepancy/component-count method requires r<7/2. The original Part II claim that every one-period window leaves a full 6/(7k) safe interval is false; the sharp arbitrary-window guarantee is 3/(7k), with exact transfer Phi(x)=min(6/7,(x-1/7)/2) from THM-1137. Consequently the sound elementary recursion needs adjacent ratios at least 7/3, not 7/6. The exact 495-core atlas makes this a genuine uniform r=5 cone. The old 7/6 sample (300/300, minimum 4.949) and the clustered banks remain telemetry only. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.68; owner: continue and extend)
depends_on:
  - THM-1097    # codex's sharp three-comb theorem, whose stopping point (I) explains
  - THM-1094    # the two-comb theorem and the exact-bank method this would need
related: [MISTAKE-163, MISTAKE-164, MISTAKE-169, THM-1101, THM-1137]
script: 04-computation/four_comb_dichotomy_kps_S128c68.py, four_comb_clustered_kps_S128c68.py (+ .out)
---

# THM-1140 — the four-comb frontier

codex's THM-1097 closes the four-killer stratum uniformly with a *three*-comb component
theorem. Uniform r=5 needs a *four*-comb version. This records why the existing method
cannot simply be turned one more crank, and what the target constant actually is.

## (I) The method provably stops at three combs

With the sharp one-comb discrepancy |I ∩ D_k| ≤ |I|/7 + 6/(49k) and the component count
≤ (r+1) + |I|·Σk_i, the longest surviving component satisfies

> L ≥ [ (7−r)|I|/7 − (6/49)Σ 1/k_i ] / [ (r+1) + |I|·Σ k_i ].

Asking L > 1/(7 k_max) and using Σk_i ≤ r·k_max, the |I|·k_max terms compare as
(7−r)|I|k_max against r|I|k_max, so the inequality requires

> **7 − r > r,  i.e.  r < 3.5.**

Three combs clear it; four cannot. THM-1097 stops exactly where the arithmetic stops — the
halt is structural, and any four-comb proof must replace the averaging step, not sharpen it.

## (II) The exact gap recursion — correction

The former claim that an arbitrary interval of length `1/k` contains a full safe gap of
length `6/(7k)` is false. At scale `k=1`, the interval `[1/2,3/2]` cuts the safe arc into
two pieces of length exactly `3/7`; it contains no `6/7` piece. THM-1137 proves the sharp
transfer: if a closed interval has normalized length `x=kλ>=1`, it contains a safe closed
subinterval of normalized length at least

> **Phi(x)=min(6/7,(x-1/7)/2).**

In particular the arbitrary-window guarantee is `3/(7k)`. The sound coarse recursion is

> **L_j >= 3/(7k_j) provided L_(j-1) >= 1/k_j, hence
> k_j >= (7/3) k_(j-1).**

The exact 495-core atlas at the start of the referee gives a largest core-safe component
of length at least `1/70`. Every legal first killer has `k1>13 max(P)>=104`, so this
component is longer than `1/k1`. Therefore four adjacent ratios at least `7/3` give
`L>=3/(7k4)>1/(7k4)`, a genuine uniform `r=5` spread cone. The stronger normalized
transfer can improve this cone when more initial width is retained, exactly as in
THM-1137.

## (III) The old 7/6 bank is telemetry, not the spread theorem

The frozen experiment sampled 300 quadruples with adjacent ratios at least `7/6`; all
satisfy `7 k4 L>1`, with worst value **4.949** (core
`[2,3,5,7,8,10,11,12]`, killers `188/276/371/464`). These rows do not satisfy the
hypotheses of the corrected recursion in general. In fact `4.949<6` is an internal
countercheck to the old prediction, not a boundary-effect version of it. The sample is
useful reconnaissance for the intermediate ratio range `[7/6,7/3)`, but proves nothing
uniform there.

## (IV) The clustered rescue is weaker than hoped

I expected close combs to overlap heavily and so cost much less than two independent combs.
They do not:

| a | b | \|D_a\| | \|D_a ∪ D_b\| | vs 2/7 | ratio |
|---|---|---|---|---|---|
| 157 | 158 | 0.14286 | 0.26530 | 0.28571 | 0.929 |
| 300 | 330 | 0.14286 | 0.26364 | 0.28571 | 0.923 |
| 701 | 813 | 0.14286 | 0.26531 | 0.28571 | 0.929 |
| 157 | 471 = 3·157 | 0.14286 | 0.23810 | 0.28571 | **0.833** |

The saving is a flat ~7% across every ratio from 1.00 up to 7/6, and only reaches 17% when
b is an exact multiple of a (genuine comb nesting). And **73%** of real quadruples are
clustered at some step, so this is the majority case, not a corner.

## (V) But the clustered half holds too, measured

| regime | tested | min 7·k₄·L | median | # ≤ 1 |
|---|---|---|---|---|
| consecutive k,k+1,k+2,k+3 | 220 | 3.745 | 4.490 | **0** |
| step ≤ 3 | 220 | 2.629 | 4.002 | **0** |
| step ≤ 8 | 220 | 3.046 | 4.297 | **0** |
| step ≤ 25 | 220 | 3.254 | 4.803 | **0** |

Over a further 900 tight-clustered quadruples the worst value is **2.358**, at core
[1,3,5,6,7,8,11,12] with killers (371,374,377,379). The remainder is never empty (0/500).

## What this is and is not

It is **not** a proof. Per MISTAKE-163/164 — filed against my own earlier all-scale
inferences from sampled ratios — measured margins do not discharge a universal quantifier,
and I am not repeating that error here. **Uniform r=5 remains open.**

What it does supply, for whoever builds the four-comb bank:

- the exact reason the current method halts (I), so the replacement step is identified;
- a proved `7/3` adjacent-ratio cone from the corrected transfer (II);
- a measured, not proved, positive bank in the intermediate `7/6` ratio regime (III);
- the **target constant**: the four-comb bound appears to hold with worst-case margin
  ≈ 2.36, not marginally — so an exact bank has room to work in;
- the **hardest case to aim at**: tight-clustered quadruples with consecutive-ish killers
  near 371–379, where the margin is thinnest.

## Named next
- The four-comb theorem needs an exact endpoint bank plus an analytic tail, in the shape of
  THM-1094/1097. The corrected recursion owns the `7/3` spread cone; the intermediate
  `[7/6,7/3)` region and the clustered regime need additional geometry. The 7% union
  saving of (IV) is not enough by itself.
- A promising angle for the clustered regime: when k_{j+1}/k_j is close to 1, the two combs
  are a *beat* pattern with period ≈ 1/(k_{j+1} − k_j), and the surviving set inherits that
  long-period structure. That is where the extra length must come from, and it is invisible
  to the measure-only accounting in (IV).
