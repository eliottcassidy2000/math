---
id: THM-1140
title: WHY THE THREE-COMB METHOD STOPS AT THREE (proved), AND RECONNAISSANCE FOR THE FOUR-COMB THEOREM. (I) PROVED — an exact reason THM-1097 halts where it does: with the sharp one-comb bound |I∩D_k| ≤ |I|/7 + 6/(49k) and the component count ≤ (r+1) + |I|·Σk_i, the longest surviving component obeys L ≥ [(7−r)|I|/7 − (6/49)Σ1/k_i] / [(r+1) + |I|Σk_i], and demanding L > 1/(7k_max) forces **7 − r > r, i.e. r < 3.5**. Three combs pass, four do not — the halt is structural, not a gap in effort. (II) THE PROPOSED REPLACEMENT, a gap recursion: inside a component of length λ, removing D_k leaves a full gap of length 6/(7k) whenever λ ≥ 1/k, so L_j ≥ 6/(7k_j) provided **k_j ≥ (7/6)k_{j−1}** at every step. (III) MEASURED — the spread half (all ratios ≥ 7/6) satisfies the four-comb target 7·k₄·L > 1 on 300/300 quadruples with worst value **4.949**. (IV) BUT 73% of real quadruples are CLUSTERED at some step, and the hoped-for rescue there is weak: |D_a ∪ D_b| = 0.2653 against 2/7 = 0.2857 for close a,b, a saving of only **7%** — and that 0.929 ratio is flat across ratios 1.00 to 7/6, dropping to 0.833 only when b is an exact multiple of a. (V) MEASURED ANYWAY, the clustered half also holds: over ~1600 clustered quadruples the worst 7·k₄·L is **2.358** (core [1,3,5,6,7,8,11,12], killers 371/374/377/379), with ZERO cases below 1 in any regime and the remainder never empty (0/500). So the four-comb theorem looks TRUE with worst-case margin ≈ 2.36, and the thinnest margin sits at tight-clustered quadruples
status: (I) PROVED — elementary arithmetic on the sharp bound and the component count; it explains THM-1097's stopping point exactly. (II) PROVED as an implication (the full-period gap argument), but it only covers the spread regime. (III),(IV),(V) MEASURED on samples — **this is reconnaissance, NOT a proof of the four-comb theorem**, and per MISTAKE-163/164 sampled margins do not discharge an all-scale quantifier. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.68; owner: continue and extend)
depends_on:
  - THM-1097    # codex's sharp three-comb theorem, whose stopping point (I) explains
  - THM-1094    # the two-comb theorem and the exact-bank method this would need
related: [MISTAKE-163, MISTAKE-164, THM-1101]
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

## (II) The gap recursion

Inside a component of length λ, the teeth of D_k are spaced 1/k with width 1/(7k), so if
λ ≥ 1/k the component contains a full gap of length 6/(7k). Iterating,

> **L_j ≥ 6/(7 k_j)  provided  L_{j−1} ≥ 1/k_j, i.e. k_j ≥ (7/6)·k_{j−1}.**

Four spread steps therefore give L ≥ 6/(7k₄), a factor 6 above the 1/(7k₄) needed.

## (III) The spread half holds

300/300 spread quadruples satisfy 7·k₄·L > 1, worst value **4.949** (core
[2,3,5,7,8,10,11,12], killers 188/276/371/464). The recursion predicts 6; boundary effects
cost some, and 4.9 is what survives.

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
- a proved implication covering the spread regime (II), which is 27% of quadruples;
- the **target constant**: the four-comb bound appears to hold with worst-case margin
  ≈ 2.36, not marginally — so an exact bank has room to work in;
- the **hardest case to aim at**: tight-clustered quadruples with consecutive-ish killers
  near 371–379, where the margin is thinnest.

## Named next
- The four-comb theorem needs an exact endpoint bank plus an analytic tail, in the shape of
  THM-1094/1097. The tail should use the gap recursion (II) for the spread regime and a
  genuinely new argument for the clustered one — the 7% union saving of (IV) is not it.
- A promising angle for the clustered regime: when k_{j+1}/k_j is close to 1, the two combs
  are a *beat* pattern with period ≈ 1/(k_{j+1} − k_j), and the surviving set inherits that
  long-period structure. That is where the extra length must come from, and it is invisible
  to the measure-only accounting in (IV).
