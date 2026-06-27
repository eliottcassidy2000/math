---
id: THM-574
title: The c-lift sieve family — |H_c| >= 14-c (c<=7) is lonely; c=7 is uniquely optimal
status: PROVED modulo the accepted LRC(<=13) input (same standing as THM-573); VERIFIED c=5,6,7
author: kind-pasteur-2026-06-27-S31ag
depends_on:
  - LRC<=13     # accepted below-frontier theorem (the paper arXiv:2604.23906), phase margin
  - THM-573     # the c=7 instance; this is its one-parameter generalization
related:
  - HYP-3087    # the polynomial-method/CRT bridge (this is the sieve-side of the apex-7 wall)
  - THM-571
results:
  - 04-computation/lrc14_clift_sieve_family_kps.py
  - 05-knowledge/results/lrc14_clift_sieve_family_kps.out
---

# THM-574 — the c-lift sieve family (and why c=7 is forced)

## Statement
Let `S` be a primitive 13-set (gcd 1) and `H_c = {v in S : c | v}`. **If `c <= 7` and
`|H_c| >= 14 - c`, then `M(S) > 1/14`.**

The threshold `14 - c` is minimized at `c = 7` (threshold `7`), recovering THM-573. THM-573 is
the `c=7` instance; THM-574 places it in a one-parameter family and explains its optimality.

## Proof (same shape as THM-573)
`H_c = c·P` with `|P| = |H_c| <= 12` (primitive ⟹ not all 13 divisible by `c`). By LRC(<=13)
there is a `P`-safe phase `v*` with `‖p v*‖ >= 1/13` for all `p in P`. Take the `c` lift points
`t_j = (v* + j)/c`, `j = 0,...,c-1`. For each `h = c p in H_c`, `‖h t_j‖ = ‖p v*‖ >= 1/13 > 1/14`,
so all of `H_c` stays safe at every lift. A speed `w` coprime to `c` meets the `c` lift points at
`c` equally-spaced phases (spacing `1/c`); its forbidden arc `{‖w t‖ < 1/14}` has length
`2/14 = 1/7`, so it catches `⌈(1/7)/(1/c)⌉ = ⌈c/7⌉` lifts. For `c <= 7`, `⌈c/7⌉ = 1`, so the
`13 - |H_c|` coprime speeds forbid at most `13 - |H_c|` of the `c` lifts. If `13 - |H_c| < c`,
i.e. `|H_c| >= 14 - c`, some lift survives and is a `1/14`-witness. ∎

## Why c=7 is uniquely optimal (the apex prime, sieve side)
`c = 7` is the **largest** `c` with `⌈c/7⌉ = 1` — the largest `c` whose lift-spacing `1/c` is
still `>=` the forbidden-arc width `1/7`. It is the unique fixed point **lift-spacing = arc-width**
(`1/7 = 1/7`). For `c > 7` the spacing is finer than the arc, one coprime speed can kill `>= 2`
lifts, and the survivor count argument fails; for `c < 7` the threshold `14 - c > 7` demands more
multiples. So `7` is forced as the optimal single-prime lift, threshold `7` — the same apex prime
that the polynomial method (arXiv:2604.23906 Prop 4.1) needs, where `φ(14) = 6` units fail to fill
the `13` indices (HYP-3087). Two mechanisms, one prime.

## Status / verification
PROVED (the argument is rigorous given LRC(<=13)). VERIFIED computationally for the firing cases
`c = 5,6,7` (`lrc14_clift_sieve_family_kps.py`: 0 failures over ~1500 trials each, `minM ≈ 0.10 ≫
1/14`). The `c = 7` case is THM-573 (independently verified 1500/1500). The residual after THM-573
remains `<= 6` multiples of 7 — no single `c`-lift closes it (the `c=2` threshold `12` is far), which
is exactly why the CRT-combined `c=2,7` lift + the analytic measure floor (CRUX 1) are needed.

→ THM-573, HYP-3087, the-polynomial-method-mod-14-why-7-is-forced.md, arXiv:2604.23906.
