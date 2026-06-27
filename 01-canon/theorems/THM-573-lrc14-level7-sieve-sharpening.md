---
id: THM-573
title: The level-7 lift sieve — every 13-speed set with >=7 multiples of 7 is lonely (sharpens THM-571)
status: PROVED modulo the accepted LRC(<=13) input (same standing as THM-571)
author: kind-pasteur-2026-06-27-S31af
depends_on:
  - LRC<=13     # accepted below-frontier theorem, used for the |P|<=12 phase margin
  - THM-571     # this generalizes & gives a single-argument proof of THM-571's domain
related:
  - THM-570     # 14-phase shift guard (subsumed by the single level-7 argument)
  - THM-568     # apex-shell divisibility (orthogonal: locates the optimum)
  - THM-523     # q=14 witness (the 14-free trivial half)
  - HYP-3084    # refuted margin conjecture + the dilation/aliasing correction
results:
  - 04-computation/lrc_level7_sieve_kps.py
  - 05-knowledge/results/lrc_level7_sieve_kps.out
---

# THM-573 — the level-7 lift sieve

## Statement
Let `S` be 13 distinct positive integers and let `H = {v in S : 7 | v}` be the speeds
divisible by 7. If `|H| >= 7`, then `M(S) > 1/14` (S is lonely at threshold `1/14`).

This **sharpens and subsumes THM-571** (which assumed `|M14| >= 7`, i.e. >= 7 multiples of
`14`). Since `14 | v => 7 | v`, we have `|H| >= |M14|`, so `{|M14| >= 7} ⊂ {|H| >= 7}`:
THM-573 covers a strictly larger family — every set with `7..12` multiples of `7`, even when
fewer than `7` of them are multiples of `14` — by a **single argument with no case split**.

## Proof
Write `H = 7P` with `P = {v/7 : v in H}`, `|P| = |H|`. `S` is taken primitive (dilation
invariance, `M(dS)=M(S)`), so not all 13 speeds are divisible by 7 and `|P| <= 12`.

**(1) A `P`-safe phase.** By the accepted LRC theorem below 14 runners (`|P| <= 12`), there is
a phase `v*` with `||p v*|| >= 1/(|P|+1) >= 1/13 > 1/14` for every `p in P`.

**(2) `H` is safe at all seven lifts.** Put `t_j = (v* + j)/7`, `j = 0..6`. For `h = 7p in H`,
```
||h t_j|| = ||7p (v*+j)/7|| = ||p(v*+j)|| = ||p v* + p j|| = ||p v*|| >= 1/13 > 1/14,
```
since `p j` is an integer. So every speed of `H` is safe at **all** seven lifts.

**(3) Each 7-coprime speed forbids at most one lift.** For `b in S` with `7 ∤ b`, the seven
values `b t_j = (b v* + b j)/7 (mod 1)`, `j=0..6`, are **seven equally spaced points** of
spacing `1/7` (because `gcd(b,7)=1` makes `{bj mod 7} = {0,..,6}`), with offset `b v*/7`. The
forbidden zone `||·|| < 1/14` is an open arc of length `2/14 = 1/7`. An open arc of length `1/7`
contains **at most one** of seven points spaced `1/7` (the next point is exactly `1/7` away, i.e.
at/after the far endpoint). So `b` forbids `<= 1` lift.

**(4) Pigeonhole.** The number of 7-coprime speeds is `13 - |H| <= 13 - 7 = 6 < 7`. They forbid
`<= 6` of the seven lifts, so some lift `t_j` is forbidden by none: every speed is safe there,
`M(S) >= 1/14`. A generic perturbation of `v*` inside the (open) `P`-safe region makes the finitely
many threshold equalities strict, giving `M(S) > 1/14`. ∎

## Verification
`lrc_level7_sieve_kps.py` (`05-knowledge/results/lrc_level7_sieve_kps.out`):
- **V1:** the level-7 certificate fired on **1500/1500** random sets with `|H| >= 7`; exact `M`
  cross-check found **0** counterexamples to `M > 1/14`.
- **V2:** over 20000 trials, a single 7-coprime speed forbids **at most 1** lift — the crux bound (3).
- **V3:** at `|H| = 6` the worst-case count is `13 - 6 = 7 = #lifts`, so the pigeonhole gives **no
  guarantee** — `|H| = 6` is exactly the boundary, and the residual begins there. (Individual
  `|H|=6` instances remain lonely, but by a *witness search*, not by this sieve as a proof.)

## What this does and does not do
- **Does:** closes every 13-set with `>= 7` multiples of 7, with one clean sieve. Folds THM-570 +
  THM-571 (Case 1 + Case 2) into a single step and **enlarges** the proved domain.
- **Relocates the residual:** LRC(14) now reduces to `[14-free: t=1/14, trivial]` ∪
  `[covering with <= 6 multiples of 7]`. The residual is **`<= 6` multiples of 7** — strictly
  smaller than THM-571's `|M14| <= 6`.
- **Does NOT** close LRC(14). The residual (`<= 6` multiples of 7, equivalently a bounded core with
  `>= 7` speeds coprime to 7) is the genuine hard part and contains the tight cases — e.g. the
  primitive `{1..12,182}` (`M = 14/183`) and the non-primitive `2·{1..13}` (`M = 1/14`); both have
  `|H| = 1`. See [[HYP-3084]] for the dilation/aliasing correction to the "covering has a margin" idea.

## The descent principle (why 7, and why it stops)
The apex prime of `14 = 2·7` is `7`. The lift sieve at level `d` proves `M > 1/14` once `>= 7`
speeds are divisible by `d` *and* the `d`-coprime residual forbids `<= 1` lift each — which needs
`d = 7` (so the forbidden arc `1/7` equals the lift spacing `1/d`). Level `2` is useless (needs
`<= 1` odd speed); there is no other guaranteed common prime. So the descent is **single-step**:
level 7 down to the `<= 6`-multiple-of-7 core, where no further arithmetic descent is available and
the bounded-core / equidistribution machinery (Node-3, S31v) must take over.

→ THM-571, THM-570, THM-568, THM-523, LRC<=13, HYP-3084, [[lrc14-thread]].
