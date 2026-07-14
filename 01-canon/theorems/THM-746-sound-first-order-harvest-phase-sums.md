---
id: THM-746
title: The sound first-order harvest -- F quadratic makes the dF_ext expansion EXACT in three terms: Phi(W) = (Xi_sv - Z(W)) - S(W)/W - T(W)/(2W^2), so C1 collapses to 2(#comp + |Xi_sv|) = 2.11 SOUNDLY (6.9x below THM-743's C1 at large W) at the price of a phase-sum C2 ballast; best bound = min(THM-743, THM-746) per W (W0 unchanged at 339/513); the phase sum S(W) = the arrangement's VERTICES AS RUNNERS observed at time W -- the perspective tower inverts, and |S(W)| observed 2-50 vs bound 516 is the next level's prize
status: PROVED + VERIFIED-EXACT (the three-term expansion is a Fraction identity at every tested W, both shapes; wedge equality re-confirms the pairing theorem; assembly zero violations; honest W0 trade-off stated)
source: opus-2026-07-14-S281 (owner prompt: bound the signed dF_ext sum with the perspective frame)
depends_on:
  - THM-745 (Phi(W) IS the first-order content -- the pairing theorem, unconditional)
  - THM-742/743 (the assembly and the standing mid-range constants)
related:
  - MISTAKE-142 (this is its sound correction: each order is an exact identity, not a crude charge)
  - the Mode-A / renormalization threads (the tower: each harvest converts a first-order constant into a symmetry-exact constant + a phase ballast one order down)
---

# THM-746 -- the sound first-order harvest (phase sums)

## The exact expansion

F(h) = h(1-h)/2 is quadratic, so F(x+d) = F(x) + psi(x) d - d^2/2 EXACTLY.  With
th1 = ceil(u1 W) - u1 W and th2 = u2 W - floor(u2 W) (the endpoints' GRID PHASES, in [0,1)),
the first-order wedge content (THM-745: exactly Phi(W) = Sum_seg orient(-dF_ext)) satisfies,
as a Fraction identity (verified at W in {90, 97, 250, 800}, both shapes):

>  Phi(W) = (Xi_sv - Z(W)) - S(W)/W - T(W)/(2W^2),
>    Xi_sv = Sum_seg orient [F(x_start) - F(x_end)]     (W-independent; +0.056122 / +0.178571)
>    S(W)  = Sum_seg orient j [psi(x_start) th1 + psi(x_end) th2]
>    T(W)  = Sum_seg orient j^2 [th1^2 - th2^2]
>    Z(W)  = the zero-crossing segments' vertex terms (du < 2/W each).

## The sound bounds and the honest trade

|Phi(W)| <= |Xi_sv| + [S1 + Z1-part]/W + S2/(2W^2) with S1 = Sum_seg j(|psi(x_s)|+|psi(x_e)|)
= 516.35 / 357.70, Z1 = Sum_seg j = 1078 / 792, S2 = Sum_seg j^2 = 9018 / 5272 (all exact).
Assembly: C1 = 2(#comp + |Xi_sv|) = **2.1122 / 8.3571**; C2 = C2(743) + 2(S1+Z1) + S2 =
15897 / 9542.  Zero violations over the W-spread.

| shape | C1: 743 -> 746 | large-W gain | C2: 743 -> 746 | W0: 743 vs 746 |
|---|---|---|---|---|
| 1 | 14.49 -> **2.11** | **6.9x** | 3690 -> 15897 | **339** vs 475 |
| 2 | 19.14 -> **8.36** | 2.3x | 1971 -> 9542 | **513** vs 564 |

(The .out's per-run label printed 19.14 for both shapes; shape 1's THM-743 C1 is 14.49 -- the
factors here are the correct ones.)  **The best sound bound is min(THM-743, THM-746) per W**:
743 below W ~ 600, 746 asymptotically; W0 unchanged at 339/513.  This is the SOUND version of
what MISTAKE-142's unsound charge attempted: each order is now an exact identity -- the
first-order 50x prize IS harvested (C1 = 2.11), and its true cost (the phase ballast) is
explicit instead of hidden.

## The vertices-as-runners inversion (the perspective punchline)

S(W) = Sum_e c_e {u_e W} + const-part: the segment endpoints u_e -- the arrangement's vertices,
rationals on the difference-runner grids (k +- 1/7)/delta -- act as RUNNERS with speeds u_e,
observed at integer TIME W, with exact coupling constants c_e = orient j psi(x_e).  Bounding
the LRC error spawned a lonely-runner-type system one level down whose runners are the pair
events of the level above.  S(W) is periodic in W modulo the lcm of the vertex denominators.
Measured: |S(W)| = 2.6 - 49 vs the absolute bound S1 = 516 -- ANOTHER 10-100x of signed
cancellation lives at this level.  The tower is now explicit: each harvest converts the current
first-order constant into a symmetry-exact tiny constant (mirror pairing) plus a phase-sum
ballast one order down; the perspective frame (origin band, pair grids, mirror = time reversal)
repeats at every level.  Diminishing W0 returns at the current level; the structure, not the
constant, is the deliverable.

## Files

04-computation/lrc14_phase_sum_harvest_thm746_opus_S281.py (+.out): exact identity checks,
wedge equality, constants, assembly, W0.
