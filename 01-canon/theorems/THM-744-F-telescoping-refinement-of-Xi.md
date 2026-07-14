---
id: THM-744
title: The F-telescoping refinement of Xi (THM-742's partial-wrap caps) via exact segment-end heights -- per-segment cost min(j du/2, 1/4, |F(x_s)-F(x_e)|), F(x)=x(1-x)/2; Xi 6.25->2.85 / 5.57->2.49; NET W0 ~unchanged (339->336 / 513->462): the per-line program's floor, reached; the STRUCTURAL yield is the loop-telescoping (signed total Xi_signed = 0.056 / 0.179, ~50x below Xi) localizing the cross-line cancellation to per-segment discretization residuals
status: CORRECTED (MISTAKE-142, opus-S278): the C2-inflation bookkeeping was UNSOUND (first-order charge Sum_seg 2j/W placed at second order); the headline W0 = 336/462 is NOT established -- THM-743's W0 = 339/513 remain the standing sound bounds. The STRUCTURAL content stands: the F-telescoping identity, exact end heights, loop chaining, and the Xi_signed diagnostic (0.056/0.179); superseded analytically by THM-745's exact identities.
source: opus-2026-07-13-S277 (owner prompt: attack the Xi partial-wrap caps with the exact segment-end heights)
depends_on:
  - THM-742/743 (the bound being refined; all other terms per THM-743)
related:
  - klein-S294/S295 (HYP-6570/6580) -- the pair-event grids and LRC(13)-localization; the residual sums below are the covering side's concrete hand-off to the Q_s class
---

# THM-744 -- the F-telescoping refinement of Xi

## Statement

With F(x) = x(1-x)/2 (an antiderivative of psi({h}) = 1/2-{h} along the height march,
CONTINUOUS across wraps since F(0) = F(1) = 0), each maximal exposed segment of a slope -j
boundary line satisfies

>  (j/W) Sum_crossings psi(s0)  =  orient * [F(x_start) - F(x_end)]  +  rho_seg,
>  |rho_seg| <= 2j/W   (head/tail rounding + per-wrap Raabe),

with x_start, x_end the segment's EXACT end heights (rationals from the arrangement: the
meeting heights of the terminating pair events).  Hence the per-segment first-order cost in
Xi may take  min( j du/2 , 1/4 , |F(x_s) - F(x_e)| )  at the price of +2j to C2 for each
segment using the third form.  All other THM-742/743 terms unchanged.

## Verified numbers

| shape | Xi old -> new | C1 old -> new | C2 (743) -> inflated | W0 (743) -> new |
|---|---|---|---|---|
| 1 | 6.246 -> **2.853** | 14.49 -> 7.71 | 3690 -> 5846 | 339 -> **336** |
| 2 | 5.569 -> **2.490** | 19.14 -> 12.98 | 1971 -> 3555 | 513 -> **462** |

Zero violations (W in {10..800}, exact).  **Honest verdict: the net W0 gain is marginal** --
the C2 inflation from the discretization residuals eats the C1 gain at the current W-scale.
The per-line constant program (S275 -> S276 -> S277: crude 1948 -> 452 -> 339 -> 336) has
reached its floor with the three terms balanced, as predicted in S276.

## The structural yield: loop telescoping and the measured cancellation

Along each closed boundary loop of R the segments CHAIN: at a same-orientation vertex
(r = 0 swap events) the F(x*) terms of the joined segments CANCEL exactly; only run
birth/death vertices (the r = +-1/7 pair events) survive, weighted +-2F(x*).  The SIGNED
telescoped total is therefore tiny:

>  **Xi_signed = |Sum_seg orient (F(x_s) - F(x_e))| = 0.0561 (shape 1), 0.1786 (shape 2)**

-- ~50x below Xi_new.  If the signed form were usable, C1 ~ 2(#comp + 0.06) ~ 2.1 and the
first-order floor would drop to W ~ 28.  What blocks it is now LOCALIZED AND MEASURED: the
per-segment residuals rho_seg (short Dedekind-type sums, exact and enumerable from the
arrangement) do not telescope; bounding their SIGNED total across segments is a joint
(cross-line) equidistribution estimate -- precisely the Q_s class (THM-729).  This is the
covering side's concrete hand-off object to the density machinery: finitely many explicit
short sums, known slopes, known end heights, measured target (~0.06 vs the naive 2.85).

## Files

04-computation/lrc14_F_telescoping_xi_thm744_opus_S277.py (+.out): exact end-height segment
walk, the three-way min with inflation bookkeeping, Xi_signed diagnostic, zero-violation
re-verification, new W0 solves.
