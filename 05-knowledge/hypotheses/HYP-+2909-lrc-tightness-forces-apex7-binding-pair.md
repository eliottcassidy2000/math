---
id: HYP-2909
title: (star) FORWARD direction PROVED -- M(S)=1/14 forces an apex-7 BINDING PAIR (14 | s_i+s_j); tightness sits at the apex-7 antipodal point
status: PROOF FRAGMENT (rigorous local-max argument + verified); the forward half of crux (star)
source: mac-mini-2026-06-22-S49 (owner: finish via the THM-079 template; prove tightness forces the apex-7 point)
related:
  - HYP-2906   # codex one-large-speed peeler (Node 3 / Move A)
  - HYP-2908   # codex forbidden-H7 state-lift guardrail
  - HYP-2900   # the dichotomy / bounded core
  - THM-029    # H=7 forbidden (the apex)
---

# HYP-2909: M(S)=1/14 FORCES the apex-7 binding pair -- the forward half of (star)

## Statement (rigorous)
Let M(S) = max_t min_{s in S} ||s t|| and suppose M(S) = 1/14, attained at a local maximum t*.
Since 1/14 < 1/2, the minimum at t* is a CROSSING of two runners, not a single sawtooth peak, so there
is an active runner on the INCREASING side, s_i with frac(s_i t*) = 1/14, and an active runner on the
DECREASING side, s_j with frac(s_j t*) = 13/14 = -1/14. (A local max of a min-of-sawtooths needs both:
an increasing active runner to cap f as t decreases, and a decreasing one as t increases. If only one
side were active, f would not be at a local max.)

Then (s_i + s_j) t* = 1/14 + (-1/14) = 0 (mod 1), so t* = m/(s_i+s_j) for an integer m; and substituting
into s_i t* = 1/14 (mod 1) gives 14 (s_i+s_j) | (14 s_i m - (s_i+s_j)), whence 14 | (s_i + s_j).

> THEOREM (star-forward): M(S)=1/14  ==>  S contains a BINDING PAIR  s_i + s_j ≡ 0 (mod 14),
> i.e. two runners in ANTIPODAL residues mod 14 (the apex-7 antipode), with the optimum at
> t* = m/(s_i+s_j) placing them at ±1/14.

## Verified (lrc tightness scan, S49)
- AP {1..13}: M=1/14 at t*=1/14, binding pair {1,13} (1+13=14). Other denom-14 optima: {5,9} at
  t*=3/14, {3,11} at t*=5/14 -- ALL sum to 14.
- GW {1..11,13,24}: M=1/14 at t*=1/14, binding pair {1,13}. (Non-tight control {1..11,13,38}: M=1/12,
  NO binding pair -- consistent with the theorem.)

## Why this is the apex-7, and the (star) connection
s_i + s_j ≡ 0 mod 14 is exactly the speed-space ANTIPODE on the 14 = 2*7 grid -- the same apex-7 as the
7 diameter-ties (S48), as H=7 = I(K_3,2) forbidden (THM-029/200), and as codex's H7 state-lift
(HYP-2908). This PROVES the forward direction of (star): M=1/14 ==> the optimum is the apex-7 antipodal
(binding-pair) point. The binding residues are the 7 antipodal pairs {1,13},{2,12},{3,11},{4,10},{5,9},
{6,8},{7,7}.

## What remains of (star) (honest)
- (ii) s_i + s_j = 14 EXACTLY (not 28, 42, ...): the L=1 / denom-exactly-14 refinement.
- (iii) COVERING ==> M != 1/14: a multiple-of-14 runner sits on the observer at t*=a/14 (apex-7 floor),
  so covering is excluded from the denom-14 optima -- needs (ii).
- (iv) the non-covering tight locus = {AP, GW} (kps/codex census).
- (v) covering ==> M > 1/14 (not <): Node 2 bounded min = 1/12 (verified) + Node 3 equidistribution.
Given the forward theorem here + (ii)-(v), (star) closes and LRC(14) follows in the THM-079 template:
Move A = R1/R2/R3 peel to the single bounded atom; Move B = this apex-7 binding-pair forcing = the
Moon/cycle-forcing step.
