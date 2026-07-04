---
id: HYP-4075
title: The covering-min is a 2-term Ostrowski (continued-fraction) ladder — M_k=[0;n-1,k]=k/(k(n-1)+1) with the AP at rung k=1 (1/n) and the deep well at rung k=n (n/Phi6); the three-gap theorem is the rigidity, the killer's unit gap IS the +1 in Phi6, and Zeckendorf=Ostrowski for the golden CF is the special case
status: SYNTHESIS (exact + verified) — unifies the covering-min, the three-gap core (g<=3), and the repo's Fibonacci/Ostrowski/Zeckendorf/Sturmian threads. Clarifies the open core (both pieces reduce to 'prove the {k alpha} structure'), does NOT close it.
source: mac-mini-2026-07-04-S38 (owner: work the core + mine Fibonacci relations)
related:
  - HYP-3739  # covering-min M-uniqueness = Zeckendorf/Ostrowski canonical (mac-mini+klein); the seed of this
  - HYP-4070  # GAP-A (non-covering tight={AP,GW}); g<=3 = 'tight => {k alpha}' = this
  - HYP-2913  # three-gap/Steinhaus for the tight locus (now identified as the CLASSICAL theorem via {k alpha})
  - HYP-1902  # codex: LRC Zeckendorf boundary normal form
  - THM-536   # Sturmian reframe; THM-486 Pisano period (Fibonacci mod n)
results:
  - 04-computation/covering_min_continued_fraction_macmini_20260704.py
  - 05-knowledge/results/covering_min_continued_fraction_macmini_20260704.out
  - 07-reflections/the-covering-min-is-an-ostrowski-ladder-and-the-ap-and-deep-well-are-its-ends.md
external: three-gap/Steinhaus theorem; Ostrowski representation; Zeckendorf (Fibonacci); LRC(14).
---

# HYP-4075 — the covering-min is an Ostrowski ladder; three-gap is the rigidity

## The ladder (verified)
`M_k = [0; n-1, k] = k/(k(n-1)+1)`. **AP = rung k=1** (`[0;13,1]=1/14`, LRC bound, non-covering); **deep well
= rung k=n** (`[0;13,14]=14/183=n/Phi6`, covering-min). The whole margin `[1/14, 14/183]` is this ladder.
`cf(14/183)=[13,14]` confirmed. Ostrowski = CF-generalization of Zeckendorf (= Ostrowski for golden `[1;1,..]`).

## Three-gap is the rigidity (verified)
At `t*=n/Phi6` the deep well's phases are `{n*k mod Phi6}` (the multiples of `alpha=n/Phi6`, a `{k alpha}` set)
plus the killer image, with **2 gaps `{1/183,14/183}`** — the CLASSICAL Steinhaus three-gap. A generic covering
family has `g=5` (not `{k alpha}`). So `g(14)<=3` `<=>` extremal config is a `{k alpha}` progression `=` "tight
=> AP-like" `=` the finiteness. Steinhaus gives the gap count FOR FREE once `{k alpha}` is known; the open part
is proving the structure, never counting gaps.

## The mechanism (verified): covering => killer => unit gap = the +1 in Phi6
Covering forces a killer runner `(n-1)n` (=182) for the moduli `n-1,n`. At `a=n`: core `{1..n-2}` -> AP
`{n*k mod Phi6}` (difference `n`); killer -> `(n-2)n+1=169`, one unit above the top core point `168`, splitting
the wrap gap into `{1, 2n}`. The solitary UNIT gap is exactly the `+1` in `Phi6=(n-1)n+1`. No killer (AP) =>
`D=n, M=1/n`; killer (covering) => `D=Phi6, M=n/Phi6`. `Phi6`'s `+1` = a three-gap defect from the covering
constraint.

## Net
Unifies covering-min + three-gap + all Fibonacci/Ostrowski threads into one 2-term-Ostrowski-ladder picture,
and localizes the open core: BOTH pieces (`covering => M>=M_n` and `g<=3`) reduce to "the extremal config is a
`{k alpha}` progression" — klein's M-minimization/budget on the Ostrowski rungs. Not a proof; a sharp map.
Golden ratio is the wrong constant (ladder is `[0;13,k]`, gap-ratio 14), but the framework is the Fibonacci one
at a different CF. -> HYP-3739, HYP-4070, HYP-2913, klein Ostrowski binding. Script: covering_min_continued_fraction.
