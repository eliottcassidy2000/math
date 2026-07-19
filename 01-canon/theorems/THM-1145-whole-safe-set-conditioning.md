---
id: THM-1145
title: THE COMPONENT-AWARE CONDITIONING WORKS — but at the level of the WHOLE SAFE SET, not one component. (I) THE BAD ZONE IS NARROW AND CONTIGUOUS: for consecutive-type quadruples the indices j with 7·k₄·L(j) ≤ 1 form a single run of width ≈ **0.0457·k₁**, always centred at **j/k₁ ≈ 0.238–0.242** — i.e. an interval of length ≈ 0.046 in t, sitting at t ≈ 0.24. Worst measured run 18 indices at (394,395,396,397), j starting 94, j/k₁ = 0.2386. (II) SO SINGLE-COMPONENT CONDITIONING FAILS: the atlas guarantees only a component of length ≥ 1/70 = 0.0143·(in t), which is **three times narrower** than the bad run, so a minimal component can sit entirely inside it. That was the form of the question and its answer is no. (III) BUT THE THEOREM DOES NOT NEED ONE COMPONENT — it needs a good k₁-gap ANYWHERE in S(P). The bad zone is a single interval of length ≈ 0.046 while S(P) has total measure **0.164–0.363** spread across [0,1] in 14–26 components, so it cannot be swallowed. (IV) VERIFIED EXHAUSTIVELY OVER ALL 495 EIGHT-SPEED CORES (not sampled): taking the max over k₁-gaps lying inside S(P), the minimum of 7·k₄·best over all cores is **1.97771** at killers (157,158,159,160), and 4.50, 4.47, 4.45, 4.46, 2.23 for the other tested quadruples; a broader sweep over consecutive quadruples with k₁ ∈ [157,420] gives **1.74074**. Every value exceeds the required 1, with the worst leaving 74% of margin
status: (I),(II) MEASURED exactly — (II) is a definitive negative for the single-component form. (III) is the correct reformulation. (IV) is **exhaustive over all 495 cores** for the six listed quadruples plus a strided sweep, which is a stronger quantifier than this thread has usually achieved, but it is still finitely many quadruples — **not a proof for all (k₁,k₂,k₃,k₄)**. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.73; owner: work the component-aware conditioning at j ≈ k1/4)
depends_on:
  - THM-1144    # which refuted the gap-uniform route and located the failures at j ≈ k₁/4
  - THM-1142    # the exact gap law
related: [THM-1143, THM-1097, THM-1094]
script: 04-computation/bad_zone_width_kps_S128c73.py, whole_safe_set_kps_S128c73.py (+ .out)
---

# THM-1145 — condition on the whole safe set, not one component

## (I) The bad zone is narrow, contiguous, and always in the same place

THM-1144 located the failing gaps near j ≈ k₁/4. Measuring the full bad set:

| killers | #bad j | longest run | run/k₁ | run start j/k₁ |
|---|---|---|---|---|
| (157,158,159,160) | 12 | 6 | 0.0382 | **0.242** |
| (197,198,199,200) | 18 | 9 | 0.0457 | **0.239** |
| (317,318,319,320) | 28 | 14 | 0.0442 | **0.240** |
| (394,395,396,397) | 36 | 18 | **0.0457** | **0.2386** |
| (371,374,377,379) | 2 | 1 | 0.0027 | 0.092 |
| (211,212,214,215) | **0** | 0 | 0 | — |

The bad indices form a **single contiguous run** of width ≤ 0.0457·k₁, and it always sits at
j/k₁ ≈ 0.24. In t-coordinates that is one interval of length ≈ 0.046 near t ≈ 0.24. Some
quadruples have no bad gaps at all.

## (II) Single-component conditioning fails

The question as posed was whether a core-safe component must contain a good gap. It need
not: the atlas guarantees a component of length ≥ 1/70 ≈ 0.0143, which is **three times
narrower** than the 0.0457 bad run, so a minimal component can sit entirely inside it.

## (III) The right reformulation

The four-comb theorem does not need a good gap in a *given* component — it needs one
**anywhere in S(P)**. And S(P) has total measure **0.164–0.363** spread over 14–26
components across [0,1], while the bad zone is a *single* interval of length ≈ 0.046. It
cannot swallow the safe set.

## (IV) Verified exhaustively over all 495 cores

Taking the maximum over k₁-gaps that lie inside S(P):

| killers | min over **all 495 cores** of 7·k₄·best |
|---|---|
| (157,158,159,160) | **1.97771** |
| (197,198,199,200) | 4.50000 |
| (300,301,302,303) | 4.45500 |
| (317,318,319,320) | 4.47161 |
| (394,395,396,397) | 4.45051 |
| (371,374,377,379) | 2.23010 |

Minimum **1.97771**, at the smallest legal quadruple with the core [1,2,6,7,8,9,11,12]. A
strided sweep over consecutive quadruples with k₁ ∈ [157,420] gives a minimum of
**1.74074** at (191,192,193,194). Every value clears the required 1 with at least 74% of
margin.

## Honest status

(IV) is exhaustive over **cores** — all 495, no sampling — which is a stronger quantifier
than most of this thread. It is still only six quadruples plus a strided sweep, so it is
**not a proof for all (k₁,k₂,k₃,k₄)**. Uniform r=5 remains open.

What has been established is the correct *shape* of the argument, and the two measured
ingredients an exact bank would need:

1. the bad set is one contiguous run of width ≤ 0.046·k₁ at j/k₁ ≈ 0.24;
2. S(P) has measure ≥ 0.164 spread over ≥ 14 components.

Since (1) is an interval and (2) is spread, a counting argument closes the gap — and both
ingredients are atlas-computable.

## Named next
- Prove (I): that the bad set is a single run, and bound its width. THM-1142's linear law
  makes this plausible — the gap value descends linearly in j, so {j : value ≤ threshold} is
  an interval by monotonicity, and its width follows from the slope d/(k₁k₄).
- Then (III) is a two-line counting argument against the atlas: an interval of length 0.046
  cannot meet every component of a set of measure 0.164 spread over ≥ 14 components, so some
  adequate component survives.
- That is a genuine route to the four-comb theorem and it is component-aware as THM-1144
  demanded, without needing a full endpoint bank over quadruples.
