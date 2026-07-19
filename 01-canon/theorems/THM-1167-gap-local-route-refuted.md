---
id: THM-1167
title: THE GAP-LOCAL ROUTE TO THE FOUR-COMB THEOREM IS REFUTED — the minimising configuration is (k, k+1, k+2, k+3) at gap index ≈ k/4, and there the longest surviving piece falls BELOW the threshold. (I) THE CONSERVATIVE REDUCTION, precisely stated: define W(k₁,k₂,k₃,k₄) = min over ALL k₁-gaps in [0,1] of the longest piece surviving the k₂,k₃,k₄ teeth. If 7·k₄·W > 1 then the four-comb theorem follows for that quadruple *regardless of which gap the core leaves available*, since the 495-core atlas gives component length ≥ 1/70 and every legal k₁ exceeds 13·max(P) ≥ 104, so a full k₁-gap always fits. This would have decoupled the core entirely. (II) IT FAILS. Exhaustive over consecutive-type quadruples with k₁ ∈ [157,340]: **7·k₄·W = 0.76980** at (317,318,319,320), and the value is monotonically DECREASING in k₁ — 0.79013 (157), 0.78193 (197), 0.77651 (237), 0.77267 (277), 0.76980 (317). Individual rows: (300,301,302,303) gives 0.77679, (157,158,159,160) gives 0.79013, (371,374,377,379) gives 0.99015. So the answer to "does the minimising configuration exceed 1.295?" is **NO — it does not even exceed 1**. (III) THE MECHANISM, via THM-1142's law: the worst gap sits at index j ≈ k₁/4 (measured 39 for k₁=157, 49 for 200, 75 for 300), where the raw gap (a − j·d)/(a·b) with d = k₄ − k₁ = 3 evaluates to ≈ 1/(4k₄); subtracting the three tooth widths leaves ≈ 0.77/(7k₄). The linear descent of THM-1142 is exactly what drives it below threshold. (IV) THE CONSEQUENCE, which is the point: **the four-comb theorem cannot be proved gap-locally.** Some k₁-gaps genuinely fail, so any proof must use WHICH gaps the core-safe component makes available. The core cannot be decoupled, and a four-comb bank must track component location rather than quantifying over all gaps
status: (I) PROVED as a valid sufficient condition. (II) REFUTED by exhaustive exact-rational computation over the stated family — the failing configurations are witnessed, so the gap-local route is definitively closed, not merely unverified. (III) explains (II) via THM-1142. (IV) follows. **Uniform r=5 remains OPEN**, and this session narrows how it can be attacked rather than advancing it
source: kind-pasteur-2026-07-18-S128 (cont.72; owner: find the minimising configuration and confirm it exceeds 1.295)
depends_on:
  - THM-1143    # the three-tooth target this tests
  - THM-1142    # the exact gap law that explains the failure
related: [THM-1137, THM-1097, MISTAKE-169]
script: 04-computation/min_config_kps_S128c72.py (+ .out)
---

# THM-1167 — the gap-local route is refuted

THM-1143 reduced the four-comb theorem to a three-tooth spacing statement inside one
k₁-gap, with a measured margin of 3.05 against a required 1.295. The obvious next step was
to find the minimising configuration and confirm the margin survives it. It does not.

## (I) The reduction that would have sufficed

> **W(k₁,k₂,k₃,k₄) := min over all k₁-gaps in [0,1] of the longest piece surviving the
> k₂,k₃,k₄ teeth.**

If 7·k₄·W > 1, the four-comb theorem holds for that quadruple **whatever gap the core
leaves**, because the 495-core atlas gives a component of length ≥ 1/70 and every legal
k₁ > 13·max(P) ≥ 104, so at least one full k₁-gap fits inside. That would have removed the
core from the problem entirely — which is why it was worth testing first.

## (II) It fails, and not narrowly

Exhaustive in exact rationals over consecutive-type quadruples, k₁ ∈ [157,340]:

| k₁ | minimising quadruple | 7·k₄·W |
|---|---|---|
| 157 | (157,158,159,160) | 0.79013 |
| 197 | (197,198,199,200) | 0.78193 |
| 237 | (237,238,239,240) | 0.77651 |
| 277 | (277,278,279,280) | 0.77267 |
| **317** | **(317,318,319,320)** | **0.76980** |

Monotonically decreasing in k₁. Individual rows confirm the picture: (300,301,302,303)
gives 0.77679 and (157,158,159,160) gives 0.79013, while looser quadruples such as
(157,170,183,196) reach 1.26677 and (157,158,159,161) reach 1.37107.

So the answer to the question asked is **no**: the minimising configuration does not exceed
1.295, and it does not even exceed 1.

## (III) Why — THM-1142's law drives it under

The worst gap sits at index **j ≈ k₁/4** (measured: 39 for k₁ = 157, 49 for 200, 75 for
300). By THM-1142 the raw gap from tooth j of a to tooth j+1 of b is (a − j·d)/(a·b); with
consecutive killers d = k₄ − k₁ = 3 and j = k₁/4 this is

> (k₁ − 3k₁/4)/(k₁k₄) = 1/(4k₄),

and subtracting the three tooth widths (≈ 3/(7k₄)) leaves ≈ 0.77/(7k₄). The linear descent
that THM-1142 identified as the *source* of the useful nonuniformity is also what pushes
the worst gap below threshold. Both facts come from the same law.

## (IV) What this means

**The four-comb theorem cannot be proved gap-locally.** Some k₁-gaps genuinely fail — a
proof cannot quantify over all of them and must use *which* gaps the core-safe component
actually makes available. The core does not decouple.

This closes a route rather than opening one, but it closes it definitively: the failing
configurations are exhibited in exact arithmetic, not inferred from samples. Given how much
of this thread has been sampled claims overturned later, a witnessed negative is worth more
than another positive census.

## Honest status

**Uniform r=5 remains open.** This session narrows how it can be attacked: any four-comb
bank must be component-aware. The three-tooth statement of THM-1143 is still the right local
object, but it must be conditioned on the gap's position within the component, not proved
uniformly over gaps.

## Named next
- Re-pose the three-tooth statement *conditionally*: for gaps that a core-safe component can
  actually contain, is the bound restored? The failing indices cluster at j ≈ k₁/4, so the
  question is whether a component of length ≥ 1/70 must contain a gap away from that index.
- That is a statement about where core-safe components sit relative to killer teeth — which
  is exactly the coupling THM-1094 and THM-1097 handle with their exact endpoint banks. The
  four-comb version likely needs the same machinery rather than an elementary shortcut.
- Do not retry any gap-uniform formulation; (II) rules out the whole family.
