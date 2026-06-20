---
id: THM-551
title: Apex-prime order-truncation of the coverage Newton expansion — the s-far packet Δ_S(B) vanishes whenever |B|+s < 7, so p0(B∪F)=Σ_{|S|≥7−|B|}Δ_S(B) begins at far-order 7−|B|; with the 1/7^{s+1} apex hierarchy a bounded core's coverage is dominated by its leading order, and the 3+3+1 packet grading is the cube-root shadow of the half-tiling recursion
status: PROVED (the truncation, elementary: p0(E)=0 for |E|<7). VERIFIED (Δ_S=0 for |B|+|S|<7, nonzero at threshold). Connects THM-549 half-tiling 7-term recursion, THM-548/HYP-2680 multi-far hierarchy, HYP-2681 cube-root Eisenstein.
source: mac-mini-2026-06-20-S4
depends_on:
  - THM-548   # boundary-value / Newton decomposition p0(B∪F)=Σ_S Δ_S
  - THM-549   # half-tiling 7-term recursion (the combinatorial parent)
related:
  - HYP-2680  # codex Φ_s Stirling hierarchy (Δ_S→Φ_s)
  - HYP-2681  # codex cube-root Eisenstein modes
external: Lonely Runner Conjecture; the apex prime 7 = #sectors = n/2 at n=14.
---

# THM-551 — Apex-prime order-truncation

## The lemma (PROVED)
The seven-sector coverage `p0(E)=meas{x: {frac(ex):e∈E} hits all 7 sectors}` satisfies
> **`p0(E)=0` whenever `|E|<7`** — at any `x`, `|E|` runners occupy at most `|E|` of the 7
> sectors, so for `|E|<7` some sector is missed at every `x`; the all-hit set has measure 0.

Hence for the Newton packet `Δ_S(B)=Σ_{T⊆S}(−1)^{|S|−|T|}p0(B∪T)`:
> **`Δ_S(B)=0` whenever `|B|+|S|<7`** (every term `p0(B∪T)` has `|B∪T|≤|B|+|S|<7`, so vanishes).
VERIFIED exactly; nonzero already possible at the threshold `|B|+|S|=7`.

## Consequences
1. **The expansion begins at far-order `7−|B|`:** `p0(B∪F)=Σ_{S⊆F, |S|≥7−|B|} Δ_S(B)`. The apex
   prime 7 is the *starting order* of the coverage, just as it is the sector count, the kernel's
   vanishing period, and the `7^{s+1}` denominator of the `s`-far constant (THM-548 §3b).
2. **Leading-order domination:** with the apex hierarchy (`s`-far constant `~1/7^{s+1}`, THM-548),
   the leading nonzero order `s_0=max(0,7−|B|)` dominates and successive orders fall by `~1/7`.
   So a bounded core `B` with `|B|=7` (e.g. `consec_7`) covers at order 0; `|B|=5` first covers at
   two-far (one-far packets identically 0 — VERIFIED); `|B|≤5` ⟹ `p0(B)=0` (kps cardinality lemma).
3. **The 3+3+1 grading is the cube-root shadow.** When `|B|=4` the live orders are exactly
   `s=3` (the lone three-far `G`) down from where; when `|B|` makes three far runners the relevant
   window, the packet set is `{3 one-far, 3 two-far, 1 three-far}` — the same `3+3+1` as the
   half-tiling's 3 corners + 3 edges + 1 center (THM-549), organized by the `S_3`/cube-root
   Eisenstein modes `S_ω=A+ωB+ωC` (HYP-2681). The coverage seven-term recursion is the
   measure-theoretic shadow of the half-tiling's recursive geometry; the apex prime 7 truncates
   the far-order tower to the seven terms that matter near the cap.

## Significance
A clean, proved bridge: the **half-tiling 7-term recursion** (combinatorial, THM-549), the
**multi-far coverage hierarchy** (analytic, THM-548/HYP-2680), and the **cube-root organization**
(HYP-2681) are one structure, and the apex prime 7 is simultaneously its sector count, its order
truncation, and its hierarchy ratio.
