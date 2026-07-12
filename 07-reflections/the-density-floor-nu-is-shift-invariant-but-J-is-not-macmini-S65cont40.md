---
source: mac-mini-2026-07-09-S65 (cont.40, 2026-07-11)
tags: [lrc14, framing, density-floor, sector-moment, shift-invariance, THM-717, endgame]
---
# The density floor ν is shift-invariant; the sector-moment J is not — a framing clarification

Working the k=9 compact reduction (my agreed lane vs klein's wide limit, THM-717), I hit a
framing subtlety worth pinning before it becomes a MISTAKE. Two objects are BOTH called "the
k=9 base check," and they are DIFFERENT:

**(ν) the density floor** ν(E) = μ(goodSet E) — THM-661's actual reach input. goodSet depends
only on the DIFFERENCE SET of E, so ν is **SHIFT-INVARIANT**: ν({0..8}) = ν({1..9}) = ν({2..10})
= 0.84014 exact. And it clears bar₉ = m_P + 1 − cap₁₀ = 114041/252252 = 0.4521 by **+0.388** —
the raw k=9 density floor is NOT tight (huge slack).

**(J) the sector-emptiness moment** J = E[N(7−N)], N = #empty sectors — the object of
THM-711/716/717. J is **SHIFT-DEPENDENT**: J({0..8}) = 5.199, J({1..9}) = 5.062, J({2..10}) =
5.392, all with the SAME difference set. Its minimizer label is therefore framing-sensitive.

**The physical consequence.** The co-offset set is e_i = Vmax − u_i with e = 0 for Vmax
(THM-527-A), so **0 ∈ E is forced** and the physically-consecutive cluster has co-offsets
{0..8}, NOT {1..9}. Hence the physical J-minimizer is {0..8} (J = 5.199), and klein's THM-717
"consec = {1..9}, J = 4465/882" is a co-offset set with no zero — a valid strict LOWER bound on
the physical J (4465/882 < 1019/196), so THM-717's bound J ≥ 4465/882 ≥ 432/91, if proved,
implies the physical J ≥ 5.199 ≥ 432/91 a fortiori. Every framing clears the threshold; nothing
breaks. But write-ups must say WHICH object and WHETHER 0 ∈ E.

**Why this matters for the endgame.** The "binding" k=9 base check is the sector-moment J
(thin +0.315 margin at its minimizer), NOT the raw density floor ν (fat +0.388). The two are
different tightnesses of the same covering phenomenon; J is klein's THM-717 covering-moment
reformulation, ν is the direct floor. The compact reduction (min J by diameter, this session)
is on J and is clean — min at the smallest diameter, rising after — so the J-extremal is
finite and consec-anchored either way.

→ THM-711/716/717 (the J object), THM-661 (the ν object), THM-527-A (0 ∈ E), the shift-vs-
dilation distinction (ν shift-invariant AND dilation-invariant; J only dilation-invariant).
