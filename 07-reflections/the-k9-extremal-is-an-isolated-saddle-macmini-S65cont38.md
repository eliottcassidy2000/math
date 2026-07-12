---
source: mac-mini-2026-07-09-S65 (cont.38, 2026-07-11)
tags: [lrc14, extremal, saddle, variance, three-gap, endgame]
---
# The k=9 base extremal is an isolated saddle, not an endpoint

The last analytic gap of the LRC(14) wide-spread direction is one inequality:
inf over primitive 9-sets of J(E) = E[N_empty(7−N_empty)] equals 4465/882, attained at the
consecutive family {1..9} (THM-711/714/716). This session's exploration reveals WHY every
structural attack on it has failed — and what the right proof shape is.

**J = μ(7−μ) − Var.** The functional is a difference of two competing structural quantities:
the mean-parabola μ(7−μ) (maximized at the "balanced" μ = 7/2) and the emptiness variance Var
(a synchronization measure — how bunched the empty sectors are, THM-715). Minimizing J wants
SMALL μ and LARGE Var simultaneously.

**These two demands are antagonistic, and consec is the unique reconciliation.** Small μ (good
average coverage) requires the phases to spread evenly — which DE-synchronizes emptiness and
crushes Var. Large Var (bunched emptiness) is achieved by resonant families (mod-7-aligned) —
which cover badly and RAISE μ. Consec sits at the saddle: it accepts a slightly-above-minimal μ
(1.446, not the achievable 1.43) in exchange for the maximal Var compatible with that μ. The
low-μ corner and the high-Var pole are BOTH empty of competitors, verified by adversarial
hill-climbing from both directions.

**This explains the cont.29 failure of compression/monotone induction.** A saddle has no
monotone approach — you cannot walk to it downhill in any single coordinate (μ or Var or a
compression parameter). That is precisely why "compression toward packing raises Φ" failed at
k ≥ 5 (mod-7 resonance) and why "Var-max at consec" is globally false (THM-715's counterexample):
each is a projection onto ONE axis of a two-axis tradeoff. The correct proof must hold BOTH the
mean and the variance simultaneously — i.e. bound J directly, as the moment-majorant route
(THM-711) does, not either factor alone.

**The three-gap connection (this session's literature merge, klein-S253).** The consec family's
emptiness structure IS the three-gap theorem's orbit {0, x, 2x, …, 8x}: Steinhaus bunching
gives it the maximal Var-per-μ. Giri-Kravitz's accumulation hierarchy (2026) is the rigorous
external form of the decorrelation limits (THM-687/700/710) that dispatch the far/spread
families to the parabola vertex. The margin-free tight handling is FORCED — no citable gap above
1/14 exists (Fan-Sun disproved the naive spectrum). So the isolated-saddle character of the
extremal is not an accident of our coordinates; it is the analytic residue of a problem whose
extremal is genuinely a tradeoff optimum, provable only by a joint (two-moment) bound.

→ THM-711/714/715/716, THM-531 (dilation), Giri-Kravitz hierarchy, three-gap theorem, the
cont.29 compression refutation (the saddle's signature).
