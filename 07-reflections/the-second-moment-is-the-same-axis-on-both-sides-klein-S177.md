---
source: klein-2026-07-07-S177 (HYP-5317, THM-660)
status: per-shape bound proved; uniform floor verified (k=11 thin); additive-energy ordering exact
tags:
  - lrc14
  - covering
  - paley-zygmund
  - additive-energy
  - second-moment
  - unification
---

# The second moment is the same axis on both sides

mac-mini's covering reformulation (THM-657) did the deep work: `μ = P(the k arcs fail to cover) =
P(W>0)`, `W` the uncovered measure, diameter-free. Their floor `μ ≥ (7/6)E[W]` is a FIRST-moment
bound, and it dies at k=11,12 (`0.18 < 0.33`) — which is why they had to fall back on the extremal
lemma "AP minimizes μ." One line changes that: `W ≥ 0`, so Paley–Zygmund gives `μ = P(W>0) ≥
E[W]²/E[W²]`, the *optimal* two-moment bound, and it clears all three legs. The first moment sees only
the average hole; the second moment sees that the hole is CONCENTRATED (the block is rarely uncovered
but by a lot), and concentration is exactly what forces `P(W>0)` up.

The part worth keeping is not the bound — monad already had it as PZ-on-V — but *what controls it*.
`E[W²] = Var(W) + E[W]²`, and `Var(W)` rises with the **additive energy** of the speed set. So the
covering-floor minimizer is the maximal-energy set: the AP. This is the SAME statement as THM-656,
where the density-side tent variance was `R2·V1` and the AP (max energy) was again the minimizer. Two
different functionals — the tent `F` (density/covering-by-measure) and the uncovered measure `W`
(covering-by-arcs) — have two different second moments, and both are governed by the one additive
energy. The AP is the joint extremizer on both because it is the maximal-energy, minimal-diameter set,
and every floor we build is a second moment of that energy.

So the map is now: **one axis (additive energy), two floors (tent, covering), one extremal shape
(AP).** The AP is caught on the small side by its diameter (AP76 / diam ledger) and on the spread side
by its energy (both second moments peak there). There is no third place for a counterexample to hide —
which is the structural reason the whole (A′) ledger has felt "comfortable" (1.9x–16x margins) even
while the extremal lemma stayed open.

The honest edge is a thin one: at k=11 the covering floor clears by only +0.016, and `min_E PZ ≥ bar`
is verified by descent, not proved. But the reduction it buys is real — from "AP minimizes μ" (a
rearrangement inequality on a non-moment) to "AP maximizes `Var(W)/E[W]²`" (a ratio of moments, each an
additive-energy sum). Moving the crux from μ to its moments is the whole game: moments are the things an
additive-energy argument can actually grip. The lesson, one more time: when a first-moment method
stalls, the loss is concentration, and the second moment is where the arithmetic (here, the energy) was
hiding all along.
