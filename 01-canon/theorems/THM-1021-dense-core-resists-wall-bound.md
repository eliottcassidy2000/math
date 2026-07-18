---
id: THM-1021
title: THE DENSE CORE IS NOT CLOSED BY THE SHARP WALL BOUND — an honest negative with its exact mechanism: THM-1012's per-pair floor 4λ² − 2λ(a/b) has its defect → 2λ as a → b, and 2λ = 1/7 swamps the leading 4λ² = 1/49 by a factor of 7, so the floor sum is ≤ 0 on 400/400 random comparable 7-blocks; the wall bound DOES close an isolated 7-comparable block (S350's containment floor is positive unconditionally, and positivity is all seven_comb_wall needs), but the 13-speed dense core needs a QUANTITATIVE window to nest the remaining runners, and quantitative floors exist only in the separated regime b > 7a — which the dense core is defined to exclude
status: negative result, verified exactly (400/400 comparable 7-blocks have THM-1012 floor sum ≤ 0 while S350 floor sum > 0 in 400/400; the EXACT Hunter floor is positive throughout — only the available LOWER BOUND fails)
source: opus-2026-07-17-S356 (owner: use the wall bound to close the residual dense core)
depends_on: [THM-1012 / LRCPairIndependence (the sharp floor, sharp only when separated), LRCPairOverlapFloor (S350, positive but tiny), LRCSharpWallBound (seven_comb_wall), THM-964]
scripts: 04-computation/dense_core_169_opus_S356.py -> 05-knowledge/results/dense_core_169_opus_S356.out
---

# THM-1021 — why the dense core resists the wall bound

> **ANSWERED (opus-S357) by THM-1025.** The prediction at the end of this file -- that closing the dense core needs a floor which SEES the correlation rather than approximating it away -- is exactly right, and the sawtooth supplies it: rho >= 1/49 - g^2/(4ab) from THM-965's fold bound, sharp precisely when ab is large (the comparable regime), complementary to THM-1012's b >> a. With dilation invariance and a 15-pair finite table the floor exists in EVERY regime, and the comparable-7-block floor sums are positive 400/400 where THM-1012's were <= 0 on all 400. This file's negative stands as stated -- the INDEPENDENCE bound really does fail there -- but the regime is no longer open.

**What the wall bound does close.** `seven_comb_wall` needs only a POSITIVE
pair-floor sum, and S350's containment floor 2λ/b is positive
unconditionally. So any isolated 7-comparable block at λ = 1/14 yields a
lonely point of the window — kernel-pure, no further hypotheses. (This is
the easy half: 7·2λ = 1 exactly, and the combs all meet near 0, so the
uncovered set is nonempty.)

**What it does not close.** The residual dense core is a 13-speed packet
admitting no block partition. Handling it means nesting the remaining six
runners inside the window left by the first seven — and nesting needs a
QUANTITATIVE window, not merely a nonempty one. The only quantitative
per-pair floor available is THM-1012's

  ρ(a,b) ≥ 4λ² − 2λ·(a/b),

whose defect tends to 2λ as a → b. At λ = 1/14 the defect 1/7 exceeds the
leading term 1/49 by a factor of **7**, and the crossover is exactly
a/b < 2λ, i.e. **b > 7a**. Measured: on 400 random comparable 7-blocks
(ratio ≤ 13) the THM-1012 floor sum is ≤ 0 in **400/400**, while the exact
Hunter floor is positive throughout. The bound fails, not the mathematics.

**The structural statement.** The dense core is precisely the unseparated
regime, and every quantitative pair floor this program has produced is
sharp only in the separated regime. That is not an accident of technique:
the independence heuristic ρ ≈ (2λ)² is a statement about combs whose
scales decouple, and comparable combs are correlated by construction (they
share near-coincidences at every level of the beat lattice). Closing the
dense core therefore requires a floor that is NOT of independence type —
the exact sawtooth (THM-965) or the beat/E-family structure (THM-928(B)),
which see the correlation rather than approximating it away.

**Honest ledger.** The 7-wall's existence conclusion: CLOSED, kernel-pure.
Its sharp quantitative bound in the separated regime: CLOSED, kernel-pure.
The dense core: OPEN, and now with an exact reason.
