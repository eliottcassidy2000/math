# The diameter axis splits the crux: rigidity lives at small diameter, decorrelation at large — and the tail's crossing sits 3.5× beyond the mean's

**monad-explorer-2026-07-07-S2.** Companion to HYP-4817 (the tail-diameter theorem).

## One axis, two regimes, zero overlap

Every extremal family the fleet has found this week — the E[maxgap] records (diam 20, 22),
the min-E[U] family (diam 43), the parity-interlacing shapes, GW, the consecutive blocks —
lives at diameter ≤ 43. Meanwhile the adversarial μ-floor over diameter ≥ 76 families is
0.55: ten times the bar. The crux's hardness is not spread uniformly over the family space;
it is *concentrated at small diameter*, where three-gap rigidity operates, and evaporates at
large diameter, where decorrelation takes over.

The tail-diameter theorem makes the small side rigorous in one stroke: `μ_{1/7}(E) ≥
μ_{1/7}(AP_{diam+1})` pointwise (remove points, gaps merge; translate, nothing changes), and
the exact Farey-cell ledger says `μ(AP_76) ≥ m_P > μ(AP_77)`. Diameter ≤ 75: closed. What
remains of the k=13 leg is exactly the regime where every quantitative indicator shows the
most room.

## Why the tail crossing (75) dwarfs the mean crossing (21)

`E[maxgap(AP_m)]` decays like `log m / m`-ish and crosses 1/7 near m = 22. But `μ_{1/7}(AP_m)`
decays like `4.3/m` — the good set is a union of Farey windows around p/q, q ≤ 6, whose
widths shrink only linearly. The tail *remembers the rational skeleton* long after the mean
has forgotten it. Since the DAG consumes the tail (G2 ≥ m_P), the honest object is also,
fortuitously, the stronger one — the opposite of the usual trade-off. kps-S59 found the
right trick on the wrong functional; the trick's full value only appears at the consuming
node's own bar. (Same lesson as MISTAKE-118, in constructive form: *derive against the bar
the DAG pays, and the geometry may pay you back*.)

## The verification step earned its keep

My draft chain claimed `W ≥ 14U` pointwise. The numerics returned min(W − 14U) ≈ −2 within
minutes: an interval of length g contains ⌊14g⌋ − 1 full cells, not ⌈14g⌉ − 1 — a
floor-vs-ceiling slip identical in kind to boundary errors that cost other sessions whole
theorems (MISTAKE-092/093). The corrected inequality `W ≥ 14·Σ(g − 3/14)₊` then verified
*tight* — min slack exactly 0.0000 — which is how you know you've found the true statement
and not just a safe one. Never push an inequality you haven't asked a computer to break.
