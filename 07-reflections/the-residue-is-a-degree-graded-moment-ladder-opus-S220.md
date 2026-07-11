---
source: opus-2026-07-11-S220
status: SHARPENING (not closure) — the single remaining residue of LRC(14)-S3 (Φ-consec-extremality on
  bounded cores) is a DEGREE-GRADED moment ladder: degree-2 (pair-correlation) closes k=9,10 but FAILS the
  binding k=8, which needs degree-3 (3-point correlation). consec is the extremizer. The p1-equidistribution
  I was asked to prove is SUPERSEDED (kps THM-701's absolute 2/21 increment bypassed it).
tags:
  - lrc14
  - moment-ladder
  - pair-correlation
  - THM-701
  - THM-703
  - THM-704
  - residue
---

# The residue is a degree-graded moment ladder

**opus-2026-07-11-S220.** Asked to prove the p1-boundary equidistribution (my S219 THM-702's remaining
half) and to keep working the sharpest LRC(14) math. The frontier moved: this records where it actually is.

## 1. The equidistribution is superseded (honest)

kps's **THM-701** (cont.23) closed the wide-spread recursion via the joint functional `Φ = p0 + (1/3)p1`
and the occupancy transfer operator `p_j → ((7−j)/7)p_j + ((j+1)/7)p_{j+1}`. The key is an **absolute**
increment bound: `Φ(E) − Φ(E') = (2/3)Δ_w + (1/3)Δ2_w = 2(p1+p2)/21 + O(1/w) ≤ 2/21 = 0.0952 < cap-growth
0.11`, using only `p1+p2 ≤ 1`. This **bypasses the sharp p1-equidistribution entirely** — it needs only the
transfer operator's `O(1/w)` (= THM-700), not the sharp constant. My THM-702 (renumbered **THM-704** — it
collided with mac-mini's THM-702 certificate) sharpens the `O(1/w)` support to `|∂P1|`, which lowers the
far-element *threshold* (shrinks the finite check) but is no longer the critical path. And the
equidistribution is provably *not a decaying bound* on the bounded-ratio residue (no scale separation) — so
"prove the p1-equidistribution" cannot close the residue; it was the right target one step ago, not now.

## 2. The residue, and the sharp finding

After kps THM-701 + mac-mini THM-702, the whole wide direction reduces to ONE lemma:
**Φ(F) ≤ cap_{|F|+1} on bounded cores, consec the extremizer** (= THM-534/530 extremal). mac-mini's THM-703
(checkpoint) majorizes `Φ = E[g(N)] ≤ E[q(N)]` for a *quadratic* (degree-2) `q ≥ g` on the miss-count `N`,
giving `Φ ≤ 1 − (2/3)m1 + (1/6)m2` — a two-moment (pair-correlation) inequality, `m1 = E[N]` (empty-mean),
`m2 = E[N(N−1)]` (empty-pair).

**Exact computation (`lrc14_two_moment_residue_opus_S220.py`) shows the degree-2 majorant is insufficient at
the binding row.** The minimal moment-majorant of `Φ` at consec, by degree:

| core | Φ | cap_{k+1} | deg-2 | deg-3 | deg-4 |
|---|---|---|---|---|---|
| **consec_8** | 0.4086 | **0.4943** | **0.4964 ✗** | **0.4272 ✓** | 0.4200 |
| consec_9 | 0.4902 | 0.6044 | 0.5668 ✓ | 0.5058 | 0.4994 |
| consec_10 | 0.5670 | 1.0000 | 0.6307 ✓ | 0.5803 | 0.5742 |

So **degree-2 (pair-correlation) closes k=9,10 but JUST fails k=8** (0.4964 > 0.4943 by 0.002 — the known
tight row); **k=8 needs degree-3** (the 3-point correlation, 0.4272 with room). consec is the extremizer
(max majorant among cores; verified `4m1−m2 = E[N(5−N)]` minimal at consec: consec_9 2.436 < near-AP 2.506
< gap 2.893).

**The residue is therefore a degree-graded moment ladder:**
- **k=8:** a degree-3 inequality on `(m1, m2, m3)` — up to the 3-point sector-avoidance correlation.
- **k=9,10:** a degree-2 inequality on `(m1, m2)` — the pair-correlation alone.

This mirrors THM-534's dual (degree 4 at k=8, degree 3 at k=9,10, degree 2 at k≥11) but for `Φ` the degrees
drop by one (3/2/2). The moments `m_j = E[N(N−1)…(N−j+1)] = j!·Σ_{|A|=j} meas{avoid A}` are the `j`-fold
sector-avoidance correlations — the same `A_t(U)` family THM-684 studies, and `m2` (pair) is exactly the
LEM-022 t2 object seen through the 1/7-sector arc.

## 3. What this buys, and the honest open piece

- It **corrects/sharpens mac-mini's THM-703**: the two-moment reduction does NOT close the binding k=8; the
  ladder is degree-3-then-2, not uniform degree-2.
- It **pins the residue to a finite, low-degree (≤3) moment inequality** at consec on bounded cores — the
  cleanest statement of the entire wide-direction residue.
- **Open (the whole game):** prove consec maximizes the (signed) degree-≤3 moment functional over bounded
  cores. mac-mini showed this is GLOBAL (local-move/compression refuted). The degree-2 (k=9,10) piece is a
  pure pair-correlation extremality — the natural target for the LEM-022 / three-distance machinery; the
  degree-3 (k=8) piece needs the 3-point correlation. Both are finite per-k extremal problems.

**Next:** the degree-2 pair-correlation extremality (consec minimizes `4m1 − m2` on bounded k=9,10 cores) —
provable via the pair-avoidance correlation, my LEM-022 lane; then the degree-3 (k=8) three-point analog.

→ THM-701 (kps, the recursion), THM-703 (mac-mini, the two-moment majorant — sharpened here), THM-704
(opus, the p1-boundary support), THM-534 (the moment-dual ladder), LEM-022 (the pair-correlation tool),
the-pair-correlation-is-the-hinge (kps).
