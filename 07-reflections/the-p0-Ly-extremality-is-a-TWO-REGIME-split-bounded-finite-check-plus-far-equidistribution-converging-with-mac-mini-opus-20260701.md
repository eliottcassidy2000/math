# The p0/L_y extremality ("consec maximizes L_y", THM-534) is a TWO-REGIME split — STRUCTURED/bounded (consec is the max, verified) + SPREAD/far (L_y → the occupancy/equidistribution limit, FAR below the cap) — and the spread limit IS mac-mini's far-element equidistribution (HYP-3786); moreover the naive boolean-moment SDP for max_E L_y is trivially LOOSE (=1) without the offset-set achievability, and that achievability is exactly the singular series where klein-S66's rigorous O(1/w) far-element correction bites — so the sector (p0/L_y) route and the far-element (equidistribution) route are the SAME endgame

*opus-2026-07-01. Owner asked to build the p0/L_y moment relaxation (task a) for "consec maximizes L_y". I
verified the target, found the naive SDP loose (achievability is the crux), and the crux resolves into a
two-regime split that converges with the team's active far-element work. (This session I also RETRACTED my
erroneous HYP-3782 — the covering-min=14/183 stands; mac-mini-S74 caught my q-range error.)*

## Verified the target (due diligence, after today's error)
THM-534: `L_y(E)=Σ_r y_r S_r(E)=Σ_t g(t) p_t(E)`, `p_t=meas{exactly t of sectors 1..6 empty}`, and consec (the
AP) maximizes L_y. Reproduced EXACTLY: **`L_y(consec_8)=0.35823=2633/7350`**, and it is the **max over 3000
offset 8-sets** (0 exceed it; cap_8=0.3815). THM-534 is solid.

## The naive SDP is loose — achievability is the crux
The clean SDP idea — maximize `L_y=Σ_r y_r e_r(χ_1,…,χ_6)` over the boolean-moment cone of the sector-empty
indicators `χ_j` (moment matrix `M[B,B']=J(B∪B')⪰0`) — is **trivially loose**: without the offset-set
achievability, `max E[g(N)] = max_t g(t) = 1` (put all mass on a single config). So the SDP does NOT bypass the
hard part. The hard part is **which moment vectors `J(B)=meas{all runners avoid the sectors in B}` are realized
by an offset set** — i.e. the **singular-series / theta structure** (THM-515/523). The SDP relaxation is exactly
as hard as the analytic extremality; the leverage is the achievability, not the moment cone.

## The two-regime split (the synthesis)
Computing `L_y(E)` as E spreads (offsets 7→21→70→280→1050) shows it **monotonically drops from consec to a
fixed limit**:
| regime | L_y(k=8) | mechanism |
|---|---|---|
| **structured / bounded** (consec) | **0.358** = MAX (< cap 0.382) | THM-534 finite check |
| intermediate spread | 0.226 → 0.066 → 0.048 | transition |
| **far / spread limit** | **0.049** ≪ cap | i.i.d. sector **occupancy** |
The far limit is the **equidistribution/occupancy value**: k−1 runners i.i.d. over 7 sectors give
`S_r = C(6,r)·((7−r)/7)^{k−1}`, so `L_y → Σ_r y_r C(6,r)((7−r)/7)^{k−1} = 0.0493` (k=8) — matching the spread
computation. **This IS mac-mini's HYP-3786** ("huge multi-patch = equidistribution on the lonely set"): as the
offsets grow, the sector occupancy equidistributes and every moment `S_r` tends to its i.i.d. value, dragging
`L_y` far below the cap.
> So **"consec maximizes L_y" = bounded finite check (THM-534) + far equidistribution (mac-mini HYP-3786 /
> klein-S66's rigorous O(1/w) correction)** — the SAME two-regime endgame as the far-element decorrelation
> (THM-546). The dangerous region is the STRUCTURED/near-consec bounded shapes (tiny cap slack, margin 0.024 at
> k=8); the spread region is comfortable (occupancy ≈ 0.049 ≪ 0.382, huge margin) and rigorously bounded by the
> far-element Fourier decay.

## What this gives / what remains
- **Gives (opus):** the sector (p0/L_y) route and the far-element (equidistribution) route are one problem; the
  spread regime of L_y is exactly the occupancy limit klein/mac-mini are bounding rigorously; the SDP crux is
  achievability, not the moment cone.
- **The remaining piece is the SAME as the far-element endgame:** a rigorous cutoff `w*` separating the
  bounded finite check (THM-534, verified spread ≤ 16) from the far equidistribution (comfortable), with the
  intermediate offsets controlled by klein-S66's `|Δ_w| = O(1/w)` correction. No new obstruction beyond the
  team's active `w*` work.
- **Recommendation:** don't build a standalone p0/L_y SDP (loose without achievability, and the achievability =
  the singular series the far-element work already attacks). Instead, **express the far regime of L_y via the
  klein-S66 correction** (`L_y(spread) ≤ occupancy + O(1/w)`) to certify "far L_y < cap", and reserve the SDP/LP
  for the BOUNDED finite check (where it's a small robust LP over the actual moments). MOSEK would help the
  bounded check at scale.

## Status
- **Verified (opus):** consec maximizes L_y (k=8, exact 2633/7350, max/3000 sets).
- **Found (opus):** the boolean-moment SDP is loose without achievability (=1); the achievability is the
  singular series.
- **Synthesis (opus):** L_y two regimes — bounded (consec MAX) + spread (occupancy limit ≪ cap = mac-mini
  HYP-3786); the p0/L_y and far-element routes converge; the remaining crux is the shared cutoff `w*`.
- **Retraction (opus, this session):** HYP-3782 was my error (covering = q=2..n incl q=n; my counterexamples
  miss q=14); covering-min=14/183 STANDS.

Related: THM-534 (p0≤L_y, consec maximizes), HYP-3786/mac-mini-S74 (equidistribution on the lonely set), HYP-3787/klein-S66
(rigorous O(1/w) far-element correction), THM-546 (far-element decorrelation), HYP-2738 (coupled/signed
extremality — why bounding moments separately fails), HYP-3782 (retracted). HYP-3790 (this).
Scripts: 04-computation/lrc_Ly_two_regime_synthesis_opus_20260701.py, lrc_verify_consec_maximizes_Ly_opus_20260701.py.
