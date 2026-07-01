# MOSEK for LRC(14): it UNBLOCKS the per-set moment-SOS M(S) SDP (exact + flat-extension where SCS was non-monotone — the right solver for mac-mini-S89's atomic covering-min certificates), but there is NO SDP/LP shortcut to the FOR-ALL-S crux — every naive relaxation of "consec maximizes L_y" (THM-534) / the covering-min lower bound returns the TRIVIAL bound 1.0, because the offset-set achievability is the 1-parameter singular series and the extremality is coupled/signed (THM-534's alternating monotonicity = HYP-2738) — so MOSEK confirms the far-element/equidistribution route (klein-S66, mac-mini-S74/89) is the necessary path, not an SDP

*opus-2026-07-01. Owner gave a MOSEK license and asked for maximum-leverage creative computation toward finishing
LRC(14). I installed+licensed MOSEK and ran a real SDP campaign. Honest verdict: MOSEK is the right solver for
the tractable per-set SDPs, but the crux is not an SDP — it is the singular series, exactly where the team's
far-element work lives.*

## What MOSEK unblocked (positive, concrete)
The **M(S) = max_t min_v ‖vt‖ moment-SOS** (univariate Chebyshev: `M̃ = max_{c∈[−1,1]} min_v (1−T_v(c))/2`):
- **MOSEK is EXACT where SCS failed:** `{1,2,3,4}→M̃=0.345492` (gap 0.000, moment-matrix rank **2** = flat
  extension = 2 atoms), `{1,2,4}→0.75` (rank **1**). SCS gave non-monotone garbage (0.360 at R=8, 0.368 at R=12).
- **Flat extension = Curto–Fialkow certificate:** the rank of the optimal moment matrix = the number of atoms
  (optimal lonely times). This is precisely **mac-mini-S89's atomic reframe** (construction = 2 atoms, AP = 6);
  **MOSEK is the accurate solver that approach needs** — the free solvers can't do it reliably.
- **Limit (Lasserre, not MOSEK):** the AP (13 speeds) converges SLOWLY (R=16→0.0915, R=32→0.0775, true 0.0495)
  because the univariate degree = max speed and the min-envelope has many near-active constraints. The
  construction (max 182) is worse. So per-set certification is practical only for small max-speed sets.

## Where MOSEK gives NO shortcut (the honest negative)
The **p0/L_y sector crux** ("consec maximizes L_y", THM-534) — the actual open piece — is a moment-functional
extremality over **offset sets**. Every relaxation I could pose is trivially loose:
- **Boolean-moment relaxation** over the sector-empty indicators `χ_j` (moment matrix `M[a,b]=J(a∨b)⪰0`) with
  the single-runner bound `J(A)≤(7−|A|)/7` + monotonicity: **MOSEK returns 1.0** (sets all `J(A)=0`, i.e.
  "sectors never empty" — which no offset set achieves: at t≈0 all runners pile in sector 0).
- **Even pinning `S_1 ≥ S_1(consec)`** (a real universal lower bound — consec minimizes S_1): **still 1.0.** The
  reason is exactly THM-534's warning: `L_y = 1 − S_1 + S_2 − 0.9 S_3 + 0.6 S_4` is a **coupled ALTERNATING**
  functional — consec minimizes S_1 but maximizes S_2, S_4 *simultaneously*; bounding moments separately is
  circular. This is **HYP-2738** (no nonnegative/per-moment certificate; the extremality is irreducibly signed).
> **The missing structure is the 1-parameter singular series:** the runners are all dilations `e·t` of ONE t, so
> the joint sector occupancy is a single-parameter curve on `{0,…,6}^k`, NOT an arbitrary joint distribution.
> The boolean/moment cone ignores this; encoding it IS the theta/singular-series lattice (THM-515/523). No clean
> SDP captures it.

## Verdict + recommendations
- **MOSEK = the right solver for mac-mini-S89's flat-extension / atomic covering-min certificates** (per-set,
  small max speed). Recommend the team route those SDPs through MOSEK (I've wired it into cvxpy;
  license at `~/mosek/mosek.lic`).
- **No SDP/LP shortcut to the for-all-S closure.** The `consec maximizes L_y` and covering-min-lower-bound cruxes
  ARE the singular-series achievability + a coupled signed extremality (HYP-2738); they are not clean
  moment programs. This **confirms the team's far-element/equidistribution route** (klein-S66's rigorous
  `|Δ_w|=O(1/w)`, mac-mini-S74/89's equidistribution + flat extension) is the necessary path — and my earlier
  two-regime synthesis (HYP-3790: L_y = bounded finite check + far equidistribution) is where the SDP hands off
  to the analytic bound.
- **The one hybrid that could work:** feed the ANALYTIC singular-series bounds on `J(A)` (from the theta lattice)
  into a MOSEK moment-LP — a hybrid, not a pure SDP. That is the only way the moment machinery bites, and it
  reduces to computing the theta bounds (the team's active work), so MOSEK is the closer, not the opener.

## Status
- **Positive (opus, MOSEK):** M(S) moment-SOS exact + flat extension for tractable sets (SCS couldn't); MOSEK is
  the solver for mac-mini-S89's atomic certificates.
- **Negative (opus, MOSEK):** the p0/L_y (and covering-min for-all-S) crux is loose in every naive relaxation
  (=1.0) — the singular-series achievability + coupled signed extremality (HYP-2738) is irreducible; no SDP
  shortcut.
- **Confirms:** the far-element/equidistribution route (klein-S66, mac-mini-S74/89) is the necessary path; the
  SDP is the per-set closer (flat extension), not the for-all opener.

Related: THM-534 (p0≤L_y, coupled signed extremality), HYP-2738 (signed-certificate impossibility), HYP-3789/mac-mini-S76
(M(S) moment relaxation + flat extension — MOSEK is its solver), HYP-3790/opus (L_y two-regime), HYP-3787/klein-S66
(far-element O(1/w)), HYP-3786/mac-mini-S74 (equidistribution), THM-515/523 (singular series). HYP-3791 (this).
Scripts: 04-computation/lrc_mosek_{mosek_M,Ly_sdp,Ly_achievable}_opus_20260701.py.
