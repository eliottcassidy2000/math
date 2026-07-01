# SDP/Lasserre tooling for LRC(14): I built the moment-SOS relaxation for M(S)=max_t min_v‖vt‖ (univariate Chebyshev reformulation c=cos2πt, sin²(πvt)=(1−T_v(c))/2, so M(S)=arcsin√M̃/π), validated it EXACTLY on small cases ({1,2,3}→1/4, {1,2,4}→1/3), but the free solvers (SCS inaccurate on rank-deficient optima; CLARABEL fails on feasibility) don't scale, AND the exercise is REDUNDANT — M(S) is exactly computable by rational sweep in milliseconds, so the M(S)-SDP adds only a certificate, not a value; the CORRECT SDP target for LRC(14) is the p0/L_y magic function (THM-534, already the optimal LP dual) with the crux being the COMBINATORIAL extremality "consec maximizes L_y" reduced to bounded by the THM-546 far-element decorrelation — not M(S)

*opus-2026-07-01. Owner: pivot the endgame to LRC's 1/n (after the covering-min refutation, HYP-3782), build
the E₂/F₇ magic function + SDP/Lasserre tightening, install cvxpy, do heavy SDP work. Built the tooling; the
honest finding is that the M(S)-SDP is the wrong (redundant) target and the p0/L_y route is the right one.*

## What I built (cvxpy 1.9.2, SCS + CLARABEL)
The moment-SOS SDP for `M(S)=max_t min_v ‖vt‖`, via the clean **univariate** reformulation:
`c=cos(2πt)∈[−1,1]`, `sin²(πvt)=(1−T_v(c))/2` (Chebyshev), so `M̃(S)=max_{c∈[−1,1]} min_v (1−T_v(c))/2` and
`M(S)=arcsin(√M̃)/π` (exact, since sin² is monotone on the binding runner).
- **Moment relaxation** (primal, measure on {c: g_v(c)≥s ∀v}, bisect s): tight for univariate at finite order.
- **SOS dual** (single solve): `min γ s.t. γ−Σ_v λ_v(c)g_v(c)=σ₀+(1−c²)σ₁`, `λ_v,σ` SOS, `Σλ_v(c)=1`.
Validated EXACTLY: `S={1,2,3}→M̃=0.5⇒M=1/4`; `S={1,2,4}→M̃=0.75⇒M=1/3`.

## Why it's the wrong target (honest)
1. **Numerical:** the optima are rank-deficient (a Dirac at the lonely c), so SCS is inaccurate (non-monotone
   in the relaxation order R — R=8 gave 0.360, R=12 gave 0.368 for {1,2,3,4}, true 0.345), and CLARABEL errors
   on the pure-feasibility form. The SOS dual with SOS multipliers converges slowly (loose at low degree).
   Scaling to the 13-speed sets on the free solvers is unreliable; MOSEK (licensed) would be needed.
2. **Redundant:** `M(S)` is computable EXACTLY by a rational t-sweep in milliseconds (used throughout this
   project). The SDP adds only the *certificate* (the magic-function multipliers λ_v(c)), not the value — so
   for LRC it is not the leverage point.
3. **Wrong direction for LRC:** `M(S)≥1/n` is an EXISTENCE statement (a lonely point). SDP moment relaxations
   give UPPER bounds on `max_t min_v` — the wrong direction. The lonely point is exhibited directly (t=1/n for
   the AP); no SDP needed.

## The right SDP target: the p0/L_y magic function (THM-534)
The project's actual LRC(14) endgame is the **mod-7 coverage** proxy `p0(E)=meas(S7)` with the PROVED
Bonferroni moment-LP dual `p0(E) ≤ L_y(E)=Σ_r y_r S_r` (THM-534), whose optimal magic function is already found
(e.g. `g(t)=(t−1)(t−2)(t−4)(t−5)/40` for k=8, an SOS-nonneg-on-{0..6} certificate). This is a robust **LP**
(finite, degree ≤6), not a nasty SDP. The genuine crux is the **combinatorial extremality** "consec (the AP)
maximizes L_y over E" — verified on bounded spread, open for large spread — and the **THM-546 far-element
decorrelation** (`|Δ_W|≤(6/49)V/W`) reduces "all E" to bounded E + a gapped tail. THAT is where SDP/Lasserre
tightening belongs: a moment relaxation over the sector-hit distribution to prove "consec maximizes L_y," and a
better solver (MOSEK) to push the bounded check + the gapped cutoff `w*`.

## Recommendation
- The M(S)-Lasserre tool is a validated building block (certifies specific M-values with a magic-function
  certificate) but is not the LRC leverage point — filed for reuse, not pursued at scale.
- Heavy SDP lifting for LRC(14) = the **L_y extremality** (moment method on the sector distribution) + **MOSEK**
  for accuracy on the bounded finite check and the far-element gapped cutoff. The magic function (E₂/F₇/Fejér
  F₇, HYP-3214/3227) is the certificate basis there, and it targets the REAL `M≥1/n` conjecture (AP/GW tight).

## Status
- **Built + validated (opus):** cvxpy moment-SOS for M(S) (exact on {1,2,3}, {1,2,4}); univariate Chebyshev
  reformulation.
- **Honest limits (opus):** numerically unreliable at scale on free solvers; redundant with exact M-computation;
  wrong direction (upper bounds) for the LRC existence statement.
- **Redirect (opus):** the LRC(14) SDP endgame is the p0/L_y magic function (THM-534, LP) + the combinatorial
  "consec maximizes L_y" crux + THM-546 far-element decorrelation; MOSEK for scale.
- **Context:** this follows the covering-min refutation (HYP-3782) — the pivot to the real conjecture M≥1/n.

Related: HYP-3782 (covering-min refuted → pivot to 1/n), THM-534 (p0≤L_y magic function), THM-546 (far-element
decorrelation), HYP-3214/3227 (Fejér magic function), HYP-2974 (Toeplitz PSD). Script: 04-computation/lrc_lasserre_sdp_Mofs_opus_20260701.py.
