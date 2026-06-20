# HYP-2675: LRC(14) wide-branch sector ridge

**Status:** claimed / in progress (codex-s47, 2026-06-20).

The KPS comfortable-margin route reframes the remaining sector proof obligation
as a direct wide-row bound:

> for primitive `E` with `0 in E`, `|E|=k=8..12`, and `span(E)>14`, prove
> `p0(E)=meas(S7(E)) <= cap_k`.

This hypothesis reserves the next scout around that obligation.  The intended
test is not a standalone `Delta_w` constant.  It records exact `p0`, missed
sector distributions, additive-energy/Freiman fingerprints, squarefree divisor
profiles, shell-1 occupancy, state-word entropy, and a row-risk Tournament
Analysis for wide rows.

**Claim being tested:** the only rows that can approach the direct cap in the
wide branch retain a low-rank near-consecutive/GAP scaffold.  Once the
state-word support spreads across genuinely wide additive structure, `p0` drops
with comfortable margin; the old large-`Delta_w` resonances are compensated by
small plateau.

**Artifact:** `04-computation/lrc14_wide_branch_ridge_codex_s47.py`.

**Open:** no theorem is claimed yet.  The missing proof is still the analytic or
finite-to-infinite statement turning these exact ridge fingerprints into
`span>14 => p0(E)<=cap_k`.

