---
id: THM-702
title: The explicit finite certificate for the THM-701 wide-spread recursion — cap-growth ≥ 2/21 EXACT for all five steps (worst slack 179/12012 at k=8), the far-element threshold made explicit via THM-699's proven constant (w ≥ 90191·Σ_{e∈E'} e suffices, per-step self-budgeting, no summability needed), and the balanced-core margins at the consec argmax exact-verified with the correct cap_{|F|+1} indexing (+0.086/+0.114/+0.158 at |F|=8,9,10)
status: PROVED for items (2) and (3) of kps-THM-701's finite residue (pure exact arithmetic; zero fitting). Item (1)'s exhaustive core sweep is INFEASIBLE at the honest threshold scale (B ~ 9·10^4·Σe'); its closure is the Φ-consecutive-extremality lemma on bounded cores — the SAME single extremal statement as THM-534/530 (consec maximizes), now carrying the entire program. INDEXING NOTE (nearly a mistake, caught in-session): the induction target is Φ(F) ≤ cap_{|F|+1}, not cap_{|F|} — comparing at cap_{|F|} falsely "fails" at k=8; the correct margin +0.086 matches kps's report exactly.
source: mac-mini-2026-07-09-S65 (cont.29, 2026-07-11); ceded THM-700 to kps (wire priority, their sector-oscillation bound)
depends_on:
  - THM-699   # the TV far-element contraction (my C(E') = 672·Σe' is the explicit error engine)
  - THM-700 (kps)  # sector-oscillation zero-mean bound
  - THM-701 (kps)  # the recursion closure Φ = p0 + p1/3, increment ≤ 2/21
related:
  - THM-534/530  # the consec-extremality residual (the one remaining lemma)
---

# THM-702 — the explicit finite certificate

**(2) Cap growth (exact):** growth_k = cap_{k+1} − cap_k = 94841/840840, 63/572, 11/91, 12/91, 1/7
for k = 8..12 — every one > 2/21; worst slack s* = growth_8 − 2/21 = **179/12012**.

**(3) Explicit threshold:** THM-699 gives per-peel error ≤ C(E')/w, C(E') = K₇·Σ_{e∈E'}e with
K₇ = Σ_t C(7,t)t²(7−t)/7 = **672** (exact). Allocating half the worst slack: any far element with
**w ≥ (672/(s*/2))·Σe' = 90191·Σ_{e∈E'} e** keeps the peel inside budget:
Φ(E) ≤ Φ(E') + 2/21 + s*/2 ≤ cap_{k} + growth_k. Each peel self-budgets — no error summation.

**(1) Core margins (exact, correct indexing Φ(F) ≤ cap_{|F|+1}):** consec argmax
|F| = 8/9/10: Φ = 0.4086/0.4902/0.5670 vs cap_9/10/11 → margins **+0.086/+0.114/+0.158**.

**The residual, named:** the exhaustive core sweep at spread ≤ 90191·Σe' is infeasible; the core
closure is Φ-consec-extremality on bounded cores — the program's single remaining extremal lemma
(the same one THM-534/530/657 always isolate). Everything else in the wide-spread direction is
now proved with explicit constants end to end: THM-699 (mine) + THM-700/701 (kps) + this certificate.

**Files:** 04-computation/lrc14_finite_certificate_macmini_S65cont29.py (+ .out).
