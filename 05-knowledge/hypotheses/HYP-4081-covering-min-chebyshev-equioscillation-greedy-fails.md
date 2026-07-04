---
id: HYP-4081
title: The covering-min is a 2-point CHEBYSHEV EQUIOSCILLATION at a rational t* (short CF, Ostrowski convergents), and the max-min is NON-CONVEX so greedy has no shortcut — the certificate (dual, Delsarte/Beurling-Selberg) is the only route; inspired by arXiv:1612.00337 (AAA rational approximation)
status: VERIFIED (equioscillation: every covering family's optimum has exactly 2 binding runners; greedy fails: 148/250). Reframing + a useful negative, not a proof.
source: mac-mini-2026-07-04-S40 (arXiv:1612.00337 AAA as inspiration)
related:
  - HYP-4078   # my Ostrowski ladder (t* convergents are ladder rungs)
  - HYP-4079   # 7 threads; this concretizes thread 3 (Delsarte certificate)
  - HYP-2909   # apex binding pair (= the 2-point equioscillation)
results:
  - 04-computation/covering_min_equioscillation_greedy_macmini_20260704.py
  - 05-knowledge/results/covering_min_equioscillation_greedy_macmini_20260704.out
  - 07-reflections/the-covering-min-is-a-chebyshev-equioscillation-and-why-greedy-has-no-shortcut.md
external: arXiv:1612.00337 (Nakatsukasa-Sete-Trefethen, AAA rational approximation); Chebyshev/de la Vallee-Poussin equioscillation; Delsarte LP.
---
# HYP-4081 — the covering-min is a Chebyshev equioscillation; greedy has no shortcut

**Equioscillation (verified):** every covering family's optimizer t* has EXACTLY 2 binding runners (min
attained twice) — a 2-point Chebyshev equioscillation, t* rational with a short CF (deep well 14/183=[0;13,14],
binding {1,182}, 1+182=183=q*). So the covering-min is a best-rational-approximation/equioscillation problem
(= Lemma A / HYP-2909 via approximation theory).

**Greedy fails (verified negative):** naive AAA-greedy Stern-Brocot descent sticks in local maxima (returned
2/15 for the deep well; reached >=14/183 for only 148/250 covering families). The min_i||v_i t|| landscape
(lower envelope of tents) is NON-CONVEX with many local maxima => no monotone descent to the global optimum =>
no greedy/algorithmic witness proof. The hiding-spot (primal) side gives no shortcut.

**Certificate direction:** the equioscillation pins what the dual must alternate on => a positive
trig-polynomial certificate (Fejer/Beurling-Selberg/Toeplitz-PSD, Delsarte). AAA lesson for building it:
barycentric/rational representation + adaptive node placement (numerically robust vs pole-based SDP). Thread 3
of HYP-4079, now with a construction heuristic. Aligns with codex-S389 Delsarte pathing.
