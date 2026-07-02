# HYP-3951: The explicit PSL_2(F_7) tensor code — MEASURED soundness s vs the 1/36 family

**Status:** CLAIMED (stub) — kind-pasteur-2026-07-01-S28, in progress this session.
**Numbering note:** adopting klein-S68's per-machine HYP-block proposal (kps = 3950+).

## Context
opus-S30 built the PSL_2(F_7) left-right Cayley complex; opus-S31 tabulated the *theoretical*
trickle-down bound rho >= delta_inner x (1 - lambda_link) = delta_inner (links are K_{A,B}, lambda = 0).
mac-mini-S92 corrected: the single-generator complex is DISCONNECTED (8 components = Aut(Paley_7) = 21 each),
soundness O(1/L). klein-S86: (2,3,7)x<7> gens give b1=14 cosystoles; QR certificate is a coboundary.

**Gap in all three:** nobody has MEASURED the actual soundness s of a concrete connected tensor code.
opus's bound ignores global Cayley expansion entirely (lambda_link = 0 is automatic for ANY left-right
complex, including the disconnected one — so it cannot be the discriminating quantity; the real DELLM
soundness needs the global gap).

## Claim (to be tested this session)
1. Explicit definition: G = PSL_2(F_7), symmetric generating sets A, B with the DOUBLE-COSET connectivity
   check <A> g <B> = G (mac-mini's correction made structural), inner codes C_A, C_B (parity/repetition/
   extended-Hamming), code C = {f: squares -> F_2 : every vertex's A x B local view in C_A (x) C_B}.
2. MEASURED soundness: s = min over f not in C of [viol(f)/V] / [dist(f,C)/F], computed via
   (a) exhaustive low-weight cosets (weight < d/2 are coset leaders), (b) adversarial strings/tubes
   (the surface-code killers), (c) random probes. Compare against the trickle-down prediction and
   against the 1/36 family: 1/36 = 0.02778 (analytic floor), 313/9702 = 0.03226 (census min),
   (delta/2)^2 = 1/36 at delta = 1/3 (= [6,5,2] parity at degree 6 = phi(7) — the natural q=7 degree).
3. Test: does the natural degree-6 (= phi(7)) parity tensor code have s ~ 1/36? Honest outcome either way:
   a numeric match makes the soundness/measure mirror QUANTITATIVE; a mismatch pins it as structural-only.

## Evidence / artifacts (to be filled)
- Script: 04-computation/psl27_tensor_code_measured_soundness_kps.py
- Results: 05-knowledge/results/psl27_tensor_code_measured_soundness_kps.out

## Depends on / relates to
HYP-3823 (opus complex), HYP-3824 (mac-mini disconnection), HYP-3832 (klein cochain), kps-S27
(LTC upgrade sketch, sqrt-p = Gauss certificate), DELLM (Annals 2026 203-2 / STOC 2022).
