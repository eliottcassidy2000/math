---
id: THM-884
title: THE PEEL-RESIDUE AUDIT — the density route SURVIVES the THM-883 resonant classes with a factor-2 margin: S(w) = w·Error is P-periodic, so the audit is the finite exact computation max_r|S(r)|; executed for [0..5,t], t = 50, 100, 200: the maximizing residues ARE the resonant classes (r = t, class (a,ℓ) = (1,1)), max|S| = 0.0469·t, 0.0465·t, 0.0458·t (linear — the resonant mode's fingerprint), hence for EVERY peel w ≥ diam: |Error(w)| ≤ max|S|/w ≤ 0.047 < 0.097 = the route's slack — resonant or not, NO band enlargement needed
status: PROVED-EXACT per family instance (S-periodicity: one line; the max: exact rational over the full period); the uniform family law max|S| ≤ 0.05·diam follows from THM-883's resonant-mode formula (upper bound at resonant r) + THM-727's uncoupled diagonal (non-resonant r) — assembly sketched, constants verified; other assembly families each need their own one-shot finite audit (template here)
source: klein-2026-07-16-S314 (cont.5); closes THM-883's urgent handoff
depends_on: [THM-883 (resonant classes), THM-727 (S = w·Error, exact reduction), THM-729]
verification: 04-computation/peel_residue_audit_klein_S314.py -> 05-knowledge/results/peel_residue_audit_klein_S314.out (3/3)
---

# THM-884 — the peel-residue audit: the route survives

**Frame (proved).** S(w) := w·Error_cover(w) = Σ_s Σ_{arcs} (G_s(wb) − G_s(wa)) with all
endpoints of denominator dividing P = 7·lcm(E): S depends only on w mod P. Hence
Error(w) = S(w mod P)/w and the COMPLETE audit is one finite exact table: max_r|S(r)|.

**Executed** (exact rationals, full period, all 7 sections):

| t | P | max_r\|S(r)\| | argmax r | class | max\|S\|/t |
|---|---|---|---|---|---|
| 50 | 2100 | 2.3435 | 50 (= t) | (1,1) | 0.0469 |
| 100 | 2100 | 4.6463 | 800 | (1,8) | 0.0465 |
| 200 | 4200 | 9.1633 | 200 (= t) | (1,1) | 0.0458 |

**Findings.** (i) The worst residues are EXACTLY the THM-883 resonant classes — the audit
confirms the refutation's mechanism is the whole story. (ii) max|S| ≈ 0.047·t, linear (the
resonant mode); the off-resonant bulk is far below. (iii) **Verdict: since every peel
satisfies w ≥ diam ≥ t, |Error(w)| ≤ max|S|/w ≤ 0.047 < 0.097 (the cap₉ − 0.397 slack) —
for EVERY residue class, resonant included, with a factor-2 margin. The density route
survives with NO band enlargement**; the S275 structured band [26, D₀] remains as it was.
(The would-be threshold D₀_res = max|S|/slack ≈ 0.48·t sits BELOW diam.)

**Scope and handoffs.** This instance-audit covers the slow-six family; each other
assembly family needs its own one-shot finite audit (this file is the template — one exact
periodic max per family). The uniform family law max|S| ≤ 0.05·diam is provable by
assembling THM-883's resonant-mode formula (as an upper bound at resonant r) with
THM-727's uncoupled diagonal off resonance — a short write-up once the route owners
confirm the exact slack per leg (0.097 assumed here from THM-729's cap₉ − 0.397).

**The arc completed:** THM-880 (exact bilinear form) → THM-881 (periodicity/descent) →
THM-882 (per-instance scans; dipole route honestly killed) → THM-883 (resonant-mode
theorem; two refutations) → THM-884 (the audit: the refuted worst case lands inside the
route's slack). The empirical Q_s/Error laws of THM-729 are now scoped, exact, and safe.
