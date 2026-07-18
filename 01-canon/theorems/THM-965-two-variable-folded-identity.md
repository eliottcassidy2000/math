---
id: THM-965
title: THE TWO-VARIABLE FOLDED IDENTITY — μ(D_a ∩ D_b) = [4ab + fold_{14g}((a+b) mod 14g) − fold_{14g}((b−a) mod 14g)]/(196ab), fold_M(r) = r(M−r), g = gcd — the pair-overlap measure in closed form for ALL pairs (verified 400/400 exact, every gcd); subsumes boxeph's consecutive form (the b = a+1 slice); makes THM-964's floor table ANALYTIC (μ ≥ 1/49 − g²·49/(196ab); the S341 table minima re-derived from the formula); the fold r(M−r) is the staircase triangle's leg product — the project's original geometry closes its last estimate
status: PROVED + KERNEL-PURE IN LEAN (opus-S346): muNum_folded in TournamentH7.LRCFoldedIdentity builds green, axioms [propext, Classical.choice, Quot.sound], no sorry/native_decide — 14·muNum a b = 4ab + fold((a+b)%14) − fold((b−a)%14) over boxeph's muNum (g=1 form; the gcd rescale is scalar). Verified 400/400 exact (all gcds, ratios ≤ 13); table re-derivation exact (1/91, 1/77, 1/77, 1/63, 1/63 at the same minimizers); proof = the min(2a,X) = X − (X−2a)⁺ telescope of the sawtooth sum into F(a+b) − F(b−a) with 14g·F(S) = S² + fold(S mod 14g), the tail by Nat.le_induction, one linear_combination (coefficients sympy-verified residual 0)
source: opus-2026-07-17-S343 (owner: finish the remaining LRC(14) proof pieces)
depends_on: [boxeph LRCTreeHunter (muNum + the consecutive slice), THM-856 (the sawtooth), THM-964 (the floor consumer)]
scripts: 04-computation/folded_identity_general_opus_S343.py -> 05-knowledge/results/folded_identity_general_opus_S343.out
---

# THM-965 — the two-variable folded identity

**Identity.** For positive integers a ≤ b, g = gcd(a,b), M = 14g,
fold_M(r) = r·(M − r):

> **μ(D_a ∩ D_b) = [ 4ab + fold_M((a+b) mod M) − fold_M((b−a) mod M) ] / (196ab).**

*Proof sketch (three lines).* In the sawtooth sum, the cap min(2a, X)
equals X − (X − 2a)⁺ (the 2b cap never binds for a ≤ b), so the sum
telescopes to F(a+b) − F(b−a) with F(S) = Σ_m (S − M|m|)⁺; the closed
form 14g·F(S) = S² + fold_M(S mod M) is the one-variable fold. ∎
(Verified 400/400 exact over all gcds; the S341 floor-table minima —
1/91 at (2,26), 1/77, 1/63 — re-derive from the formula at the same
minimizers.)

**Consequences.** (1) The floor table is ANALYTIC: μ = 1/49 +
[fold₊ − fold₋]/(196ab) ≥ 1/49 − g²·49/(196ab) — no decide sweep; the
Lean floor-table task reduces to this identity over boxeph's muNum plus
per-class minimization of an explicit rational function. (2) boxeph's
consecutive_closed_form is the b = a+1 slice (fold₋ ≡ 13 constant).
(3) The deviation of any pair overlap from independence (1/49) is TWO
FOLDS — the beat coordinate (b−a) and the sum coordinate (a+b), each
folded at modulus 14g: THM-928(B)'s two Bézout edges, in measure form.

**The geometry.** fold_M(r) = r(M − r) is the staircase triangle's leg
product k(n−1−k) — the object the project began with. The original
tournament-parity geometry supplies the closed form that finishes the
last wall's floor: one reflection, first to last.
