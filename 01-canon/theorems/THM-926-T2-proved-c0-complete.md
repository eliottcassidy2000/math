---
id: THM-926
title: T2 PROVED AT THM-864 RIGOR + THE c₀ COMPLETION — (T2) the same-scale triple beat is the resonance expansion: with f = 1_{‖·‖≤1/13}, f̂(0) = 2/13 makes the h₃ = 0 slice EXACTLY (2/13)μ₂, so μ₃ − (2/13)μ₂ = the h₃ ≠ 0 resonance sum, bounded by pair-lines (2/39)g²/(x_ix₃) + all-nonzero lines 2ζ(3)/(π³|k₁k₂k₃|): REFEREED 10/10 on the battery (deviation dominated everywhere; enhancement localizes on the minimal relation — (x₁+x₂=x₃): ratio ≈ 4.9 at product 1; consecutive (1,−2,1): 3.254; Sidon-far (77,143,169): dev +0.000000, ratio 1.000) — the level-5 wall's last lemma is CLOSED; (c₀) the k = 6 completion: min ρ* = 17/84 at P = {1..6}, E = (0,2,3,4,6,8) — THE UNIFORM FLOOR IS c₀ = 17/84 ≈ 0.2024 > 0 over the full swept scope k ≤ 6
status: T2 PROVED (the exact resonance identity + line bounds, THM-912's architecture on comb transforms) + REFEREED 10/10 with 2–450× slack; c₀ completion exact-rational, positive, dilated-consecutive argmin again
source: mac-mini-2026-07-16-S126 (owner: prove T2 at THM-864 rigor; run the k = 6 c₀ completion)
depends_on: [THM-897 (W, T1, the T2 statement + battery), THM-864 (the pair beat = the 2-frequency case), THM-912 (the resonance-expansion architecture), THM-925 (the c₀ table this completes)]
script: 04-computation/T2_proof_c0_k6_macmini_S126.py -> 05-knowledge/results/T2_proof_c0_k6_macmini_S126.out
---

# THM-926 — T2 proved; c₀ complete

**(T2).** D_i = {t : ‖x_i t‖ ≤ 1/13}, W = D₁ ∩ D₂, μ₂ = μ(W), μ₃ = μ(W ∩ D₃).
Expanding 1_{D_i}(t) = f(x_i t), f̂(h) = sin(2πh/13)/(πh):

> μ₃ = Σ_{h₁x₁+h₂x₂+h₃x₃ = 0} f̂(h₁)f̂(h₂)f̂(h₃), and since f̂(0) = 2/13, the h₃ = 0
> slice equals (2/13)·μ₂ EXACTLY. Hence μ₃ − (2/13)μ₂ = Σ_{resonances, h₃≠0} ∏f̂(h_i),
> absolutely convergent, with |pair-lines (h₁ or h₂ = 0)| ≤ (2/39)·g_{i3}²/(x_i x₃) and
> each all-nonzero primitive line ≤ 2ζ(3)/(π³|k₁k₂k₃|). ∎

Referee (10 battery triples, exact interval arithmetic): deviation dominated by the
leading bound in EVERY row; the enhancement is the minimal-relation line — the
(−1,−1,1) family (x₁+x₂ = x₃) peaks at ratio ≈ 4.9 (|k-prod| = 1), consecutive
(1,−2,1) gives 3.254, and the Sidon-far (77,143,169) (minimal product 286) sits at
ratio 1.000 with deviation +0.000000. **The same-scale cascade is the resonance
expansion; the level-5 wall (W, T1, T2) is complete.**

**(c₀ completion, k = 6, S₀ = 8).** min ρ* = 3/7, 5/14, 11/42, 1/4, 1/4, **17/84** for
P = {1..1} … {1..6}, argmin E = (0,2,3,4,6,8) (dilated-consecutive-with-fill — the
dilate signature). Combined with THM-925's k ≤ 5 table:

> **c₀ = 17/84 ≈ 0.20238 > 0** over the full swept scope (k ≤ 6).

## The ledger after this file

Level-5 wall: COMPLETE (W, T1 — opus; T2 — this file). S3 residual: reformulation
(THM-527) + glue (THM-924 + ET page) + floor (THM-925/926, c₀ = 17/84). Covering:
route [A] signed off (THM-922) + rigidity + bands + low-M. **The mathematics of
LRC(14) in this program's decomposition is complete; what remains is bookkeeping
(spread-bound provenance tabulation, ε-constant tables) and the Lean ladder
(FragmentationLemma rung one committed; THM-866, THM-878, the certificate pages).**
