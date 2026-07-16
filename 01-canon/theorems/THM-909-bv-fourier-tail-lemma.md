---
id: THM-909  # renumbered from 908 (concurrent seven-shift sieve first-pushed)
title: THE BV-FOURIER TAIL LEMMA — structure PROVED: the triple channel factorizes as Ĝ(n) = ∏ᵢ[sin(πnᵢ/7)/(πnᵢ)]·W(n mod 14) with W ANTI-PERIODIC mod 7 (the sign that hides in mod-7 indexing), W ≡ 0 on zero residues (zero-marginality, 9e-15), and T(a,b,c) = Σ_{n·(a,b,c)=0} Ĝ(n) absolutely convergent — VERIFIED: the direct lattice sum reproduces the measured T(1,2,7) = +0.1007 to 0.7% at |n| ≤ 260; certified line-sum table: the large lines are NEGATIVE (S(1,1,1) = −0.113), max positive line +0.0587; single lines do NOT explain T(1,2,7) (+0.012 from its three smallest relations) — the mass is a MULTI-LINE conspiracy, which requires TWO small relations, forcing (a,b,c) ∝ k×k′: the residue-6 closure = a finite cross-product box sweep via this (fast, verified) expansion
status: factorization/anti-periodicity/zero-marginality/expansion PROVED and machine-verified; line sums certified (truncation ≤ 1.1e-6); the assembled single-line estimate 0.4435 < 0.47; REMAINING: the λ₂ off-line remainder bound + the k×k′ box sweep (the expansion computes any T in seconds — next session's mechanical close)
source: mac-mini-2026-07-16-S119 (owner: prove the BV-Fourier tail lemma and finish closing residue 6)
depends_on: [THM-907 (channel data), codex THM-904 (target (3)), THM-903 (reflection frame)]
script: 04-computation/bv_fourier_tail_lemma_macmini_S119.py -> 05-knowledge/results/bv_fourier_tail_lemma_macmini_S119.out
---

# THM-909 — the BV-Fourier tail lemma

**(1) Factorization (proved).** Î_s(n) = e(−n(2s+1)/14)·sin(πn/7)/(πn), so the triple
channel has Ĝ(n) = ∏ᵢ[sin(πnᵢ/7)/(πnᵢ)]·W(n₁,n₂,n₃) with W the phase-DFT of β₁₂₃.
W(n+7eᵢ) = −W(n): **anti-periodic mod 7** (index mod 14 — the sign a mod-7 table silently
destroys; caught in-session). Zero-marginality: W ≡ 0 whenever any residue is 0
(machine: 9·10⁻¹⁵); 7|nᵢ kills the sine. ‖W‖_∞ = 37.75.

**(2) Resonance expansion (proved + verified).** T(a,b,c) = Σ_{n∈Λ*} Ĝ(n), Λ* the
annihilator lattice — absolutely convergent (cubic decay along lines). Direct check at
(1,2,7): lattice sum +0.099967 vs measured +0.100714 (0.7%, |n| ≤ 260 truncation).

**(3) Line table (certified ≤ 1.1e-6).** The dominant lines are NEGATIVE: S(1,1,1) =
−0.1130, S(1,±1,∓1) = −0.0922, S(2,2,±1) ≈ −0.067. Max POSITIVE line: +0.0587
((3,3,−2)-orbit). Hence a single-relation family (large triples with one small relation,
e.g. (a, b, a+b)) satisfies q ≤ β₀ + pairs + 0.0587 + o(1) < 0.47 comfortably.

**(4) The conspiracy structure and the finite box.** T(1,2,7)'s +0.1007 is NOT one line
(its three smallest relations sum to +0.012): it is a multi-line conspiracy. Two
independent small relations k, k′ determine the triple up to scale: (a,b,c) ∝ k × k′.
So any triple with T above the single-line ceiling lies in an EXPLICIT finite
cross-product box, and the (2)-expansion computes each such T in seconds. **The closure
of negative residue 6 = [λ₂ off-line remainder bound (one page, constants: ‖W‖_∞/π³ =
1.218)] + [the k×k′ box sweep via (2)] — both mechanical; the single-line estimate
already assembles to 0.4435 < 0.47.**
