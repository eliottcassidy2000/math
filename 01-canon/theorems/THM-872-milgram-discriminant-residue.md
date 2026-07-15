---
id: THM-872
title: THE RESIDUE LAWS ARE DISCRIMINANT-FORM ARITHMETIC — odd n: d = 2e with e ∈ A_{n−1} (an EVEN root lattice), so 8 | x is literally lattice evenness; even n: d is an odd-coset vector, x ≡ n (mod 8) is the coset norm (odd² ≡ 1 mod 8); the Gauss–Milgram sum of disc(A_{n−1}) = (Z/n, q(k) = (n−1)k²/2n) equals √n·e^{2πi(n−1)/8} (signature n−1 mod 8; verified n ≤ 40) — the metagraph's 8 IS Milgram's 8; and the bridge ladder continues: at n = 16 every score half-vector lies in D₁₆⁺, the OTHER even unimodular 16-dim lattice — Milnor's isospectral partner of E8⊕E8 (which fails in the split frame: odd half-sum witnesses)
status: PROVED (coset-norm arguments: two lines each; Gauss sums: verified numerically n ≤ 40 + classical Milgram; D₁₆⁺ membership: one line as in THM-868, verified on random boxes; split-frame E8⊕E8 failure: explicit witnesses)
source: klein-2026-07-15-S313 (cont.3); answers THM-868 "Named next steps" item 1
depends_on: [THM-866/867/868, classical: Milgram's formula, Milnor 1964 (isospectral tori from E8⊕E8 vs D₁₆⁺)]
verification: 04-computation/icosian_kakeya_milgram_sedenion_klein_S313.py (section A; 24/24)
---

# THM-872 — Milgram/discriminant-form formalization of the residue laws

## 1. Odd n: the law IS lattice evenness

d = 2e with Σe = 0: e ranges over (a box in) the root lattice A_{n−1}, which is
EVEN (Σe² ≡ (Σe)² = 0 mod 2). Hence x = 4Σe² ≡ 0 (mod 8). The mod-8 level law at
odd n is nothing but "A_{n−1} is an even lattice," scaled by 4.

## 2. Even n: the law is the coset norm

All d_v odd: x = Σd² ≡ n (mod 8) since odd² ≡ 1 (mod 8). Equivalently v = d/2 lies
in the half-integer coset (Z+½)ⁿ ∩ {Σ = 0}: writing v_i = a_i + ½,
Σv² = Σ(a² + a) + n/4 ≡ n/4 (mod 2) — the coset's discriminant-form value. Both
laws are instances of one statement: **x mod 8 = 4·(norm of the relevant coset in
the discriminant group), and the 8 is Milgram's 8.**

## 3. The Gauss–Milgram signature of the score lattice

disc(A_{n−1}) = A*/A ≅ Z/n with q(k) = (n−1)k²/(2n) mod 1. Milgram:
Σ_k e^{2πi q(k)} = √n · e^{2πi σ/8}, σ = n − 1 (mod 8) — verified numerically
n ≤ 40. The signature ladder σ(n) = n−1 mod 8 is the period-8 clock that the CD
rungs n = 2^k + 1 ride (σ = 2^k mod 8 = 2, 4, 0, 0, … from rung H onward: every
rung from n = 9 has σ ≡ 0 — the octonion rung and all higher rungs sit at the
UNIMODULAR-compatible signature, which is why an E8-type bridge first exists at
n = 8/9 and not before).

## 4. The bridge ladder: n = 16 is Milnor's other drum

At n = 16, v = d/2 has all coordinates in Z+½ and Σv = 0 ∈ 2Z: **v ∈ D₁₆⁺** — the
even unimodular lattice that is E8⊕E8's isospectral partner (Milnor's famous pair).
In the split frame E8⊕E8 requires each 8-block half-sum even, which fails for
explicit tournaments (odd half-sum witnesses found): coordinate-wise, **the
16-vertex score world lives in D₁₆⁺, not E8⊕E8.** The metagraph hears one of the
two drums. Next rung of the ladder: n = 24 (Leech/Niemeier territory — OPEN which
of the 24 Niemeier lattices, if any, is score-selected; the same one-line membership
gives D₂₄⁺-coset containment, not unimodularity).
