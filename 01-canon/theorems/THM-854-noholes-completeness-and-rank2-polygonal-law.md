---
id: THM-854
title: TWO PROOFS AND A LAW — (I) NO-HOLES COMPLETENESS (the F3 exchange walk, PROVED for all n): every residue-allowed E4 level is realized. Lemma: a tournament has two equal-degree vertices iff it is non-transitive (all-distinct degrees force {0,…,n−1} = transitive); flipping the arc between an equal-degree pair gives Δc₃ = d_u − d_v − 1 = −1 EXACTLY (THM-833 atom). Hence every tournament descends by unit c₃-steps to the transitive, hitting EVERY intermediate level; with x = E4 affine in c₃ and the level laws (odd n: x ≡ 0 mod 8; even n: x ≡ n mod 8), the populated levels form the FULL step-8 progression — no holes, all n (the concurrent probe session's n ≤ 8 exact census becomes a theorem). (II) THE RANK-2 POLYGONAL LAW: the polygonal triangle (columns = 1s, naturals, triangulars, squares, pentagonals, …) is the rank-2 affine shadow of Pascal — on diagonal d the column-j entry is (j−1)·T_d + (d+1), AFFINE in j — so the differences from Pascal are Δ_d(j) = C(j+d, d) − (j−1)·T_d − (d+1): d=2 gives the triangulars (1,3,6,10,15), d=3 gives 4, 13, 28, 50, 80, with binomial-basis rows [0,1,1], [0,4,5,1], [0,10,15,6,1], [0,20,35,21,7,1] (leading coefficients = tetrahedral numbers); polygonal row sums = A000127 (Moser) and skew sums = the Fibonacci-analogue 1,1,2,3,5,8,13,21,33,51,76,111,157,218 (verified 14 terms) — polygonal is to Pascal exactly what Moser is to 2^n, a rank truncation. (III) LOCKER PARITY upgraded to mechanism level: c_odd(D_n) = 2, 4, 12, 40, 136 — EVEN at n = 5..9 (so H(D_n) ≡ 1 mod 4 holds via THM-466 at the c_odd level); the spectra ({3:3, 5:1} at n=6, both odd) PROVE any pairing involution must MIX adjacent odd cycle lengths — the divisor-pairing d ↔ x/d candidate is so constrained (proof open)
status: (I) PROVED (lemma refereed exhaustively n=4..6 + the two-line classical argument; explicit descent). (II) PROVED (one-line algebra P_k(m) = (k−2)T_{m−1} + m) + verified 14 terms / 5 diagonals. (III) VERIFIED n=5..9 at mechanism level; the length-mixing involution is the open piece
source: kind-pasteur-2026-07-15-S128 (cont.21; owner brief - locker proof, F3 no-holes, polygonal vs polyhedral diagonals)
depends_on:
  - THM-833   # the flip atom powering (I)
  - THM-466   # H ≡ 1 + 2 c_odd (mod 4) behind (III)
related:
  - THM-853 (the locker tournament), the concurrent level-law probe session (mod-8 level laws), OPEN-Q-108 A..G corners (the 3-triangular frame), A000127/Moser
---

# THM-854 — no-holes, the rank-2 polygonal law, locker parity

**(I) Proof.** If all degrees are distinct they are {0,…,n−1} and T is transitive (peel the
score-(n−1) vertex and induct); so c₃ > 0 forces an equal-degree pair u, v. Their arc, say u→v,
flipped changes c₃ by d_u − d_v − 1 = −1 (THM-833). Iterating, any T walks down to the transitive
touching every integer c₃ in [0, c₃(T)]; from a c₃-maximal start, every level in [0, c₃max] is
realized. E4 is affine in c₃, so within the proved mod-8 residue class the populated E4 levels are
the full ±8-step progression. ∎ (Lemma refereed on all tournaments n = 4, 5, 6.)

**(II) Proof.** P_k(m) = (k−2)·T_{m−1} + m. On diagonal d, column j holds P_{j+1}(d+1) =
(j−1)·T_d + (d+1) — affine in j — while Pascal holds C(j+d, d). Subtract. ∎

**(III)** c_odd(D_n) even, n = 5..9; spectra {3:2,5:0}, {3:3,5:1}, {3:6,5:6,7:0}, {3:10,5:20,7:10},
{3:16,5:51,7:69,9:0}. At n=6 both odd counts are odd, so any parity involution exchanges 3-cycles
with 5-cycles: the divisor-pairing candidate must change cycle length by ±2.

## Evidence log
- [x] all three referees (noholes_locker_polygonal_kps_S128c21.py + .out)
- [ ] (III) the length-mixing involution
