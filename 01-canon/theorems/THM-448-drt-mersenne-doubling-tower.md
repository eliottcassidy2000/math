# THM-448 — The DRT/Mersenne Doubling Tower: Paley at 7, then a NEW non-Paley DRT family with frozen automorphisms F_21

**Status:** PART (a) PROVED; parts (b)–(e) VERIFIED computationally (orders ≤ 64) with exact
certificates; conjectures flagged. Adversarially re-verified (66/66 independent checks,
`verify_B_tower_kps1.out`).
**Source:** kind-pasteur-2026-06-09-S1 (branch B of the doubling fan-out + verifier).
**Related:** THM-447 (the doubling), THM-451 (Hadamard chirality), THM-452 (tiling coset law),
HYP-2333 (now resolved into this), HYP-2351/2352 (successor conjectures, renumbered twice from 2339/2340), THM-067 (Mersenne
vanishing), S18h BIBD trichotomy (links reappear in T15's vertex links).

## The tower

Seed S_1 = [1]; iterate S → [[S, S], [S−2I, 2I−S]] (THM-447). Write W_{2^k} for the order-2^k
member and T_{2^k−1} for the core tournament of the (automatically) normalized matrix.

**(a) PROVED: every core is a doubly regular tournament on the Mersenne number 2^k − 1.**
Proof: the tower is skew-Hadamard at every order (THM-447(2), induction from order 1); its first
row is all-ones at every order — analytically, from the closed-form arc law (verified exactly at
orders 2..16, provable by induction):
```
S[i,j] = (−1)^( popcount((i AND j) >> (v+1)) + bit_v(i) ),   v = lsb(i XOR j)
```
("skew-Walsh function": the Sylvester character ⟨i,j⟩ truncated above the lowest differing bit,
plus the twin-direction bit). At i=0 the exponent is 0, so row 0 is all +1 — the tower is BORN
normalized, no renormalization needed. The core of a normalized skew-Hadamard matrix of order m
is a DRT on m−1 vertices (classical). ∎

**(b) T_7 ≅ Paley T_7** (explicit isomorphism; H = 189). But the tower then LEAVES Paley:

**(c) T_31 is NOT isomorphic to Paley T_31** — a non-Paley, non-circulant DRT(31). Three
independent proofs:
1. Triple-intersection distributions: T31 has |N⁺(u)∩N⁺(v)∩N⁺(w)| spectrum
   {0:28, 1:84, 2:252, 3:3192, 4:812, 5:84, 6:28, 7:15}; Paley31's support is only {2,3,4}
   ({2:930, 3:2015, 4:1550}).
2. |Aut(T31)| = 21 vs |Aut(Paley T_31)| = 465 (the 465 reproduced exactly by the same searcher).
3. Exhaustive backtracking isomorphism search: NO map (2,266,906 nodes, search completed).
Also: T31 is not circulant — exhaustive 2^15 enumeration shows the ONLY circulant DRT(31)
S-sets are QR and NQR (both = Paley up to op).

**(d) The automorphism tower FREEZES at the Frobenius group F_21:**
|Aut| = 21 with element orders {1:1, 3:14, 7:6} at T_7, T_15, T_31, AND T_63. Doubling adds
vertices but never adds or loses symmetry past level 7. Orbits: blocks of 7 (core labels
[8k..8k+6]); fixed points are exactly core labels ≡ 7 (mod 8), counts 0, 1, 3, 7 at levels
7, 15, 31, 63. (Conjecture HYP-2352: F_21 for all k.)

**(e) Perfect link self-similarity — the tower is its own vertex-link recursion:**
B_0(T_{2m−1}) = T_{m−1} EXACTLY as a labeled submatrix (identity map), verified:
B_0(T63) = T31, B_0(T31) = T15, B_0(T15) = Paley T_7. The out-neighborhood of vertex 0
at each level IS the previous level. (Strengthened by the verifier from "isomorphic" to
"labeled-identical".)

## T_15 — the first 15-vertex DRT analyzed in this project

- 7-regular; M² = J − 15I exactly; spectrum {0 ×1, ±i√15 ×7}; c3 = 140
  (NOTE: c3 = C(15,3) − 15·C(7,2) = 140 is forced for EVERY 7-regular 15-tournament — the
  prompt's guess 105 was wrong).
- **H(T15) = 198,335,025** = 2.4850 × (15!/2^14 = 79,814,109 = E[H]). Five random 15-tournaments
  gave 43.4M..109.8M (mean 81.6M). Paley T_11's ratio is ≈ 2.44 — T15 is the natural
  H-maximizer candidate at n = 15 (HYP-2351, conjectural).
- **NO circulant DRT(15) exists at all** (0 of 128 antisymmetric S-sets of Z_15) — at the
  Mersenne order 15, the tower provides a DRT where Z_15 provides none.
- |Aut(T15)| = 21, orbits [7,7,1]; NOT vertex-transitive, NOT self-converse — T15 has ZERO
  anti-automorphisms of any order (complete search), so its iso class is PURE BLACK (no
  grid-symmetric tiling in any Hamiltonian-path frame), unlike Paley DRTs (which are SC).
- Link-H spectrum: the vertex links of T15 have H ∈ {189, 171} — exactly two of the three rigid
  regular-7-tournament classes of S18h/THM-027 (BIBD/Paley 189 and the α₂=10 class 171). The
  S18h trichotomy reappears inside the tower.
- Link-H is a strictly stronger vertex invariant than (scores, per-vertex c3, per-vertex cycle
  counts): on T31 it splits vertices into 5 classes ({198335025 ×9, 197568027 ×7, 197147697 ×7,
  196908495 ×1, 196179375 ×7}) while all cheap per-vertex invariants are constant and EQUAL to
  Paley's (which has constant link-H 196,922,055 — a fourth non-isomorphism witness).

## Walsh comparison

S_tower and Sylvester H differ in exactly N²/2 entries at every order (N = 2..64) — the naive
"Walsh + low-weight correction" guess is wrong; the exact relation is the recursive correction
tower S = W + Corr, Corr_{2m} = [[C, C], [C−2I, 2I−C]] (see THM-451 for equivalence-class
chirality and the butterfly transform).

## Scripts

`drt_mersenne_tower_kps1.py`, `drt_tower31_paley_kps1.py`, `drt_tower_followup_kps1.py`,
`verify_B_tower_kps1.py` (+ .out files in 05-knowledge/results/).
