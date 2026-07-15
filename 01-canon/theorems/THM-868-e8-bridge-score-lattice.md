---
id: THM-868
title: THE E8 BRIDGE — at n = 8 the score-deviation half-vectors d/2 (d_v = 2s_v − 7) lie EXACTLY in the Σ=0 slice of E8's half-integer coset (D8⁺): membership is immediate (all coordinates in Z+½, sum 0 ∈ 2Z) and verified over the full tournament box (1,012,664 vectors); the near-regular shell x = 8 is precisely the 70 = C(8,4) roots of E8 with Σ = 0 (of its 240); the tie-split climb (THM-867) is an E8 SHELL WALK (+8 in x = +2 in E8 norm) visiting every shell; the Landau filter is the SHELL TRUNCATION with a 5-level corona (trivial ≤ 128; bites 10/12, 10/14, 11/15, 6/12, 1/10 at x = 136..168; empty beyond the transitive ceiling 168). Plus the Feit–Thompson face (all tournament |Aut| odd ⟹ solvable — the quintic's non-solvable symmetry CANNOT occur as tournament symmetry; A5 reaches tournaments only through the lattice: E8 = McKay(2I), 2I = SL(2,5)) and the climb-increment law climb(n) − climb(n−1) = T_⌊(n−1)/2⌋ (paired triangulars)
status: PROVED (the bridge membership: one line; the shell/corona census: exact; the climb law: exact from THM-867's formula) + REFEREED (e8_bridge_referee_opus_S318.out: full box, 40 shells, Feit–Thompson check over all 528 classes n ≤ 7)
source: opus-2026-07-15-S318 (owner: relate the +8 climb to the Cayley–Dickson tower and the geometry/topology of quintic unsolvability)
depends_on:
  - THM-867 (the tie-split walk this geometrizes)
  - HYP-6935/S316 (the residue laws — now read as discriminant-form arithmetic)
related: [the CD tower canon (Mode-B rungs n = 2,3,5,9,17), cd-tower-architecture (kps), THM-855 F3, the truncation grammar (S315/S316/S317)]
verification: 05-knowledge/results/e8_bridge_referee_opus_S318.out
---

# THM-868 — the E8 bridge

## 1. The bridge (n = 8)

For an 8-vertex tournament with scores s_v, put d_v = 2s_v − 7 (all ODD,
Σd = 0). Then v := d/2 has all coordinates in Z + ½ and Σv = 0 ∈ 2Z, i.e.

> **v ∈ E8** (the D8⁺ model: Z⁸-vectors with even sum ∪ (Z+½)⁸-vectors with
> even sum) — every 8-tournament's score deviation is (twice) an E8 vector
> in the Σ=0 slice of the half-integer coset. ∎

Verified over the whole tournament box (|d| ≤ 7): 1,012,664 vectors, all in
E8. The axis is x = Σd² = 4·|v|²_{E8}: **the metagraph's +8 level step is
+2 in E8 norm — the minimal possible step in an even lattice.**

## 2. The root shell

The floor x = 8 (near-regular scores 3,3,3,3,4,4,4,4) is |v|² = 2: the ROOT
shell. The Σ=0 half-coset roots of E8 are the (±½)⁸ vectors with four minus
signs: exactly **C(8,4) = 70 of E8's 240 roots**. Near-regular 8-tournaments
sit on 70 roots of E8; the other 170 roots (the integer-coset ones and the
unbalanced half-coset ones) are not score-realizable.

## 3. The Landau corona (the tournament shadow of the shells)

Per shell (sorted orbits, exact): the Landau filter is **TRIVIAL for
x ≤ 128** (every lattice orbit in the box is a score sequence), **bites in a
5-level corona** x = 136 (10/12), 144 (10/14), 152 (11/15), 160 (6/12),
168 (1/10 — the transitive sequence alone survives at the ceiling), and is
**empty for x > 168** (24 further shells up to 392 are pure lattice). The
tournament world is the E8-slice sharply truncated at the transitive ceiling
(n³−n)/3 with a thin combinatorial corona — the truncation grammar's fourth
instance (after Farey-14/golden, the OCF 2-adic tower, and
polygonal/polyhedral).

## 4. The climb through the tower

THM-867's tie-split walk visits every shell from the root shell to the
ceiling: an E8 shell walk in minimal steps. Across n, the climb length
((n³−n)/24 − [n even]·n/8) satisfies

> **climb(n) − climb(n−1) = T_⌊(n−1)/2⌋** — the triangulars, each twice
> (1, 3, 3, 6, 6, 10, 10, 15, 15, …; verified n = 3..12).

The Cayley–Dickson tower touches this at two off-by-one ladders: the CD
dimensions 1, 2, 4, 8 (R, C, H, O) give the +8 step its octonionic rank;
the Mode-B rungs n = 2, 3, 5, 9, 17 = (CD dim) + 1 are the VERTEX ladder
(the stationary-runner shift). At n = 9 the deviations are even, d/2 ∈ Z⁹
with Σ = 0 — the **A8 slice**; E8's glue vectors (the ±⅓-classes of
A8* / A8) are invisible to integer scores: **the octonions meet tournaments
at n = 8 through the half-coset trick, one rung before Mode B's n = 9, and
the exceptional structure degrades to the A-series exactly at the tower
rung.**

## 5. The Feit–Thompson face (why the quintic cannot enter by symmetry)

Every tournament automorphism has odd order (an even-order automorphism
would contain an involution swapping some pair u, v — impossible, since it
would have to reverse the u→v arc). Hence |Aut(T)| is odd and, **by
Feit–Thompson, every tournament symmetry group is SOLVABLE** (verified
explicitly: all 528 classes at n ≤ 7 have odd |Aut|). The quintic's
obstruction — the non-solvable A5, topologically the icosahedral monodromy,
the perfect group with H₁ = 0 whose double cover 2I = SL(2,5) bounds the
Poincaré sphere S³/2I with its E8 plumbing — therefore CANNOT act as the
symmetry of any tournament. It reaches the tournament world through the
LATTICE instead: E8 = the McKay correspondent of 2I. **The quintic's
geometry enters at n = 8 through the score geometry (§1), not through
symmetry — the same E8 seen from its other face.** The mod-8 arithmetic
powering our residue laws (odd² ≡ 1 mod 8) is the Milgram/Gauss-sum
periodicity whose minimal even unimodular witness is E8 — Bott's 8, the
octonions' 8, and the metagraph's 8 are one 8.

## Named next steps

- Formalize the residue laws as discriminant-form (Milgram) statements for
  the score lattice at general n (the A_{n−1}-dual quotient's Gauss sum).
- The sedenion rung n = 17: zero divisors appear in the algebra — what
  degeneration appears in the metagraph at the corresponding rung? (open)
- A5 as MONODROMY of class families (not symmetry) — e.g. the permutation
  action on the 5-element fibers/lines of suitable n = 8 structures. (open)
- The corona's exact combinatorics at general even n (widths, survivor
  counts — corona(8) = 5 levels; is the width ⌊n/2⌋+1?).
