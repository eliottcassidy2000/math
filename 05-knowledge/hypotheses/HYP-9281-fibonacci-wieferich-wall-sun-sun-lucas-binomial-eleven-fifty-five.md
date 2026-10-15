---
id: HYP-9281
title: "The golden eleven (interpretive): the eleven/fifty-five coincidences of the eleven-square torus, of Paley P_11 and of the biplane read as one structure, the golden ring at the Fibonacci fixed point F_5 = 5 (11 = L_5 = 2F_5 + 1, 55 = F_10 = C(11,2)); exact anchors: C(L_n,2) = F_2n iff n in {2,4,5}; PSL(2,11) = A_5 Aut(P_11) with intersection C_5 (classical); an icosahedral ANALOGY for the torus colouring; Wall-Sun-Sun primes are Wieferich primes to the base phi (CITED; orbit law THM-4656)"
status: >
  Anchors: PROVED (elementary) or FINITE-EXACT (group computations) where marked; several are the eleven-squares paper's own
  results (its Lemma 5.1, Table 2, Theorem 6.4) or classical (Galois; Klein; Luo; Wall; Freyd-Brown), and are credited as such.
  The hypothesis - that the golden fixed point F_5 = 5 is the common source, so the "golden lane" coincidences are one
  structure rather than several - is INTERPRETIVE (ANALOGY). mac-mini's THM-4656 types 55 = F_10 as a COINCIDENCE (Luo);
  no mechanism is proved here either. The icosahedron-torus statement is ANALOGY (it matches only the C_5-orbit pattern).
session: opus-2026-10-14-S36
source: 05-knowledge/results/strong_balance_pisano_golden_eleven_20261014.md
scripts:
  - 04-computation/experiments/eleven_icosahedron_paley_psl211_20261014.py
  - 04-computation/experiments/fibonacci_pisano_wall_sun_sun_20261014.py
related: [THM-4656, THM-4657, THM-4644, THM-4592, THM-135, HYP-9273, HYP-9277]
---

# HYP-9281: the golden eleven

**Sources.** The eleven-squares paper is "The Fibonacci geometry of eleven squares" (PINGYOU LTD, 2026, UNREFEREED). It works with `R_n = Z[φ]/(φ^n − 1)` and its group `G_n = ⟨x + 1, φx⟩`.

## Exact anchors

1. **Norms (the paper's Lemma 5.1 and Table 2; re-derived).** `N(φ^n − 1) = (−1)^n + 1 − L_n`, so `|R_n| = |L_n − 1 − (−1)^n|`.
   - That is `L_n` for odd `n` and `L_n − 2` for even `n`: `1, 1, 4, 5, 11, 16, 29, 45, 76, 121, …`.
   - Small cases: `R_2 = 0`, `R_4 ≅ F_5` (with `G_4 = AGL(1, 5)`, of order 20), and `R_5 ≅ F_11`.
2. **Regularity on pairs (the paper's Theorem 6.4; re-derived).** For `n >= 3`, `|G_n| = n|R_n|`, so `G_n` can be regular on the unordered pairs of `R_n` only if `|R_n| = 2n + 1`.
   - For odd `n` this reads `L_n = 2n + 1`, which holds only at `n = 5`.
   - For even `n` it reads `L_n = 2n + 3`, which never holds.

   So the paper's "exactly when `n = 5`" follows. The earlier formula `|G_n| = nL_n` was wrong for even `n`; for example `|G_4| = 20` and `|G_6| = 96`.
3. **Luo via Lucas (PROVED, elementary).** `C(L_n, 2) = F_(2n)` iff `L_n = 2F_n + 1` iff `F_(n−3) = 1` iff `n ∈ {2, 4, 5}`.
   - This gives `3 = F_4`, `21 = F_8` and `55 = F_10`, three of the four triangular Fibonacci numbers (Luo 1989).
   - They occur at `L_n = 3, 7, 11`, the Paley primes whose tournaments `P_3, P_7, P_11` maximise the number of Hamiltonian paths at `n = 3, 7, 11` (A038375).
   - A038375 is known only for `n <= 11`. Beyond that the Paley tournament need not be the maximiser: THM-135 shows `P_19` is not even the best circulant on `Z_19`. So the match is NUMEROLOGY within the known range.
4. **Fixed points (PROVED).** `F_n = n` iff `n ∈ {0, 1, 5}`. At `n = 5` the conditions `L_n = 2n + 1` and `L_n = 2F_n + 1` coincide.
5. **`PSL(2, 11)` (FINITE-EXACT; classical).**
   - `|PSL(2, 11)| = 660`.
   - `B = Aut(P_11) = {x ↦ sx + t : s ∈ QR_11}` has order 55 and acts regularly on the 55 arcs. It is the paper's group `G_5`.
   - There are two non-conjugate subgroups `A ≅ A_5` containing `x ↦ 4x`. Each is transitive on `P¹(F_11)` with stabilisers `C_5`, and `A ∩ B = C_5`.
   - `PSL(2, 11) = A·B` as a product of subsets (`60·55/5 = 660`). This is not an exact factorisation, since `A ∩ B ≠ 1`; the exact one is `A_5·C_11`. The factorisation is classical: it is Galois' 11-point action, and factorisations with an `A_5` factor are tabulated by Liebeck–Praeger–Saxl (1990).
   - The `A`-invariant valency-5 graph on the 12 points is the icosahedron. `∞` is antipodal to `0`, and the pentagons are `QR_11 = {1, 3, 4, 5, 9}` and its complement.
   - The 12 vertices carry 15 golden rectangles, one for each two of the 6 antipodal pairs. They form 5 orthogonal triples.
6. **Wall–Sun–Sun as Wieferich to the base `φ` (CITED; orbit law PROVED in THM-4656).**
   - For split `p`, Wall–Sun–Sun ⟺ `φ^(p−1) ≡ 1 (mod p²)` (Bragman–Rowland 2022, Theorem 10; "Fibonacci–Wieferich", McIntosh–Roettger 2007).
   - The circulant tournaments on `Z_(p²)` invariant under the odd part of `φ` explode exactly there (THM-4656 part 1; independently THM-4657).
   - At golden Paley primes (`11, 19, 31, …`) the invariant circulants on `Z_p` are `P_p` and its reverse (THM-4656 part 4).
7. **Pisano periods (CITED and FINITE-EXACT).**
   - `π(n) <= 6n`, with equality iff `n = 2·5^k`, `k >= 1` (Freyd–Brown 1992; K. S. Brown; OEIS A001175). Note `π(2) = 3 ≠ 12`.
   - Split primes have `π(p) | p − 1`, and the owner's "one below" is this law. Inert primes `p >= 3` have `π(p) | 2(p + 1)` but not `p + 1`. The exception is `p = 2`, where `π(2) = 3` divides `p + 1`. Finally `π(5) = 20`.

## The icosahedral picture (ANALOGY)

- **What matches.** The paper colours `H = {1, 3, 4, 5, 9}` gold, `−H` teal and `0` grey. This is the `C_5`-orbit pattern `{0}, H, −H` of `Z_11` under `x ↦ 4x`. It is also the pattern of `P¹(F_11) − {∞}` seen from the icosahedron's vertex `∞`: the south pole `0` and the two pentagons.
- **What does not match.**
  - The paper's "mark" is a boundary point of the flat torus: the canonical lift `z ∈ Z_w` that encodes the side `T` (its Proposition 3.4 and Theorem 3.5). It is not a twelfth cell.
  - The torus cell-adjacency graph is the 6-regular circulant `C_11(1, 3, 4)`. The icosahedron minus a vertex has degrees 4 and 5.
  - `A_5` does not act on the 11 cells.

So "the eleven-square torus is the icosahedron with a marked vertex" holds only as an analogy, at the level of the orbit pattern.

## Hypothesis (interpretive)

The following "golden lane" items are one structure seen through different lenses:
- the eleven-square torus;
- `P_11` maximising Hamiltonian paths among 11-vertex tournaments (A038375);
- the 11-point biplane;
- `π(11) = 10`;
- `|Aut P_11| = 55 = F_10`;
- `11 = L_5`;
- the icosahedral `A_5 ⊂ PSL(2, 11)`.

That structure would be the golden ring `Z[φ]` at its Fibonacci fixed point `n = 5`, where `φ^5 ≡ 1` in a field of `L_5` elements.

**Status of the interpretation.** No mechanism links these items beyond the exact anchors. mac-mini's THM-4656 types `55 = F_10 = F_(π(11))` as a coincidence (by Luo's theorem, 55 is the largest triangular Fibonacci number). This file agrees that nothing beyond the anchors is proved.

**Tests.** The Lucas-binomial cases `n = 2, 4` do not carry a torus structure.
- The identity `C(L_4, 2) = 21 = F_8` concerns `L_4 = 7`, but `R_4 ≅ F_5` has 5 elements, and `G_4 = AGL(1, 5)` is not regular on its 10 pairs.
- `R_2 = 0`.
- So the earlier "`n = 4` = Fano plane / `P_7`" prediction confused `L_n` with `|R_n|`. Any `n = 4` analogue of the torus must come from a different ring (for example `Z[φ]/(7)`, of order 49, since 7 is inert). That is an open item, not a prediction.
