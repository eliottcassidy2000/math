---
id: HYP-9162
title: "The Sierpinski tournament tower: H_2 = [[1,1],[-1,1]], H_(2n) = [[H,H],[-H^T,H^T]] gives skew Hadamard matrices of every order 2^k and hence doubly regular tournaments T_k of every Mersenne order 2^k - 1 (the zero-triangle sides of Gilbreath's single-seed sea); T_3 is the Paley heptagon, Aut(T_k) contains its Frobenius group F_21 for all k >= 3, and the conjecture is that Aut(T_k) = F_21 exactly for all k >= 3 (so T_k would never be Paley for k >= 4; T_5 != P_31 and T_7 != P_127 are established outright by the automorphism orders)"
status: >
  RESOLVED 2026-10-06: PROVED (THM-4557). The key fact is A. Hanaki, arXiv:2011.06141 (2020), Thm 3.4
  (doubling a doubly regular tournament on >= 7 vertices keeps its automorphism group), so the conjecture
  already followed from the literature when posed; THM-4557 gives an independent proof. Aut(T_k) = F_21 for all k >= 3. History: PROVED: the doubling preserves skewness and orthogonality (two-line check),
  double regularity is the classical skew-Hadamard correspondence, the
  recursion T_(k+1) = T_k + {0'} + T_k' is explicit, and every automorphism of
  T_k extends diagonally, so F_21 <= Aut(T_k) for k >= 3. FINITE-EXACT:
  |Aut(T_k)| = 3, 21, 21, 21, 21, 21, 21 for k = 2..8 (orders 3 to 255),
  T_3 = P_7, T_5 not isomorphic to P_31 (|Aut(P_31)| = 465) though both have
  the same 4-vertex census. The independent audit added |Aut(T_9)| = 21 (n = 511) by two methods and
  proved that the stabilizer of 0' in Aut(T_(k+1)) consists exactly of the
  diagonal extensions of Aut(T_k), so the conjecture is equivalent to: 0' is
  fixed by every automorphism of T_(k+1), k >= 3. OPEN: that equality.
source: opus-2026-09-26 session gilbreath-fermat-platonic-20260926
related: [THM-871 (Fermat-rung rigidity of rotational tournaments), HYP-3805 (Paley heptagon as extremal object), S9 note collatz_oscillation_gilbreath_20260926.md (zero triangles of the sea)]
verification: 04-computation/experiments/gilbreath_fermat_tower_20260926_tournaments.py -> .out
---

# HYP-9162 -- the Sierpinski tournament tower keeps exactly the heptagon's 21 symmetries

**Construction.** `H_2 = [[1, 1], [-1, 1]]`, `H_(2n) = [[H, H], [-H^T, H^T]]`.
If `H = S + I` with `S` antisymmetric then `H_(2n) - I = [[S, S+I], [S-I, -S]]`
is antisymmetric and `H_(2n) H_(2n)^T = 2n I`; normalization is inherited.
Deleting row and column `0` gives the tournament `T_k` on `2^k - 1` vertices,
doubly regular (checked to `n = 255`), with the recursion

```text
T_(k+1) = T_k + {0'} + T_k' (arcs reversed),  i -> i',  i -> j' iff i -> j,  i' -> j iff i -> j,  0' -> T_k,  T_k' -> 0'.
```

**Computed.** `|Aut(T_2)| = 3`; `T_3 = P_7`, `|Aut| = 21`; `|Aut(T_k)| = 21` for `k = 4, 5, 6, 7, 8` (`n = 15, 31, 63, 127, 255`), with `Aut(T_5)` having element orders
`{3: 14, 7: 6}` (the Frobenius group `Z_7 x| Z_3`) and orbits of sizes
`7, 1, 7, 1, 7, 1, 7`. `T_5` is not isomorphic to `P_31`, and `T_7` is not `P_127`
(`|Aut(P_127)| = 8001`).

**Conjecture.** `Aut(T_k) = F_21` for every `k >= 3`; equivalently (audit:
the stabilizer of `0'` is exactly the set of diagonal extensions, proved via
the constant number `(n+1)/4` of `3`-cycles through every arc), no
automorphism of `T_(k+1)` moves `0'` once `k >= 3`. Verified for
`k = 3..9` (orders `7` to `511`). If true, `T_k` is
not vertex-transitive for `k >= 4` and never Paley (a Paley tournament on a
Mersenne prime `2^k - 1` has `|Aut| = (2^k - 1)(2^(k-1) - 1) > 21`).

**Why it matters here.** The sides of the zero triangles of Gilbreath's
single-seed sea are exactly the Mersenne numbers `2^m - 1` (session note
`gilbreath_fermat_platonic_20260926.md`, section 2), and this tower is the
doubling of the sea itself read with signs: it makes the owner's "zeros of
tournament size edged by 2s" precise (each such size carries a doubly regular
tournament) and shows what symmetry the tower keeps: the heptagon's.

**Update 2026-10-06 (mac-mini, [mod 18/19/7/63 note](../results/mod18_mod19_seven_sixtythree_fractal_20261006.md), section 5; independently audited, corrections in MISTAKE-569).**

* **Arc rule (PROVED; checked `k <= 9`).** `H_(2^k)(x,y) = (-1)^q(x,y)`, where `q(x,y) = sum over bits l with x_l = 1 of (1 + y_l + [x, y differ below bit l])` over `F_2`.
  * Proof: induction on the top bit. `H` is skew off the diagonal (`H + H^T = 2I`, by induction from the block form).
  * If `x`'s top bit is 0, the entry is `H_n(x', y')`. If it is 1, the entry is `∓H_n(y', x') = ∓(-1)^[x' ≠ y'] H_n(x', y')`, which is the new top-bit term of `q`.
* **Self-similarity (PROVED).** `q(8x, 8y) = q(x, y)`, so for every `k` the multiples of 8 induce exactly `T_(k-3)`.
* **Frobenius action (PROVED).** `q(x,y)` is `q_3(x mod 8, y mod 8)` plus terms that depend only on the high bits and on `[x ≢ y (mod 8)]`. Hence every automorphism of `T_3 = P_7`, acting on the low three bits and fixing the multiples of 8, is an automorphism of every `T_k`. So `F_21 <= Aut(T_k)`, with `2^(k-3)` orbits `{8m+1, ..., 8m+7}` of size 7 and `2^(k-3) - 1` fixed points.
* **Equality `|Aut(T_k)| = 21` (FINITE-EXACT `k <= 8`, nauty).** At 63 vertices this is `63 = 8·7 + 7`: eight heptagon orbits around a fixed Paley heptagon.
* **Reduction (PROVED).** Let `D` be the set of vertices with a doubly regular out-neighbourhood. The conjecture follows from one lemma: *no vertex of the second copy `T_k'` lies in `D`* (FINITE-EXACT `k <= 9`).
  * `N^+(0') = T_k` is doubly regular, so `0' ∈ D`.
  * The lemma gives `D ⊆ T_k ∪ {0'}`.
  * Since `0' -> T_k`, `0'` is the unique source of the `Aut`-invariant set `D`.
  * The stabiliser theorem of this hypothesis's audit (`Stab(0')` = the diagonal extensions of `Aut(T_k)`) and induction from `Aut(P_7) = F_21` finish the proof.
  * Separately, `D` is exactly the base heptagon plus the apex chain for `k <= 9` (FINITE-EXACT). The reduction does not use this.
* **The lemma's failures.** They occur only at pairs `(x, y')` with `x ∈ N^+(i)` and `y ∈ N^-(i)`. Pairs of other types provably have count `λ`, and the failing type's average is exactly `λ` (double counting). So the lemma is a variance statement about triple intersections.
* **Status.** OPEN at the time of this update; RESOLVED later the same day (see below).

**Resolution 2026-10-06 (mac-mini, [THM-4557](../../01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md)).** The second-copy lemma is PROVED, for every doubly regular tournament `T` on `n = 4t+3 >= 7` vertices.
* `N^+(x')` is `T` with the arcs inside `N^−(x)` reversed.
* Double regularity would force `S_R 1_(O_j) = 0` for the skew matrix `S_R` of `R = T[N^−(x)]`, which has odd order.
* `S_R ≡ J − I (mod 2)` has `F_2`-rank `m − 1`, so `ker S_R = span(1)`. But `|O_j| = t + 1`.

Hence `Aut(D(T)) = Aut(T)`, and `Aut(T_k) = F_21` for all `k >= 3`.

**Prior art (found by the independent audit; MISTAKE-571).** The theorem is Theorem 3.4 of A. Hanaki, arXiv:2011.06141 (2020), proved there with triple intersection numbers. The intransitivity is already in Faradžev–Klin–Muzichuk (1994), Theorem 2.6.6. This hypothesis was a corollary of the literature when it was posed.
