---
id: HYP-9162
title: "The Sierpinski tournament tower: H_2 = [[1,1],[-1,1]], H_(2n) = [[H,H],[-H^T,H^T]] gives skew Hadamard matrices of every order 2^k and hence doubly regular tournaments T_k of every Mersenne order 2^k - 1 (the zero-triangle sides of Gilbreath's single-seed sea); T_3 is the Paley heptagon, Aut(T_k) contains its Frobenius group F_21 for all k >= 3, and the conjecture is that Aut(T_k) = F_21 exactly for all k >= 3 (so T_k would never be Paley for k >= 4; T_5 != P_31 and T_7 != P_127 are established outright by the automorphism orders)"
status: >
  PROVED: the doubling preserves skewness and orthogonality (two-line check),
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

**Update 2026-10-06 (mac-mini, [mod 18/19/7/63 note](../results/mod18_mod19_seven_sixtythree_fractal_20261006.md), section 5).**

* **Arc rule (PROVED, checked `k <= 9`).** `H_(2^k)(x,y) = (-1)^q(x,y)`, where `q(x,y) = sum over bits l with x_l = 1 of (1 + y_l + [x, y differ below bit l])` over `F_2`.
* **Orbit structure (FINITE-EXACT `k <= 8`).**
  * `Aut(T_k)` has `2^(k-3)` orbits of size 7 and `2^(k-3) - 1` fixed points (H-index `= 0 mod 8`).
  * The fixed points induce exactly `T_(k-3)`.
  * At 63 vertices: `63 = 8·7 + 7`, with a fixed Paley heptagon.
* **Reduction (PROVED).** The conjecture follows from one lemma: *no vertex of the second copy `T_k'` has a doubly regular out-neighbourhood* (FINITE-EXACT `k <= 9`). The argument:
  * The vertices with doubly regular out-neighbourhood are exactly the base heptagon plus the apex chain (`k <= 9`).
  * The top apex dominates all of them, so it is the unique source of that `Aut`-invariant set and is fixed by every automorphism.
  * The audit's stabiliser result and induction then finish the proof.
* **The lemma's failures.** They occur only at pairs `(x, y')` with `x ∈ N^+(i)` and `y ∈ N^-(i)`. Their common-out-neighbour count has exact average `λ`, so the lemma is a variance statement about triple intersections.
* **Status.** OPEN.
