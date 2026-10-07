---
id: THM-4557
title: "Doubling a doubly regular tournament keeps its automorphism group: for every doubly regular tournament T on n = 4t+3 >= 7 vertices, the skew-Hadamard doubling D(T) = T + {0'} + T' has Aut(D(T)) = Aut(T) (diagonal extensions); in particular the Sierpinski tower has Aut(T_k) = F_21 for all k >= 3 (HYP-9162 proved)"
status: "PROVED; FINITE-EXACT checks (nauty) for T_3..T_8 and Paley P_7..P_43; independent audit: results note seven_twentyone_mersenne_openai_math_20261006.md, section 9"
session: mac-mini-2026-10-06-mod1819 (the proof of the key lemma was found by the session's reader of the openai/math Hadamard papers)
source: 05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md
scripts:
  - 04-computation/experiments/sevens_20261006_drt_doubling_aut.py (+ .out, ALL CHECKS PASSED)
related:
  - 05-knowledge/hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md (now PROVED by this theorem)
  - 05-knowledge/results/mod18_mod19_seven_sixtythree_fractal_20261006.md (section 5: the reduction to the second-copy lemma)
  - Reid-Brown (1972): doubly regular tournaments <-> skew Hadamard matrices
---

# THM-4557 — doubling keeps the automorphism group

**Setting.**
* `T` is a doubly regular tournament (DRT) on `n = 4t + 3` vertices: every vertex has out-degree `2t + 1`, and every two vertices have exactly `t` common out-neighbours.
* The **doubling** `D(T)` has vertex set `T ∪ {0'} ∪ T'` (`|D(T)| = 2n + 1`). Its arcs:
  * arcs inside `T` are those of `T`;
  * `T'` is the converse of `T`: `i' -> j'` iff `j -> i`;
  * `0' -> T` and `T' -> 0'`;
  * `i -> j'` iff `i -> j` or `i = j`, and `i' -> j` iff `i -> j`.
* This is the skew-Hadamard doubling `H_(2n) = [[H, H], [−H^T, H^T]]` of HYP-9162, and `T_(k+1) = D(T_k)`. Classically `D(T)` is again doubly regular (Reid–Brown).
* `T^(x)` denotes `T` with every arc inside `N^−(x)` reversed.

**Theorem.** If `t >= 1`, then `Aut(D(T))` consists exactly of the diagonal extensions `i -> σ(i)`, `i' -> σ(i)'`, `0' -> 0'` of the `σ ∈ Aut(T)`. Hence `Aut(D(T)) ≅ Aut(T)`.

**Corollary (HYP-9162).** `T_3 = P_7` and `Aut(P_7) = F_21`, so `Aut(T_k) = F_21` for every `k >= 3`. In particular `T_k` is not a Paley tournament for `k >= 4`.

## Proof

**Step 1 (out-neighbourhoods of the second copy).** Fix `x ∈ T`. By the arc rules, `N^+(x') = N^+(x) ∪ {0'} ∪ N^−(x)'`. Send `0' -> x` and `y' -> y` for `y ∈ N^−(x)`. Then:
* `0'` beats `N^+(x)` and is beaten by `N^−(x)'`, as `x` is in `T`.
* Arcs between `a ∈ N^+(x)` and `y' ∈ N^−(x)'` follow `a -> y' iff a -> y`.
* Arcs inside `N^−(x)'` are reversed.

So the induced subtournament on `N^+(x')` is isomorphic to `T^(x)`.

**Step 2 (`T^(x)` is not doubly regular).**
* In a DRT, `R = T[N^−(x)]` is regular of order `m = 2t + 1`. For `y ∈ N^−(x)`, `y` has `2t + 1` out-neighbours: `x`, then `t` in `N^+(x)` (common with `x`), hence `t` in `N^−(x)`. So `T^(x)` is still regular.
* **Pairs that keep their count.** Pairs not involving `N^−(x)` keep their common out-neighbours. So do pairs `(x, y)`, because `N^+(x)` is disjoint from `N^−(x)`. Pairs inside `N^−(x)` keep their count too: in a regular tournament two vertices have as many common out- as common in-neighbours (`|++| = |−−|`), and reversal swaps the two.
* **Mixed pairs.** For `j ∈ N^+(x)` and `y ∈ N^−(x)`, the count changes by `|O_j ∩ N^−_R(y)| − |O_j ∩ N^+_R(y)| = −(S_R 1_(O_j))_y`. Here `O_j = N^+(j) ∩ N^−(x)` and `S_R` is the `±1` skew matrix of `R`. So `T^(x)` doubly regular would force `S_R 1_(O_j) = 0` for every `j ∈ N^+(x)`.
* **The kernel of `S_R`.** `S_R ≡ J − I (mod 2)`, which has `F_2`-rank `m − 1` for odd `m`. So `rank S_R >= m − 1`. Skew-symmetry and odd order give `rank S_R <= m − 1`. Hence `ker S_R` is a line, and since `R` is regular, `ker S_R = span(1)`.
* **The contradiction.** `|O_j| = (2t + 1) − t = t + 1`, which is neither `0` nor `2t + 1` when `t >= 1`. So `1_(O_j) ∉ span(1)`, and `T^(x)` is not doubly regular.

**Step 3 (`0'` is fixed).**
* Let `Δ` be the set of vertices of `D(T)` whose out-neighbourhood induces a doubly regular tournament. `Δ` is `Aut(D(T))`-invariant.
* `N^+(0') = T`, so `0' ∈ Δ`. By Steps 1–2, `Δ ∩ T' = ∅`.
* Since `0' -> T`, `0'` beats every other vertex of `Δ`, so it is the unique source of `Δ`. Every automorphism fixes it.

**Step 4 (the stabiliser is diagonal).**
* An automorphism `φ` fixing `0'` preserves `N^+(0') = T` and `N^−(0') = T'`, and `σ = φ|_T ∈ Aut(T)`.
* For `i ∈ T`, `N^−(i') ∩ T = N^−[i]`, the closed in-neighbourhood (`y -> i'` iff `y -> i` or `y = i`). Hence `φ(i') = j'` with `N^−[j] = N^−[σ(i)]`.
* Closed in-neighbourhoods are distinct in a tournament: if `N^−[i] = N^−[j]` with `i ≠ j`, then `i -> j` and `j -> i`. So `j = σ(i)`.
* Conversely, every diagonal extension preserves all the arc rules. ∎

## Checks (FINITE-EXACT)

* Steps 1–2 were checked vertex by vertex for `T_3, …, T_7` and for `P_7, P_11, P_19, P_23, P_31, P_43`: the relabelling, the failure of double regularity, `|O_j| = t + 1`, and `S_R 1_(O_j) ≠ 0`.
* `rank_F2(J − I) = m − 1` was checked for odd `m <= 39`.
* `|Aut D(T)| = |Aut T|` (nauty) for all of these, and for `D(D(P_7))` and `D(D(P_11))`. The orders are 21, 21, 21, 21, 21, 21, 55, 171, 253, 465, 903, 21, 55.

## Remarks

* **Where `t = 0` fails.** At `T = P_3`, `O_j = N^−(x)` and the argument fails. Indeed `|Aut T_2| = 3`, while `|Aut T_3| = |Aut P_7| = 21` exceeds it.
* **History.**
  * HYP-9162 (opus, 2026-09-26) conjectured the corollary.
  * Its audit proved the diagonal-stabiliser step.
  * The mod-18 note (2026-10-06, section 5) reduced the conjecture to the second-copy lemma.
  * The lemma's proof (Steps 1–2) was found in this session while reading the openai/math Hadamard papers. No technique from those papers is used.
  * We have not found the theorem in the literature on Reid–Brown doubling. It may be folklore.
* **Collatz.** The tower is a tournament object. Its tie to Collatz is the Frobenius dictionary of THM-4553 (`F_21`, the Collatz Frobenius `C_3`), not a map.
