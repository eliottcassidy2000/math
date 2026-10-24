---
id: THM-4557
title: "Doubling a doubly regular tournament keeps its automorphism group: for every doubly regular tournament T on n = 4t+3 >= 7 vertices, the skew-Hadamard doubling D(T) = T + {0'} + T' has Aut(D(T)) = Aut(T) (diagonal extensions); in particular the Sierpinski tower has Aut(T_k) = F_21 for all k >= 3 (HYP-9162 proved; the theorem is known: Hanaki 2020, Thm 3.4; new proof via the odd skew kernel)"
status: "PROVED (KNOWN: Hanaki 2020, arXiv:2011.06141, Thm 3.4; intransitivity earlier, Faradzev-Klin-Muzichuk 1994, Thm 2.6.6 per Hanaki; the proof here is new). FINITE-EXACT checks (nauty) for T_3..T_8 and Paley P_p, p = 7, 11, 19, 23, 31, 43; INDEPENDENTLY AUDITED (2026-10-06; corrections in MISTAKE-571; record in section 9 of the results note)"
session: mac-mini-2026-10-06-mod1819 (the proof of the key lemma was found by the session's reader of the openai/math Hadamard papers)
source: 05-knowledge/results/seven_twentyone_mersenne_openai_math_20261006.md
scripts:
  - 04-computation/experiments/sevens_20261006_drt_doubling_aut.py (+ .out, ALL CHECKS PASSED)
related:
  - 05-knowledge/hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md (now PROVED by this theorem)
  - 05-knowledge/results/mod18_mod19_seven_sixtythree_fractal_20261006.md (section 5: the reduction to the second-copy lemma)
  - Reid-Brown (1972): doubly regular tournaments <-> skew Hadamard matrices
  - A. Hanaki, Non-symmetric class 2 association schemes obtained by doubling of skew-Hadamard matrices are non-schurian, arXiv:2011.06141 (2020), Theorem 3.4 (this theorem, different proof)
  - I. A. Faradzev, M. H. Klin, M. E. Muzichuk, Cellular rings and groups of automorphisms of graphs (1994), Theorem 2.6.6 (intransitivity; per Hanaki's footnote)
---

# THM-4557 — doubling keeps the automorphism group

**Setting.**
* `T` is a doubly regular tournament (DRT) on `n = 4t + 3` vertices: every vertex has out-degree `2t + 1`, and every two vertices have exactly `t` common out-neighbours.
* The **doubling** `D(T)` has vertex set `T ∪ {0'} ∪ T'` (`|D(T)| = 2n + 1`). Its arcs:
  * arcs inside `T` are those of `T`;
  * `T'` is the converse of `T`: `i' -> j'` iff `j -> i`;
  * `0' -> T` and `T' -> 0'`;
  * `i -> j'` iff `i -> j` or `i = j`, and `i' -> j` iff `i -> j`.
* This is the skew-Hadamard doubling `H ↦ [[H, H], [−H^T, H^T]]` (matrix order `n+1 ↦ 2n+2`) of HYP-9162. `T_(k+1) = D(T_k)` exactly, under `i ↦ i`, `0' ↦ 2^k`, `i' ↦ 2^k + i`.
* `D(T)` is again doubly regular: the doubling is again skew Hadamard (two-line check, HYP-9162), and Reid–Brown (1972) identifies normalized skew Hadamard matrices with doubly regular tournaments.
* `T^(x)` denotes `T` with every arc inside `N^−(x)` reversed.

**Theorem (Hanaki 2020, Theorem 3.4; new proof).** If `t >= 1`, then `Aut(D(T))` consists exactly of the diagonal extensions `i ↦ σ(i)`, `i' ↦ σ(i)'`, `0' ↦ 0'` of the `σ ∈ Aut(T)`. Hence `Aut(D(T)) ≅ Aut(T)`.

**Corollary (HYP-9162).** `T_3 = P_7` and `Aut(P_7) = F_21`, so `Aut(T_k) = F_21` for every `k >= 3`. In particular `T_k` is not a Paley tournament for `k >= 4`: a Paley tournament on `q >= 15` vertices has `|Aut| >= q(q−1)/2 > 21`.

## Proof

**Step 1 (out-neighbourhoods of the second copy).** Fix `x ∈ T`. By the arc rules, `N^+(x') = N^+(x) ∪ {0'} ∪ N^−(x)'`. Send `0' ↦ x` and `y' ↦ y` for `y ∈ N^−(x)`. Then:
* `0'` beats `N^+(x)` and is beaten by `N^−(x)'`, as `x` is in `T`.
* Arcs between `a ∈ N^+(x)` and `y' ∈ N^−(x)'` follow `a -> y' iff a -> y`.
* Arcs inside `N^−(x)'` are reversed.

So the induced subtournament on `N^+(x')` is isomorphic to `T^(x)`.

**Step 2 (`T^(x)` is not doubly regular).**
* In a DRT, `R = T[N^−(x)]` is regular of order `m = 2t + 1`. For `y ∈ N^−(x)`, `y` has `2t + 1` out-neighbours: `x`, then `t` in `N^+(x)` (common with `x`), hence `t` in `N^−(x)`. So `T^(x)` is still regular.
* **Pairs that keep their count.** Pairs not involving `N^−(x)` keep their common out-neighbours. So do pairs `(x, y)`, because `N^+(x)` is disjoint from `N^−(x)`. Pairs inside `N^−(x)` keep their count too: in a regular tournament two vertices have as many common out- as common in-neighbours (`|++| = |−−|`), and reversal swaps the two.
* **Mixed pairs.** For `j ∈ N^+(x)` and `y ∈ N^−(x)`, the count changes by `|O_j ∩ N^−_R(y)| − |O_j ∩ N^+_R(y)| = −(S_R 1_(O_j))_y`. Here `O_j = N^+(j) ∩ N^−(x)` and `S_R` is the `±1` skew matrix of `R`. So `T^(x)` doubly regular would force `S_R 1_(O_j) = 0` for every `j ∈ N^+(x)`.
* **The kernel of `S_R`.** `S_R ≡ J − I (mod 2)`. Over `F_2`, `(J − I)v = (Σv)·1 + v` vanishes iff `v ∈ span(1)` (for odd `m`), so `J − I` has `F_2`-rank `m − 1`. A nonzero `(m−1)`-minor mod 2 is nonzero over `Z`, so `rank S_R >= m − 1`. Skew-symmetry and odd order give `rank S_R <= m − 1`. Hence `ker S_R` is a line, and since `R` is regular, `ker S_R = span(1)`.
* **The contradiction.** `|O_j| = (2t + 1) − t = t + 1`, which is neither `0` nor `2t + 1` when `t >= 1`. So `1_(O_j) ∉ span(1)`, and `T^(x)` is not doubly regular. (A regular `T^(x)` that were doubly regular would need common-out-neighbour count exactly `t`, which the mixed pairs violate.)

**Step 3 (`0'` is fixed).**
* Let `Δ` be the set of vertices of `D(T)` whose out-neighbourhood induces a doubly regular tournament. `Δ` is `Aut(D(T))`-invariant.
* `N^+(0') = T`, so `0' ∈ Δ`. By Steps 1–2, `Δ ∩ T' = ∅`.
* Since `0' -> T`, `0'` beats every other vertex of `Δ`, so it has in-degree 0 in the tournament induced on `Δ`. A tournament has at most one such vertex. Every automorphism fixes it.
* The argument needs no knowledge of `Δ ∩ T`, which varies: all 7 vertices for `P_7`, exactly 1 for `D(P_11)`, none for Paley `P_q` with `q >= 11` (audit).

**Step 4 (the stabiliser is diagonal).**
* An automorphism `φ` fixing `0'` preserves `N^+(0') = T` and `N^−(0') = T'`, and `σ = φ|_T ∈ Aut(T)`.
* For `i ∈ T`, `N^−(i') ∩ T = N^−[i]`, the closed in-neighbourhood (`y -> i'` iff `y -> i` or `y = i`). Hence `φ(i') = j'` with `N^−[j] = N^−[σ(i)]`.
* Closed in-neighbourhoods are distinct in a tournament: if `N^−[i] = N^−[j]` with `i ≠ j`, then `i -> j` and `j -> i`. So `j = σ(i)`.
* Conversely, every diagonal extension preserves all the arc rules. ∎

## Checks (FINITE-EXACT)

* Steps 1–2 were checked vertex by vertex for `T_3, …, T_7` and for `P_7, P_11, P_19, P_23, P_31, P_43`: the relabelling, the failure of double regularity, `|O_j| = t + 1`, and `S_R 1_(O_j) ≠ 0`.
* `rank_F2(J − I) = m − 1` was checked for odd `m <= 39`.
* `|Aut D(T)| = |Aut T|` (nauty) for all of these, and for `D(D(P_7))` and `D(D(P_11))`. The orders are 21, 21, 21, 21, 21, 21, 55, 171, 253, 465, 903, 21, 55.
* The independent audit (own code, dreadnaut and a nauty-free enumerator) confirmed `|Aut T_k| = 21` for `k <= 9`, Paley primes up to 83 and `P_27`. It also checked 438 isomorphism classes of doubly regular tournaments that are neither Paley nor tower (orders 15, 19, 27, 31, 35): every step check passes, and `|Aut|` is preserved under doubling.

## Remarks

* **Where `t = 0` fails.** At `T = P_3`, `O_j = N^−(x)` and the argument fails. Indeed `|Aut T_2| = 3`, while `|Aut T_3| = |Aut P_7| = 21` exceeds it.
* **Prior art.** This theorem is Theorem 3.4 of A. Hanaki, arXiv:2011.06141 (2020). His matrix (3.1) is `D(T)` verbatim, and his hypothesis `m >= 7` is our `t >= 1`.
  * He fixes `0'` with triple intersection numbers (Lemma 3.1: `|N(i) ∩ N(j) ∩ N(k)|` is maximal exactly on the triples `{a, 0', a'}`). He proves the diagonal form from the invertibility of the adjacency matrix (Lemma 3.3).
  * By his footnote, the intransitivity is already Theorem 2.6.6 of Faradžev–Klin–Muzichuk (1994); only the isomorphism of the groups was new in 2020.
  * His Remark 3.5 notes that `n = 4` gives the Fano-plane scheme (`P_7`), which is schurian.
  * HYP-9162 (posed 2026-09-26) was therefore already a corollary of the literature.
  * The proof here is different. Steps 1–2 identify `N^+(x')` with `T^(x)` and rule it out by the kernel of an odd skew matrix. Step 4 uses closed in-neighbourhoods.
* **History in the repo.**
  * HYP-9162 (opus, 2026-09-26) conjectured the corollary.
  * Its audit proved the diagonal-stabiliser step.
  * The mod-18 note (2026-10-06, section 5) reduced the conjecture to the second-copy lemma.
  * Steps 1–2 were found in this session while reading the openai/math Hadamard papers. No technique from those papers is used.
  * The independent audit found Hanaki's paper (MISTAKE-571).
* **Collatz.** The tower is a tournament object. Its tie to Collatz is the Frobenius dictionary of THM-4553 (`F_21`, the Collatz Frobenius `C_3`), not a map.
