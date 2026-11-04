# The sole-exit tower is tribonacci — proof

**Setting.** `T` is a tournament on the vertex set `V`, `t ∈ V`. Write `h(T)` for the number of Hamiltonian paths, `e_t(T)` for the number ending at `t`, `s_t(T)` for the number starting at `t`, `h(T[∅]) = 1`, and

`B_t(T) = Σ_{A ⊆ V∖{t}} h(T[A]) · s_t(T[V∖A])`   (the *root response*; `A = ∅` contributes `s_t(T)`, `A = V∖{t}` contributes `h(T − t)`).

This is the lineage's `S = [P A_t]` (nonzero_preparation_20261103) and the near-join's `B(T, t)` (HYP-9453): it counts the pairs `(P₁, t P₂)` of a Hamiltonian path `P₁` of some `T[A]` (possibly empty) and a Hamiltonian path of the complement starting at `t`.

**The one-exit extension.** `T + x` is `T` with a new vertex `x` whose only out-neighbour is `t` (`N⁺(x) = {t}`, `N⁻(x) = V∖{t}`).

*Lemma 1 (one-exit law; THM-4862 (a)).* `h(T + x) = h(T) − e_t(T) + B_t(T)`.
*Proof.* A Hamiltonian path of `T + x` either ends at `x` — its predecessor is any vertex other than `t`, so these are the `h(T) − e_t(T)` Hamiltonian paths of `T` not ending at `t`, each extended by one arc — or contains `x → t`; then it is `P₁ x t P₂` with `P₁` a Hamiltonian path of `T[A]` for some `A ⊆ V∖{t}` (every vertex beats `x`, so any `P₁` may precede `x`; `A = ∅` is the case `x` first) and `t P₂` a Hamiltonian path of `T[V∖A]` starting at `t`: `B_t(T)` of them. ∎

*Lemma 2.* `e_x(T + x) = h(T) − e_t(T)`. *Proof.* The first case above. ∎

*Lemma 3.* `B_x(T + x) = B_t(T) + h(T)`.
*Proof.* `B_x(T + x) = Σ_{A ⊆ V} h((T+x)[A]) · s_x((T+x)[(V∖A) ∪ {x}])`. For a set `W ∋ x`, a Hamiltonian path of `(T+x)[W]` starting at `x` continues to `t` (the only out-neighbour of `x`), so `s_x((T+x)[W]) = s_t(T[W∖{x}])` when `t ∈ W`, and `= [W = {x}]` when `t ∉ W`. Hence the terms with `t ∈ A` vanish except `A = V`, which contributes `h(T) · 1`, and the terms with `A ⊆ V∖{t}` contribute `Σ_{A ⊆ V∖{t}} h(T[A]) s_t(T[V∖A]) = B_t(T)`. ∎

**The tower.** `T_0 = T`, `x_0 := t`, and `T_k = T_{k−1} + x_k` with `N⁺(x_k) = {x_{k−1}}`. Put `h_k = h(T_k)`, `e_k = e_{x_k}(T_k)`, `B_k = B_{x_k}(T_k)` (so `e_0 = e_t(T)`, `B_0 = B_t(T)`). Lemmas 1–3 at each step give, for `k ≥ 1`,

`h_k = h_{k−1} − e_{k−1} + B_{k−1}`,  `e_k = h_{k−1} − e_{k−1}`,  `B_k = B_{k−1} + h_{k−1}`.   (R)

**Theorem.** `h_1 = h − e + B`, `h_2 = h + 2B`, and `h_k = h_{k−1} + h_{k−2} + h_{k−3}` for all `k ≥ 3`, where `h = h(T)`, `e = e_t(T)`, `B = B_t(T)`. Consequently `h_k = h·A_k − e·C_k + B·D_k`, where `(A_k), (C_k), (D_k)` are the tribonacci-recurrence sequences with initial values `(1, 1, 1)`, `(0, 1, 0)`, `(0, 1, 2)`: `A = A000213 = 1, 1, 1, 3, 5, 9, 17, 31, 57, 105, 193, …`, `C = A001590 = 0, 1, 0, 1, 2, 3, 6, 11, 20, 37, 68, …`, `D_k = C_{k+2} = 0, 1, 2, 3, 6, 11, 20, 37, 68, 125, 230, …`.

*Proof.* `h_1 = h − e + B` is Lemma 1. By (R), `h_2 = h_1 − e_1 + B_1 = (h − e + B) − (h − e) + (B + h) = h + 2B`. For `k ≥ 1`, subtracting consecutive instances of the first relation in (R) and using the other two,

`h_{k+1} − h_k = (h_k − h_{k−1}) − (e_k − e_{k−1}) + (B_k − B_{k−1}) = (h_k − h_{k−1}) − (h_{k−1} − 2e_{k−1}) + h_{k−1}`,

i.e. `h_{k+1} = 2h_k − h_{k−1} + 2e_{k−1}` for `k ≥ 1`.   (★)

For `k ≥ 2`, (★) at `k − 1` gives `2e_{k−2} = h_k − 2h_{k−1} + h_{k−2}`, and `e_{k−1} = h_{k−2} − e_{k−2}` gives `2e_{k−1} = 2h_{k−2} − 2e_{k−2} = 2h_{k−2} − h_k + 2h_{k−1} − h_{k−2}`. Substituting into (★):

`h_{k+1} = 2h_k − h_{k−1} + 2h_{k−2} − h_k + 2h_{k−1} − h_{k−2} = h_k + h_{k−1} + h_{k−2}`   (`k ≥ 2`),

which is the claim for `k + 1 ≥ 3`. The closed form follows because `(h, e, B) ↦ (h_k)_{k ≥ 0}` is linear (by (R)) and the initial triples `(h_0, h_1, h_2)` are `(1, 1, 1)` for `(h, e, B) = (1, 0, 0)`, `(0, −1, 0)` for `(0, 1, 0)` and `(0, 1, 2)` for `(0, 0, 1)`; every sequence obeying the tribonacci recurrence is determined by its first three terms, which gives `A_k`, `−C_k` and `D_k`. `D_k = C_{k+2}` holds because both sides satisfy the recurrence and agree at `k = 0, 1, 2` (`C_2, C_3, C_4 = 0, 1, 2`). ∎

**Corollaries.**

1. *(The lineage's C2 law at any vertex.)* The lineage's two births — `x` with `x → a` and every other old vertex `→ x`, then `y` with `y → x` and every old vertex `→ y` — are the tower of height two over `(T, a)`, so the child has `h + 2B_a(T) = 3h + 2σ` with `σ = B_a(T) − h(T) = S − h`, for **every** vertex `a`, which is the lineage's law (10) in `multiport_defect_grammar_20261103/construction/RESULT.md` (`σ = S − h ≥ 0` there); the sole-exit condition is its *joint* guard (for the root functional `F`, `τ = 0`), not a condition for the `H`-law.
2. *(The digit identity, THM-4862 (b), in one line.)* If `a` is a sole exit `a → r` of `T = U + a`, then `B_a(T) = B_r(U) + h(U)` and `h(T) = h(U) − e_r(U) + B_r(U)` (Lemmas 3 and 1), so `σ = B_a(T) − h(T) = e_r(U) = e_r(T − a)`.
3. *(The tribonacci partner.)* The tower over a single vertex is the full tiling tournament `T*_n` (`i → i+1`, `j → i` for `j ≥ i+2`) with `(h, e, B) = (1, 1, 1)`: `h_k = A_k − C_k + D_k = A000213(k+1)`, which is THM-4731 (6) / kind-pasteur-2026-03-05-S11 (`04-computation/verify_tribonacci_proof.py`); its `e(n) = h(n−1) − e(n−1)` is the second relation of (R).
4. *(Near-join.)* HYP-9453's family `T → M_k → z` (with `z → t` reversed) has `h + 2^k B_t(T)` with the same `B`; at `k = 1` it has the same value as the height-two tower.
5. *(Re-rooted tower.)* For the digit-1 class whose entry attachment is the pure backward cycle (`F = ∅` in the chord picture of `data/chords_all.txt`), `T` is the tower over a point of height `n − 2` with one more sole-exit vertex pointing at the base vertex `x_0`, so `H = h_{n−2} − e_{x_0}(T_{n−2}) + B_{x_0}(T_{n−2})`: `5, 13, 31, 77, 193, 489, 1255, 3253, 8493, 22293` at orders 4–13 (`scripts/fempty_check.py`; not in OEIS on 2026-11-04).
