# Proofs and verified laws

klein-2026-11-04. Terms: for a tournament `T` on `V`, `h(T)` is its number of Hamiltonian paths, `s_v(T)` and `e_v(T)` those starting and ending at `v`, and for a vertex `a`, `[P A_a](T) = Σ_{S ⊆ V − a} h(T[S]) · s_a(T[V − S])` (ordered two-slot covers with the second slot starting at `a`; `h(∅) = 1`), the quantity `S` of `multiport_defect_grammar_20261103/construction/RESULT.md` §1. A vertex `a` is a **sole exit into `r`** if `N⁺(a) = {r}`. The inherited C2 law (RESULT.md §2, JOINT_CLASSIFICATION.md): for a strong input with sole exit `a → r`, the two births give `H′ = 3h + 2σ` with `σ = [P A_a] − h = e_r(T − a)`.

## 1. The one-exit extension law (PROVED)

**Lemma 1.** Let `T` be a tournament on `V`, `a ∈ V`, and let `T + x` be `T` with a new vertex `x` such that `x → a` and `v → x` for every `v ≠ a` (the first C2 birth). Then

    h(T + x) = h(T) − e_a(T) + [P A_a](T),

equivalently `h(T + x) = 2h(T) + σ − e_a(T)` with `σ = [P A_a] − h`.

*Proof.* Classify the Hamiltonian paths of `T + x` by the position of `x`. If `x` is last, its predecessor `v` satisfies `v → x`, i.e. `v ≠ a`, and the rest is a Hamiltonian path of `T` ending at `v`: there are `h(T) − e_a(T)` such paths. If `x` is not last, the vertex after `x` is its only out-neighbour `a`; the vertices before `x` form a Hamiltonian path of `T[S]` for some `S ⊆ V − a` (possibly empty; its last vertex is automatically `≠ a`, hence beats `x`), and the vertices from `a` on form a Hamiltonian path of `T[V − S]` starting at `a`. Conversely every such pair gives one path. The count is `Σ_S h(T[S]) s_a(T[V − S]) = [P A_a](T)`. ∎

(Checked exhaustively on every tournament of orders 3–7 and every `a`: [`scripts/one_exit_lemma.py`](scripts/one_exit_lemma.py), 1,794 cases, no mismatch.)

## 2. The digit is a cycle count of the flipped attachment (PROVED)

**Lemma 2.** Let `U` be a tournament, `r ∈ U`, and attach a new vertex in the two possible one-arc ways: `U + a_exit` with `a_exit → r` and `v → a_exit` for all `v ≠ r`; `U + a_entry` with `r → a_entry` and `a_entry → v` for all `v ≠ r`. Then the number of directed Hamiltonian cycles satisfies

    C(U + a_exit) = s_r(U),      C(U + a_entry) = e_r(U).

Hence the readiness digit of `(U + a_exit, a_exit)`, which is `σ = e_r(U)` by the inherited classification, equals `C(U + a_entry)`.

*Proof.* A Hamiltonian cycle of `U + a_exit` passes through `a_exit` as `(…, v, a_exit, r, …)`; deleting `a_exit` leaves a Hamiltonian path of `U` from `r` to `v`, and every Hamiltonian path of `U` starting at `r` closes up since its end `v ≠ r` beats `a_exit`. The entry case is the converse statement. ∎

**Corollary 2.1 (the Fibonacci coincidence).** The map `Φ: (T, a) ↦ ((T − a)^op + a_exit, a)` is a bijection from strong tournaments with a sole-exit vertex `a → r` and digit `σ` onto strong tournaments with a sole-exit vertex and exactly `σ` Hamiltonian cycles. (`C((T − a)^op + a_exit) = s_r((T − a)^op) = e_r(T − a) = σ`; strongness of the image follows because `e_r(T − a) > 0` forces `T − a` strong, as a Hamiltonian path of a non-strong tournament ends in its last component while strongness of `T` puts `r` in the first.) This is why the numbers of `σ = 1`-ready classes, `1, 3, 8, 21, 55, 144, 377` at orders 4–10, are the Fibonacci numbers `F(2n − 6)` of Garey's theorem on tournaments with a unique Hamiltonian cycle (J. Combin. Theory Ser. B 13 (1972) 266–269, CITED; reproduced in THM-4846), and why the `σ = 2` counts `1, 5, 18, 60, 191, 589` agree with the two-cycle classes. `Φ` does not preserve `H`, so the ready spectra `R_σ` are not the path spectra of the unique-cycle classes; they are their images under the functional of Lemma 1 applied to the converse core.

## 3. Readiness of the near-join (PROVED)

**Lemma 3.** For `T′ = T → M_k → z` with the single arc `z → t₀` reversed (THM-4860's packet, `spectrum_analogies_nearjoin_20261102`), the last chain vertex `m_k` is a sole exit into `z`, and its digit is `e_z(T′ − m_k) = h(T)` for `k ≥ 2` and `h(T) − e_{t₀}(T)` for `k = 1`.

*Proof.* `m_k` loses to all of `T` and to `m_1, …, m_{k−1}` and beats only `z`. A Hamiltonian path of `T′ − m_k` ending at `z` cannot use `z → t₀`, so it is a Hamiltonian path of `T → M_{k−1}` followed by `z`; for `k ≥ 2` every such path ends in the chain and beats `z`, giving `h(T)`; for `k = 1` it is a Hamiltonian path of `T` ending at a vertex that beats `z`, i.e. not `t₀`. ∎ (Verified on 246 cases of orders 3–5.)

Consequently every tournament `T` and every attained value `h(T)` occurs as a digit: the C2 child of `(T′, m_k)` has `H = 3(h + 2^k B(T, t₀)) + 2h = 5h + 3·2^k B(T, t₀)`.

## 4. The sandwich law (VERIFIED, orders 3–6)

For a tournament `X`, let `S(X)` be the tournament on `X ∪ {o, m, r, a}` with `o → X → m → o`, `o → r`, `r → X`, `r → m`, `a → r`, and every other vertex beating `a`. Then `a` is a sole exit into `r`, the digit is `σ = e_r(S(X) − a) = h(X)`, and

    H(S(X)) = −h(X) + 4K(X) + 2K₃(X),

where `K = [P²]` and `K₃ = [P³]` are the ordered two- and three-slot cover counts of `X` (empty slots allowed). Verified on all 74 tournaments of orders 3–6 ([`scripts/sandwich.py`](scripts/sandwich.py)); the digit identity follows from Lemma 2 since the only in-neighbour of `r` is `o` and the only Hamiltonian paths of `S(X) − a` ending at `r` are `(path of X), m, o, r`. The law shows that the C2 sweep's digit can be any attained value and that `H` of a ready tournament can be a linear form in slot responses; it is not a dense family by itself (see README §4).

## 5. Prefix joins scale the count and the digit together (PROVED)

**Lemma 4.** Let `T` have a sole exit `a → r` with digit `σ`, and let `J ▷ T` be the join in which every vertex of `J` beats every vertex of `T`. Then `a` is still a sole exit into `r` in `J ▷ T`, `H(J ▷ T) = h(J) · H(T)`, and the digit of `(J ▷ T, a)` is `h(J) · σ`.

*Proof.* Every Hamiltonian path of a join traverses `J` first (THM-1862), so paths of `J ▷ T` are pairs of paths, and paths of `(J ▷ T) − a` ending at `r` are pairs with the second path a Hamiltonian path of `T − a` ending at `r`. `a` loses to all of `J`, so its out-neighbourhood is unchanged. ∎ (Verified on 98 cases: all sole-exit 5-tournaments with all prefixes of orders 2–4.)

Consequently the lineage's padding by a transitive prefix (`h(J) = 1`) preserves both coordinates, while a nontrivial prefix multiplies the digit: digits `σ ≥ 3` are produced in bulk by joins, and the C2 child of `J ▷ T` is `h(J)(3H + 2σ)`, the join of `J` with the child of `T`. Readiness with small digits cannot be manufactured by joins; it lives in the strong core.

## 6. The block lemma behind the pruned generators (PROVED)

**Lemma 5.** Let `Y` be a tournament with a Hamiltonian cycle `C = (c₀, c₁, …, c_{n−1})` and let `Y + w` be an extension by one vertex. Write `b_i = 1` if `c_i → w` and `b_i = 0` if `w → c_i`. Every cyclic position `i` with `b_i = 1`, `b_{i+1} = 0` gives a Hamiltonian cycle of `Y + w` (insert `w` between `c_i` and `c_{i+1}`), and distinct positions give distinct cycles. Hence if `Y + w` has at most `k` Hamiltonian cycles, the cyclic word `b` has at most `k` blocks of 1s; if `Y + w` has exactly one Hamiltonian cycle, `b` is exactly one nonempty, proper cyclic block (a constant word makes `w` a source or a sink, which has no Hamiltonian cycle at all).

*Proof.* Immediate from the definitions: a `1 → 0` transition at `i` means `c_i → w → c_{i+1}`, so replacing the arc `c_i → c_{i+1}` of `C` by the path `c_i, w, c_{i+1}` is a Hamiltonian cycle; the resulting cycles differ in the neighbours of `w`. ∎

**Application.** For a strong tournament `T` with a sole exit `a → r` of digit σ, THM-4862 (b) says `Y = Φ(T) = (T − a)^op + a_exit` has exactly σ Hamiltonian cycles. An extension `T + w` keeps `a` a sole exit only if `w → a`, and its digit is the cycle count of `Y + w'` where `w'` carries the conversed pattern on the core. By Lemma 5, digit 1 forces the "beats `w'`" positions along the unique cycle of `Y` to form one cyclic block not containing `a`; digit ≤ 2 forces at most two blocks along any Hamiltonian cycle of `Y`. The generators [`scripts/gen1.c`](scripts/gen1.c) and [`scripts/gen2.c`](scripts/gen2.c) enumerate exactly these patterns plus the `n` patterns in which `w` itself becomes the sole exit, test each candidate exactly, and are therefore complete for all extensions in which the old sole exit survives or the new vertex is the sole exit. Together with heredity (verified exhaustively through order 11: every ready class of order `n + 1` is such an extension of a ready class of order `n`) and Garey's count `F(2n − 6)`, matching counts certify that no digit-1 class is missed.
