# The path splice: HYP-9862 proved, and the order-statistics bound for any overlap

Session klein-2026-11-16-splice. Script `scripts/shuffle_check.py` (exhaustive, orders 4–7); receipts
`shuffle_check_4..7.out`.

## 1. Statement

Let `U` be a tournament and `A`, `B` induced subtournaments with `W = V(A) ∩ V(B)`. For a linear order `σ` of `W`
write `h_A(σ)` for the number of Hamiltonian paths of `A` whose restriction to `W` (as a subsequence) is `σ`, and
`h_B(σ)` likewise.

**Theorem (path splice).** The number of Hamiltonian paths `Z` of `U` whose restrictions `Z|_A` and `Z|_B` are
Hamiltonian paths of `A` and of `B` is at least `Σ_σ h_A(σ) h_B(σ)`. Hence

    H(U) ≥ Σ_σ h_A(σ) h_B(σ),

and when `|W| ≤ 1` (one order `σ`, `h_A(σ) = H(A)`, `h_B(σ) = H(B)`)

    H(U) ≥ H(A) · H(B)        — HYP-9862, PROVED.

## 2. Proof

**Merge lemma.** Let `X = x_1 → … → x_p` and `Y = y_1 → … → y_q` be vertex-disjoint directed paths of a tournament.
Then the tournament on `V(X) ∪ V(Y)` has a Hamiltonian path that is a shuffle of `X` and `Y` (its restriction to
`V(X)` is `X` and to `V(Y)` is `Y`).

*Proof.* Induction on `p`, the case `p = 0` being `Y` itself. Let `Z′` be a shuffle of `x_1 → … → x_(p−1)` and `Y`.
Since `x_(p−1)` is the last vertex of `X` in `Z′`, everything after it in `Z′` is a suffix `y_i → … → y_q` of `Y`
(possibly empty). Insert `x_p` into the segment `x_(p−1), y_i, …, y_q`: if the segment ends at `x_(p−1)`, or if
`x_p → y_i`, place `x_p` right after `x_(p−1)` (the arc `x_(p−1) → x_p` is an arc of `X`); otherwise `y_i → x_p`,
and with `j ≥ i` the largest index such that `y_j → x_p`, place `x_p` after `y_j` (then `x_p → y_(j+1)` by
maximality, or `y_j` is last). The result is a directed path, and the relative orders of `X` and of `Y` are kept. ∎

**Proof of the theorem.** Let `P_A`, `P_B` be Hamiltonian paths of `A` and `B` inducing the same order
`σ = (w_1, …, w_r)` on `W`. Write `P_A = α_0 w_1 α_1 w_2 … w_r α_r` and `P_B = β_0 w_1 β_1 … w_r β_r`, where the
`α_i` are (possibly empty) paths in `A − W` and the `β_i` paths in `B − W`. For each `i` let `μ_i` be a merge of `α_i`
and `β_i` given by the lemma (empty if both are empty), and set

    Z = μ_0 w_1 μ_1 w_2 … w_r μ_r.

`Z` is a Hamiltonian path of `U`: inside each `μ_i` the consecutive vertices are joined by arcs of `U` (lemma); the
first vertex of a nonempty `μ_i` is the first vertex of `α_i` or of `β_i`, which `w_i` beats in `P_A` or in `P_B`; the
last vertex of a nonempty `μ_i` is the last vertex of `α_i` or of `β_i`, which beats `w_(i+1)`; and if `μ_i` is empty
then `w_i → w_(i+1)` is an arc of `P_A`. Since the merges are shuffles, `Z|_A = P_A` and `Z|_B = P_B`. The map
`(P_A, P_B) ↦ Z` is therefore injective — its image determines the pair by restriction — and lands in the set of
Hamiltonian paths of `U` restricting to Hamiltonian paths of both. Pairs inducing different orders on `W` are
incompatible (no common extension), so the count is `Σ_σ h_A(σ) h_B(σ)` and not `H(A) H(B)` in general. ∎

The whole argument is the standard proof that every tournament has a Hamiltonian path (insert one vertex at a time),
applied segment by segment between the shared vertices; what HYP-9862's odd-cycle splice lacked — a decodable image —
the path splice has for free, because restriction is the decoder.

## 3. Corollaries

- **HYP-9862** (one shared vertex buys nothing): `H(A ∪ B) ≥ H(A)H(B)` for `|A ∩ B| ≤ 1`, every order. The census of
  klein-2026-11-16-join (218,925,648 instances at order 9, no violation) is now a corollary.
- **THM-5312 (A)** is a corollary: for odd primes `p, q` not both 3, `pq ≥ p + q + 7` (as `(p−1)(q−1) ≥ 8`), with
  equality at `(3, 5)`; and the one-vertex excess over `p + q` is at least `(p − 1)(q − 1) − 1`, the bound HYP-9862
  predicted.
- **Two shared vertices** `u, v`: `H(U) ≥ h_A(uv) h_B(uv) + h_A(vu) h_B(vu)`. For two triangles on a common arc `u → v`
  this gives `2·2 + 1·1 = 5 = H(U)`: the bound is tight at the very example that refutes the two-vertex product bound
  (`5 < 9`). So the second shared vertex buys exactly the order statistics.
- **Excess formula.** `H(U) − Σ_σ h_A(σ) h_B(σ) = Σ_(compatible pairs) (m(P_A, P_B) − 1) + #{Hamiltonian paths of U
  whose restriction to A or to B is not Hamiltonian}`, where `m(P_A, P_B) ≥ 1` is the number of common shuffles.
  Equality in HYP-9862 holds iff every Hamiltonian path of `U` restricts to Hamiltonian paths of `A` and `B` and every
  pair has a unique merge (orders 4–7: 114 / 1,476 / 22,398 / 479,878 equalities among the `|W| ≤ 1` instances; at
  order 7 that is 19.8 % of 2,419,536).
- **The odd-cycle form (the splice injection in OCF language).** With `Mix(U)` the weight of the odd-cycle packings of
  `U` (THM-002) that contain a cycle meeting both `A − v` and `B − v`, the packings with no such cycle are the direct
  unions of packings of `A` and of `B` in which `v` lies on a nontrivial cycle in at most one of them, so

      H(U) = H(A) H(B − v) + H(A − v) H(B) − H(A − v) H(B − v) + Mix(U),

  and HYP-9862 is equivalent to `Mix(U) ≥ (H(A) − H(A − v)) (H(B) − H(B − v))`: the mixed packings outweigh the pairs
  of `v`-cycles. The path splice proves this inequality without exhibiting the injection on packings that the
  join session looked for; the triangle splices that "forget the cycles" are the reason that injection is hard, and
  the inequality can be strict only through mixed packings with singletons (example: two triangles at `v` with every
  arc from `A − v` to `B − v` has `Mix = 4 = (3 − 1)(3 − 1)`, equality).

## 4. Exhaustive check (FINITE-EXACT)

`shuffle_check.py n` enumerates, for every tournament of order `n` (gentourng) and every pair of vertex subsets
`S_A, S_B` (sizes ≥ 2, neither contained in the other, `|S_A ∩ S_B| ≤ 2`), the Hamiltonian paths of `A`, `B`, `U`, the
order statistics, the bound and the shuffle count.

| order | instances (`|W| ≤ 2`) | violations of the bound | shuffle count below the bound | equalities, `|W| ≤ 1` |
|---|---|---|---|---|
| 4 | 264 | 0 | 0 | 114 |
| 5 | 5,040 | 0 | 0 | 1,476 |
| 6 | 117,040 | 0 | 0 | 22,398 |
| 7 | 4,143,216 | 0 | 0 | 479,878 |

Pre-registered P-1 (0 violations at order 7; shuffle count never below the bound) **confirmed**; P-2 (equality in
more than 15 % of the `|W| ≤ 1` instances at order 7) **confirmed** (19.8 %).
