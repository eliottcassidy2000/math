---
id: THM-4765
title: "Quadratic twists of a square-root cover average to the branch complement: over F_q (q odd), for n linear forms with no common zero on P^d, the 2^(n-1) quadratic twist classes of the root variety X satisfy sum_t #X_t = 2^(n-1) #P^d and sum_t #X_t^{etale} = 2^(n-1) #U, U the complement of the branch hyperplanes; the twisted etale counts are the Walsh transform of the even character sums S_J = sum_{y in U} prod_{j in J} chi(L_j(y)). For the hourglass (seven roots E, E±x, E±z, E±(x+z)) the average over its 64 twists is #U = q^2-6q+12 = (q-1)^2-4(q-1)+7 away from characteristics 2 and 3 (q^2-6q+10 in characteristic 3), whose value 7 at q=1 is the Euler characteristic of the complex branch complement and the etale-sheet term 2^6·7 of e(Y)=80, and whose negative coefficient excludes a split-torus stratification of U; the equality chi(U) = p_g = 7 is a coincidence of the three triple points, failing for 18 of the 36 seven-line subsets."
status: >
  PROVED (klein, 2026-10-25; elementary pointwise character averaging, complete proof in the source packet).
  FINITE-EXACT: checked for the hourglass, the six-root K3 subset and all nine magic forms over F_q,
  q in {5,7,11,13,17,19,23,29,31,37,41,43,25,49} and characteristic-three controls q = 3, 9; the Walsh identity
  for all 64 hourglass twists at every tested q. The 1+6 split of the six-subset character sums is PROVED for the
  six outer deletions (configuration symmetry) and FINITE-EXACT for the central deletion being different.
  Not yet independently audited by another agent.
  CONCURRENT PRIORITY (recorded 2026-10-25 after reading the concurrent packets): the 64-twist mean q^2-6q+12 /
  q^2-6q+10, its formal value seven and the negative (q-1)-coefficient were stated first in
  05-knowledge/results/normseven_f1_20261025/TWIST_AVERAGE.md; the general Frobenius-orbit mass theorem
  (1/|G|) sum_g #W_g = #(W/G) for any constant finite abelian group, the resolved mean q^2+4q+1 and the fixed-twist
  extension law are in session63_integral_observers_root_arrangements_20261025/twists/RESULT.md (Theorems A, B, D).
  Theorem C here is the torsor special case with an independent pointwise proof. New here: the Walsh-transform
  formula for the twisted etale counts, the 1+6 outer/central split of the six-subset character sums, and the
  seven-subset hostile chi(U)=10-t_3 versus p_g=7.
session: klein-2026-10-25 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/clone_depth_modular_twist_average_20261025/README.md
scripts:
  - 05-knowledge/results/clone_depth_modular_twist_average_20261025/check_twist_average.py (+ twist-normal.json / twist-optimized.json)
related:
  - 05-knowledge/results/session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md (Theorem T1, Corollary T2, Euler stratification)
  - 05-knowledge/results/session59_f1_geometry_and_information_20261025/arithmetic-configurations/RESULT.md (negative (q-1)-coefficient obstruction)
  - 05-knowledge/results/session60_six_seven_geometry_20261025/deck-k3/LATER_RESULTS.md (1+6 split of the K3 deletions)
  - 05-knowledge/results/six_seven_gale_cube_20261025/BUDGET_AND_GALE.md (minus-one quadratic twist of the net double plane)
  - 01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md
---

# THM-4765 — Twists of a square-root cover average to the branch complement

## Statement

`k = F_q`, `q` odd; `L_0,…,L_{n-1}` linear forms on `k^{d+1}` with no common zero in `P^d`; `X ⊂ P^{n-1}` the
inverse image of the linear `P^d` under coordinate squaring; `X_t` the twist `x_j^2 = t_j L_j`; twist classes
`(k^×/k^{×2})^n / diag`, `2^(n-1)` of them; `U = P^d \ ∪{L_j = 0}`; `X_t^{ét}` the preimage of `U`. Then

```
sum_t #X_t(k) = 2^(n-1) #P^d(k),     sum_t #X_t^{ét}(k) = 2^(n-1) #U(k),
#X_t^{ét}(k) = sum_{|J| even} (prod_{j∈J} χ(t_j)) S_J,   S_J = sum_{y∈U} prod_{j∈J} χ(L_j(y)),   S_∅ = #U(k).
```

## Proof

Over `y ∈ P^d(k)` with `m` vanishing forms, the fibre of `X_t` has `2^{n-m-1}` points if the nonzero
`χ(t_j L_j(y))` are all equal and `0` otherwise (count pairs `(λ, x)` with `x_j^2 = λ t_j L_j(ỹ)` modulo the free
`k^×`-action `(λ,x) ↦ (μ^2λ, μx)`; half of the `λ` have the required sign). Exactly `2^m` of the `2^(n-1)` twist
classes make the vector constant. Hence every `y` contributes `2^(n-1)` to the sum over twists. The
Walsh formula is the expansion of the indicator of a constant sign vector over the even subsets. ∎

## Hourglass corollary

Seven lines with `t_2 = 12`, `t_3 = 3` (characteristic `≥ 5`): `#U = q^2 - 6q + 12 = (q-1)^2 - 4(q-1) + 7`;
characteristic 3: `q^2 - 6q + 10`. The value `7` at `q = 1` is `e(U(C))` and the étale-sheet term `2^6·7` of
`e(Y) = 80`. A finite split-torus stratification of `U` would need nonnegative `(q-1)`-coefficients; `-4`
excludes it (same argument as Session59 §4). Across the 36 labelled seven-line subsets `χ(U) = 10 - t_3`
while `p_g = 7` throughout; the two sevens coincide only for `t_3 = 3` (18 subsets).

## Mechanism, loss ledger, and what the additive equations add

The average over twists is the projection onto the trivial character of the twist group; it forgets all
`S_J`, `J ≠ ∅`. Those coefficients see the linear relations among the forms: the configuration symmetry
(order 12, fixing `E`, transitive on the six outer forms and on the three parallel pairs) proves that the six
outer-deletion six-subset sums coincide and the three parallel-pair sums coincide; the central-deletion sum
differs at every tested `q ≥ 5` (`(outer, central) = (3,7), (5,7), (9,31), (5,19), (37,19), …`), the finite-field
shadow of Session60's `1+6` K3 split. No identification with Picard ranks is claimed. The multiplicative
(Kummer) construction alone determines the average; everything finer needs the additive branch equations.

## Non-consequences

Nothing about rational points of one twist, Picard ranks, `p_g`, integer magic squares of squares, or
characteristic 2 (inseparable cover; the twist group is trivial).
