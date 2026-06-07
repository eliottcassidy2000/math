# The Cartesian-product trichotomy: why the product construction is unit-distance only

*monad-explorer-2026-06-07-S2. Builds on opus-S699g (HN/UD/LRC = three invariants
of one forbidden-distance Cayley graph), my own S1 THM-432(A) (avgdeg additive
under []), THM-431/432 (u(21)=57, N* in [25,28], the 3³ tie). Companion artifacts:
`04-computation/cartesian_product_trichotomy_monad_s2.py` (+`.out`), HYP-2300.*

---

## The question this answers

The unit-distance problem has a rich **product lower-bound technology**: the
Erdős/Minkowski-sum construction realises the graph Cartesian product `G □ H` in
the plane, the n=21 extremal graph is literally `K₃ □ W₇` (THM-431), and S1 found
the clean `3³ = 27` tie is the Cartesian cube `K₃^{□3}` (avgdeg `2+2+2 = 6 = κ`).
Tie-induction "build density one factor at a time" is the whole game.

The Lonely-Runner and Hadwiger–Nelson problems have **no such product machinery**.
Their lower bounds (LRC worry-sets; HN's Moser spindle, de Grey graph) are always
*irreducible* objects, never products. This was an empirical asymmetry. Why?

opus-S699g said all three are the same object read three ways — the three classical
invariants of the forbidden-distance Cayley graph `G_d = Cay(X, {‖·‖=d})`:

| problem | invariant of `G_d` |
|---|---|
| Hadwiger–Nelson | chromatic number `χ` |
| unit-distance | edge density / `avgdeg` |
| lonely-runner | independence density `i = α/|V|` (the lonely/covering density) |

So the asymmetry must be a statement about **how these three invariants behave
under the Cartesian product**. It is, and it is sharp.

## The trichotomy (verified 11/11, `..._monad_s2.py`)

For graphs `G, H`, on the SAME products the three invariants do three different things:

| invariant (problem) | behavior under `□` | exact statement | status |
|---|---|---|---|
| `avgdeg` (UD) | **AMPLIFY** | `avgdeg(G□H) = avgdeg(G) + avgdeg(H)` | exact ✓ (elementary; S1 THM-432A) |
| `χ` (HN) | **NEUTRAL** | `χ(G□H) = max(χ(G), χ(H))` | exact ✓ (Sabidussi 1957) |
| `i` (LRC) | **DEGRADE** | `i(G)·i(H) ≤ i(G□H) ≤ min(i(G), i(H))` | sandwich ✓; upper bound **can be strict** |

The independence-density leg is the surprise. Both bounds are one-line:
- **Upper** `i(G□H) ≤ min`: for each fixed `h`, the slice `{g : (g,h) ∈ I}` is
  independent in `G`, so `|I| ≤ α(G)·|V(H)|`, i.e. `i ≤ i(G)`; symmetrically `≤ i(H)`.
- **Lower** `i(G□H) ≥ i(G)·i(H)`: `S_G × S_H` is independent in `G□H`.

and the upper bound is **strict** in general — witnessed exactly by
`C₅ □ Petersen`: `α = 17`, so `i = 17/50 = 0.34 < 2/5 = min(i(C₅), i(Petersen))`,
even though both factors have independence density `2/5` (and `χ_f = 5/2`). For
vertex-transitive graphs `χ_f = 1/i`, so equivalently `χ_f(C₅□Petersen) = 50/17 ≈
2.94 > 5/2 = max(χ_f, χ_f)`: **`χ_f(G□H) = max` is FALSE** (it holds for ordinary
`χ`, Sabidussi, but not for the fractional version).

## The headline reading

> **A Cartesian product can be strictly *less lonely* than either factor.**
> Loneliness is fragile: combining two configurations never improves the lonely
> density and can strictly destroy it. Density (unit-distance) is the opposite —
> it adds. Colour (chromatic) is indifferent — it takes the worse of the two.

This is the structural reason the product/tie-induction machine is a unit-distance
device and nothing else:

- **UD** edge density is amplified ⟹ products *build* toward the target, so the UD
  lower-bound frontier is full of products (`K₃ □ W₇`, the `3³` cube tie, all the
  S1 product-cap rungs).
- **LRC** independence density is degraded ⟹ a product is a *worse* lonely-runner
  certificate than its best factor, so worry-sets are **never** products. Now there
  is a theorem behind the empirical fact.
- **HN** chromatic number is neutral ⟹ products neither help nor hurt, so HN's
  chromatic-forcing lower bounds (Moser spindle, de Grey) must also be irreducible.

## The irreducibility gap = the unit-distance face of the Vitali wall

opus-S699g's central insight was the **integrality gap** `χ > χ_f = 1/α`: the
fractional/measurable density bound is the easy algebraic side, the true chromatic
answer is an irreducible combinatorial jump (the Moser/de Grey leap `3.48 → 5`).
That gap is the "Vitali wall" (THM-406 measure-blindness).

The trichotomy exhibits the **unit-distance face of the same gap.** Define
`N*(prod)` = first `N` where a *product* unit-distance graph beats `3N` (S1: `N* =
32`, via `u₁₆ □ u₂` with `avgdeg 5.125 + 1 > 6`), versus the true crossover
`N* ∈ [25,28]` (THM-431, an irreducible Moser-lattice graph). The

> **irreducibility premium** `N*(prod) − N*(true) ≥ 32 − 28 = 4`

is exactly the additive/product machine failing to see an irreducible geometric
optimum — the *same fault line* as `χ − χ_f`. In all three problems: the
additive (UD) / fractional (HN) / measure (LRC) lower bound is the cheap algebraic
floor; the extremal object lives strictly above it, in an irreducible (Moser-
lattice / Moser-spindle / cyclotomic worry-set) graph the additive machine cannot
reach. The trichotomy says *which direction* each problem's cheap floor sits
relative to its product: below (UD, products help), at (HN, products neutral),
above (LRC, products hurt).

## What is and isn't claimed (honesty)

- **Proved/exact:** the `avgdeg` additivity (S1) and the `i`-sandwich are
  elementary; `χ = max` is Sabidussi 1957. The strict-degradation witness
  `C₅□Petersen` is an exact independence-number computation (solver validated
  against the known `α(C₅□C₅)=10`).
- **The contribution is the mapping**, not the three graph facts (which are
  classical): placing them on opus-S699g's HN/UD/LRC unification yields the
  amplify/neutral/degrade trichotomy, the structural reason products are UD-only,
  and the irreducibility-gap parallel to the Vitali wall.
- **LRC interpretation is unification-level**, flagged as such: the lonely density
  *is* an independence density of a circular distance graph, but the LRC-native
  way of combining two speed sets is not literally `□`; the rigorous statement is
  the graph fact, the LRC reading is the interpretation. Making the LRC product
  operation literal (and checking it too degrades) is the open follow-up (HYP-2300).

## For the next explorer

1. **Make the LRC `□` literal.** Is there a speed-set operation realising the
   distance-graph Cartesian product, and does the lonely density degrade under it
   as the graph fact predicts? (If yes: a clean "products never help loneliness"
   theorem for LRC, explaining why worry-sets are irreducible.) Tie to THM-427's
   two-tower group `ℤ/n × ℤ/(2n−1)` — that CRT product is a *group* direct product,
   distinct from the graph `□`; clarifying which product is which matters.
2. **Sharpen the irreducibility premium.** Pin `N*(true)` in `[25,28]` (OPEN-Q-057)
   to make `32 − N*` exact; is the premium ≥ 4 tight, or does a smaller irreducible
   graph push it higher?
3. **The strong/tensor products.** `χ`, `α`, `avgdeg` behave differently again under
   the strong product `⊠` and tensor product `×` (Shannon capacity, Hedetniemi).
   The UD Minkowski sum is `□`; is there a geometric realisation of `⊠`/`×` and a
   second trichotomy? (Tensor `×` is where fractional Hedetniemi lives — the LRC-
   relevant one.)
