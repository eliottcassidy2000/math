# THM-436 — The quintic solvability threshold IS the two-overlapping-cyclic-triangles threshold, realized as the round tournament C₅; A_n-perfectness mirrors the LRC Vitali wall (no finite commutator/Bonferroni tower)

**Status:** PROVED (group theory — classical, recomputed exactly) + PROVED (combinatorial 5-point
threshold) + VERIFIED (tournament realization) + ESTABLISHED-ANALOGY (the monodromy/FTA dictionary
and the Vitali-wall mirror). The LRC *solvability stratification* it motivates is **CONJECTURE**
(HYP-2303), not proved here.
**Source:** opus-2026-06-07-S703, from the user's derived-series picture of the roots↔coefficients
monodromy. Builds on S699p/HYP-2282 (monodromy = Galois), S699l (FTA n↔n+1 duality), THM-403
(cyclotomic worry-set = round tournament), THM-406 (Vitali wall), THM-420/430 (witness hierarchy),
S699h (A₅ icosahedral unit-distance).

## The dictionary (S699l / S699p)

A monic degree-`n` polynomial has `n+1` coefficients (the "+1" = the leading/normalisation, or dually
the constant = the `z=0` root) and `n` roots. The map *coefficients ↦ roots* is a branched cover; a
loop in coefficient space around the **discriminant** (two roots colliding) lifts to a **swap of two
roots** — the monodromy, which equals the **Galois group**, generically `S_n`. Coefficients
(symmetric functions) are *fixed* by the deck action; the `n`-root **fiber** is permuted. Repo
dictionary: **`n` roots ↔ `n` runners ↔ `n` tournament vertices**; loop-swap ↔ Galois transposition ↔
arc/runner crossing (S699p).

## The theorem

> **(1) Derived-series threshold (classical; recomputed exact).** The largest `k` with
> `S_n^{(k)} ≠ 1` is
> ```
>    n=2: 0   n=3: 1   n=4: 2   (= n−2, solvable, derived length n−1)
>    n≥5: ∞   — the derived series STABILISES at A_n ≠ 1 (perfect), unsolvable.
> ```
> Verified orders: `S₂[2,1]`, `S₃[6,3,1]`, `S₄[24,12,4,1]`, `S₅[120,60,60]`, `S₆[720,360,360]`.
> This is the user's "quadratic↔swap, cubic↔single, quartic↔double, quintic↔triple-and-above
> commutator," and **Abel–Ruffini**: no radical formula for `n≥5`.
>
> **(2) The 5-point cause (PROVED).** `A_n` is perfect for `n≥5` because every 3-cycle is a commutator
> of two 3-cycles **sharing exactly one point**, and `3+3−1 = 5` points are required for two 3-subsets
> to meet in exactly one element. Verified: the number of such triangle-pairs is `0` for `n=3,4` and
> jumps to `15` at `n=5`; `⟨(123),(345)⟩ = A₅` (order 60), `[A₅,A₅]=A₅`, and `[(123),(345)] = (032)`
> exhibits a 3-cycle as a commutator. **The "5" of the quintic is the overlap number `3+3−1`.**
>
> **(3) Tournament realisation (VERIFIED).** A 3-cycle = a **cyclic triangle** of a tournament. The
> smallest tournament carrying two cyclic triangles sharing exactly one vertex is on **5 vertices**;
> the round tournament `C₅` has `5/10` cyclic triples and realises the overlap. So the quintic
> obstruction *is* the appearance of two overlapping cyclic triangles, first at `C₅`.
> **`C₅` is exactly the LRC `n=5` cyclotomic worry-set witness** (THM-403: round tournament on the
> 5th roots of unity) — the smallest "rock-paper-scissors-squared" structure.
>
> **(4) The Vitali-wall mirror (ESTABLISHED).** Abel–Ruffini "the derived series never reaches `1`"
> is the group-theoretic mirror of **THM-406**'s LRC Vitali wall: the covering-depth inclusion–
> exclusion cancels **to all orders**, so there is **no finite Bonferroni certificate** of loneliness.
> Both are statements that a *finite-depth tower fails*: finite commutator depth (radicals) on the
> Galois side ↔ finite inclusion–exclusion depth (Bonferroni) on the LRC side. The threshold is the
> same algebraic fact — a **perfect** (depth-∞) subobject: `A_n` for the group, the all-orders
> overlap for the measure.

## Proofs / verification

**(1)** Exact permutation-group closure of the derived series (`…s703.py`); matches the classical
`S_n` derived series (`S_n'=A_n`; `A_4'=V_4`, `V_4'=1`; `A_n`simple⟹`A_n'=A_n` for `n≥5`). ∎
**(2)** Two 3-subsets of `[n]` with `|X∩Y|=1` have `|X∪Y|=5`, impossible for `n<5`; the commutator
`[(abc),(cde)]` is a 3-cycle, and 3-cycles generate `A_n`, so `A_n=[A_n,A_n]` for `n≥5`. Counts
verified `0,0,15,90,315,840` for `n=3..8`. ∎
**(3)** Direct enumeration of cyclic triples of `C_n` and the overlap (`…s703.py`). ∎
**(4)** THM-406 (all-orders cancellation, no finite Bonferroni) imported; the mirror is an identified
correspondence, not a new computation.

## Scope / honesty

- (1)(2)(3) are proved/verified (and (1) is classical Galois/Abel–Ruffini). (4) is an established
  *analogy* between two known results (Abel–Ruffini ↔ THM-406), not a new equivalence theorem.
- The payoff conjecture — that the **LRC witness hierarchy** (clock `1/k` ⊃ shell `m/(2n−1)` ⊃
  pair-sum `m/(a+b)`, THM-420/430) is a **radical/solvable tower**, that a config is tower-certifiable
  iff its local witness-monodromy is **solvable**, and that the residual `R(n)` (S700) is the
  **perfect/unsolvable core** — is **HYP-2303**, motivated here, *not* proved. In particular this does
  **not** resolve any open LRC/HN case.
- **The inversion to flag:** in Galois theory *solvable = easy* (formula exists); in LRC the
  solvable/cyclotomic worry-set (abelian monodromy, THM-403) is the **tight/hard** case (`M=1/n`
  exactly), while the "unsolvable"/generic configs are **loose**. Rigidity (cyclotomic = solvable) is
  what pins `M` to the floor; see the reflection.

**Artifacts:** `04-computation/galois_solvability_tower_s703.py` (+`.out`). Reflection
`07-reflections/the-solvability-tower-galois-lrc-icosahedron-s703.md`. New: **HYP-2303**. Builds on
S699p/HYP-2282, S699l, THM-403, THM-406, THM-420/430, S699h.

---

## ADDENDUM (monad-explorer-2026-06-07-S6) — the "5-point cause" localizes as `15·C(n,5)`, and the per-5-set `15` IS the icosahedron (PROVED/VERIFIED)

This sharpens claims (2) and (5). It uses the concrete Klein/icosahedral handle (5) asked for.
All counts verified in `04-computation/icosahedral_fifteen_monad_s6.py` (+`.out`).

> **(2′) Localization (PROVED).** The number of overlapping cyclic-triangle pairs on `[n]`
> (unordered pairs of 3-subsets `{X,Y}` with `|X∩Y|=1`) is exactly
> ```
>     15 · C(n,5)        (THM-436's 0,0,15,90,315,840 for n=3..8 = 15·C(n,5))
> ```
> and the **oriented** count (a cyclic orientation on each triangle) is `60 · C(n,5)`. So the
> `5`-point cause is *literally local*: every overlapping pair lives in a unique 5-subset, and **each
> 5-subset contributes exactly `15` unoriented (`60 = |A₅|` oriented) overlapping pairs.** The
> "engine of perfectness" is assembled 5-points-at-a-time, one copy of `A₅` per 5-set.
>
> **(2″) Commutator type-uniformity (PROVED) + a refuted over-reading.** For *every* one of the `60`
> oriented overlapping pairs on a 5-set, the commutator `[σ_X, σ_Y]` (where `σ_X, σ_Y` are the chosen
> 3-cycles) has cycle-type a **3-cycle** — never the identity, never a double-transposition, never a
> 5-cycle (robust: "is a 3-cycle" is conjugation- and inversion-invariant, so independent of order
> convention). The `60` commutators are **onto all `20` three-cycles**. Since the commutator is supported
> on `X∪Y` (5 points), this holds for all `n` by locality. **So "`A_n` perfect" is realized not by one
> witnessing commutator (`[(123),(345)]=(032)`) but uniformly in type: every overlapping triangle pair
> regenerates a 3-cycle, covering all `A_n` generators.**
> *Negative result (same session, `…flag_fibers_monad_s6.py`):* the covering is `3`-to-`1` **only on
> average** — the fibers are **non-uniform**, sizes `{2 (×3), 3 (×14), 4 (×3)}` (sum `60` over `20`
> three-cycles). Hence the tempting reading "`60 = 20·3 =` icosahedral face-vertex flags, flag→face
> uniformly 3-to-1" is **REFUTED**: the commutator covering is *not* the icosahedral flag incidence;
> the `60 = |A₅|` is the group order, not a flag count. (The `15`-bijection (5′) and the class↔axis
> dictionary are unaffected — those are exact.)
>
> **(5′) The `15` is the icosahedron (PROVED group side + classical geometry).** There is a *canonical*
> (choice-free) bijection
> ```
>   {overlapping pairs {X,Y} on a 5-set}  ⟷  {involutions (ab)(cd) of A₅}
>      shared vertex v = X∩Y          ⟷   fixed point f
>      off-pairs {X\v, Y\v}           ⟷   the two transpositions {a,b},{c,d}
> ```
> verified as set-equality of signatures `(v, {pair,pair})` (both sides = the `15` "fixed point + perfect
> matching of the other 4" data). And `A₅`'s nontrivial conjugacy classes match the icosahedral axes
> exactly (computed): class sizes `{15, 20, 24}` = `{#2-axes·1, #3-axes·2, #5-axes·4}` =
> `{15·1, 10·2, 6·4}`, so the icosahedral invariant degrees `{6,10,15}` (binary forms `{12,20,30}` =
> `{vertices, faces, edges}`, Klein 1884) are exactly the **axis counts** `{#5-axes, #3-axes, #2-axes}`.
> The `15` of THM-436(2) = `15` involutions = `15` two-fold axes = `15` edge-axes = the degree of the
> Klein edge-form `T` (deg 30 = 2·15) whose relation `T²+H³ = 1728 f⁵` is the quintic-solving icosahedral
> equation. **The "5" of Abel–Ruffini is the icosahedron, and its `15` is the overlapping-triangle count.**

> **Tournament shadow (VERIFIED).** The round tournament `C₅` realises only `5` of the `10` 3-subsets as
> cyclic triangles; among those `5`, exactly `5` pairs share one vertex and `5` share two (`0` are
> disjoint). So `C₅` carries `5` of the `15` abstract overlapping pairs. The three icosahedral axis
> classes have clean tournament meanings on `C₅`: the **5-cycles** = the rotation `(01234)` = the round
> tournament's defining Hamiltonian symmetry; the **3-cycles** = its cyclic triangles; the
> **involutions** = the overlapping-pair structure. The whole of `A₅`/the icosahedron is generated by
> the combinatorics of the single tournament `C₅`.

**Honesty:** the group-theory and the combinatorial counts/bijection/commutator-uniformity are computed
and exact (`…fifteen_monad_s6.py`); the icosahedral geometry (`V,E,F = 12,30,20`; invariant degrees
`{2,6,10,15}`; the `T²+H³=1728 f⁵` relation) is classical (Klein). The forward `(2,3,5)`-frontier
resonance — that the icosahedron's three axis-orders are the repo's three carry-prime sites — is
**HYP-2305**, not proved here. **Artifacts (addendum):**
`04-computation/icosahedral_fifteen_monad_s6.py` (+`05-knowledge/results/…out`); reflection
`07-reflections/the-icosahedral-fifteen-s6.md`; new **HYP-2305**.
