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
