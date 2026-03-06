# THM-021: Elementary proof of real-rootedness via Turán's theorem

**Status:** PROVED for n ≤ 8
**Certainty:** 5
**Author:** opus-2026-03-06-S15
**Dependencies:** Turán's theorem (classical), AM-GM inequality
**Relation to THM-020:** Alternative proof of the same result. THM-020 uses claw-freeness + Chudnovsky-Seymour (2007); this proof uses only classical combinatorics.

---

## Statement

For any tournament T on n ≤ 8 vertices, the independence polynomial I(Ω(T), x) has all real (non-positive) roots.

## Proof

Write I(Ω(T), x) = 1 + α₁x + α₂x² + ··· + α_d x^d where d = α(Ω(T)) ≤ ⌊n/3⌋.

**Case n ≤ 5:** d ≤ 1. A polynomial 1 + α₁x trivially has one real root at x = -1/α₁. ∎

**Case n ∈ {6, 7}:** d ≤ 2. Write I = 1 + α₁x + α₂x². This has all real roots iff the discriminant α₁² − 4α₂ ≥ 0.

At n ≤ 7, the only vertex-disjoint odd-cycle pairs are pairs of 3-cycles (since 3 + 5 = 8 > n for any other combination). Let G be the graph whose vertices are the c₃ directed 3-cycles of T, with an edge between two cycles iff they are vertex-disjoint.

**Key observation:** Three pairwise disjoint triples require 3 × 3 = 9 > n vertices. Therefore ω(G) ≤ 2, i.e., G is **triangle-free**.

By Turán's theorem for triangle-free graphs: |E(G)| ≤ ⌊c₃²/4⌋.

Since α₂ = |E(G)| (disjoint cycle pairs) and α₁ ≥ c₃ (α₁ also counts 5- and 7-cycles):

**α₁² − 4α₂ ≥ c₃² − 4⌊c₃²/4⌋ ≥ 0.** ∎

**Case n = 8:** d ≤ 2. Now vertex-disjoint pairs include:
- **(3,3)-pairs:** Two disjoint 3-cycles (using 6 of 8 vertices). Count: α₂^{33}.
- **(3,5)-pairs:** A 3-cycle and a disjoint 5-cycle (using all 8 vertices). Count: B.
- Other combinations (5+5, 3+7, etc.) require > 8 vertices: impossible.

So α₂ = α₂^{33} + B.

**Bound on α₂^{33}:** Same Turán argument (ω ≤ 2 since 9 > 8). So 4α₂^{33} ≤ c₃².

**Bound on B:** Each (3,5)-pair consists of a directed 3-cycle S and a directed 5-cycle on V\S. Each 5-cycle has a unique complement triple, so B ≤ c₅.

**Combining:** 4α₂ = 4α₂^{33} + 4B ≤ c₃² + 4c₅.

We verify α₁² ≥ c₃² + 4c₅ where α₁ = c₃ + c₅ + c₇:

- If c₃ = 0: B = 0 (no cyclic triple), so α₂ = 0 and the bound is trivial.
- If c₃ = 1: α₂^{33} = 0. Then 4α₂ = 4B ≤ 4c₅. And α₁² = (1 + c₅ + c₇)² ≥ (1 + c₅)² = 1 + 2c₅ + c₅² ≥ 4c₅ since (c₅ − 1)² ≥ 0. ∎
- If c₃ ≥ 2: α₁² ≥ (c₃ + c₅)² = c₃² + 2c₃c₅ + c₅² ≥ c₃² + 4c₅ + c₅² ≥ c₃² + 4c₅ since 2c₃ ≥ 4. ∎

---

## Tightness

The bound is **tight** at n = 6 and n = 7: tournaments with exactly two directed 3-cycles on complementary vertex sets (and no other odd cycles) achieve discriminant = 0, giving a double root at x = −1.

- **n = 6:** 80 tournaments (1 isomorphism class), score sequence (4,4,4,1,1,1). I(Ω,x) = (1+x)².
- **n = 7:** 1680 tournaments, score sequence (6,4,4,4,1,1,1). I(Ω,x) = (1+x)².
- **n = 8:** Tightest observed ratio α₁²/(4α₂) = 9/8 = 1.125, at c₃=3, α₂=2.

## Verification record

| n | Method | Tournaments | Failures | Min ratio |
|---|--------|-------------|----------|-----------|
| 5 | Exhaustive | 1,024 | 0 | ∞ (degree ≤ 1) |
| 6 | Exhaustive | 32,768 | 0 | 1.0000 (tight) |
| 7 | Exhaustive | 2,097,152 | 0 | 1.0000 (tight) |
| 8 | Sampled | 2,000 | 0 | 1.1250 |

## Significance

1. **Elementary:** Uses only Turán's theorem and AM-GM, versus the deep structure theory of claw-free graphs (Chudnovsky-Seymour).
2. **Self-contained:** Does not require computing or classifying the graph structure of Ω(T).
3. **Tight:** Achieves equality at the discriminant boundary, showing the proof is optimal for degree-2 polynomials.
4. **Structural insight:** The proof reveals WHY real-rootedness holds at small n: vertex-disjoint odd-cycle pairs are constrained by the pigeonhole principle (few vertices available), and Turán's triangle-free bound matches the quadratic discriminant exactly.

## Limitation

For n ≥ 9, α(Ω(T)) ≥ 3 and the independence polynomial has degree ≥ 3. The discriminant approach no longer suffices — one needs to control higher-order terms. The Newton inequalities (α_k² ≥ (k+1)/(k) · α_{k-1} · α_{k+1}) are necessary conditions for real roots, and ULC (ultra-log-concavity, T071) is computationally verified but not proved.

## See also

- THM-020 (Chudnovsky-Seymour proof of same result)
- THM-019 (Omega perfectness fails at n=8)
- T071 (ultra-log-concavity conjecture)
- OPEN-Q-015 (real roots for general n)
