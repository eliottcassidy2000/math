# THM-208: H = 1 + 2t₃ for n ≤ 4

**Status:** VERIFIED (exhaustive, n=3,4), PROOF NEEDED
**Found by:** opus-2026-03-14-S89b
**Verified in:** `04-computation/H_6t3_regression_89b.py`

## Statement

For tournaments on n ≤ 4 vertices:

**H(T) = 1 + 2·t₃(T)**

where t₃(T) is the number of directed 3-cycles in T.

This is **exact** with zero residual for every tournament.

## Key Facts

- **n=3**: t₃ ∈ {0, 1}. H ∈ {1, 3}. Formula gives 1, 3. ✓
- **n=4**: t₃ ∈ {0, 1, 2}. H ∈ {1, 3, 5}. Formula gives 1, 3, 5. ✓
  - Note: t₃ = 3 cannot occur at n=4 (no tournament on 4 vertices has exactly 3 of its 4 triangles as 3-cycles).
  - Note: t₃ = 4 cannot occur at n=4 either (would require all C(4,3)=4 triangles to be 3-cycles, which requires scores (1,1,2,2) — but these exist with t₃=2, not 4).
- **n=5**: The formula breaks — t₃=4 gives H ∈ {11, 13, 15}. The **5-cycle count** t₅ discriminates: t₅=1→H=11, t₅=2→H=13, t₅=3→H=15.
- **n≥5**: The residual ε = H - 1 - 2t₃ is always **even** and **non-negative**.

## The Residual Structure

For n ≥ 5, define ε(T) = H(T) - 1 - 2·t₃(T).

| n | ε range | Mean(ε) | ε=0 fraction |
|---|---------|---------|--------------|
| 3 | [0,0]   | 0       | 1.0000       |
| 4 | [0,0]   | 0       | 1.0000       |
| 5 | [0,6]   | 1.5     | 0.4688       |
| 6 | [0,28]  | 11.5    | 0.1172       |

The residual ε captures the contribution of **higher-order cycles** (5-cycles, 7-cycles, etc.) to the Hamiltonian path count.

## Regression Slopes

The exact regression slopes b = Cov(H, t₃)/Var(t₃):

| n | Slope b | R²    |
|---|---------|-------|
| 3 | 2       | 1     |
| 4 | 2       | 1     |
| 5 | 3       | 18/19 |
| 6 | 6       | 12/13 |

## Connection

- The slope 2 for n ≤ 4 means each 3-cycle contributes exactly 2 additional Hamiltonian paths.
- The slope 6 = |S₃| for n=6 connects to the triangle symmetry group.
- The slope 3 = Φ₃(1) for n=5 connects to the cyclotomic polynomial.
