# THM-302: Same-End Quadratic Coefficient Formula

**Status:** CONJECTURED, verified n=5,6,7 (68/68 matches)
**Discovered by:** opus-2026-04-03-S28, refined opus-2026-04-04-S3

## Statement

For two tiles (x₁,y₁) and (x₂,y₂) that share a same-end vertex (x₁=x₂ or y₁=y₂), with skips s₁ = x₁-y₁ and s₂ = x₂-y₂:

$$c_{ij} = -2^{\max(1, |s_1 - s_2| - 1)}$$

In particular:
- Adjacent tiles (|s₁-s₂| ≤ 1): c_{ij} = -2
- Tiles with skip difference d = |s₁-s₂| ≥ 2: c_{ij} = -2^{d-1}

## Proof Sketch (from OCF, THM-287)

Same-end tile pairs have exactly one directed 3-cycle passing through both tiles (the cycle that traverses the shared vertex and the two "distant" endpoints). In this cycle, the two tiles are used in **opposite** directions (one forward, one backward), giving chi coefficient = -1 and weight contribution 2·(-1) = -2.

When |s₁-s₂| ≥ 2, additional odd cycles (5-cycles, 7-cycles, etc.) pass through both tiles, each contributing -2. The total count of such cycles is 2^{|s₁-s₂|-2}, giving total: -2·2^{|s₁-s₂|-2} = -2^{|s₁-s₂|-1}.

## Verification

| n | Tested pairs | Formula matches |
|---|-------------|-----------------|
| 5 | 8 | 8/8 |
| 6 | 20 | 20/20 |
| 7 | 40 | 40/40 |
