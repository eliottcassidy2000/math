# POS Class Synthesis — The Self-Complementary Score as OCR Residual Source

**Session:** kind-pasteur-2026-03-21-S16b

## The Core Discovery

At n=5, the ENTIRE OCR residual (1/19 of Var(H)) comes from a SINGLE score class: **(1,2,2,2,3)**. This class is:

1. **Self-complementary** (POS): s_i + s_{n-1-i} = n-1 for all i.
2. The **only ambiguous class** at n=5: H ∈ {11, 13, 15} within a single score class.
3. The class with the **lowest nonzero S₂ = 2** (closest to regular).
4. Contains the **"cone over C₃"** tournament (the canonical non-regular maximizer).

## What Determines H Within the Class

H = 9 + 2·c5_dir, where c5_dir ∈ {1, 2, 3} is the directed Hamiltonian cycle count.

The distinction is controlled by a **single binary choice**: does the score-1 vertex beat the score-3 vertex?

| c5_dir | H | Count | Weak beats Strong? |
|--------|---|-------|--------------------|
| 3 | 15 | 40 | **Yes** |
| 2 | 13 | 120 | No |
| 1 | 11 | 120 | No |

When the weakest vertex (score 1) beats the strongest (score 3), the tournament has maximum Hamiltonian cycle count (c5_dir = 3) and maximum H = 15 within the class. This is the **"upset creates cycles"** principle: the single "wrong-way" arc forces additional directed 5-cycles.

The c5_dir = {1, 2} distinction within the "No" group is determined by the cycle structure among the three middle vertices (all score 2).

## The ANOVA vs Regression Gap

| Measure | Value | Interpretation |
|---------|-------|----------------|
| R²(S₂ regression) | 18/19 ≈ 0.947 | Linear prediction from S₂ |
| R²(ANOVA, score class) | 129/133 ≈ 0.970 | Optimal within-class prediction |
| Gap | 76/2527 ≈ 0.030 | Nonlinearity of H vs S₂ |

The S₂ regression loses 3% to nonlinearity (H is not exactly linear in S₂). But the S₂ collisions (three score classes all with S₂ = 8) lose NOTHING because those classes all have the same mean H.

## Generalization to n=6

At n=6:
- 6/22 score classes are self-complementary
- 3/6 SC classes are ambiguous (multiple H values)
- 6/16 non-SC classes are ambiguous
- **SC classes contribute 73.8% of the within-class variance** (despite being only 27% of classes)
- The worst class is SC (1,2,2,3,3,4): 6 distinct H values, spread 14

## The Universal Principle

**Self-complementary score sequences have the MOST internal freedom.** Their symmetry constraint (s_i + s_{n-1-i} = n-1) ensures that the "top half" and "bottom half" of the score distribution are perfectly balanced, but the CONNECTIONS between these halves are not score-determined. The number of directed 5-cycles (and hence H) depends on HOW the mid-ranked vertices connect to each other.

Non-SC score sequences have an asymmetry that CONSTRAINS cycle counts more tightly. If you know who wins a lot and who loses a lot, and the asymmetry is pronounced, there's less room for variation in cycle structure.

This is why **regularity (S₂ = 0) determines H less precisely than irregularity (S₂ >> 0)**:
- Regular scores give maximum freedom for cycle arrangement → multiple H values
- Highly irregular scores constrain cycle structure → unique H value

The OCR residual measures this freedom. It comes primarily from the most symmetric (SC) score classes.

## Connection to the Paley Maximizer

Paley tournaments are simultaneously:
- Regular (S₂ = 0)
- Self-complementary
- Maximizers of H within their score class

The SC Maximizer Theorem (THM-255) says the maximum H within each SC score class is achieved by an SC tournament. Combined with the above analysis: the SC score classes are the ones with the most within-class variation, and the SC tournaments within those classes achieve the maximum.

**The Paley tournament is the maximizer because it occupies the maximally symmetric position within the maximally symmetric score class.** Its H = 189 at n=7 comes from having the most directed cycles (c3 = 14, c5_dir and c7_dir all maximized among regular tournaments).

## Open Questions

1. **Does the SC residual fraction (73.8% at n=6) increase or decrease with n?**
2. **Is the "weak beats strong" binary at n=5 the simplest example of a general pattern?**
3. **Can we compute the exact within-class Var(H) for the regular class at each n?**
   At n=7: regular H ∈ {171, 175, 189} with counts {1680, 720, 240}.
   What determines the distinction?
4. **Is there a "self-complementary OCR" — the R² restricted to SC score classes?**
