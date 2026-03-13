# THM-168: Lambda Pair-Coverage Completeness

**Status:** VERIFIED (exhaustive n=5,6; n=7 regular; n=7 general pending)
**Session:** kind-pasteur-2026-03-13-S61

## Statement

### Definition

For a tournament T on n vertices, define the **pair-coverage function**:

    lambda_{uv}(T) = #{3-cycle vertex sets containing both u and v}

The **lambda histogram** is the sorted multiset {lambda_{uv} : {u,v} in C(n,2)}.

### Main result

For the following score classes, the lambda histogram uniquely determines H(T):

- **n=5, ALL score classes:** Lambda determines H. (Exhaustive, 1024 tournaments)
- **n=6, score (2,2,2,3,3,3):** 4 lambda classes map to 4 distinct H values.
  Lambda determines H AND c5_dir independently.
- **n=7, regular (3,3,3,3,3,3,3):** 3 lambda classes map to 3 distinct H values.
  Lambda determines H via quadratic formula (THM-167).

### Full n=7 results (exhaustive, 2097152 tournaments)

- **288 distinct (score, lambda) classes**
- **233 complete** (lambda determines H): 80.9%
- **55 ambiguous** (lambda does NOT determine H): 19.1%
- **39/59 score classes** are lambda-complete (66.1%)
- **20/59 score classes** are lambda-incomplete (33.9%)

Complete score classes include: transitive, near-transitive, regular, and many others.

### Failure cases

- **n=6:** 3 ambiguous (score, lambda) classes. Extended lambda (3+5 cycle coverage) still ambiguous for 2.
- **n=7:** 55 ambiguous lambda classes in 20 score classes. Even c5 does NOT resolve (same lambda + c5 can give different H).
- First ambiguous: score (0,2,3,3,4,4,5), lambda with specific pattern: H in {25, 29} with same c5=6.

## The Measurability Hierarchy

    Level 0 (Score): Determines c3. Local information (out-degrees).
    Level 1 (Lambda): Determines disj_33 via sum C(lambda,2) - correction.
                      For special score classes, determines H.
    Level 2 (Full):   Contains c5, c7, etc. Always determines H.

At n=7 regular, the constraint c5_dir + 2*disj_33 = 56 makes the Level 2
information (c5) redundant — it's determined by Level 1 (lambda).

At n=6 for most score classes, Level 1 doesn't determine Level 2,
but for the near-regular class (2,2,2,3,3,3) it does.

## Verified Cases (n=6, score (2,2,2,3,3,3))

| Lambda histogram (sorted) | disj | H  | c5 | Count |
|---------------------------|------|----|----|-------|
| (0,0,0,2^12)              | 4    | 45 | 6  | 240   |
| (0,1^6,2^6,3^2)           | 2    | 41 | 8  | 720   |
| (1^6,2^6,3)               | 1    | 43 | 11 | 1440  |
| (1^6,2^9)                 | 1    | 45 | 12 | 240   |

## Connection to Vitali sets

The lambda completeness boundary is analogous to measurability:
- Score = "Lebesgue measurable" (constructive from local data)
- Lambda = "Borel measurable" (requires pairwise information)
- Full tournament = "arbitrary set" (may contain non-constructive choices)

When lambda suffices to determine H, the Vitali partition collapses.
When it doesn't, the 5-cycle count acts as a "non-measurable" choice
that cannot be deduced from the 3-cycle structure alone.

## Verification Scripts

- quadratic_formula_exploration.py: Lambda completeness at n=6, n=7 regular
- vitali_overlap_hidden_structure.py: Overlap weight analysis
- lambda_completeness_boundary.py: Full boundary determination (running)
