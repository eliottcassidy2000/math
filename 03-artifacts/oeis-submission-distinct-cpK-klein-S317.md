# OEIS submission draft — distinct skew-adjacency characteristic polynomials of tournaments

**Status:** ready to submit (checked NOT in OEIS 2026-07-16, query `1,2,2,6,11,50`)
**Instance:** klein-2026-07-16-S317 (data from the cospectral census, THM-924/925/931 arc)

## Sequence

**Name:** Number of distinct characteristic polynomials of the skew-adjacency matrix K = Aᵀ − A over all tournaments on n vertices (equivalently, the number of distinct symmetrized adjacency spectra).

**Data (offset 1):**
`1, 1, 1, 2, 2, 6, 11, 50`

(All eight values recomputed fresh 2026-07-16 by exact integer Faddeev–LeVerrier over all iso classes — n = 1..7 by full orbit census, n = 8 over the 6880 extension-census reps; scripts in repo, scratch runs `distinct_cpk.py`/`cpk8.py`.)

**Comments:**
- K = Aᵀ − A where A is the adjacency (dominance) matrix of a tournament; cpK is even/odd according to n, and its roots are purely imaginary.
- THE WALK RECIPROCITY (this project, THM-924): cpK is an explicit linear functional of cpA — cpK(y) = 2^{n−1}·[cpA((y−1)/2) + (−1)ⁿ·cpA((−y−1)/2)] — so a(n) also counts the distinct images of the tournament A-spectra under the symmetrization x ↦ (roots of the reciprocal average). In particular a(n) ≤ #distinct cpA(n), with dramatic collapse: distinct cpA counts are 1, 1, 2, 3, 9, 28, 168, 1523 (n = 1..8, verified same run) vs 1, 1, 1, 2, 2, 6, 11, 50 for cpK — at n = 8 the symmetrization collapses 1523 A-spectra onto 50 K-spectra (30×).
- a(n) < A000568(n) (tournaments): the skew spectrum is far from complete as an invariant; at n = 7 the 456 classes fall into just 11 skew-spectral classes, and at n = 8 the 6880 into 50.

**Crossrefs:** A000568 (tournaments), A002854 (even graphs, this project's dual object).

**Program:** Python via exact Faddeev–LeVerrier over ℚ; see the repo script `cospectral_tie_census_klein_S315.py` (census function + charpoly of K per iso class; count distinct tuples).

**Keyword suggestions:** `nonn,more,hard` (n = 9 computable from the 191,536-class census machinery in `n9_wild_hunt_klein_S317.py` with one extra pass).

## Companion sequence (also OEIS-absent as of 2026-07-16): distinct cpA

Number of distinct adjacency characteristic polynomials of tournaments on n vertices:
`1, 1, 2, 3, 9, 28, 168, 1523` (offset 1, same runs). Query `2,3,9,28,168,1523` returned no results — submit as a pair with a(n) ≤ this ≤ A000568(n) and the THM-924 reciprocity linking them.

## Submission checklist
- [x] Recompute n = 1..8 fresh (done 2026-07-16, exact; see data note). Optional: n = 9 term via one census pass.
- [ ] Register OEIS account under owner's email; paste name/data/comments/program.
- [ ] After acceptance, backlink the A-number in METAGRAPH-ATLAS section 8 and THM-924.
