# THM-167: Quadratic H Formula for Regular n=7

**Status:** VERIFIED (exhaustive, 2640 regular tournaments)
**Session:** kind-pasteur-2026-03-13-S61

## Statement

For every regular tournament T on n=7 vertices:

    H(T) = disj_33^2 - 23*disj_33 + 301

where disj_33 is the number of vertex-disjoint pairs of 3-cycle vertex sets.

### Equivalent formulations

Using the constraint c5_dir + 2*disj_33 = 56:

    c7_dir = (disj_33^2 - 23*disj_33 + 160) / 2

And in terms of pair-coverage variance:

    H = (441/4)*var_lambda^2 - (189/2)*var_lambda + 189

where var_lambda is the variance of the pair-coverage function lambda_{uv} = #{3-cycles containing both u and v}.

### Achieved values

| disj_33 | H   | c7_dir | c5_dir | var_lambda | Count | Description   |
|---------|-----|--------|--------|------------|-------|---------------|
| 7       | 189 | 24     | 42     | 0          | 240   | Paley (BIBD)  |
| 10      | 171 | 15     | 36     | 6/21       | 1680  |               |
| 14      | 175 | 17     | 28     | 14/21      | 720   |               |

### Key properties

1. **Parabola with minimum at disj=11.5:** H_min = 168.75 (not achievable by any tournament).
   Paley (disj=7) is FARTHEST from minimum, giving maximum H.

2. **BIBD characterization:** Paley T_7 has lambda=2 for ALL vertex pairs (uniform pair-coverage).
   This is the unique (7,3,2)-BIBD. The uniformity minimizes disj via Jensen's inequality.

3. **Disjointness formula:** disj = sum_{edges} C(lambda_e, 2) - 14 = (sum lambda^2 - 42)/2 - 14.

4. **Three rigid constraints:**
   - c3_dir = 14 (constant for all regular n=7)
   - c5_dir + 2*disj_33 = 56 (constraint)
   - c7_dir = (disj_33^2 - 23*disj_33 + 160)/2 (quadratic constraint)

## Proof sketch

From THM-166 (OCF decomposition): H = 1 + 2*alpha_1 + 4*alpha_2

For regular n=7:
- alpha_1 = c3_dir + c5_dir + c7_dir (total directed odd cycles)
- alpha_2 = disj_33 (vertex-disjoint cycle pairs; 5+3 > 7 so only 3-3 pairs)
- c3_dir = 14 (constant)

Substituting:
- H = 29 + 2*c5_dir + 2*c7_dir + 4*disj_33
- Using c5_dir = 56 - 2*disj_33: H = 141 + 2*c7_dir

So the formula reduces to showing c7_dir = (disj_33^2 - 23*disj_33 + 160)/2.

This second constraint (c7 quadratic in disj) appears to follow from the
rigid combinatorial structure of regular 7-vertex tournaments, where the
3-cycle pair-coverage lambda completely determines the tournament's cycle
spectrum (up to permutation).

## Connections

- **THM-165** (Coefficient-2): The linear term 2*c7_dir comes from Hamiltonian cycles each contributing 2 to H.
- **THM-166** (OCF Decomposition): The full OCF reduces to 3 free parameters at n=7 regular.
- **THM-027** (BIBD Maximization): Paley/BIBD minimizes disj but maximizes total cycles, giving max H.
- **HYP-746** (Paley Uniformity): Paley = unique BIBD with lambda=2 uniform.

## Verification

- overlap_mechanism_deep.py: Full derivation and verification of quadratic formula
- quadratic_formula_exploration.py: Pair-coverage analysis, Vitali connection
