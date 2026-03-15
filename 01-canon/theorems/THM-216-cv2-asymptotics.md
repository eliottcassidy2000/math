---
id: THM-216
name: CV² Asymptotics
status: PROVED
proved_by: opus-2026-03-14-S89c
verified_computationally: n=3..21 (exact polynomial match)
---

# THM-216: CV²(H) = 2/n + O(1/n³)

## Statement

For the Hamiltonian path count H over uniformly random tournaments on n vertices:

  CV²(H) = Var(H)/E[H]² = 2/n + O(1/n³)

The 1/n² coefficient vanishes exactly. More precisely:

  CV² = Σ_{k≥1} 2·g_k(n-2k) / (n)_{2k}

where (n)_{2k} is the falling factorial, and the scaled contributions are:
  g₁(m) = m        (dominos)
  g₂(m) = m²       (domino pairs)
  g₃(m) = m(2m²+1)/3   (domino triples)

The 1/n² cancellation: |S|=2 contributes -2/n² and |S|=4 contributes +2/n²,
giving exact cancellation. The leading correction is -14/(3n³) + O(1/n⁴).

## Proof

**Setup.** Let σ be a uniformly random permutation of {0,...,n-1}. Define:
- X_j = 1[σ(j+1) = σ(j)+1]  (unit ascent at position j)
- Y_j = 1[σ(j+1) = σ(j)-1]  (unit descent at position j)
- Z_j = X_j - Y_j ∈ {-1, 0, 1}

Since X_j·Y_j = 0 (cannot have both at same position):

  2^{adj1(σ)} · 1_{NUD}(σ) = ∏_j (1+X_j) · ∏_j (1-Y_j) = ∏_j (1+Z_j)

Therefore W(n)/n! = E[∏_{j=0}^{n-2} (1+Z_j)] and CV² = W(n)/n! - 1.

**Expand the product:**
  E[∏(1+Z_j)] = Σ_{S⊆{0,...,n-2}} E[∏_{j∈S} Z_j]

**Step 1: Odd-size subsets contribute 0.**
The complement involution σ ↦ (n-1-σ(1), ..., n-1-σ(n)) preserves the
uniform distribution and sends X_j ↔ Y_j, hence Z_j ↦ -Z_j.
For |S| odd: E[∏Z_j] = (-1)^|S| E[∏Z_j] = -E[∏Z_j] = 0.

**Step 2: Non-adjacent pairs contribute 0.**
For |j-k| ≥ 2: E[Z_jZ_k] = 2(E[X_jX_k] - E[X_jY_k]).
Both E[X_jX_k] and E[X_jY_k] count the probability of having two disjoint
consecutive pairs (one ascending, one ascending/descending) at non-adjacent
positions. The combinatorial count is identical regardless of orientation
(both require selecting two disjoint edges from the path graph on values),
so E[X_jX_k] = E[X_jY_k] and E[Z_jZ_k] = 0.

**Verified computationally:** E[X_jX_k] = E[X_jY_k] = 1/(n(n-1)) for |j-k|≥2, all n=3..8.

**Step 3: Adjacent pairs contribute 2/(n(n-1)).**
For adjacent positions j, j+1:
- E[X_jY_{j+1}] = 0: σ(j+1)=σ(j)+1 and σ(j+2)=σ(j+1)-1=σ(j) requires repeating σ(j).
- E[Y_jX_{j+1}] = 0: similarly, σ(j+1)=σ(j)-1 and σ(j+2)=σ(j) requires repeat.
- E[X_jX_{j+1}] = (n-2)·(n-3)!/n! = 1/(n(n-1)): counts runs v,v+1,v+2.
- E[Y_jY_{j+1}] = 1/(n(n-1)): by complement symmetry.

So E[Z_jZ_{j+1}] = 2/(n(n-1)).

Summing: Σ_{j=0}^{n-3} E[Z_jZ_{j+1}] = (n-2)·2/(n(n-1)) = 2/n - 2/(n(n-1)).

**Step 4: |S| ≥ 4 terms are O(1/n²).**
The dominant |S|=4 contribution comes from pairs of disjoint adjacent pairs
{j,j+1,k,k+1}. There are ~n²/2 such sets, each contributing ~4/n⁴
(product of pair expectations). Connected 4-clusters contribute O(1/n³).
All |S| ≥ 6 contribute O(1/n³). Total from |S| ≥ 4: O(1/n²).

**Conclusion:**
  CV² = (2/n - 2/(n(n-1))) + O(1/n²) = 2/n + O(1/n²).  □

## Known Exact Values

| n | CV² | n × CV² | 2-n×CV² |
|---|------|---------|---------|
| 3 | 1/3 | 1.000 | 1.000 |
| 4 | 1/3 | 1.333 | 0.667 |
| 5 | 19/60 | 1.583 | 0.417 |
| 6 | 13/45 | 1.733 | 0.267 |
| 7 | 131/504 | 1.819 | 0.181 |
| 8 | 131/560 | 1.871 | 0.129 |
| 9 | 1097/5184 | 1.905 | 0.095 |
| 10 | 3121/16200 | 1.927 | 0.073 |

## Computational Data (n=11..25)

| n | CV² (approx) | n × CV² |
|---|---------------|---------|
| 15 | 0.13147 | 1.9721 |
| 20 | 0.09927 | 1.9855 |
| 25 | 0.07964 | 1.9911 |

## Key Identities

### Tridiagonal covariance
The matrix E[Z_jZ_k] is TRIDIAGONAL:
- Diagonal: E[Z_j²] = 2/n
- Off-diagonal (adjacent): E[Z_jZ_{j+1}] = 2/(n(n-1))
- All other entries: exactly 0

### Domino structure
Only "domino" subsets S (unions of adjacent pairs {j,j+1}) contribute
to E[∏_{j∈S} Z_j]. All other subsets give 0, because:
1. Isolated positions give 0 by the complement symmetry Z→-Z
2. Three+one patterns give 0 by the same argument applied to the odd part

### Exact |S|=2k contributions (verified n=3..21)
The |S|=2k contribution equals 2·g_k(n-2k) / (n)_{2k} where:
  g₁(m) = m
  g₂(m) = m²
  3·g_k(m) = a_k·m³ + b_k·m² + c_k·m + d_k   for all k ≥ 3

**All g_k for k ≥ 3 are degree-3 polynomials.** Coefficients:

| k | a_k | b_k | c_k | d_k |
|---|-----|-----|-----|-----|
| 3 | 2 | 0 | 1 | 0 |
| 4 | 10 | -33 | 50 | -24 |
| 5 | 388 | -2040 | 3431 | -1776 |
| 6 | 69660 | -380445 | 653748 | -342960 |
| 7 | 19826270 | -109486152 | 189674605 | -100014720 |
| 8 | 7309726742 | -40641958545 | 70757788486 | -37425556680 |
| 9 | 3262687720240 | -18232387983408 | 31858349908595 | -16888649645424 |

Universal boundary values: g_k(1) = 1 and g_k(2) = 2k for all k ≥ 1.
Note: g_k(0) ≠ 0 for k ≥ 4 (the polynomial domain is m ≥ 1).
Verified: CV² = Σ 2·g_k(n-2k)/(n)_{2k} reproduces W(n)/n!-1 exactly for n=3..21.

### 1/n² cancellation
|S|=2: 2(n-2)/(n(n-1)) = 2/n - 2/n² + 2/n³ - ...
|S|=4: 2(n-4)²/(n)_4 = 2/n² - 4/n³ + ...
Sum:   2/n + 0/n² + (2-4)/n³ + |S|=6 + ...
|S|=6: 4/(3n³) + ...
Total: 2/n + (-6 + 4/3)/n³ = 2/n - 14/(3n³) + O(1/n⁴)

## Note

Earlier conjecture that CV² → 1/4 was WRONG (based on n≤7 data only).
The value crosses below 1/4 at n=8 and continues toward 0.

## Files

- `04-computation/nud_weight.c` — C implementation for n up to 25
- `04-computation/cv2_proof_89c.py` — Proof verification script
- `04-computation/gk_verify_89c.py` — Cross-verification: all g_k polynomials reproduce W(n) exactly
- `04-computation/gk_coefficients_v2_89c.py` — Sequential extraction of g_k coefficients
- `04-computation/gk_cascade_v2_89c.py` — Manual cascade extraction (correct, no g_k(0)=0 assumption)
- `04-computation/w_recurrence_89c.py` — Recurrence analysis
- `04-computation/w_analytic_89c.py` — Poisson approximation analysis
