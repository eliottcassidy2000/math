# THM-059: Universal Master Polynomial Theorem

**Type:** Theorem (computational discovery, verified n=4..9)
**Certainty:** 3 -- VERIFIED computationally (22/22 cases across 6 values of n, both parities)
**Status:** VERIFIED, algebraic proof pending
**Added by:** opus-2026-03-06-S30
**Tags:** #W-polynomial #hierarchy #central-factorial #master-sequence #universal

---

## Statement

**Theorem.** For a tournament T on n vertices, let W(r) = sum_P prod_{i=0}^{n-2} (r + s_i) be the weighted path polynomial (s_i = A[p_i, p_{i+1}] - 1/2). Decompose:

W(r) = C_0(r) + sum_I C_I(r) * I(T)

where the sum is over all independent set types I in the odd-cycle conflict graph Omega(T), and C_I(r) is the per-invariant r-polynomial.

Then:

**(a)** C_I(r, n) = 2^{parts(I)} * F_f(r)

where parts(I) is the number of disjoint odd cycles in independent set type I, and f = (n-1) - 2*|pi_I| is the **free position count** (|pi_I| = sum of block lengths in the position partition).

**(b)** The universal master polynomials F_j(r) are determined by the **central factorial number triangle** (OEIS A036969) via:

F_{2k}(r) = sum_{j=0}^k b_{k,j} * (2j+1)! * (r^2 - 1/4)^j

F_{2k+1}(r) = r * sum_{j=0}^k b_{k,j} * (2j+2)! * (r^2 - 1/4)^j

where b_{k,j} satisfies the recurrence:

b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j},  with b_{k,0} = 1, b_{0,0} = 1.

**(c)** Properties of F_j:
- F_j(1/2) = 1 for all j (from OCF: H = W(1/2) = I(Omega(T), 2))
- Leading coefficient of F_j is (j+1)!
- F_j has parity j: only r^{2m} terms for even j, only r^{2m+1} for odd j
- At odd n, free positions are always even; at even n, always odd

---

## The Central Factorial Number Triangle

| k\j | 0 | 1 | 2 | 3 | 4 |
|-----|---|---|---|---|---|
| 0   | 1 |   |   |   |   |
| 1   | 1 | 1 |   |   |   |
| 2   | 1 | 5 | 1 |   |   |
| 3   | 1 | 21| 14| 1 |   |
| 4   | 1 | 85|147| 30| 1 |

Recurrence: b_{k,j} = b_{k-1,j-1} + j^2 * b_{k-1,j}

Column properties: b_{k,0} = 1 always, b_{k,k} = 1 always, b_{k,1} = (4^k - 1)/3.

---

## The Master Sequence

| j | F_j(r) | Leading | Evaluation at r=1/2 |
|---|--------|---------|---------------------|
| 0 | 1 | 1 = 1! | 1 |
| 1 | 2r | 2 = 2! | 1 |
| 2 | -1/2 + 6r^2 | 6 = 3! | 1 |
| 3 | -4r + 24r^3 | 24 = 4! | 1 |
| 4 | 1 - 30r^2 + 120r^4 | 120 = 5! | 1 |
| 5 | 17r - 240r^3 + 720r^5 | 720 = 6! | 1 |
| 6 | -17/4 + 231r^2 - 2100r^4 + 5040r^6 | 5040 = 7! | 1 |
| 7 | (predicted) | 40320 = 8! | 1 |
| 8 | 31 - 2640r^2 + 40320r^4 - 211680r^6 + 362880r^8 | 362880 = 9! | 1 |

---

## Verification Table

All 22 cases verified with zero error:

| Invariant | n | parts | |pi| | free | OCF wt | Predicted = Actual |
|-----------|---|-------|------|------|--------|-------------------|
| t3 | 4 | 1 | 1 | 1 | 2 | [4] |
| t3 | 5 | 1 | 1 | 2 | 2 | [-1, 12] |
| t3 | 6 | 1 | 1 | 3 | 2 | [-8, 48] |
| t3 | 7 | 1 | 1 | 4 | 2 | [2, -60, 240] |
| t3 | 8 | 1 | 1 | 5 | 2 | [34, -480, 1440] |
| t3 | 9 | 1 | 1 | 6 | 2 | [-17/2, 462, -4200, 10080] |
| t5 | 5 | 1 | 2 | 0 | 2 | [2] |
| t5 | 6 | 1 | 2 | 1 | 2 | [4] |
| t5 | 7 | 1 | 2 | 2 | 2 | [-1, 12] |
| t5 | 8 | 1 | 2 | 3 | 2 | [-8, 48] |
| t5 | 9 | 1 | 2 | 4 | 2 | [2, -60, 240] |
| t7 | 7 | 1 | 3 | 0 | 2 | [2] |
| t7 | 8 | 1 | 3 | 1 | 2 | [4] |
| t7 | 9 | 1 | 3 | 2 | 2 | [-1, 12] |
| t9 | 9 | 1 | 4 | 0 | 2 | [2] |
| bc | 6 | 2 | 2 | 1 | 4 | [8] |
| bc | 7 | 2 | 2 | 2 | 4 | [-2, 24] |
| bc | 8 | 2 | 2 | 3 | 4 | [-16, 96] |
| bc | 9 | 2 | 2 | 4 | 4 | [4, -120, 480] |
| bc35 | 8 | 2 | 3 | 1 | 4 | [8] |
| bc35 | 9 | 2 | 3 | 2 | 4 | [-2, 24] |
| a3 | 9 | 3 | 3 | 2 | 8 | [-4, 48] |

---

## Consequences

### 1. The Shift Principle (Corollary)

C_{t_{2j+1}}(r) at n = C_{t_{2j-1}}(r) at n-2.

This follows because both have the same free position count f = (n-1) - 2j, and the same OCF weight 2.

### 2. Complete Predictability

The W-coefficient of ANY invariant at ANY n is determined by exactly one entry in the central factorial triangle. No computation needed beyond knowing the partition type.

### 3. The Penalty Formula

H - w_0 = W(1/2) - W(0) = sum_I I(T) * 2^{parts(I)} * [F_f(1/2) - F_f(0)]
        = sum_I I(T) * 2^{parts(I)} * [1 - F_f(0)]

The penalty for each invariant is proportional to 1 - F_f(0), which is the "deviation from evaluation" of the master polynomial.

### 4. Even-n Hierarchy

At even n, the same structure applies with odd powers of r. The complete coefficient table at n=8:

| Coeff | Formula | Error |
|-------|---------|-------|
| w_7 = 40320 = n! | universal | 0 |
| w_5 = -20160 + 1440*t_3 | t_3 | 0 |
| w_3 = 3024 - 480*t_3 + 48*t_5 + 96*bc | t_3, t_5, bc | 0 |
| w_1 = -124 + 34*t_3 - 8*t_5 + 4*t_7 - 16*bc + 8*bc35 | all | 0 |

---

## Predictions (testable at n=11)

From b_{4,j} = [1, 85, 147, 30, 1]:

F_8(r) = 31 - 2640r^2 + 40320r^4 - 211680r^6 + 362880r^8

Therefore: C_{t3}(r) at n=11 = 2*F_8(r) = 62 - 5280r^2 + 80640r^4 - 423360r^6 + 725760r^8

In particular:
- w_9 at n=11 = 11! (universal)
- w_7 at n=11 has t_3 coefficient = 725760 = 2*9! (THM-058)
- w_5 at n=11 has t_3 coefficient = -423360

---

## Scripts

- `04-computation/universal_master_polynomial.py` (22-case verification)
- `04-computation/w_generating_function_test.py` (shift principle discovery)
- `04-computation/w_even_n_hierarchy.py` (even-n computations)
- `04-computation/w1_n8_complete.py` (complete n=8 table)

---

## Open Questions

1. **Algebraic proof:** Can the central factorial number recurrence be derived from the Newton identity / position pattern decomposition?

2. **Constant polynomial C_0(r):** The "background" polynomial (for the transitive tournament) also follows from the master sequence but with a different structure. What determines it?

3. **Connection to central factorials:** The central factorial numbers arise in the expansion x^n = sum T(n,k) * x^[k] where x^[k] is the central factorial. Is there a direct interpretation where r^j expands in central factorials of (r^2 - 1/4)?
