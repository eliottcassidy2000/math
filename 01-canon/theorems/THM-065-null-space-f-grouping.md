# THM-065: Null Space of Forward-Edge Distribution = f-Level Grouping

**Type:** Theorem (PROVED)
**Certainty:** 5 -- PROVED (algebraic, verified computationally n=7,9)
**Status:** PROVED
**Added by:** kind-pasteur-2026-03-07-S28
**Tags:** #forward-edge #null-space #f-grouping #reduced-polynomial

---

## Statement

The forward-edge distribution a_k(T) of a tournament T determines the OCF invariants only up to f-level equivalence.

### (i) f-level grouping

Each OCF invariant I has a free-position count f_I = n - 1 - (sum of cycle sizes - parts(I)). Two invariants I, J with f_I = f_J have proportional coefficient vectors in the a_k formula:

**coeff_k(I) = 2^{parts(I)} * c_k^{(f_I, n-1)}**

Since c_k^{(f,d)} depends only on f, the ratio coeff_k(I)/coeff_k(J) = 2^{parts(I)}/2^{parts(J)} is k-independent.

### (ii) Null space dimension formula

The null space of the map (invariants) -> (a_k) has dimension:

**null_dim(n) = Q(n) - (n+1)/2**

where Q(n) = A000009(n) = number of partitions of n into distinct parts, and (n+1)/2 = number of distinct f-values (one per even f in {0, 2, ..., n-3}).

The number of OCF invariants (= partitions into odd parts >= 3 with sum <= n) equals Q(n) - 1.

| n  | #types (Q(n)-1) | #f-values ((n-1)/2) | null_dim |
|----|-----------------|---------------------|----------|
| 5  | 2               | 2                   | 0        |
| 7  | 4               | 3                   | 1        |
| 9  | 7               | 4                   | 3        |
| 11 | 11              | 5                   | 6        |
| 13 | 17              | 6                   | 11       |

At each f-level with invariants I_1, ..., I_r having parts p_1, ..., p_r, there are r-1 null vectors of the form:

**(2^{p_j} * e_{I_i} - 2^{p_i} * e_{I_j})** for i < j

### (iii) Explicit null vectors

**n=7** (4 invariants, 3 distinct f-values = {6, 2, 0}, null dim = 1):
- f=2 group: t5 (p=1), bc (p=2). Null vector: (dt5, dbc) = (-2, 1).

**n=9** (7 invariants, 4 distinct f-values = {6, 4, 2, 0}, null dim = 3):
- f=4 group: t5 (p=1), bc33 (p=2). Null vector: (-2, 1) in (t5, bc33).
- f=2 group: t7 (p=1), bc35 (p=2), a3 (p=3). Two null vectors:
  - (-2, 1, 0) in (t7, bc35, a3)
  - (-4, 0, 1) in (t7, bc35, a3)

### (iv) What a_k CAN see

The forward-edge distribution determines exactly the **f-level weighted sums**:

**S_f = sum_{I: f_I = f} 2^{parts(I)} * I(T)**

These are the coefficients in the expansion of the independence polynomial on a coarser basis.

---

## Proof

Part (i): From THM-062, a_k(T) = A(n,k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T). The coefficient of I(T) in a_k is 2^{parts(I)} * c_k^{(f_I, n-1)}, which factors as (I-dependent scalar) * (f-dependent function of k).

Part (ii): The coefficient matrix has columns proportional within each f-group. The rank equals the number of distinct f-values, since columns from different f-groups are linearly independent (different polynomial degrees).

---

## Significance

This explains the "null space discovery" of opus-S33 at n=7. The null vector (0, -2, 0, 1) in (t3, t5, t7, bc) exists because t5 (f=2, p=1) and bc (f=2, p=2) have proportional coefficient vectors.

The key insight: **The forward-edge distribution cannot distinguish a 5-cycle from a disjoint 3-cycle pair** because both have the same number of free positions (f=2). The "internal structure" of the cycle collection (whether it's one 5-cycle or two 3-cycles) is invisible to the permutation statistics.

However, **OCF CAN distinguish them** via the independence polynomial, since alpha_1 counts individual cycles while alpha_2 counts disjoint pairs. The information lost by a_k is precisely the intra-f-level decomposition.

---

## Scripts

- `04-computation/regular_n7_classes.py` -- Verified at n=7
- `04-computation/null_space_n9.py` -- n=9 null space computation
