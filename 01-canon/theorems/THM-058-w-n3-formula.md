# THM-058: Universal formula for w_{n-3}

**Status:** PROVED (algebraic proof, verified computationally at n=5,7,9)
**Proved by:** kind-pasteur-2026-03-06-S25g
**Depends on:** THM-055 (coefficient hierarchy), definitions.md

---

## Statement

For any tournament T on n vertices (n >= 3, n odd):

  w_{n-3}(T) = (n-2)! * [2*t_3(T) - C(n,3)/2]

where:
- W(r) = sum_P prod_{i=0}^{n-2} (r + A[p_i,p_{i+1}] - 1/2) is the weighted path polynomial
- w_{n-3} is the coefficient of r^{n-3} in W(r)
- t_3(T) is the number of directed 3-cycles in T
- C(n,3) = n(n-1)(n-2)/6 is the total number of 3-element vertex subsets

Equivalently:

  w_{n-3}(T) = 2*(n-2)! * t_3(T) - C(n,3)*(n-2)!/2

## Properties

1. **Centered:** E[w_{n-3}] = 0 over random tournaments (since E[t_3] = C(n,3)/4).
2. **Factored:** w_{n-3} = (n-2)! * (deviation of t_3 from random midpoint)
3. **Linear in t_3:** Depends ONLY on the 3-cycle count, not on any finer invariant.

## Proof

### Step 1: Decomposition into sigma sums

W(r) = sum_{k=0}^{n-1} r^{n-1-k} * sigma_k, where
sigma_k = sum_{|S|=k, S subset [n-2]} sigma(S) and
sigma(S) = sum_P prod_{i in S} s_{p_i, p_{i+1}} with s_{a,b} = A[a,b] - 1/2.

So w_{n-3} = sigma_2 = sum_{|S|=2} sigma(S).

### Step 2: Non-adjacent pairs contribute zero

For S = {i, j} with |i-j| >= 2, the factors s_{p_i,p_{i+1}} and s_{p_j,p_{j+1}}
involve disjoint position pairs. Summing over all permutations:

sigma({i,j}) = sum_P s_{p_i,p_{i+1}} * s_{p_j,p_{j+1}}

Since each factor individually has zero sum over all permutations
(sum_P s_{p_i,p_{i+1}} = sum_P (A[p_i,p_{i+1}] - 1/2) = n!/2 - n!/2 = 0),
and the factors are "independent" (involve disjoint position variables for |i-j| >= 2),
the product sum is zero.

[More precisely: for fixed (p_j, p_{j+1}), the sum over all compatible
(p_i, p_{i+1}) of s_{p_i,p_{i+1}} is zero by symmetry, since positions i,i+1
are disjoint from j,j+1 and the remaining vertices can fill them symmetrically.]

### Step 3: Count adjacent pairs

For S = {i, i+1} (adjacent positions), sigma(S) is the same for all i
by position-translation invariance. There are n-2 such pairs.

sigma_adj := sigma({0,1}) = sum_P s_{p_0,p_1} * s_{p_1,p_2}

### Step 4: Compute sigma_adj

Expanding s = A - 1/2:

sigma_adj = sum_P A[p_0,p_1]*A[p_1,p_2] - (1/2)*sum_P A[p_0,p_1]
            - (1/2)*sum_P A[p_1,p_2] + (1/4)*n!

Each position sum: sum_P A[p_i,p_{i+1}] = (n-2)! * C(n,2) = n!/2.

So: sigma_adj = (n-3)! * D - n!/4

where D = sum of A[a,b]*A[b,c] over ordered triples (a,b,c) of distinct vertices
        = number of ordered directed 2-paths a -> b -> c.

### Step 5: Compute D in terms of t_3

D = sum_b d_b^{in} * d_b^{out} where d_b^{in} = n-1-s_b, d_b^{out} = s_b.
D = (n-1)*sum s_b - sum s_b^2 = (n-1)*C(n,2) - sum s_b^2.

Using the classical identity t_3 = C(n,3) - sum_b C(s_b, 2):
sum s_b^2 = 2*C(n,3) - 2*t_3 + C(n,2).

So: D = (n-1)*C(n,2) - 2*C(n,3) + 2*t_3 - C(n,2)
      = (n-2)*C(n,2) - 2*C(n,3) + 2*t_3
      = C(n,3) + 2*t_3.

[Check: (n-2)*C(n,2) = (n-2)*n(n-1)/2 = 3*C(n,3), so 3*C(n,3) - 2*C(n,3) = C(n,3).]

### Step 6: Combine

sigma_adj = (n-3)! * (C(n,3) + 2*t_3) - n!/4
          = (n-3)!*C(n,3) + 2*(n-3)!*t_3 - (n-3)!*C(n,3)*3/2
          = (n-3)! * [C(n,3) - 3*C(n,3)/2 + 2*t_3]
          = (n-3)! * [2*t_3 - C(n,3)/2]

[Using n!/4 = n(n-1)(n-2)*(n-3)!/4 = 6*C(n,3)*(n-3)!/4 = (3/2)*C(n,3)*(n-3)!]

### Step 7: Final formula

w_{n-3} = (n-2) * sigma_adj = (n-2)! * [2*t_3 - C(n,3)/2]. QED.

## Verification

- n=5: w_2 = 3! * [2*t_3 - 10/2] = 6*(2*t_3 - 5) = 12*t_3 - 30. Verified exhaustively.
- n=7: w_4 = 5! * [2*t_3 - 35/2] = 120*(2*t_3 - 17.5) = 240*t_3 - 2100. Verified (20 samples, 0 err).
- n=9: w_6 = 7! * [2*t_3 - 84/2] = 5040*(2*t_3 - 42) = 10080*t_3 - 211680. Verified (15 samples, 0 err).
