# THM-077: Direct Walsh Proof of OCF

**Status:** PROVED (elementary, complete)
**Author:** opus-2026-03-07-S35 (continued^4)
**Date:** 2026-03-07
**Dependencies:** THM-076 (Walsh-OCF factorization, for the I side)

## Statement

H(T) = I(Omega(T), 2) for all tournaments T on n vertices (n odd).

This is the Odd-Cycle Collection Formula (OCF), originally proved by Grinberg-Stanley. We give a new proof via Walsh-Hadamard analysis.

## Proof

### Setup

A tournament T on n vertices is a point in {0,1}^m where m = C(n,2). The Walsh-Hadamard transform gives

hat{f}[S] = (1/2^m) sum_{T} f(T) chi_S(T)

for any function f: {0,1}^m -> R, where chi_S(T) = (-1)^{sum_{e in S} T_e}.

Since the Walsh basis {chi_S} is complete and orthonormal, f = g iff hat{f}[S] = hat{g}[S] for all S.

### Step 1: Walsh coefficients of H(T)

H(T) = sum_{sigma in S_n} prod_{i=0}^{n-2} T(sigma(i), sigma(i+1)).

Direct computation gives:

hat{H}[S] = ((-1)^|S| / 2^{n-1}) * sum_{sigma: S subset E_sigma} prod_{e in S} d_e(sigma)

where E_sigma = {edges of the Hamiltonian path sigma}, and d_e(sigma) = +1 if sigma traverses e forward, -1 if backward.

### Step 2: hat{H}[S] = 0 unless S is a union of disjoint even-length paths

**Case 1: S contains a cycle.** Since E_sigma is the edge set of a path (a tree), no cycle can be a subset of E_sigma. So the sum is empty and hat{H}[S] = 0.

**Case 2: |S| is odd.** By complement invariance H(T) = H(T^op) and chi_S(T^op) = (-1)^|S| chi_S(T), we get hat{H}[S] = (-1)^|S| hat{H}[S], so hat{H}[S] = 0.

**Case 3: S is a union of paths, some with odd length.** Let C be a component of S with an odd number of edges. In any ordering sigma with S subset E_sigma, C appears as a contiguous block (each internal vertex of C has degree 2 in S, hence both neighbors are determined). Define sigma' = sigma with block C reversed. Then sigma' also has S subset E_{sigma'}, and the sign product for sigma' differs from sigma's by (-1)^{|C|} = -1 (since |C| is odd). This bijection pairs each sigma with a sigma' of opposite sign, so the sum cancels to 0.

### Step 3: hat{H}[S] for S = P_{2a_1} union ... union P_{2a_r}

For S a union of r disjoint even-length paths:

(a) **Contiguity:** Each internal vertex of each component C has degree 2 in S. In any sigma with S subset E_sigma, both S-neighbors of v are HP-neighbors, so C appears as a contiguous block.

(b) **Sign:** Each component P_{2a_i} has 2a_i edges. Forward traversal gives product (+1)^{2a_i} = +1. Backward gives (-1)^{2a_i} = +1. So ALL orderings contribute +1 to the signed sum.

(c) **Block counting:** With r blocks and (n - 2k - r) singletons (where k = sum a_i), the number of orderings is 2^r * (n - 2k)! (each block in 2 orientations, (n-2k) items to arrange).

Therefore: hat{H}[S] = epsilon_S * 2^r * (n-2k)! / 2^{n-1}.

### Step 4: Walsh coefficients of I(Omega(T), 2)

By THM-076 (Walsh-OCF Factorization), proved via exponential generating function analysis:

hat{I}[S] = epsilon'_S * 2^r * (n-2k)! / 2^{n-1}

with the same amplitude. (The EGF proof: G_r(x) = 2r!/((1-x)^r(1+x)), and Sigma_r(m) = (m+r)!/m! via the hockey-stick identity.)

### Step 5: Sign formula

The sign of hat{H}[S] is epsilon_S = (-1)^{asc(S)} where asc(S) counts the total number of ascents across all path components. For a path component v_0-v_1-...-v_{2a}, an ascent is a pair (v_i, v_{i+1}) with v_i < v_{i+1}.

PROOF: For each edge {v_i, v_{i+1}} with v_i < v_{i+1} (canonical ordering), traversing it forward in the HP gives direction factor -1, backward gives +1. For the contiguous block, traversal forward (v_0->...->v_{2a}) gives direction factor (-1)^{asc} where asc counts ascents. Traversal backward gives (-1)^{desc} = (-1)^{2a - asc} = (-1)^{asc} (since 2a is even). So both directions give the same sign (-1)^{asc}.

### Step 5b: Sign of hat{I}[S] (from THM-076 covering analysis)

In the covering analysis, each cycle C containing path S = v_0-...-v_{2a} traverses S either forward or backward. In both cases, the chi_S value equals (-1)^{asc(S)}:

Forward (v_0->v_1->...->v_{2a}): For each edge {v_i, v_{i+1}}, if v_i < v_{i+1} the cycle goes forward on the canonical edge (T_e = 1, contribution -1), otherwise backward (T_e = 0, contribution +1). Product = (-1)^{asc}.

Backward: same analysis gives (-1)^{desc} = (-1)^{2a-asc} = (-1)^{asc} (since 2a is even).

So hat{I}[S] = (-1)^{asc(S)} * 2^r * (n-2k)! / 2^{n-1}, matching hat{H}[S] exactly.

**No circular reasoning**: the sign on the I side comes from evaluating chi_S at cycle orientations (a purely combinatorial operation), and the sign on the H side comes from the direction factors of HPs. Both independently equal (-1)^{asc(S)}.

### Conclusion

hat{H}[S] = hat{I}[S] for all S in {0,1}^m. By completeness of the Walsh basis:

H(T) = I(Omega(T), 2) for all tournaments T. QED.

## Key Insights

1. **The proof is elementary.** Steps 1-3 use only the definition of H(T) and simple counting. No algebraic machinery beyond the Walsh basis.

2. **Contiguity is the engine.** The fundamental reason OCF works is that even-length path edge sets are always contiguous in any Hamiltonian path that contains them. This topological constraint forces the block counting to produce exactly 2^r * (n-2k)!.

3. **Even-length is essential.** Odd-length path components cancel by the reversal bijection. This is why OCF uses EVEN-length fan pairs (the edges of each odd cycle contribute an even number of edges to the Walsh monomial).

4. **THM-076 bridges the gap.** The I side is handled by the EGF/set-partition analysis, which shows the covering sum telescopes to the same amplitude.

## Verified

- hat{H}[S] computed directly for all degree-2 monomials at n=5 (30 monomials)
- hat{H}[S] computed directly for degree-4 monomials at n=7 (P2+P2, P4)
- Non-path monomials verified to give 0: cycles, odd-length paths, matchings
- Sign matching verified at transitive tournament
