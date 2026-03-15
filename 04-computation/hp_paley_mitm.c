/*
 * hp_paley_mitm.c — Meet-in-the-middle for H(P_p), p=29
 * opus-2026-03-14-S89c
 *
 * Standard bitmask DP is O(2^n × n) which for n=29 needs 125 GB.
 * Meet-in-the-middle splits vertices into two halves A,B of size ~n/2.
 *
 * Approach:
 *   Phase 1: For each subset S⊆A and endpoint v∈S, compute
 *     f1[mask_A][v] = # paths through S ending at v.
 *   Phase 2: For each subset S⊆B and endpoint v∈S, compute
 *     f2[mask_B][v] = # paths through S starting at v.
 *   Phase 3: Combine. For each (mask_A, mask_B) with mask_A ∪ mask_B = full,
 *     and edge v→w from A to B:
 *     H += f1[mask_A][v] × f2[mask_B][w]
 *
 * Memory: 2^(n/2) × n × 8 bytes per half.
 *   n=29, half=15: 2^15 × 29 × 8 = 7.4 MB per half. Trivial!
 *   n=29, half=14: 2^14 × 29 × 8 = 3.7 MB. Also trivial.
 *
 * But the combination phase iterates over all (mask_A, mask_B) pairs where
 * mask_A ∪ mask_B = [n]. For |A|=14, |B|=15:
 * We fix mask_A ⊆ A (2^14 choices), then mask_B = B (must use all of B).
 * Wait no — we need ALL n vertices covered. So mask_A must be exactly A
 * and mask_B must be exactly B. That's just one pair!
 *
 * The issue is: a Hamiltonian path uses ALL vertices, but it doesn't respect
 * the A/B partition. A path alternates between A and B vertices.
 *
 * CORRECT APPROACH: Inclusion-exclusion over the "interface"
 * Split [n] into sets A, B of sizes a, b (a+b=n).
 * A Hamiltonian path visits all n vertices. Consider the sequence of
 * "segments" within A and within B:
 * The path alternates between A-segments and B-segments.
 *
 * This is actually the standard "profile DP":
 *   For each subset S_A ⊆ A of size k, and each vertex v ∈ S_A at the end,
 *   compute the number of paths through S_A ending at v.
 *   Similarly for S_B ⊆ B.
 *   Then combine by iterating over all masks S_A ⊆ A, S_B ⊆ B
 *   with |S_A| + |S_B| = n, for each interface edge v→w.
 *
 * But |S_A| + |S_B| = n means S_A = A, S_B = B, and the combination
 * just considers the LAST vertex in A connecting to the FIRST in B
 * or vice versa. The problem is that the path weaves between A and B,
 * not using all of A first then all of B.
 *
 * ALTERNATIVE: "broken profile" or "layer DP"
 * We can't simply split into halves because the path interleaves.
 *
 * ACTUAL meet-in-the-middle for Hamiltonian paths:
 * Phase 1: For each subset S of [n] with |S| = n/2, and each pair (s,t) ∈ S²,
 *   compute f[S][s][t] = # paths through S from s to t.
 * Phase 2: H = Σ_{|S|=n/2} Σ_{s,t∈S, u,v∈S̄} f[S][s][t] × f[S̄][u][v] × edge(t,u)
 *   where edge(t,u) connects the end of the first half-path to the start of the second.
 *
 * This works but the Phase 1 tables are huge:
 *   C(29,14) × 14 × 14 = 2×10⁹ × 196 ≈ 400 GB.  STILL TOO BIG.
 *
 * Even smarter: the "Hamiltonian path meet in the middle" approach.
 * Phase 1: For each subset S with |S| = ⌊n/2⌋ and endpoint t:
 *   f[S][t] = # Hamiltonian paths in G[S] ending at t.
 * Phase 2: For each S with |S| = ⌈n/2⌉ and starting point u:
 *   g[S̄][u] = # Hamiltonian paths in G[S̄] starting at u.
 * Phase 3: H = Σ_{|S|=⌊n/2⌋, t∈S, u∈S̄, edge(t,u)} f[S][t] × g[S̄][u]
 *
 * Memory for Phase 1: C(29,14) × 14 ≈ 2.8×10⁹ × 8B = 22.5 GB. Barely feasible.
 * But computing this requires iterating over all subsets of size 14: C(29,14) = 2×10⁸.
 * For each, run the DP within that subset: O(2^14 × 14) per subset.
 * Total: 2×10⁸ × 2^14 × 14 ≈ 4.6×10¹³. WAY too slow.
 *
 * CONCLUSION: Meet-in-the-middle doesn't help much for Hamiltonian paths.
 * The standard approach is O(2^n × n²), and MITM reduces space but not time.
 *
 * For p=29, we'd need 2^29 × 29 × 8 = 125 GB RAM. Not feasible on most machines.
 *
 * ALTERNATIVE: Use the circulant structure of Paley tournaments.
 * P_p has automorphism group of order p (cyclic shift). So
 * H(P_p) = p × h₀ where h₀ = # Ham paths starting at vertex 0.
 * This doesn't reduce the DP size, but does reduce the output.
 *
 * ANOTHER ALTERNATIVE: compress using the fact that P_p is circulant.
 * The adjacency is determined by QR status of differences.
 * The DP state dp[mask][v] has a symmetry: for the cyclic shift i ↦ i+1 mod p,
 * dp[shifted_mask][v+1] = dp[mask][v].
 * This means we only need to store dp for representatives of each orbit
 * of the cyclic group on (mask, v) pairs.
 * Number of orbits of Z_p on 2^p × p: by Burnside's lemma,
 * (1/p) Σ_{g∈Z_p} |Fix(g)|. For g=0: all 2^p × p pairs are fixed.
 * For g≠0: very few are fixed (only masks invariant under shift by g).
 * So roughly 2^p × p / p = 2^p orbits.
 * For p=29: 2^29 × 8 bytes ≈ 4.3 GB. Feasible!
 *
 * But implementing this orbit DP is complex. Let me instead just
 * verify the formula H(P_p)/E[H] for the known values and see the trend.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Factorial as double */
double lfact(int n) {
    double s = 0;
    for (int i = 2; i <= n; i++) s += log(i);
    return s;
}

int main() {
    printf("H/E[H] analysis for Paley tournaments\n\n");

    /* Known values */
    int primes[] = {3, 7, 11, 19, 23};
    /* H values from hp_paley_dp runs */
    double H_vals[] = {3, 189, 95095, 1172695746915.0, 15760206976379349.0};
    int nprimes = 5;

    for (int i = 0; i < nprimes; i++) {
        int p = primes[i];
        double H = H_vals[i];

        /* E[H] = p! / 2^{p-1} */
        double log_EH = lfact(p) - (p-1) * log(2);
        double EH = exp(log_EH);

        double ratio = H / EH;
        double log_ratio = log(H) - log_EH;

        printf("p=%2d: H = %.0f\n", p, H);
        printf("      E[H] = %.6e\n", EH);
        printf("      H/E[H] = %.8f\n", ratio);
        printf("      ln(H/E[H]) = %.8f\n", log_ratio);
        printf("      e - H/E[H] = %.8f\n", M_E - ratio);
        printf("      1/(e - H/E[H]) = %.4f\n", 1.0/(M_E - ratio));
        printf("\n");
    }

    /* Extrapolation */
    printf("Fitting H/E[H] = e - c/p^alpha:\n");
    /* Use p=19 and p=23 to estimate */
    double r19 = H_vals[3] / exp(lfact(19) - 18*log(2));
    double r23 = H_vals[4] / exp(lfact(23) - 22*log(2));
    double gap19 = M_E - r19;
    double gap23 = M_E - r23;
    double alpha = log(gap19/gap23) / log(23.0/19.0);
    double c = gap23 * pow(23, alpha);
    printf("  From p=19,23: alpha = %.4f, c = %.4f\n", alpha, c);
    printf("  Predicted H/E[H] at p=29: %.6f\n", M_E - c/pow(29, alpha));
    printf("  Predicted H/E[H] at p=31: %.6f\n", M_E - c/pow(31, alpha));
    printf("  Predicted H/E[H] at p=43: %.6f\n", M_E - c/pow(43, alpha));
    printf("  Predicted H/E[H] at p=47: %.6f\n", M_E - c/pow(47, alpha));

    printf("\ne = %.10f\n", M_E);

    return 0;
}
