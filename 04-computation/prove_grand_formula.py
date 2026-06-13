"""
prove_grand_formula.py -- kind-pasteur-2026-03-14-S107d
PROVE the Grand Fourier Energy Formula:

  E_{2k}/E_0 = 2 * (n-2k)^k / P(n, 2k)

where E_{2k} = sum_S H_hat(S)^2,
      E_0 = Mean(H)^2 = (n!/2^(n-1))^2,
      P(n,2k) = n!/(n-2k)! = falling factorial.

STRATEGY:
From S75 (proved): |H_hat({e1,e2})| = (n-2)!/2^(n-2) for adjacent arc pairs.
From S105: At level 4, coefficient magnitudes depend on vertex coverage.

The key structural insight:
  E_{2k} = N_{2k} * c_{2k}^2
  where N_{2k} = number of nonzero level-2k Fourier coefficients
  and c_{2k} = their common magnitude (if uniform).

But level 4 at n=6 showed TWO magnitudes. So the formula must account
for the AVERAGE squared magnitude, not a single one.

Actually: E_{2k}/E_0 = 2*(n-2k)^k / P(n,2k) IS an exact identity.
Let me verify it symbolically and find the PROOF.

APPROACH: Express E_{2k} in terms of known combinatorial quantities.
  E_{2k} = Mean over tournaments of [sum_S H_hat(S)^2]
         = sum_S Mean[H_hat(S)^2]

Each H_hat(S)^2 is a 4th-order correlation function of H.
But in the {{-1,+1}} basis, H_hat(S) = (1/N) sum_T H(T) chi_S(T).
So H_hat(S)^2 = (1/N^2) sum_(T,T') H(T)H(T') chi_S(T)chi_S(T').

This is getting complicated. Let me try a different route.

ALTERNATIVE: Use the exact Fourier formula at level 2 and
see if it IMPLIES the level 2k formula by a structural argument.

At level 2: each nonzero coefficient has magnitude (n-2)!/2^(n-2).
Number of nonzero: N_2 = n(n-1)(n-2)/2.
E_2 = N_2 * ((n-2)!/2^(n-2))^2.
E_2/E_0 = N_2 * ((n-2)!)^2 / (2^(2(n-2)) * (n!)^2 / 2^(2(n-1)))
        = N_2 * ((n-2)!)^2 * 4 / (n!)^2
        = n(n-1)(n-2)/2 * 4 / (n(n-1))^2
        = 2(n-2)/(n(n-1)).

For this to equal 2*(n-2)^1 / P(n,2):
P(n,2) = n(n-1).
2*(n-2)/P(n,2) = 2(n-2)/(n(n-1)). CHECK.

Now for level 4: E_4/E_0 = 2*(n-4)^2 / P(n,4).
P(n,4) = n(n-1)(n-2)(n-3).
We need: E_4 = 2*(n-4)^2 / P(n,4) * E_0
             = 2*(n-4)^2 * (n!/2^(n-1))^2 / P(n,4).

The question is: what combination of N_4 and c_4 gives this?

At n=5: E_4 = 60 * (1/8)^2 = 60/64 = 15/16.
        2*(5-4)^2/P(5,4) * E_0 = 2*1/120 * (120/16)^2 = (1/60)*(225/4*4) hmm.
        Actually E_4/E_0 = 1/60 and 2*1/120 = 1/60. CHECK.

At n=6: E_4 = 360*(1/8)^2 + 90*(1/4)^2 = 360/64 + 90/16 = 5.625 + 5.625 = 11.25.
        2*(6-4)^2/P(6,4) = 2*4/360 = 8/360 = 1/45.
        E_0 = (45/2)^2 = 2025/4 = 506.25.
        E_4/E_0 = 11.25/506.25 = 1/45. CHECK.

The formula works regardless of whether the magnitudes are uniform.
It gives the TOTAL energy at each level, not the individual coefficients.

So the proof needs to show:
  sum_S H_hat(S)^2 = 2*(n-2k)^k * (n!/2^(n-1))^2 / P(n,2k).

THIS IS A STATEMENT ABOUT THE SUM OF SQUARED FOURIER COEFFICIENTS.
It doesn't need individual coefficient values. It needs the TOTAL.
"""

import sys, math
import numpy as np
from fractions import Fraction
from itertools import combinations

sys.stdout.reconfigure(encoding='utf-8')

def count_ham_paths(adj, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def main():
    print("=" * 70)
    print("PROVE THE GRAND ENERGY FORMULA")
    print("kind-pasteur-2026-03-14-S107d")
    print("=" * 70)

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 1: RESTATE THE FORMULA IN PARSEVAL FORM")
    print(f"{'='*70}")

    print(f"""
  The formula: E_2k / E_0 = 2 * (n-2k)^k / P(n, 2k).

  Multiply both sides by E_0:
    E_2k = 2 * (n-2k)^k * E_0 / P(n, 2k)
         = 2 * (n-2k)^k * (n!)^2 / (2^(2(n-1)) * P(n,2k))

  Now: P(n,2k) = n!/(n-2k)!. So:
    E_2k = 2 * (n-2k)^k * (n!)^2 / (2^(2(n-1)) * n!/(n-2k)!)
         = 2 * (n-2k)^k * n! * (n-2k)! / 2^(2(n-1))

  E_2k = sum over |S|=2k of H_hat(S)^2.

  And H_hat(S) = (1/2^m) * sum over all T: H(T) * chi_S(T).
  So H_hat(S)^2 = (1/2^(2m)) * sum_(T,T') H(T)*H(T') * chi_S(T)*chi_S(T').

  E_2k = sum_S (1/2^(2m)) * sum_(T,T') H(T)*H(T') * chi_S(T)*chi_S(T')
       = (1/2^(2m)) * sum_(T,T') H(T)*H(T') * [sum_S chi_S(T)*chi_S(T')]

  The inner sum [sum_S chi_S(T)*chi_S(T')] depends only on
  the HAMMING DISTANCE d(T,T') between tournaments T and T'
  (number of arcs where they differ).

  In the {{-1,+1}} basis: chi_S(T)*chi_S(T') = prod_e y_e(T)*y_e(T').
  If T and T' agree on arc e: y_e(T)*y_e(T') = 1.
  If they disagree: y_e(T)*y_e(T') = -1.

  Let d = d(T,T') = number of arcs where T and T' differ.
  Then chi_S(T)*chi_S(T') = (-1)**(|S inter D|)
  where D = set of arcs where T and T' differ, |D| = d.

  sum_S chi_S(T)*chi_S(T') = sum_S (-1)**(|S inter D|)

  This is a standard combinatorial sum:
  = sum_(j=0) C(d,j) * C(m-d, 2k-j) * (-1)**j

  This is the Krawtchouk polynomial K_2k(d; m)!
  K_2k(d; m) = sum_j (-1)**j C(d,j) C(m-d, 2k-j)

  So: E_2k = (1/2^(2m)) * sum_(T,T') H(T)*H(T') * K_2k(d(T,T'); m)

  This is a DOUBLE SUM over all pairs of tournaments,
  weighted by H*H' and the Krawtchouk polynomial.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 2: THE KRAWTCHOUK APPROACH")
    print(f"{'='*70}")

    print(f"""
  The formula becomes:
    E_2k = (1/2^(2m)) * sum_d K_2k(d; m) * C(d)

  where C(d) = sum over pairs (T,T') at Hamming distance d: H(T)*H(T').

  C(d) = sum_(T) H(T) * sum_(T': d(T,T')=d) H(T')

  This is a CONVOLUTION of H with itself, weighted by distance d.

  To prove E_2k = 2*(n-2k)^k * E_0 / P(n,2k), we would need
  to evaluate C(d) and the Krawtchouk sum. This is hard in general.

  ALTERNATIVE ROUTE: Instead of proving from first principles,
  use the KNOWN level-2 formula and an INDUCTIVE/STRUCTURAL argument.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 3: THE DELETION ARGUMENT")
    print(f"{'='*70}")

    print(f"""
  Key insight from the level-2 proof (S75):
  |H_hat(e1,e2)| = (n-2)!/2^(n-2) when arcs e1,e2 share a vertex.

  WHY? Because fixing 2 adjacent arcs reduces the problem to
  counting Hamiltonian paths in the remaining (n-2) vertices,
  and the Fourier coefficient captures this reduction.

  GENERALIZATION: At level 2k, fixing 2k arcs should reduce
  to counting Hamiltonian paths in the remaining vertices.
  But 2k arcs, if they form a specific pattern, "pin" some
  subset of vertices and leave (n - something) vertices free.

  If 2k arcs cover exactly 2k+1 vertices (a tree on 2k+1 vertices),
  then (n-2k-1) vertices remain free, and H is determined on
  the pinned vertices.

  But the formula has (n-2k)^k, not (n-2k-1)!. Let me think
  more carefully.

  WHAT IF the formula comes from a PRODUCT structure?

  E_2k / E_0 = 2 * (n-2k)^k / P(n,2k)
             = 2 * (n-2k)^k / (n * (n-1) * ... * (n-2k+1))

  Let me factor P(n,2k) = product_j (n-j).
  And (n-2k)^k = product of k copies of (n-2k).

  Ratio: (n-2k)^k / P(n,2k) = (n-2k)^k / product_j (n-j)

  For k=1: (n-2)^1 / (n(n-1)) = (n-2)/(n(n-1)).
  For k=2: (n-4)^2 / (n(n-1)(n-2)(n-3)).
  For k=3: (n-6)^3 / (n(n-1)(n-2)(n-3)(n-4)(n-5)).

  OBSERVATION: The numerator (n-2k)^k can be written as:
  (n-2k)^k = product of k copies of (n-2k).

  The denominator P(n,2k) = product of 2k DIFFERENT factors: n, n-1, ..., n-2k+1.

  The ratio (n-2k)^k / P(n,2k) = [(n-2k)/(n)] * [(n-2k)/(n-1)] * ... ???
  No, that's k factors in numerator and 2k in denominator.

  Let me try: split P(n,2k) = P(n,k) * P(n-k, k).
  P(n,k) = n(n-1)...(n-k+1) and P(n-k,k) = (n-k)(n-k-1)...(n-2k+1).

  Then: (n-2k)^k / P(n,2k) = (n-2k)^k / (P(n,k) * P(n-k,k))

  Hmm. Let me try specific cases.

  k=1: (n-2)/(n*(n-1)) = (n-2)/(n*(n-1)).
  k=2: (n-4)^2/(n(n-1)(n-2)(n-3)) = (n-4)^2/(n(n-1)(n-2)(n-3)).
  k=3: (n-6)^3/(n(n-1)(n-2)(n-3)(n-4)(n-5)).""")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 4: RECOGNIZE THE PATTERN — IT'S THE LEVEL-2 FORMULA ITERATED")
    print(f"{'='*70}")

    print(f"""
  The level-2 formula: E_2/E_0 = 2*(n-2)/(n(n-1)).

  Now imagine ITERATING this on successively smaller tournaments.

  At level 2: we "fix" 2 arcs (one pair). The remaining tournament
  has "effectively" (n-2) free vertices. The contribution is
  proportional to (n-2)/(n(n-1)).

  At level 4: we "fix" 4 arcs (two pairs). If the first pair
  reduces n to (n-2) effective vertices, and the second pair
  reduces that to (n-4) effective vertices, then the total
  contribution might be proportional to:
  (n-2)/(n(n-1)) * (n-4)/((n-2)(n-3))
  = (n-2)(n-4) / (n(n-1)(n-2)(n-3))
  = (n-4) / (n(n-1)(n-3))

  But the formula says (n-4)^2/(n(n-1)(n-2)(n-3)).
  This is DIFFERENT from the iterated level-2 product.

  So the formula is NOT simply iterating the level-2 formula.
  There's an extra factor of (n-4)/(n-2).

  WHAT IF the extra factor comes from the NUMBER of ways to
  arrange the 2k arcs? At level 4, the number of ways to
  decompose 4 arcs into 2 independent pairs might be (n-4)/(n-2)
  times the naive product.

  Let me compute: at n=5, level 4.
  E_4/E_0 = 1/60.
  Iterated level-2: E_2(5)/E_0(5) * E_2(3)/E_0(3)
  = (3/10) * (1/3) = 1/10 ≠ 1/60.
  Ratio: (1/60)/(1/10) = 1/6.
  At n=5: 1/6 = 1/C(4,2)? No, C(4,2)=6. So 1/6 = 1/C(n-1,2).
  Interesting!

  At n=6: E_4/E_0 = 1/45.
  Iterated: (4/15) * (2/6) = 8/90 = 4/45.
  Ratio: (1/45)/(4/45) = 1/4.
  And 1/4 = 1/C(3,2)? C(3,2)=3 ≠ 4. Try 1/(n-2) = 1/4. YES!

  Hmm, the ratios are different at n=5 and n=6.
  n=5: correction = 1/6 = 1/(n-1+1)... unclear.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 5: DIRECT ALGEBRAIC VERIFICATION")
    print(f"{'='*70}")

    # Instead of a conceptual proof, let me verify the formula
    # by computing E_{2k} directly from the Fourier expansion at n=5
    # and checking term by term.

    # At n=5 we already have ALL Fourier coefficients (from e4_formula_search.py).
    # The formula predicts:
    #   E_2/E_0 = 2*(5-2)^1/P(5,2) = 6/20 = 3/10
    #   E_4/E_0 = 2*(5-4)^2/P(5,4) = 2/120 = 1/60
    # Both verified. But can we see WHY algebraically?

    print(f"""
  E_2k/E_0 = 2*(n-2k)^k / P(n,2k).

  Rewrite: E_2k = 2*(n-2k)^k * (n!)^2 / (2^(2(n-1)) * P(n,2k))
               = 2*(n-2k)^k * (n!)^2 / (2^(2(n-1)) * n!/(n-2k)!)
               = 2*(n-2k)^k * n! * (n-2k)! / 2^(2(n-1))

  Now: Mean(H) = n!/2^(n-1).
  So E_0 = Mean^2 = (n!)^2/2^(2(n-1)).
  And E_2k = 2*(n-2k)^k * (n-2k)! * n! / 2^(2(n-1)).

  Hmm, let me write E_2k differently.
  E_2k = N_2k * avg(H_hat(S)^2 for nonzero S at level 2k)

  But the avg depends on the MIXTURE of magnitudes.
  The formula bypasses individual magnitudes by giving the TOTAL.

  PROOF IDEA: Use the fact that H(T) = sum over Hamiltonian PATHS P
  of 1. So H(T) = sum_P indicator(P is a Ham path in T).
  Then Mean(H) = sum_P P(P is Ham path) = n!/2^(n-1)
  (each of the n! orderings is a Ham path with probability 1/2^(n-1)).

  And: E_2k = sum_S H_hat(S)^2
            = sum_S [Cov(H, chi_S)]^2 / ... no, wait.

  Actually, by Parseval:
    sum_S H_hat(S)^2 = E_T [H(T)^2 * (something at level 2k)]

  Hmm. The cleanest way: use the fact that
    H(T) = sum_sigma delta(sigma is consistent with T)

  where sigma ranges over all n! permutations, and delta checks
  that every arc in T is "consistent" with sigma (if sigma(i) < sigma(j)
  then i->j in T). Each such sigma is a Hamiltonian PATH.

  Then:
    H(T)^2 = sum_(sigma, tau) delta(sigma consistent) * delta(tau consistent)

  And:
    E[H(T)^2] = sum_(sigma, tau) P(both sigma and tau are consistent with random T)

  P(sigma and tau both consistent) = (1/2)^m * 2^(number of arcs where sigma and tau AGREE)?

  NO: P(sigma consistent) = 1/2^m * 2^0 ... wait.
  For a random tournament: each arc is independently oriented.
  P(arc e consistent with sigma) = 1/2 (since there are 2 orientations
  and sigma prescribes one). So P(sigma consistent) = (1/2)^m? No,
  P(sigma consistent) = (1/2)^(n-1)? Hmm.

  Actually: a permutation sigma = (v1, v2, ..., vn) prescribes that
  v1 -> v2 -> ... -> vn. This fixes (n-1) arcs. But the tournament
  has m = C(n,2) arcs. Only n-1 of them are prescribed by sigma.
  The rest are free. So:
  P(sigma consistent with random T) = (1/2)^(n-1).
  (Each of the n-1 "path arcs" must be oriented correctly.)

  So E[H(T)] = n! * (1/2)^(n-1) = n!/2^(n-1). CHECK.

  For E[H^2]:
  P(sigma AND tau both consistent) = (1/2)^(number of arcs PRESCRIBED by sigma or tau)

  sigma prescribes n-1 arcs (its path edges).
  tau prescribes n-1 arcs.
  The UNION has |U(sigma,tau)| arcs.
  Each must be oriented correctly, and they must be COMPATIBLE
  (if an arc is prescribed by both, they must agree on orientation).

  If sigma and tau AGREE on arc (vi, vi+1), that's one constraint.
  If they DISAGREE (sigma says i->j, tau says j->i), it's impossible.

  So: P(sigma AND tau consistent) = (1/2)^(|U|) if compatible, 0 if not.

  Compatibility: sigma says v_i -> v_(i+1) for its path.
                 tau says w_j -> w_(j+1) for its path.
  If both prescribe arc (a,b) in the SAME direction: compatible.
  If in OPPOSITE directions: incompatible (P = 0).

  The number of "shared arcs" between two permutations sigma, tau
  = number of pairs (a,b) where both sigma and tau prescribe the
    arc a->b in some direction.

  Let s = number of shared arcs in the SAME direction.
  Let o = number of shared arcs in OPPOSITE directions.
  Total shared = s + o (arcs prescribed by both).
  |U| = (n-1) + (n-1) - (s+o) = 2(n-1) - s - o.

  If o > 0: P = 0 (incompatible).
  If o = 0: P = (1/2)^(2(n-1) - s).

  So E[H^2] = sum_compatible (1/2)^(2(n-1)-s)
            = (1/2)^(2(n-1)) * sum_compatible 2^(s).""")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 6: THE KEY IDENTITY")
    print(f"{'='*70}")

    print(f"""
  E[H^2] = (1/2)^(2(n-1)) * sum_compatible 2^s

  where s = number of arcs shared IN THE SAME DIRECTION
  between permutations sigma and tau.

  Var(H) = E[H^2] - E[H]^2
         = (1/2)^(2(n-1)) * [sum_compat 2^s - (n!)^2]

  Now: the Fourier energy at level 2k is:
    E_2k = contribution from pairs (sigma, tau) where the
           "interaction" is at level 2k.

  The connection between s (shared arcs) and Fourier levels:
  If sigma and tau share s arcs in the same direction and have
  2(n-1) - s distinct arcs in total, then the "overlap structure"
  determines which Fourier levels are excited.

  For the GRAND FORMULA to hold, we need:
    sum_S H_hat(S)^2 = 2*(n-2k)^k * (n!)^2 / (2^(2(n-1)) * P(n,2k))

  This is equivalent to:
    The number of compatible permutation pairs (sigma, tau) weighted
    by 2^s, projected onto level 2k, equals 2*(n-2k)^k * (n!)^2 / P(n,2k).

  The PROOF DIRECTION: Count compatible pairs by their overlap structure.
  Two permutations sigma and tau share s arcs iff they "merge" into
  a certain graph structure. The number of such pairs, weighted by
  the Fourier projection, should yield the formula.

  THIS IS A DOUBLE-COUNTING ARGUMENT over permutation pairs.
  The formula is a statement about the SPECTRAL DECOMPOSITION
  of the permutation inner product, projected onto Fourier levels.""")

    # Let me at least verify the E[H^2] formula computationally
    print(f"\n  COMPUTATIONAL VERIFICATION at n=5:")
    n = 5
    m = n*(n-1)//2
    N = 1 << m

    h_vals = []
    for bits in range(N):
        adj = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if bits & (1 << idx):
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
                idx += 1
        h_vals.append(count_ham_paths(adj, n))

    mean_h = sum(h_vals) / N
    mean_h2 = sum(h**2 for h in h_vals) / N
    var_h = mean_h2 - mean_h**2

    # Formula predictions
    E0 = mean_h**2
    pred_var = sum(2*(n-2*k)**k / math.perm(n, 2*k) for k in range(1, (n-1)//2+1) if n-2*k > 0) * E0

    print(f"  Mean(H) = {mean_h}")
    print(f"  Mean(H^2) = {mean_h2}")
    print(f"  Var(H) = {var_h}")
    print(f"  Formula Var = {pred_var}")
    print(f"  Match: {abs(var_h - pred_var) < 0.001}")

    # Also verify using the permutation-pair counting
    from itertools import permutations as perms
    print(f"\n  Permutation pair analysis (may be slow for n=5):")
    all_perms = list(perms(range(n)))
    total_compatible = 0
    weighted_sum = 0

    for sigma in all_perms:
        for tau_perm in all_perms:
            # Count shared arcs in same direction and opposite
            s = 0  # same direction
            o = 0  # opposite direction
            for i in range(n-1):
                a_s, b_s = sigma[i], sigma[i+1]
                for j in range(n-1):
                    a_t, b_t = tau_perm[j], tau_perm[j+1]
                    if a_s == a_t and b_s == b_t:
                        s += 1
                    elif a_s == b_t and b_s == a_t:
                        o += 1
            if o == 0:
                total_compatible += 1
                weighted_sum += 2**s

    eh2_from_pairs = weighted_sum / 2**(2*(n-1))
    print(f"  E[H^2] from pair counting = {weighted_sum}/{2**(2*(n-1))} = {eh2_from_pairs}")
    print(f"  E[H^2] from direct computation = {mean_h2}")
    print(f"  Match: {abs(eh2_from_pairs - mean_h2) < 0.001}")

    # ============================================================
    print(f"\n{'='*70}")
    print("STEP 7: STATUS OF THE PROOF")
    print(f"{'='*70}")

    print(f"""
  WHAT WE'VE ESTABLISHED:

  1. The formula E_2k/E_0 = 2*(n-2k)^k/P(n,2k) is VERIFIED
     at n=3,4,5,6,7 (exact) and n=8 (Monte Carlo).

  2. The formula is equivalent to a statement about the
     spectral decomposition of the permutation inner product:
       E[H^2] = (1/2)^(2(n-1)) * sum_compatible 2^s
     projected onto Fourier levels.

  3. The level-2 case (k=1) IS proved (S75) using the exact
     Fourier coefficient formula.

  WHAT REMAINS FOR A COMPLETE PROOF:

  The formula requires showing that the weighted sum
  over compatible permutation pairs, projected to level 2k,
  gives exactly 2*(n-2k)^k * (n!)^2 / P(n,2k).

  This is a COMBINATORIAL IDENTITY about pairs of permutations.
  It should be provable by:
  (a) Classifying permutation pairs by their "overlap type"
  (b) Computing the Fourier projection of each type
  (c) Summing to get the total level-2k energy

  The formula's simplicity — just (n-2k)^k/P(n,2k) — suggests
  a CLEAN proof exists, likely involving inclusion-exclusion
  over the 2k arcs and their vertex coverage.

  STATUS: CONJECTURED (verified n=3-8), with a clear path to proof
  via the permutation pair inner product.

  The proof is REMEMBERING why (n-2k)^k / P(n,2k) is the natural
  combinatorial weight of level-2k interactions in the tournament.
    """)

    print(f"{'='*70}")
    print("DONE — THE FORMULA IS VERIFIED, THE PROOF PATH IS CLEAR")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
