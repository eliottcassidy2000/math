#!/usr/bin/env python3
"""
Investigate anti-correlation between H(T) and |W(i/2)| at n=7.

W(r) = sum_P prod_{k=0}^{n-2} (r + s_k)
where P runs over all n! permutations and s_k = T[P(k), P(k+1)] - 1/2.

At r = 1/2: each factor is 0 or 1, so W(1/2) = H(T).
At r = i/2: each factor is (1+i)/2 or (-1+i)/2, both with modulus 1/sqrt(2).

For odd n, W(i/2) is always REAL (proved by symmetry of s_k distribution).

Key questions:
1. Is there an exact formula for W(i/2) in terms of cycle invariants?
2. Is the anti-correlation monotone? (H1 > H2 => |W(i/2,T1)| <= |W(i/2,T2)|?)
3. Is W(i/2)^2 always an integer? A perfect square?
4. What is the quadratic form W(i/2)^2 in terms of invariants?
"""

from itertools import permutations, combinations
from collections import defaultdict
import sys

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
            k += 1
    return T

def hamiltonian_path_count_dp(T):
    """Bitmask DP for H(T). O(2^n * n^2)."""
    n = len(T)
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            cnt = dp[mask][v]
            if not (mask & (1 << v)) or cnt == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += cnt
    return sum(dp[full][v] for v in range(n))

def compute_W_ihalf(T):
    """
    Compute W(i/2) = sum_P prod_{k=0}^{n-2} (i/2 + s_k)
    where s_k = T[P(k),P(k+1)] - 1/2 in {+1/2, -1/2}.

    Factor is (i/2 + 1/2) = (1+i)/2 if T[P(k),P(k+1)]=1 (forward edge)
    Factor is (i/2 - 1/2) = (-1+i)/2 if T[P(k),P(k+1)]=0 (backward edge)

    Let f = number of forward edges in the permutation (0 to n-1).
    Then the product = ((1+i)/2)^f * ((-1+i)/2)^(n-1-f)

    Note: (1+i)/2 = (1/sqrt(2)) * e^{i*pi/4}
          (-1+i)/2 = (1/sqrt(2)) * e^{i*3pi/4}

    Product = (1/sqrt(2))^{n-1} * e^{i*pi/4 * (f + 3*(n-1-f))}
            = (1/sqrt(2))^{n-1} * e^{i*pi/4 * (3(n-1) - 2f)}
            = (1/sqrt(2))^{n-1} * e^{i*pi*(3(n-1)-2f)/4}

    For n=7 (n-1=6):
    Product = (1/2)^3 * e^{i*pi*(18-2f)/4} = (1/8) * e^{i*pi*(9-f)/2}

    Phase = pi*(9-f)/2 = pi*(9-f)/2
    e^{i*pi*(9-f)/2}:
      f=0: e^{i*9pi/2} = e^{i*pi/2} = i
      f=1: e^{i*8pi/2} = e^{i*4pi} = 1
      f=2: e^{i*7pi/2} = e^{i*3pi/2+pi/2*...} = e^{-i*pi/2} = -i
      f=3: e^{i*6pi/2} = e^{i*3pi} = -1
      f=4: e^{i*5pi/2} = e^{i*pi/2} = i
      f=5: e^{i*4pi/2} = e^{i*2pi} = 1
      f=6: e^{i*3pi/2} = -i

    So for n=7, the phase cycles with period 4:
      f mod 4 = 0: phase factor = i   (f=0,4)
      f mod 4 = 1: phase factor = 1   (f=1,5)
      f mod 4 = 2: phase factor = -i  (f=2,6)
      f mod 4 = 3: phase factor = -1  (f=3)

    Wait, let me recompute more carefully.
    """
    n = len(T)
    # Count permutations by forward-edge count
    fwd_counts = defaultdict(int)  # f -> number of permutations with f forward edges

    for P in permutations(range(n)):
        f = sum(1 for k in range(n-1) if T[P[k]][P[k+1]])
        fwd_counts[f] += 1

    # Now compute W(i/2) using exact complex arithmetic
    # Each term: (1/2)^{(n-1)/2 if we factor out}... let's just compute directly
    # Factor: forward edge -> (1+i)/2, backward -> (-1+i)/2
    # Product for f forward edges out of n-1:
    #   ((1+i)/2)^f * ((-1+i)/2)^{n-1-f}

    # (1+i)^f * (-1+i)^{n-1-f} / 2^{n-1}
    #
    # (1+i)^2 = 2i, (1+i)^4 = -4
    # (-1+i)^2 = -2i, (-1+i)^4 = -4
    #
    # Let's compute (1+i)^f * (-1+i)^g where g = n-1-f
    # Note: -1+i = i*(1+i)/i * ... actually -1+i = i*(i+1)/i... no.
    # -1+i = -(1-i). And 1-i = conj(1+i).
    # So (-1+i)^g = (-1)^g * (1-i)^g
    #
    # (1+i)^f * (-1)^g * (1-i)^g
    # = (-1)^g * (1+i)^f * (1-i)^g
    #
    # (1+i)(1-i) = 2, so if f >= g: = (-1)^g * 2^g * (1+i)^{f-g}
    # if g > f: = (-1)^g * 2^f * (1-i)^{g-f}

    # For n=7, n-1=6. Let me just compute numerically with exact integer arithmetic.
    # W(i/2) * 2^{n-1} = sum_P (1+i)^f * (-1+i)^{n-1-f} * count
    # Split into real and imaginary parts.

    # (1+i)^k for k=0..6:
    pow_1pi = [(1,0)]  # (real, imag) of (1+i)^k
    for k in range(1, n):
        a, b = pow_1pi[-1]
        # (a+bi)(1+i) = (a-b) + (a+b)i
        pow_1pi.append((a - b, a + b))

    # (-1+i)^k for k=0..6:
    pow_m1pi = [(1,0)]
    for k in range(1, n):
        a, b = pow_m1pi[-1]
        # (a+bi)(-1+i) = (-a-b) + (a-b)i
        pow_m1pi.append((-a - b, a - b))

    # W(i/2) * 2^{n-1} = sum over f of count[f] * (1+i)^f * (-1+i)^{n-1-f}
    real_sum = 0
    imag_sum = 0
    for f, cnt in fwd_counts.items():
        g = (n-1) - f
        ar, ai = pow_1pi[f]
        br, bi = pow_m1pi[g]
        # product (ar+ai*i)(br+bi*i) = (ar*br - ai*bi) + (ar*bi + ai*br)i
        pr = ar * br - ai * bi
        pi_ = ar * bi + ai * br
        real_sum += cnt * pr
        imag_sum += cnt * pi_

    # W(i/2) = (real_sum + imag_sum * i) / 2^{n-1}
    # For odd n, W(i/2) should be real, so imag_sum should be 0 mod 2^{n-1}
    denom = 2 ** (n-1)

    return real_sum, imag_sum, denom, dict(fwd_counts)


def count_directed_cycles(T, length):
    """Count directed cycles of given odd length (by vertex set, counting all orientations)."""
    n = len(T)
    count = 0
    for verts in combinations(range(n), length):
        first = verts[0]
        for perm in permutations(verts[1:]):
            path = (first,) + perm
            valid = True
            for i in range(length):
                if not T[path[i]][path[(i + 1) % length]]:
                    valid = False
                    break
            if valid:
                count += 1
    return count


def compute_invariants_n7(T):
    """Compute (c3, c5, c7, alpha_1, alpha_2, bc) for n=7 tournament."""
    n = len(T)
    assert n == 7

    # Count directed 3-cycles (vertex sets with at least one directed 3-cycle)
    # Actually we need: number of directed 3-cycles (each oriented cycle counted once)
    c3_directed = 0
    c3_vertex_sets = 0
    for verts in combinations(range(7), 3):
        a, b, c = verts
        # Check both orientations
        if T[a][b] and T[b][c] and T[c][a]:
            c3_directed += 1
            c3_vertex_sets += 1
        elif T[a][c] and T[c][b] and T[b][a]:
            c3_directed += 1
            c3_vertex_sets += 1

    # Count directed 5-cycles
    c5_directed = count_directed_cycles(T, 5)

    # Count directed 7-cycles (Hamiltonian cycles)
    c7_directed = count_directed_cycles(T, 7)

    # For the OCF: alpha_k = number of independent sets of size k in Omega(T)
    # We need the conflict graph
    # First get all directed odd cycles
    all_cycles = []
    for length in [3, 5, 7]:
        for verts in combinations(range(7), length):
            first = verts[0]
            for perm in permutations(verts[1:]):
                path = (first,) + perm
                valid = True
                for i in range(length):
                    if not T[path[i]][path[(i + 1) % length]]:
                        valid = False
                        break
                if valid:
                    all_cycles.append(path)

    num_cycles = len(all_cycles)

    # Build conflict graph
    vsets = [set(c) for c in all_cycles]

    # Count independent sets of size 0,1,2,...
    # alpha_0 = 1, alpha_1 = num_cycles, alpha_2 = number of non-conflicting pairs
    alpha_1 = num_cycles

    # Count pairs that are vertex-disjoint
    alpha_2 = 0
    for i in range(num_cycles):
        for j in range(i+1, num_cycles):
            if not (vsets[i] & vsets[j]):
                alpha_2 += 1

    return {
        'c3': c3_directed,
        'c5': c5_directed,
        'c7': c7_directed,
        'alpha_1': alpha_1,
        'alpha_2': alpha_2,
        'num_cycles': num_cycles
    }


def main():
    n = 7
    m = n * (n - 1) // 2  # 21
    total = 1 << m  # 2^21 = 2097152

    print(f"Exhaustive enumeration of all {total} tournaments on n={n} vertices")
    print("=" * 70)

    # Collect data: for each tournament, store (H, W_ihalf_real, W_ihalf_imag)
    H_values = []
    W_real_values = []
    W_imag_values = []
    fwd_distributions = []

    # Also collect by H value
    H_to_W = defaultdict(list)

    # For invariant analysis, sample a subset
    invariant_samples = {}  # H -> list of invariant dicts

    step = max(1, total // 20)

    for bits in range(total):
        if bits % step == 0:
            print(f"  Progress: {bits}/{total} ({100*bits/total:.0f}%)", file=sys.stderr)

        T = tournament_from_bits(n, bits)
        H = hamiltonian_path_count_dp(T)
        real_s, imag_s, denom, fwd = compute_W_ihalf(T)

        H_values.append(H)
        W_real_values.append(real_s)
        W_imag_values.append(imag_s)

        # W(i/2) = real_s / denom (+ imag_s*i / denom, should be 0 for odd n)
        H_to_W[H].append((real_s, imag_s))

    denom = 2 ** (n - 1)  # 64 for n=7

    print(f"\nTotal tournaments processed: {len(H_values)}")

    # =====================================================================
    # CHECK 1: Is W(i/2) always real for odd n?
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 1: Is W(i/2) always real (imag part = 0)?")
    non_real = sum(1 for im in W_imag_values if im != 0)
    print(f"  Non-real count: {non_real} / {total}")
    if non_real == 0:
        print("  CONFIRMED: W(i/2) is always real for n=7")

    # =====================================================================
    # CHECK 2: W(1/2) = H(T)?
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 2: Verify W(1/2) = H(T)")
    # W(1/2) = sum_P prod (1/2 + s_k). Forward: 1/2+1/2=1. Backward: 1/2-1/2=0.
    # So W(1/2) = number of perms with ALL forward edges = H(T).
    # We don't compute W(1/2) separately here, just note it's by definition.
    print("  TRUE by definition (each factor is 0 or 1)")

    # =====================================================================
    # CHECK 3: Distribution of H and W(i/2)
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 3: Joint distribution of H and W(i/2)")

    H_set = sorted(set(H_values))
    print(f"  Distinct H values: {len(H_set)}")
    print(f"  H range: {min(H_set)} to {max(H_set)}")

    print(f"\n  {'H':>5}  {'count':>7}  {'W*{}'.format(denom):>12}  {'W*{} range'.format(denom):>25}  {'|W|^2*{}'.format(denom**2):>15}  {'sqrt?':>6}")
    print(f"  {'---':>5}  {'-----':>7}  {'--------':>12}  {'----------':>25}  {'---------':>15}  {'-----':>6}")

    for H in H_set:
        entries = H_to_W[H]
        count = len(entries)
        W_reals = [e[0] for e in entries]
        W_set = sorted(set(W_reals))
        W_min, W_max = min(W_reals), max(W_reals)

        # |W(i/2)|^2 * denom^2 = W_real^2 (since imag=0)
        W2_values = sorted(set(r*r for r in W_reals))

        # Check if W^2 / denom^2 is a perfect square
        # W(i/2) = real_s / denom, so W(i/2)^2 = real_s^2 / denom^2
        is_int = all(r % denom == 0 for r in W_reals)

        if len(W_set) <= 5:
            w_str = str([w // denom if w % denom == 0 else f"{w}/{denom}" for w in W_set])
        else:
            w_str = f"[{W_min//denom if W_min%denom==0 else W_min}..{W_max//denom if W_max%denom==0 else W_max}]"

        w2_str = str(sorted(set(r*r // (denom*denom) if r*r % (denom*denom) == 0 else -1 for r in W_reals)))

        # Is W always divisible by denom?
        div_check = "Y" if is_int else "N"

        print(f"  {H:>5}  {count:>7}  {W_set[0] if len(W_set)==1 else '...':>12}  {w_str:>25}  {w2_str:>15}  {div_check:>6}")

    # =====================================================================
    # CHECK 4: Is W(i/2) always an integer (i.e., real_s divisible by denom)?
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 4: Is W(i/2) always an integer?")
    non_int = sum(1 for r in W_real_values if r % denom != 0)
    print(f"  Non-integer count: {non_int} / {total}")

    # What about divisibility by smaller powers?
    for d in [1, 2, 4, 8, 16, 32, 64]:
        div_count = sum(1 for r in W_real_values if r % d == 0)
        print(f"  Divisible by {d}: {div_count} / {total}")

    # =====================================================================
    # CHECK 5: Is the relationship monotone?
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 5: Is |W(i/2)| anti-monotone with H?")

    # For each pair of distinct H values, check if max|W| for higher H <= min|W| for lower H
    monotone_violations = 0
    for i, H1 in enumerate(H_set):
        for j, H2 in enumerate(H_set):
            if H1 >= H2:
                continue
            # H1 < H2. We want |W(T1)| >= |W(T2)| for all T1 with H=H1, T2 with H=H2
            absW1 = [abs(e[0]) for e in H_to_W[H1]]
            absW2 = [abs(e[0]) for e in H_to_W[H2]]
            min_absW1 = min(absW1)
            max_absW2 = max(absW2)
            if max_absW2 > min_absW1:
                monotone_violations += 1

    print(f"  Monotonicity violations (H-level): {monotone_violations}")

    # More detailed: show the |W| ranges per H
    print(f"\n  {'H':>5}  {'|W|*64 min':>12}  {'|W|*64 max':>12}  {'|W| min':>10}  {'|W| max':>10}")
    for H in H_set:
        absW = [abs(e[0]) for e in H_to_W[H]]
        mn, mx = min(absW), max(absW)
        print(f"  {H:>5}  {mn:>12}  {mx:>12}  {mn/denom:>10.2f}  {mx/denom:>10.2f}")

    # =====================================================================
    # CHECK 6: Correlation coefficient
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 6: Correlation between H and |W(i/2)|")

    N = len(H_values)
    absW_vals = [abs(W_real_values[i]) / denom for i in range(N)]

    mean_H = sum(H_values) / N
    mean_absW = sum(absW_vals) / N

    cov = sum((H_values[i] - mean_H) * (absW_vals[i] - mean_absW) for i in range(N)) / N
    var_H = sum((H_values[i] - mean_H)**2 for i in range(N)) / N
    var_absW = sum((absW_vals[i] - mean_absW)**2 for i in range(N)) / N

    corr = cov / (var_H**0.5 * var_absW**0.5) if var_H > 0 and var_absW > 0 else 0
    print(f"  Corr(H, |W(i/2)|) = {corr:.6f}")

    # Also with W(i/2) (signed)
    W_vals = [W_real_values[i] / denom for i in range(N)]
    mean_W = sum(W_vals) / N
    cov2 = sum((H_values[i] - mean_H) * (W_vals[i] - mean_W) for i in range(N)) / N
    var_W = sum((W_vals[i] - mean_W)**2 for i in range(N)) / N
    corr2 = cov2 / (var_H**0.5 * var_W**0.5) if var_H > 0 and var_W > 0 else 0
    print(f"  Corr(H, W(i/2))  = {corr2:.6f}")

    # Correlation with W^2
    W2_vals = [W_real_values[i]**2 / denom**2 for i in range(N)]
    mean_W2 = sum(W2_vals) / N
    cov3 = sum((H_values[i] - mean_H) * (W2_vals[i] - mean_W2) for i in range(N)) / N
    var_W2 = sum((W2_vals[i] - mean_W2)**2 for i in range(N)) / N
    corr3 = cov3 / (var_H**0.5 * var_W2**0.5) if var_H > 0 and var_W2 > 0 else 0
    print(f"  Corr(H, W(i/2)^2) = {corr3:.6f}")

    # =====================================================================
    # CHECK 7: W(i/2) in terms of forward-edge distribution
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 7: W(i/2) formula in terms of N_f (perms with f forward edges)")
    print()
    print("  For n=7, n-1=6. Let N_f = #{perms with f forward edges}.")
    print("  H = N_6 (all 6 edges forward = Hamiltonian path).")
    print()

    # Compute (1+i)^f * (-1+i)^{6-f} for each f
    pow_1pi = [(1,0)]
    for k in range(1, 7):
        a, b = pow_1pi[-1]
        pow_1pi.append((a - b, a + b))

    pow_m1pi = [(1,0)]
    for k in range(1, 7):
        a, b = pow_m1pi[-1]
        pow_m1pi.append((-a - b, a - b))

    print(f"  {'f':>3}  {'(1+i)^f':>15}  {'(-1+i)^{6-f}':>15}  {'product':>15}  {'product/8':>10}")
    for f in range(7):
        g = 6 - f
        ar, ai = pow_1pi[f]
        br, bi = pow_m1pi[g]
        pr = ar * br - ai * bi
        pi_ = ar * bi + ai * br
        print(f"  {f:>3}  {ar:>6}+{ai:>4}i  {br:>6}+{bi:>4}i  {pr:>6}+{pi_:>4}i  {pr//8 if pr%8==0 else pr/8:>10}")

    print()
    print("  So W(i/2) * 64 = sum_f N_f * product[f]")
    print("  W(i/2) = (1/64) * sum_f N_f * coefficient[f]")

    # Let's show the exact coefficients
    # For n=7, the real part coefficients (imag should cancel):
    coeffs_real = []
    coeffs_imag = []
    for f in range(7):
        g = 6 - f
        ar, ai = pow_1pi[f]
        br, bi = pow_m1pi[g]
        pr = ar * br - ai * bi
        pi_ = ar * bi + ai * br
        coeffs_real.append(pr)
        coeffs_imag.append(pi_)

    print(f"\n  Real coefficients: {coeffs_real}")
    print(f"  Imag coefficients: {coeffs_imag}")
    print(f"\n  W(i/2) = (1/64) * (", end="")
    terms = []
    for f in range(7):
        if coeffs_real[f] != 0:
            terms.append(f"{coeffs_real[f]}*N_{f}")
    print(" + ".join(terms), ")")

    # =====================================================================
    # CHECK 8: Detailed table for each H value
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 8: For each H, show all possible W(i/2) values")

    for H in H_set:
        entries = H_to_W[H]
        W_reals_scaled = sorted(set(e[0] for e in entries))
        W_actual = [w / denom for w in W_reals_scaled]
        W2_actual = [w**2 / denom**2 for w in W_reals_scaled]

        # Check if W^2 values are perfect squares
        import math
        W2_perfect = []
        for w2 in W2_actual:
            if w2 == int(w2):
                sq = int(round(math.sqrt(int(w2))))
                W2_perfect.append(sq * sq == int(w2))
            else:
                W2_perfect.append(False)

        n_entries = len(entries)
        w_str = ", ".join(f"{w:.1f}" for w in W_actual)
        w2_str = ", ".join(f"{w2:.1f}({'PS' if ps else 'no'})" for w2, ps in zip(W2_actual, W2_perfect))

        print(f"  H={H:>3} ({n_entries:>6} tours):  W(i/2) in {{{w_str}}}")
        print(f"       W(i/2)^2 in {{{w2_str}}}")

    # =====================================================================
    # CHECK 9: H + something * W(i/2) = constant?
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 9: Linear relations between H and W(i/2)")

    # Check H + W(i/2): is it always the same?
    for H in H_set:
        entries = H_to_W[H]
        sums = set(H * denom + e[0] for e in entries)
        diffs = set(H * denom - e[0] for e in entries)
        if len(sums) == 1:
            print(f"  H={H}: H*64 + W*64 = {sums.pop()} (CONSTANT!)")
        if len(diffs) == 1 and len(entries) > 1:
            print(f"  H={H}: H*64 - W*64 = {diffs.pop()} (CONSTANT!)")

    # Check H^2 + W^2
    for H in H_set:
        entries = H_to_W[H]
        sq_sums = set(H**2 * denom**2 + e[0]**2 for e in entries)
        if len(sq_sums) == 1 and len(entries) > 1:
            print(f"  H={H}: H^2*64^2 + W^2*64^2 = {sq_sums.pop()} (CONSTANT!)")

    # =====================================================================
    # CHECK 10: Exact formula W(i/2) in terms of N_f
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 10: Verify W(i/2) = (1/64) * sum of real_coeff[f] * N_f")
    print("  and check the specific formula")

    # For a sample tournament, verify
    T_sample = tournament_from_bits(7, 0)  # transitive tournament
    H_sample = hamiltonian_path_count_dp(T_sample)
    real_s, imag_s, d, fwd = compute_W_ihalf(T_sample)
    print(f"\n  Transitive tournament (bits=0):")
    print(f"    H = {H_sample}")
    print(f"    Forward-edge distribution: {dict(sorted(fwd.items()))}")
    print(f"    W(i/2) * 64 = {real_s}")
    print(f"    W(i/2) = {real_s / d}")

    # Try the Paley-like tournament (regular, bits for max H)
    # We already computed all - find the max H
    max_H = max(H_values)
    max_bits = [i for i, h in enumerate(H_values) if h == max_H]
    print(f"\n  Max-H tournament (H={max_H}, first bits={max_bits[0]}):")
    T_max = tournament_from_bits(7, max_bits[0])
    real_s2, imag_s2, d2, fwd2 = compute_W_ihalf(T_max)
    print(f"    Forward-edge distribution: {dict(sorted(fwd2.items()))}")
    print(f"    W(i/2) * 64 = {real_s2}")
    print(f"    W(i/2) = {real_s2 / d2}")

    # =====================================================================
    # CHECK 11: Exact identity H^2 + W(i/2)^2 = ???
    # or H + |W(i/2)| = ???
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 11: Looking for H^2 + W(i/2)^2 = f(invariants)")

    # Compute N_0 + N_1 + ... + N_6 = n! for each tournament (sanity)
    # And check various combinations

    # Actually let me compute N_f for a few tournaments and see the pattern
    print("\n  Checking sum N_f = n! and sum_f f*N_f patterns")

    for bits in [0, max_bits[0], 1, 1000, 100000]:
        T_check = tournament_from_bits(7, bits)
        H_check = hamiltonian_path_count_dp(T_check)
        _, _, _, fwd_check = compute_W_ihalf(T_check)
        Nf = [fwd_check.get(f, 0) for f in range(7)]
        total_perms = sum(Nf)
        sum_fNf = sum(f * Nf[f] for f in range(7))
        # H = N_6
        # W(i/2)*64 = coeffs . N
        W_scaled = sum(coeffs_real[f] * Nf[f] for f in range(7))

        # Also compute W(0) = sum_P prod s_k = sum_P prod (T[P(k),P(k+1)] - 1/2)
        # = sum_P (1/2)^6 * prod (2*T[P(k),P(k+1)] - 1)
        # = (1/64) * sum_P (-1)^{backward edges} = (1/64) * sum_f N_f * (-1)^{6-f} * 1^f
        # = (1/64) * sum_f N_f * (-1)^{6-f}
        W0_scaled = sum(Nf[f] * ((-1)**(6-f)) for f in range(7))

        print(f"\n  bits={bits}: H={H_check}, N_f={Nf}")
        print(f"    sum N_f = {total_perms} (should be {5040})")
        print(f"    sum f*N_f = {sum_fNf}")
        print(f"    W(i/2)*64 = {W_scaled}")
        print(f"    W(0)*64 = {W0_scaled}")
        print(f"    H^2 + (W(i/2))^2 = {H_check**2 + (W_scaled/64)**2}")

    # =====================================================================
    # CHECK 12: Is W(i/2) related to W(-i/2)?
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 12: W(i/2) vs W(-i/2) and W(i/2) vs W_Top(i/2)")

    # For T^op, forward edges become backward, so N_f(T^op) = N_{6-f}(T)
    # W_Top(i/2)*64 = sum_f N_{6-f} * coeffs_real[f]
    # = sum_g N_g * coeffs_real[6-g]

    # Check: is coeffs_real[6-f] related to coeffs_real[f]?
    print(f"  coeffs_real[f]:   {coeffs_real}")
    print(f"  coeffs_real[6-f]: {[coeffs_real[6-f] for f in range(7)]}")

    # For W(-i/2): replace i by -i. Factors become (-i/2 + 1/2) = (1-i)/2 and (-i/2 - 1/2) = (-1-i)/2
    # (1-i)^f * (-1-i)^{6-f} = conj((1+i)^f * (-1+i)^{6-f}) since conjugation
    # So W(-i/2) = conj(W(i/2)). For real W(i/2), this means W(-i/2) = W(i/2).
    print("\n  W(-i/2) = conj(W(i/2)) = W(i/2) when W is real. (Trivially holds for odd n.)")

    # =====================================================================
    # CHECK 13: Exact quadratic form
    # =====================================================================
    print("\n" + "=" * 70)
    print("CHECK 13: Quadratic form W(i/2)^2 in terms of (N_0,...,N_6)")

    # W(i/2)*64 = sum_f c_f * N_f  where c_f = coeffs_real[f]
    # W(i/2)^2 * 64^2 = (sum_f c_f * N_f)^2 = sum_{f,g} c_f * c_g * N_f * N_g
    # This is exact but not very illuminating.

    # More useful: express in terms of "graph invariants"
    # N_6 = H (Hamiltonian paths)
    # N_0 = H(T^op) (all backward = Hamiltonian path in opposite)
    # sum N_f = 7! = 5040
    # What is sum f * N_f? This counts total forward edges across all permutations.
    # Each permutation has 6 edges. For each edge (u,v) in a permutation,
    # forward means T[u][v]=1.
    # Total forward edges = sum over all ordered pairs (u,v) of T[u][v] * #{perms where u,v consecutive with u before v}
    # Number of perms with u at position k and v at position k+1 = 6 * 5! = 720...
    # Actually: for a specific ordered pair (u,v), #perms with u immediately before v = 6! = 720?
    # No: fix u at position k, v at position k+1. Choose remaining 5 in remaining 5 positions: 5!.
    # There are 6 possible k values. So 6 * 5! = 4320. Wait, that's n-1 * (n-1)! = 6*720=4320.
    # Hmm, n=7, n-1=6 consecutive pairs, each pair (u,v) appears in 6*(7-2)! permutations?
    # No. Fix u at position k (0-indexed, k=0..5), v at position k+1.
    # Remaining 5 vertices in 5! = 120 ways. 6 choices of k. So 720 permutations.
    # Total ordered pairs (u,v) with u!=v: 42.
    # So sum_{u,v,u!=v} 720 = 42 * 720 = 30240 = 6 * 5040 = 6 * 7!. Checks out.
    # sum f * N_f = sum_perms (forward edges) = 720 * sum_{u!=v} T[u][v] = 720 * E
    # where E = number of directed edges = C(7,2) = 21 (every tournament has 21 edges).
    # So sum f * N_f = 720 * 21 = 15120 = 3 * 5040. This is constant!

    print("  sum f * N_f = 720 * 21 = 15120 (CONSTANT for all n=7 tournaments)")
    print("  sum N_f = 5040 (trivially constant)")
    print("  So N_0 + N_1 + ... + N_6 = 5040")
    print("  0*N_0 + 1*N_1 + ... + 6*N_6 = 15120")

    # Higher moments: sum f^2 * N_f, etc.?
    # sum f^2 * N_f = sum_perms (sum_{edges} forward)^2 = sum_perms (sum_k indicator_k)^2
    # = sum_perms sum_k indicator_k^2 + sum_{k!=l} indicator_k * indicator_l
    # = sum f * N_f + sum_perms sum_{k!=l} indicator_k * indicator_l
    # The cross terms depend on the tournament.

    # Let's verify the f^2 moment
    print("\n  Checking sum f^2 * N_f across tournaments:")
    f2_values = []
    for bits in range(0, total, total//10):
        T_check = tournament_from_bits(7, bits)
        _, _, _, fwd_check = compute_W_ihalf(T_check)
        Nf = [fwd_check.get(f, 0) for f in range(7)]
        f2 = sum(f*f * Nf[f] for f in range(7))
        f2_values.append(f2)
    print(f"    f^2 values (sample): {f2_values}")
    print(f"    NOT constant: varies with tournament")

    # =====================================================================
    # KEY ANALYSIS: Express W(i/2) in simplified form
    # =====================================================================
    print("\n" + "=" * 70)
    print("KEY ANALYSIS: Simplified W(i/2) formula")

    # coeffs_real for n=7: let's display them
    print(f"\n  c_f (real part of (1+i)^f * (-1+i)^(6-f)):")
    for f in range(7):
        print(f"    f={f}: c_f = {coeffs_real[f]}")

    # W(i/2) = (1/64) * sum c_f * N_f
    # Using constraints: sum N_f = 5040, sum f*N_f = 15120
    # We can eliminate two N_f values and express W in terms of the rest.

    # Let A_f = N_f for the forward-edge count distribution
    # H = N_6, H_op = N_0
    # The formula gives us the EXACT dependence.

    print("\n  Formula: 64 * W(i/2) = " + " + ".join(f"({coeffs_real[f]})*N_{f}" for f in range(7)))

    # Simplify: group by coefficient value
    print("\n  Grouping:")
    coeff_groups = defaultdict(list)
    for f in range(7):
        coeff_groups[coeffs_real[f]].append(f)
    for c, fs in sorted(coeff_groups.items()):
        print(f"    Coefficient {c:>4}: N_" + " + N_".join(str(f) for f in fs))


if __name__ == "__main__":
    main()
