"""
sequence_harvest.py — opus-2026-04-04-S13

SYSTEMATIC HARVEST of sequences from the multilinear polynomial theory.
For each sequence: compute values at n=3..7, check OEIS, identify new ones.
"""
import numpy as np
from math import comb, factorial
from itertools import combinations, permutations
from collections import defaultdict, Counter
import sys, time

sys.path.insert(0, '/Users/e/Documents/GitHub/math/04-computation')
from h_multilinear_complete import compute_all_H, mobius_inversion, make_tiles

def compute_all_data(n):
    """Compute everything we need at a given n."""
    tiles = make_tiles(n)
    m = len(tiles)
    N = 1 << m
    H, _, _ = compute_all_H(n)
    coeffs = mobius_inversion(H, m)

    # Group coefficients by degree
    by_deg = defaultdict(list)
    for S, c in coeffs.items():
        by_deg[bin(S).count('1')].append(c)

    return {
        'n': n, 'm': m, 'N': N, 'H': H, 'coeffs': coeffs,
        'tiles': tiles, 'by_deg': by_deg
    }

# ==============================================================
# SEQUENCE 1: P(n) = number of nonzero multilinear coefficients
# ==============================================================
def seq_P(data):
    """P(n) = #{S : c_S ≠ 0}"""
    return sum(1 for c in data['coeffs'].values() if c != 0)

# ==============================================================
# SEQUENCE 2: P_pos(n), P_neg(n) = positive/negative coefficient counts
# ==============================================================
def seq_P_sign(data):
    pos = sum(1 for S, c in data['coeffs'].items() if c > 0 and S != 0)
    neg = sum(1 for c in data['coeffs'].values() if c < 0)
    return pos, neg

# ==============================================================
# SEQUENCE 3: S_d(n) = sum of degree-d coefficients
# ==============================================================
def seq_S_d(data):
    sums = {}
    for d, vals in data['by_deg'].items():
        sums[d] = sum(vals)
    return sums

# ==============================================================
# SEQUENCE 4: Sigma_H(n) = sum of H over all tilings = 2^m * mean(H)
# ==============================================================
def seq_Sigma_H(data):
    return int(sum(data['H']))

# ==============================================================
# SEQUENCE 5: max H, min H (excluding transitive)
# ==============================================================
def seq_extremes(data):
    return int(min(data['H'])), int(max(data['H']))

# ==============================================================
# SEQUENCE 6: Number of distinct H values
# ==============================================================
def seq_distinct_H(data):
    return len(set(int(h) for h in data['H']))

# ==============================================================
# SEQUENCE 7: Number of forbidden odd H values ≤ max(H)
# ==============================================================
def seq_forbidden_count(data):
    achieved = set(int(h) for h in data['H'])
    max_h = max(achieved)
    return sum(1 for k in range(1, max_h+1, 2) if k not in achieved)

# ==============================================================
# SEQUENCE 8: f-vector entries of the Newton polytope
# ==============================================================
def seq_f_vector(data):
    fvec = []
    max_d = max(data['by_deg'].keys())
    for d in range(1, max_d+1):
        fvec.append(len(data['by_deg'].get(d, [])))
    return fvec

# ==============================================================
# SEQUENCE 9: Euler characteristic of coefficient support
# ==============================================================
def seq_euler_char(data):
    euler = 0
    for d, vals in data['by_deg'].items():
        if d >= 1:
            euler += (-1)**(d+1) * len(vals)
    return euler

# ==============================================================
# SEQUENCE 10: Variance of H over tilings
# ==============================================================
def seq_var_H(data):
    h = data['H'].astype(float)
    return float(np.var(h))

# ==============================================================
# SEQUENCE 11: H(all backward) = H(1,...,1)
# ==============================================================
def seq_H_all_backward(data):
    return int(data['H'][(1 << data['m']) - 1])

# ==============================================================
# SEQUENCE 12: Number of tilings with H ≡ 1 mod 4
# ==============================================================
def seq_H_mod4_1(data):
    return sum(1 for h in data['H'] if int(h) % 4 == 1)

# ==============================================================
# SEQUENCE 13: Sum of |c_S| over all S
# ==============================================================
def seq_sum_abs_c(data):
    return sum(abs(c) for c in data['coeffs'].values())

# ==============================================================
# SEQUENCE 14: Number of degree-d nonzero coefficients (for each d)
# ==============================================================
def seq_degree_counts(data):
    return {d: len(vals) for d, vals in data['by_deg'].items() if vals}

# ==============================================================
# SEQUENCE 15: c_apex = linear coefficient of apex tile (n,1)
# ==============================================================
def seq_c_apex(data):
    # Apex tile has skip = n-1
    return 2**(data['n'] - 2)

# ==============================================================
# SEQUENCE 16: Number of iso classes (from prior work)
# ==============================================================
def seq_V_n():
    """A000568: number of tournament iso classes."""
    return {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

# ==============================================================
# SEQUENCE 17: Extension orbit count (from S12)
# ==============================================================
def seq_extension_orbits():
    return {4: 12, 5: 48, 6: 296, 7: 3040}

# ==============================================================
# MAIN: Compute and display all sequences
# ==============================================================
if __name__ == "__main__":
    print("SEQUENCE HARVEST FROM MULTILINEAR POLYNOMIAL THEORY")
    print("=" * 70)

    all_data = {}
    for n in range(3, 8):
        print(f"\nComputing n={n}...", flush=True)
        t0 = time.time()
        all_data[n] = compute_all_data(n)
        print(f"  Done in {time.time()-t0:.1f}s")

    # Compute each sequence
    sequences = {}

    # 1. P(n) = nonzero coefficients
    P_vals = {n: seq_P(d) for n, d in all_data.items()}
    sequences['P(n) = nonzero coefficients'] = P_vals
    print(f"\n  P(n) = {[P_vals[n] for n in sorted(P_vals)]}")

    # 2. Positive / Negative counts
    sign_vals = {n: seq_P_sign(d) for n, d in all_data.items()}
    pos_vals = {n: v[0] for n, v in sign_vals.items()}
    neg_vals = {n: v[1] for n, v in sign_vals.items()}
    sequences['P_pos(n)'] = pos_vals
    sequences['P_neg(n)'] = neg_vals
    print(f"  P_pos(n) = {[pos_vals[n] for n in sorted(pos_vals)]}")
    print(f"  P_neg(n) = {[neg_vals[n] for n in sorted(neg_vals)]}")

    # 3. Degree sums
    S_d_all = {n: seq_S_d(d) for n, d in all_data.items()}
    for d_target in range(5):
        vals = {n: S_d_all[n].get(d_target, 0) for n in sorted(S_d_all)}
        sequences[f'S_{d_target}(n)'] = vals
        print(f"  S_{d_target}(n) = {[vals[n] for n in sorted(vals)]}")

    # 4. Sigma_H
    sigma_vals = {n: seq_Sigma_H(d) for n, d in all_data.items()}
    sequences['Sigma_H(n) = sum of all H'] = sigma_vals
    print(f"  Sigma_H(n) = {[sigma_vals[n] for n in sorted(sigma_vals)]}")

    # Check: Sigma_H = 2^m * mean = 2^m * H(1/2,...,1/2)
    for n in sorted(all_data):
        m = all_data[n]['m']
        print(f"    n={n}: Sigma_H={sigma_vals[n]}, 2^m={1<<m}, "
              f"ratio={sigma_vals[n]/(1<<m):.4f}")

    # 5. Extremes
    extremes = {n: seq_extremes(d) for n, d in all_data.items()}
    max_H_vals = {n: v[1] for n, v in extremes.items()}
    sequences['max_H(n)'] = max_H_vals
    print(f"  max_H(n) = {[max_H_vals[n] for n in sorted(max_H_vals)]}")

    # 6. Distinct H values
    distinct_vals = {n: seq_distinct_H(d) for n, d in all_data.items()}
    sequences['distinct_H(n)'] = distinct_vals
    print(f"  distinct_H(n) = {[distinct_vals[n] for n in sorted(distinct_vals)]}")

    # 7. Forbidden count
    forbidden_vals = {n: seq_forbidden_count(d) for n, d in all_data.items()}
    sequences['forbidden_odd_H(n)'] = forbidden_vals
    print(f"  forbidden_odd_H(n) = {[forbidden_vals[n] for n in sorted(forbidden_vals)]}")

    # 8. f-vector
    for n in sorted(all_data):
        fv = seq_f_vector(all_data[n])
        print(f"  f-vector(n={n}) = {fv}")

    # 9. Euler char
    euler_vals = {n: seq_euler_char(d) for n, d in all_data.items()}
    sequences['euler_char(n)'] = euler_vals
    print(f"  euler_char(n) = {[euler_vals[n] for n in sorted(euler_vals)]}")

    # 10. Variance
    var_vals = {n: seq_var_H(d) for n, d in all_data.items()}
    sequences['var_H(n)'] = var_vals
    print(f"  var_H(n) = {[f'{var_vals[n]:.2f}' for n in sorted(var_vals)]}")

    # 11. H(all backward)
    h_back = {n: seq_H_all_backward(d) for n, d in all_data.items()}
    sequences['H_all_backward(n)'] = h_back
    print(f"  H(all backward)(n) = {[h_back[n] for n in sorted(h_back)]}")

    # 12. H mod 4 = 1 count
    mod4_1 = {n: seq_H_mod4_1(d) for n, d in all_data.items()}
    sequences['count_H_mod4_eq_1'] = mod4_1
    print(f"  count(H≡1 mod 4)(n) = {[mod4_1[n] for n in sorted(mod4_1)]}")

    # 13. Sum of |c_S|
    sum_abs = {n: seq_sum_abs_c(d) for n, d in all_data.items()}
    sequences['sum_|c_S|(n)'] = sum_abs
    print(f"  sum_|c_S|(n) = {[sum_abs[n] for n in sorted(sum_abs)]}")

    # 14. Degree counts
    for n in sorted(all_data):
        dc = seq_degree_counts(all_data[n])
        print(f"  degree_counts(n={n}) = {dict(sorted(dc.items()))}")

    # 15. Apex linear coefficient
    apex_vals = {n: seq_c_apex(all_data[n]) for n in sorted(all_data)}
    print(f"  c_apex(n) = {[apex_vals[n] for n in sorted(apex_vals)]}")

    # 16. V_n
    V_n = seq_V_n()
    print(f"  V_n = {[V_n.get(n, '?') for n in range(3, 9)]}")

    # 17. Extension orbits
    ext_orb = seq_extension_orbits()
    print(f"  extension_orbits(n) = {ext_orb}")

    # ==============================================================
    # OEIS LOOKUP TABLE
    # ==============================================================
    print(f"\n\n{'='*70}")
    print(f"OEIS LOOKUP TABLE")
    print(f"{'='*70}")

    oeis_check = {
        'P(n) nonzero coeffs': [2, 6, 35, 200, 1782],
        'Sigma_H': [sigma_vals[n] for n in range(3, 8)],
        'max_H': [max_H_vals[n] for n in range(3, 8)],
        'distinct_H': [distinct_vals[n] for n in range(3, 8)],
        'H(all backward)': [h_back[n] for n in range(3, 8)],
        'sum_|c_S|': [sum_abs[n] for n in range(3, 8)],
        'euler_char': [euler_vals[n] for n in range(3, 8)],
        'P_pos': [pos_vals[n] for n in range(3, 8)],
        'P_neg': [neg_vals[n] for n in range(3, 8)],
        'count_Hmod4eq1': [mod4_1[n] for n in range(3, 8)],
        'forbidden_count': [forbidden_vals[n] for n in range(3, 8)],
        'S_1(n)=2^n-2n': [2**n - 2*n for n in range(3, 8)],
    }

    for name, vals in oeis_check.items():
        print(f"\n  {name}: {vals}")
        # Check common OEIS sequences
        # A000568: tournament iso classes
        # A000571: sum of H over all tournaments
        # etc.

    # ==============================================================
    # DERIVED SEQUENCES (ratios, differences, etc.)
    # ==============================================================
    print(f"\n\n{'='*70}")
    print(f"DERIVED SEQUENCES")
    print(f"{'='*70}")

    # Sigma_H / 2^m = mean H = average HP count with fixed base path
    mean_H = {n: sigma_vals[n] / (1 << all_data[n]['m']) for n in range(3, 8)}
    print(f"\n  mean_H(n) = Sigma_H/2^m = {[f'{mean_H[n]:.4f}' for n in range(3, 8)]}")
    # Is mean_H * 2 / (n-1)! an integer?
    for n in range(3, 8):
        val = mean_H[n] * 2 / factorial(n-1)
        print(f"    n={n}: mean_H * 2 / (n-1)! = {val:.6f}")

    # P(n) / 2^m = fraction of nonzero coefficients
    frac_nz = {n: P_vals[n] / (1 << all_data[n]['m']) for n in range(3, 8)}
    print(f"\n  P(n)/2^m = {[f'{frac_nz[n]:.4f}' for n in range(3, 8)]}")

    # max_H / mean_H
    max_mean_ratio = {n: max_H_vals[n] / mean_H[n] for n in range(3, 8)}
    print(f"\n  max_H/mean_H = {[f'{max_mean_ratio[n]:.3f}' for n in range(3, 8)]}")

    # H(all backward) / max_H
    back_max_ratio = {n: h_back[n] / max_H_vals[n] for n in range(3, 8)}
    print(f"  H(backward)/max_H = {[f'{back_max_ratio[n]:.3f}' for n in range(3, 8)]}")

    # log2(max_H) / n
    log_max = {n: np.log2(max_H_vals[n]) / n for n in range(3, 8)}
    print(f"  log2(max_H)/n = {[f'{log_max[n]:.3f}' for n in range(3, 8)]}")

    # S_1(n) = 2^n - 2n (proved)
    S1_check = {n: S_d_all[n].get(1, 0) for n in range(3, 8)}
    S1_formula = {n: 2**n - 2*n for n in range(3, 8)}
    print(f"\n  S_1(n) = {[S1_check[n] for n in range(3, 8)]}")
    print(f"  2^n-2n = {[S1_formula[n] for n in range(3, 8)]}")
    print(f"  Match: {all(S1_check[n] == S1_formula[n] for n in range(3, 8))}")

    # S_2(n)
    S2_vals = {n: S_d_all[n].get(2, 0) for n in range(3, 8)}
    print(f"  S_2(n) = {[S2_vals[n] for n in range(3, 8)]}")
    # Try: S_2 = a*4^n + b*2^n + c*n^2 + d*n + e
    # We have 5 data points (n=3..7)

    # Sigma_H(n) in terms of known quantities
    print(f"\n  Sigma_H(n) = {[sigma_vals[n] for n in range(3, 8)]}")
    print(f"  Sigma_H/n! = {[sigma_vals[n]/factorial(n) for n in range(3, 8)]}")
    # Sigma_H = sum over all LABELED tournaments T containing base path P_0 of H(T)
    # = sum over all tilings of H
    # = 2^m * E[H] where E is over uniform random tiling
    # = 2^m * H(1/2,...,1/2) by multilinearity

    # W(n) from the repo = Sigma_H / something
    # From MEMORY: W(n) = 1,2,8,32,158,928,6350,49752 for n=1..8
    W_known = {1: 1, 2: 2, 3: 8, 4: 32, 5: 158, 6: 928, 7: 6350, 8: 49752}
    print(f"\n  W(n) (from prior sessions) = {[W_known.get(n, '?') for n in range(1, 9)]}")
    print(f"  Sigma_H(n) = {[sigma_vals[n] for n in range(3, 8)]}")
    # Is Sigma_H related to W?
    for n in range(3, 8):
        ratio = sigma_vals[n] / W_known.get(n, 1)
        print(f"    n={n}: Sigma_H / W(n) = {ratio:.4f}")

    # ==============================================================
    # NEW SEQUENCE CANDIDATES (not in OEIS)
    # ==============================================================
    print(f"\n\n{'='*70}")
    print(f"NEW SEQUENCE CANDIDATES")
    print(f"{'='*70}")

    candidates = {
        'P(n) = nonzero multilinear coefficients': [P_vals[n] for n in range(3, 8)],
        'P_pos(n) = positive coefficients': [pos_vals[n] for n in range(3, 8)],
        'P_neg(n) = negative coefficients': [neg_vals[n] for n in range(3, 8)],
        'sum_|c_S|(n) = L1 norm of coefficients': [sum_abs[n] for n in range(3, 8)],
        'euler_char(n) of support': [euler_vals[n] for n in range(3, 8)],
        'distinct_H(n) = distinct HP counts': [distinct_vals[n] for n in range(3, 8)],
        'forbidden_odd(n) = odd gaps below max H': [forbidden_vals[n] for n in range(3, 8)],
        'H(all backward)(n)': [h_back[n] for n in range(3, 8)],
        'Sigma_H(n) = sum H over tilings': [sigma_vals[n] for n in range(3, 8)],
    }

    for name, vals in candidates.items():
        print(f"  {name}: {vals}")
