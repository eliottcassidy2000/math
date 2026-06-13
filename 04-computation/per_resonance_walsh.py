"""
per_resonance_walsh.py -- kind-pasteur-2026-03-13-S60

Decompose H Walsh coefficients into per-resonance-level contributions.

The cascade analysis (resonance_cascade_p7_p11_p19_p23.out) showed that
at each resonance level q, the sign follows: sign*chi(ab) = chi(q).

This script verifies that the TOTAL Walsh coefficient is a weighted sum
of per-resonance contributions, each with the chi(q) sign law.

Key formula to verify:
  h_hat_H[{a,b}] = sum_q R_q * chi(q) * chi(ab)

where R_q depends on q but NOT on the specific pair (a,b) within that
resonance class.
"""

import numpy as np
from itertools import combinations, permutations
from math import gcd

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def build_circulant(p, S):
    n = p
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for s in S:
            A[i, (i + s) % n] = 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n + 1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u, v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_k_cycles(A, n, k):
    """Count directed k-cycles using traces."""
    A_np = np.array(A, dtype=np.float64)
    Ak = np.linalg.matrix_power(A_np, k)
    return int(round(np.trace(Ak))) // k

def compute_all_data(p):
    """Compute H, c_k for all 2^m orientations."""
    m = (p - 1) // 2
    N = 1 << m
    data = []

    for bits in range(N):
        S = []
        for k in range(m):
            if bits & (1 << k):
                S.append(k + 1)
            else:
                S.append(p - (k + 1))
        A = build_circulant(p, S)
        H = count_ham_paths(A, p)

        # Count cycles of each odd length
        cycles = {}
        for k in range(3, p + 1, 2):
            cycles[k] = count_directed_k_cycles(A, p, k)

        data.append({'bits': bits, 'S': S, 'H': H, 'cycles': cycles})

    return data

def walsh_coeff(values, a, b, m):
    """Compute Walsh coefficient h_hat[{a,b}]."""
    N = 1 << m
    total = 0.0
    for bits in range(N):
        sa = 1 if (bits & (1 << a)) else -1
        sb = 1 if (bits & (1 << b)) else -1
        total += values[bits] * sa * sb
    return total / N

def find_all_resonances(a, b, p):
    """Find ALL resonances q for chord pair (a+1, b+1).

    A resonance at q means: q*(a+1) = +/-(b+1) mod p
    or equivalently, (b+1)*(a+1)^{-1} = +/-q mod p.
    """
    ca, cb = a + 1, b + 1
    ca_inv = pow(ca, p - 2, p)

    q_pos = (cb * ca_inv) % p  # q such that q*ca = cb
    q_neg = ((-cb) * ca_inv) % p  # q such that q*ca = -cb

    # Reduce to range [1, m]
    m = (p - 1) // 2
    resonances = set()
    for q in [q_pos, q_neg]:
        q_red = q if q <= m else p - q
        if 1 <= q_red <= m:
            resonances.add(q_red)

    return sorted(resonances)

def per_cycle_walsh_decomposition(p, data):
    """Decompose the Walsh coefficient of H into per-cycle-length contributions.

    H = 1 + sum_k 2^{alpha_k equivalent} contributions from c_k cycle counts.
    But more directly: H = I(Omega(T), 2) = sum_{j>=0} 2^j * alpha_j.

    For Walsh analysis, we decompose h_hat_H into contributions from
    each cycle length k = 3, 5, 7, ..., p.
    """
    m = (p - 1) // 2
    N = 1 << m

    H_vals = [d['H'] for d in data]

    print(f"\n{'='*70}")
    print(f"PER-CYCLE-LENGTH WALSH DECOMPOSITION at p={p}")
    print(f"{'='*70}")

    # For each cycle length k, compute h_hat of c_k
    for a in range(min(m, 3)):
        for b in range(a + 1, min(m, 4)):
            ca, cb = a + 1, b + 1
            chi_ab = legendre(ca * cb, p)
            resonances = find_all_resonances(a, b, p)

            h_H = walsh_coeff(H_vals, a, b, m)

            print(f"\n  Pair ({ca},{cb}), chi(ab)={'+' if chi_ab>0 else '-'}, "
                  f"resonances q={resonances}")
            print(f"    h_hat_H = {h_H:.4f}")

            # Decompose by cycle length
            total_check = 0
            for k in range(3, p + 1, 2):
                ck_vals = [d['cycles'][k] for d in data]
                h_ck = walsh_coeff(ck_vals, a, b, m)
                if abs(h_ck) > 1e-6:
                    sign_ck = 1 if h_ck > 0 else -1
                    print(f"    c_{k}: h_hat = {h_ck:12.4f}, sign = {'+' if sign_ck>0 else '-'}, "
                          f"sign*chi(ab) = {'+' if sign_ck*chi_ab>0 else '-'}")

def overlap_from_cycle_counts(p, data):
    """Use THM-155 to compute overlap counts from cycle counts.

    c5 + 2*ov2 = p(p^2-1)(p^2-9)/160 for regular tournaments.
    So ov2 = (const - c5) / 2.
    Also: c3 is constant, ov1 = C(c3, 2) - ov2 - disj3.
    And K = c5 - 2*ov1 - 2*ov2 = const.
    """
    m = (p - 1) // 2
    c3 = p * (p**2 - 1) // 24  # constant for regular
    target = p * (p**2 - 1) * (p**2 - 9) // 160

    print(f"\n{'='*70}")
    print(f"OVERLAP ANALYSIS via THM-155 at p={p}")
    print(f"{'='*70}")
    print(f"  c3 = {c3} (constant for all regular)")
    print(f"  c5 + 2*ov2 = {target}")
    print(f"  K = c5 - 2*ov1 - 2*ov2 = {-3*p*(p**2-1)*(p**2-9)//320}")

    for d in data[:6]:
        c5 = d['cycles'][5]
        ov2 = (target - c5) // 2
        K = -3 * p * (p**2 - 1) * (p**2 - 9) // 320
        ov1 = (c5 - K - 2 * ov2) // 2
        disj3 = c3 * (c3 - 1) // 2 - ov1 - ov2
        print(f"  bits={d['bits']:05b}: c5={c5}, ov2={ov2}, ov1={ov1}, disj3={disj3}, "
              f"H={d['H']}")

def magnitude_formula_test(p, data):
    """Test if magnitudes follow |h_hat| = p * f(q) for some function f."""
    m = (p - 1) // 2
    N = 1 << m

    H_vals = [d['H'] for d in data]

    print(f"\n{'='*70}")
    print(f"MAGNITUDE FORMULA TEST at p={p}")
    print(f"{'='*70}")

    by_q = {}
    for a in range(m):
        for b in range(a + 1, m):
            h = walsh_coeff(H_vals, a, b, m)
            if abs(h) < 1e-6:
                continue
            resonances = find_all_resonances(a, b, p)
            q_key = tuple(resonances)
            if q_key not in by_q:
                by_q[q_key] = []
            by_q[q_key].append({'a': a, 'b': b, 'h': h})

    for q_key in sorted(by_q.keys()):
        entries = by_q[q_key]
        mags = [abs(e['h']) for e in entries]
        mag = mags[0]
        all_eq = all(abs(v - mag) < 1e-6 for v in mags)

        # Check if magnitude is p * integer / power_of_2
        ratio = mag / p
        # Try ratio * 2^k for k=0..10
        found = False
        for k in range(11):
            val = ratio * (2**k)
            if abs(val - round(val)) < 1e-4:
                print(f"  resonances q={q_key}: count={len(entries)}, "
                      f"|h|/p = {ratio:.4f} = {int(round(val))}/2^{k}, "
                      f"all_equal={all_eq}")
                found = True
                break
        if not found:
            print(f"  resonances q={q_key}: count={len(entries)}, "
                  f"|h|/p = {ratio:.6f}, all_equal={all_eq}")

        # Also check signs vs chi(q) law
        for e in entries[:2]:
            ca, cb = e['a'] + 1, e['b'] + 1
            chi_ab = legendre(ca * cb, p)
            sign_h = 1 if e['h'] > 0 else -1
            product = sign_h * chi_ab
            chi_qs = [legendre(q, p) for q in q_key]
            print(f"    ({ca},{cb}): sign={'+' if sign_h>0 else '-'}, chi(ab)={chi_ab:+d}, "
                  f"sign*chi(ab)={product:+d}, chi(q)={chi_qs}")


# Run analysis
print("=" * 70)
print("PER-RESONANCE WALSH ANALYSIS")
print("=" * 70)

# p=7
print("\n### p=7 ###")
data7 = compute_all_data(7)
per_cycle_walsh_decomposition(7, data7)
overlap_from_cycle_counts(7, data7)
magnitude_formula_test(7, data7)

# p=11
print("\n### p=11 ###")
data11 = compute_all_data(11)
per_cycle_walsh_decomposition(11, data11)
overlap_from_cycle_counts(11, data11)
magnitude_formula_test(11, data11)
