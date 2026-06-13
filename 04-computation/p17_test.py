"""
p17_test.py — Test the Hessian sign at p=17 (mod 4 = 1).

Prediction: lambda_H > 0 (center = local min), consistent with p=13.
This would confirm the mod 4 dichotomy pattern.

Strategy: Only compute H for orbit representatives (not all 256 tournaments).
Z_17^* acts on connection sets; orbits share the same H.

Author: kind-pasteur-2026-03-12-S56c
"""

import sys
import cmath
import math
import time
from itertools import combinations
from collections import defaultdict

sys.path.insert(0, '04-computation')


def all_circulant_tournaments(n):
    pairs, used = [], set()
    for a in range(1, n):
        if a not in used:
            b = n - a
            if a == b: return []
            pairs.append((a, b)); used.add(a); used.add(b)
    results = []
    for bits in range(2 ** len(pairs)):
        S = [a if (bits >> i) & 1 else b for i, (a, b) in enumerate(pairs)]
        results.append(tuple(sorted(S)))
    return results


def ham_count_dp(n, S):
    """Held-Karp DP for Hamiltonian path count."""
    S_set = set(S)
    adj = [0] * n  # bitmask adjacency
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in S_set:
                adj[i] |= (1 << j)

    full_mask = (1 << n) - 1
    # Use dict for sparse DP
    dp = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        popcount = bin(mask).count('1')
        if popcount >= n:
            continue
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            # Neighbors of v not in mask
            candidates = adj[v] & ~mask
            w = 0
            while candidates:
                if candidates & 1:
                    dp[(mask | (1 << w), w)] += cnt
                candidates >>= 1
                w += 1

    return sum(dp.get((full_mask, v), 0) for v in range(n))


def circulant_eigenvalues(n, S):
    omega = cmath.exp(2j * cmath.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def elementary_symmetric(xs, k):
    if k == 0: return 1.0
    if k > len(xs): return 0.0
    total = 0.0
    for subset in combinations(range(len(xs)), k):
        prod = 1.0
        for i in subset:
            prod *= xs[i]
        total += prod
    return total


def multiplicative_orbit(p, S):
    """Compute the Z_p^* orbit of connection set S."""
    S_set = frozenset(S)
    orbit = {S_set}
    for a in range(2, p):
        if math.gcd(a, p) != 1:
            continue
        new_S = frozenset((a * s) % p for s in S)
        orbit.add(new_S)
    return orbit


def find_orbit_representatives(p):
    """Find one representative per Z_p^* orbit of connection sets."""
    all_S = all_circulant_tournaments(p)
    seen = set()
    reps = []
    for S in all_S:
        S_frozen = frozenset(S)
        if S_frozen in seen:
            continue
        orbit = multiplicative_orbit(p, S)
        reps.append(tuple(sorted(S)))
        for member in orbit:
            seen.add(member)
    return reps


def solve_linear_system(A, b):
    n = len(b)
    M = [row[:] + [bi] for row, bi in zip(A, b)]
    for col in range(n):
        max_row = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[max_row][col]) < 1e-12: return None
        M[col], M[max_row] = M[max_row], M[col]
        for row in range(col + 1, n):
            f = M[row][col] / M[col][col]
            for j in range(col, n + 1): M[row][j] -= f * M[col][j]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][j] * x[j] for j in range(i + 1, n))) / M[i][i]
    return x


def main():
    p = 17
    m = (p - 1) // 2  # = 8
    print(f"=" * 70)
    print(f"p = {p} (mod 4 = {p % 4}), m = {m}")
    print(f"PREDICTION: lambda_H > 0 (center = local min)")
    print(f"=" * 70)

    # Find orbit representatives
    t0 = time.time()
    reps = find_orbit_representatives(p)
    print(f"\n  {len(reps)} orbit representatives (of {2**m} total circulants)")
    print(f"  Found in {time.time() - t0:.2f}s")

    # Compute H and spectral data for each representative
    data = []
    for i, S in enumerate(reps):
        t1 = time.time()
        H = ham_count_dp(p, S)
        elapsed = time.time() - t1
        eigs = circulant_eigenvalues(p, S)
        y2 = [eigs[k].imag ** 2 for k in range(1, m + 1)]
        esyms = [elementary_symmetric(y2, k) for k in range(m + 1)]
        orbit_size = len(multiplicative_orbit(p, S))

        data.append({
            'S': S, 'H': H, 'y2': y2, 'esyms': esyms,
            'orbit_size': orbit_size
        })
        print(f"  [{i + 1}/{len(reps)}] S={list(S)}: H={H}, "
              f"orbit={orbit_size}, time={elapsed:.1f}s")

    data.sort(key=lambda d: d['H'], reverse=True)
    H_vals = sorted(set(d['H'] for d in data), reverse=True)
    n_H = len(H_vals)

    print(f"\n  {n_H} distinct H values")
    print(f"  H range: [{H_vals[-1]}, {H_vals[0]}]")

    # Table
    print(f"\n  {'H':>15} {'orbit':>5} {'e_1':>10} {'e_2':>12} {'e_3':>12} {'e_m':>12}")
    for H_val in H_vals:
        d = next(x for x in data if x['H'] == H_val)
        print(f"  {H_val:>15} {d['orbit_size']:>5} {d['esyms'][1]:>10.2f} "
              f"{d['esyms'][2]:>12.2f} {d['esyms'][3]:>12.2f} "
              f"{d['esyms'][m]:>12.6f}")

    # Fit H = c_0 + c_2*e_2 + ... + c_m*e_m
    reps_per_H = [next(d for d in data if d['H'] == H_val) for H_val in H_vals]

    # Try increasing terms — use ALL e_2 through e_m
    for n_terms in range(1, m):
        k_start = 2
        k_end = k_start + n_terms
        n_vars = 1 + n_terms
        if n_H < n_vars:
            continue

        X = [[1.0] + [d['esyms'][k] for k in range(k_start, k_end)] for d in reps_per_H]
        y = [float(d['H']) for d in reps_per_H]

        XtX = [[sum(X[i][r] * X[i][c] for i in range(n_H)) for c in range(n_vars)] for r in range(n_vars)]
        Xty = [sum(X[i][r] * y[i] for i in range(n_H)) for r in range(n_vars)]

        coeffs = solve_linear_system(XtX, Xty)
        if coeffs is None:
            print(f"  n_terms={n_terms}: SINGULAR")
            continue

        residuals = [d['H'] - (coeffs[0] + sum(coeffs[j + 1] * d['esyms'][k_start + j]
                                           for j in range(n_terms)))
                     for d in reps_per_H]
        max_res = max(abs(r) for r in residuals)
        ss_res = sum(r**2 for r in residuals)
        ss_tot = sum((d['H'] - sum(d2['H'] for d2 in reps_per_H)/n_H)**2 for d in reps_per_H)
        r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 1
        print(f"\n  n_terms={n_terms} (e_2..e_{k_end-1}): R^2={r2:.8f}, max|res|={max_res:.2f}")

        if max_res < 1.0:
            formula = f"H = {coeffs[0]:.2f}"
            signs = []
            for j, k in enumerate(range(k_start, k_end)):
                sign = "+" if coeffs[j + 1] >= 0 else ""
                formula += f" {sign}{coeffs[j + 1]:.4f}*e_{k}"
                signs.append("+" if coeffs[j + 1] >= 0 else "-")
            print(f"\n  FIT: {formula}")
            print(f"  max |residual| = {max_res:.4f}")
            print(f"  Sign pattern: {signs}")

            # Hessian
            v = p / 4.0
            lambda_H = 0
            for j, k in enumerate(range(k_start, k_end)):
                if k >= 2:
                    binom = math.comb(m - 2, k - 2)
                    term = -coeffs[j + 1] * binom * v ** (k - 2)
                    lambda_H += term

            print(f"  lambda_H = {lambda_H:.4f}")
            if lambda_H > 0:
                print(f"  *** POSITIVE HESSIAN CONFIRMED — center is LOCAL MIN ***")
                print(f"  Consistent with p = 1 mod 4 pattern")
            else:
                print(f"  *** NEGATIVE HESSIAN — center is LOCAL MAX ***")
                print(f"  UNEXPECTED for p = 1 mod 4!")

            # Check highest coefficient sign
            c_highest = coeffs[-1]
            k_highest = k_end - 1
            print(f"\n  Highest coefficient: c_{k_highest} = {c_highest:.4f} "
                  f"({'POSITIVE' if c_highest > 0 else 'NEGATIVE'})")
            break

    # QR analysis
    qr = sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})
    qr_valid = (p - 1) not in qr
    print(f"\n  QR_{p} = {qr}, valid = {qr_valid}")
    if qr_valid:
        qr_result = next((d for d in data if set(d['S']) == set(qr)), None)
        if qr_result:
            print(f"  QR: H = {qr_result['H']}, max: {qr_result['H'] == H_vals[0]}")

    # Cyclic interval
    ci = tuple(range(p - m, p))
    ci_result = next((d for d in data if d['S'] == ci), None)
    if ci_result:
        print(f"  CI: H = {ci_result['H']}, max: {ci_result['H'] == H_vals[0]}")

    # Flattest
    flattest = min(data, key=lambda d: max(d['y2']) - min(d['y2']))
    print(f"  Flattest: H = {flattest['H']}, max: {flattest['H'] == H_vals[0]}, "
          f"spread = {max(flattest['y2']) - min(flattest['y2']):.4f}")


if __name__ == '__main__':
    main()
