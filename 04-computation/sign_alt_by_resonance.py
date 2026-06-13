"""
sign_alt_by_resonance.py -- kind-pasteur-2026-03-13-S60

The sign alternation formula sign(h_hat_{c_k}) = (-1)^{(k-3)/2} * chi(ab)
holds perfectly at p=11 (35/35) but fails at p=13 (29/48).

This script checks whether the alternation depends on the RESONANCE CLASS q
of the pair, not just the cycle length k. Specifically, test:

    sign(h_hat_{c_k}[{a,b}]) = (-1)^{(k-q_min)/2} * chi(ab)  ?

where q_min is the minimum resonance of pair (a,b).

Or perhaps: the sign depends on both k and q modularly:
    sign = (-1)^{(k-q)/2} * chi(q) * chi(ab) ?

Also test: does the alternation hold WITHIN each resonance class
even if the global formula fails?
"""

import numpy as np
from math import gcd
from itertools import combinations

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

def find_all_resonances(a, b, p):
    ca, cb = a + 1, b + 1
    m = (p - 1) // 2
    resonances = set()
    for q in range(1, m + 1):
        if (q * ca - cb) % p == 0 or (q * ca + cb) % p == 0:
            resonances.add(q)
    return sorted(resonances)

def compute_all_cycle_walsh(p, max_k=None):
    m = (p - 1) // 2
    N = 1 << m
    if max_k is None:
        max_k = p

    all_ck = {}
    for k in range(3, max_k + 1, 2):
        all_ck[k] = [0] * N

    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))
        A = build_circulant(p, S)
        A_np = np.array(A, dtype=np.float64)
        for k in range(3, max_k + 1, 2):
            Ak = np.linalg.matrix_power(A_np, k)
            all_ck[k][bits] = int(round(np.trace(Ak))) // k

    walsh = {}
    for k in all_ck:
        walsh[k] = {}
        for a in range(m):
            for b in range(a + 1, m):
                total = 0.0
                for bits in range(N):
                    sa = 1 if (bits & (1 << a)) else -1
                    sb = 1 if (bits & (1 << b)) else -1
                    total += all_ck[k][bits] * sa * sb
                walsh[k][(a, b)] = total / N

    return walsh


def analyze_by_resonance(p, walsh):
    """Analyze sign patterns grouped by resonance class."""
    m = (p - 1) // 2

    print(f"\n{'='*70}")
    print(f"SIGN ANALYSIS BY RESONANCE CLASS at p={p}")
    print(f"{'='*70}")

    # Group pairs by their resonance class
    pair_info = {}
    for a in range(m):
        for b in range(a + 1, m):
            resonances = find_all_resonances(a, b, p)
            pair_info[(a, b)] = resonances

    # For each resonance class, track sign patterns
    res_classes = {}
    for (a, b), res in pair_info.items():
        key = tuple(res)
        if key not in res_classes:
            res_classes[key] = []
        res_classes[key].append((a, b))

    print(f"\n  Resonance classes:")
    for key in sorted(res_classes.keys()):
        pairs = res_classes[key]
        print(f"    q={key}: {len(pairs)} pairs")

    # For each class, show the sign*chi(ab) pattern by k
    print(f"\n  Sign*chi(ab) pattern by cycle length within each resonance class:")
    for key in sorted(res_classes.keys()):
        pairs = res_classes[key]
        print(f"\n    Resonance q={key}:")

        for k in sorted(walsh.keys()):
            signs = []
            for (a, b) in pairs:
                val = walsh[k][(a, b)]
                if abs(val) < 1e-6:
                    signs.append(0)
                else:
                    ca, cb = a + 1, b + 1
                    chi_ab = legendre(ca * cb, p)
                    sign_val = 1 if val > 0 else -1
                    signs.append(sign_val * chi_ab)

            nonzero = [s for s in signs if s != 0]
            if nonzero:
                all_same = all(s == nonzero[0] for s in nonzero)
                sign_str = '+' if nonzero[0] > 0 else '-'
                n_zero = signs.count(0)

                # Test formula: sign*chi(ab) = (-1)^{(k-3)/2}
                expected = (-1) ** ((k - 3) // 2)
                match_expected = (nonzero[0] == expected)

                # Test formula: sign*chi(ab) = (-1)^{(k-q_min)/2} * chi(q_min)
                q_min = key[0]
                if (k - q_min) % 2 == 0:
                    alt_expected = (-1) ** ((k - q_min) // 2) * legendre(q_min, p)
                    match_alt = (nonzero[0] == alt_expected)
                else:
                    alt_expected = None
                    match_alt = None

                print(f"      c_{k:2d}: sign*chi = {sign_str}, consistent={all_same}, "
                      f"n_nonzero={len(nonzero)}/{len(pairs)}, "
                      f"(-1)^(k-3)/2 {'MATCH' if match_expected else 'FAIL'}", end='')
                if alt_expected is not None:
                    print(f", alt(-1)^(k-q)/2*chi(q) {'MATCH' if match_alt else 'FAIL'}", end='')
                print()

def test_gauss_sum_sign(p, walsh):
    """Test: sign(h_hat_{c_k}[{a,b}]) = chi(q) * (-1)^{(k-q_min)/2} * chi(ab)

    This formula combines the resonance q, cycle length k, and pair character chi(ab).
    """
    m = (p - 1) // 2
    total = 0
    matches = 0

    print(f"\n  Testing Gauss sum sign formula at p={p}:")
    print(f"  sign = chi(q_min) * (-1)^((k-q_min)/2) * chi(ab)")

    for k in sorted(walsh.keys()):
        k_total = 0
        k_match = 0
        for (a, b), val in sorted(walsh[k].items()):
            if abs(val) < 1e-6:
                continue

            ca, cb = a + 1, b + 1
            chi_ab = legendre(ca * cb, p)
            sign_val = 1 if val > 0 else -1
            resonances = find_all_resonances(a, b, p)
            q_min = resonances[0]

            if (k - q_min) % 2 != 0:
                # k and q_min have different parity — this formula doesn't apply
                continue

            expected = legendre(q_min, p) * ((-1) ** ((k - q_min) // 2)) * chi_ab

            total += 1
            k_total += 1
            if sign_val == expected:
                matches += 1
                k_match += 1

        if k_total > 0:
            print(f"    c_{k:2d}: {k_match}/{k_total}")

    print(f"\n  Overall: {matches}/{total}")
    return matches, total


def test_alternating_sign_patterns(p, walsh):
    """Test various sign formulas to find what works."""
    m = (p - 1) // 2
    print(f"\n{'='*70}")
    print(f"COMPREHENSIVE SIGN FORMULA SEARCH at p={p}")
    print(f"{'='*70}")

    formulas = {}

    for k in sorted(walsh.keys()):
        for (a, b), val in sorted(walsh[k].items()):
            if abs(val) < 1e-6:
                continue

            ca, cb = a + 1, b + 1
            chi_ab = legendre(ca * cb, p)
            sign_val = 1 if val > 0 else -1
            observed = sign_val * chi_ab

            resonances = find_all_resonances(a, b, p)
            q_min = resonances[0]

            # Formula 1: (-1)^{(k-3)/2}
            f1 = (-1) ** ((k - 3) // 2)
            # Formula 2: chi(q_min) * (-1)^{(k-q_min)/2} if same parity
            f2 = legendre(q_min, p) * ((-1) ** ((k - q_min) // 2)) if (k - q_min) % 2 == 0 else None
            # Formula 3: (-1)^{(k-1)/2} * legendre(q_min, p)
            f3 = (-1) ** ((k - 1) // 2) * legendre(q_min, p)
            # Formula 4: (-1)^{(k+1)/2} * legendre(q_min, p)
            f4 = (-1) ** ((k + 1) // 2) * legendre(q_min, p)
            # Formula 5: (-1)^{(k-3)/2} * legendre(q_min, p)
            f5 = (-1) ** ((k - 3) // 2) * legendre(q_min, p)

            for name, f in [("(-1)^(k-3)/2", f1), ("chi(q)*(-1)^(k-q)/2", f2),
                           ("(-1)^(k-1)/2*chi(q)", f3), ("(-1)^(k+1)/2*chi(q)", f4),
                           ("(-1)^(k-3)/2*chi(q)", f5)]:
                if f is None:
                    continue
                if name not in formulas:
                    formulas[name] = [0, 0]
                formulas[name][1] += 1
                if observed == f:
                    formulas[name][0] += 1

    for name in sorted(formulas.keys(), key=lambda n: -formulas[n][0]):
        m_count, t_count = formulas[name]
        print(f"  {name:30s}: {m_count}/{t_count} = {100*m_count/t_count:.1f}%")


# Main
for p in [7, 11, 13]:
    print(f"\n{'#'*70}")
    print(f"# p = {p}")
    print(f"{'#'*70}")

    walsh = compute_all_cycle_walsh(p)
    analyze_by_resonance(p, walsh)
    test_gauss_sum_sign(p, walsh)
    test_alternating_sign_patterns(p, walsh)
