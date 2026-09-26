"""Finite positive-kernel calculus for the exact global P2 pairing tree.

Exact integer/Fraction controls; standard library only. Checks remain active -O.
Companion proof: 05-knowledge/results/crossroads10_20260926_flow.md.
"""
from array import array
from fractions import Fraction
from itertools import product
import json
from math import fsum, inf, isqrt, log, sqrt


def require(test, message):
    if not test:
        raise AssertionError(message)


def clauses(X):
    eq, le = [], []
    for i in range(2, X + 1):
        if i % 2 == 0:
            if 3 * i // 2 <= X:
                eq.append((i, 3 * i // 2))
        else:
            j, k = (3 * i - 1) // 2, (3 * i + 1) // 2
            if j <= X:
                le.append((j, i))
            if k <= X:
                le.append((i, k))
    return eq, le


def cost_from_cutoffs(X, terms):
    """terms is a list of (integer cutoff, integer coefficient)."""
    zero, one = [0] * (X + 1), [0] * (X + 1)
    for i in range(X, 1, -1):
        a, b = 0, sum(w for cutoff, w in terms if i <= cutoff)
        if i % 2 == 0:
            j = 3 * i // 2
            if j <= X:
                a += one[j]
                b += zero[j]
        else:
            j, k = (3 * i - 1) // 2, (3 * i + 1) // 2
            if j <= X:
                a += zero[j]
                b += min(zero[j], one[j])
            if k <= X:
                a += min(zero[k], one[k])
                b += one[k]
        zero[i], one[i] = a, b
    return min(zero[2], one[2]) if X >= 2 else 0


def optimum(X, a):
    return cost_from_cutoffs(X, [(2**k * X // 3**k, w) for k, w in enumerate(a)])


def convolution(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def normalizer(a):
    return sum((Fraction(w) * Fraction(2, 3)**k for k, w in enumerate(a)), Fraction(0))


def full_cost_certificates(a, H):
    """Direct fixed-root arrays; imports no Delta/toll recurrence."""
    require(all(w >= 0 for w in a) and sum(a) > 0, 'kernel domain')
    lam = normalizer(a)
    zero, one = array('q', [0]), array('q', [a[0]])
    values = [Fraction(0)]
    totals = [0]
    for h in range(1, H + 1):
        period, w = len(zero), sum(a[:h + 1])
        newzero, newone = array('q'), array('q')
        total = 0
        for i in range(2 * period):
            if i % 2 == 0:
                j = (3 * i // 2) % period
                x, y = one[j], w + zero[j]
            else:
                j, k = ((3 * i - 1) // 2) % period, ((3 * i + 1) // 2) % period
                x = zero[j] + min(zero[k], one[k])
                y = w + min(zero[j], one[j]) + one[k]
            newzero.append(x)
            newone.append(y)
            total += min(x, y)
        zero, one = newzero, newone
        totals.append(total)
        values.append(Fraction(total, 3**(h + 1)) / lam)
    return values, totals


def tv(p, q):
    n = max(len(p), len(q))
    p, q = p + [Fraction(0)] * (n - len(p)), q + [Fraction(0)] * (n - len(q))
    return sum((abs(x - y) for x, y in zip(p, q)), Fraction(0)) / 2


def harmonic_optimum(X, prefix=None):
    """Exact Fraction costs; infinity only marks a forbidden root choice."""
    zero, one = [Fraction(0)] * (X + 1), [Fraction(0)] * (X + 1)
    for i in range(X, 1, -1):
        a, b = Fraction(0), Fraction(1, i)
        if i % 2 == 0:
            j = 3 * i // 2
            if j <= X:
                a += one[j]
                b += zero[j]
        else:
            j, k = (3 * i - 1) // 2, (3 * i + 1) // 2
            if j <= X:
                a += zero[j]
                b += min(zero[j], one[j])
            if k <= X:
                a += min(zero[k], one[k])
                b += one[k]
        if prefix is not None and i < len(prefix):
            if prefix[i] == 0:
                b = inf
            else:
                a = inf
        zero[i], one[i] = a, b
        if prefix is None:
            require(abs(b-a) <= Fraction(3, i-1), 'harmonic root-bit gap')
    answer = min(zero[2], one[2]) if X >= 2 else Fraction(0)
    require(isinstance(answer, Fraction), 'finite exact harmonic optimum')
    return answer


def phase_zero_bands(X):
    """Exact membership in [R^j,R^(j+1/2)), R=3/2, using integer squares."""
    bands = []
    j = 0
    while (3**j + 2**j - 1) // 2**j <= X:
        lo = (3**j + 2**j - 1) // 2**j
        hi = min(X, isqrt((3**(2*j+1) - 1) // 2**(2*j+1)))
        if lo <= hi:
            bands.append((lo, hi))
        j += 1
    return bands


def main():
    # Ten vertices are pair indices 2,...,11, with genuine P2 legal edges.
    X = 11
    eq, le = clauses(X)
    feasible = []
    for choice in product([0, 1], repeat=10):
        bits = [0, 0] + list(choice)
        if all(bits[i] + bits[j] == 1 for i, j in eq) and all(bits[i] <= bits[j] for i, j in le):
            feasible.append(bits)
    require(len(eq) + len(le) == 9 and len(feasible) == 16, 'ten-vertex carrier')
    objectives = [lambda b: sum(b) + sum(b[:5]),
                  lambda b: sum(b[:5]) + sum(b[:2]),
                  lambda b: sum(b) + 2 * sum(b[:5]) + sum(b[:2]),
                  lambda b: sum(b) + 2 * sum(b[:5]) + sum(b[:3])]
    costs = [min(f(b) for b in feasible) for f in objectives]
    require(costs == [3, 1, 4, 5], 'rounding hostile control')
    print(json.dumps({'ten_vertices': [2, 11], 'complement_edges': eq, 'order_edges': le,
                      'all_assignments': 1024, 'feasible_assignments': len(feasible),
                      'costs_separate1_separate2_nested_joint_fused_joint': costs}))

    a = [1, 0, 1]
    aa = convolution(a, a)
    # Monotone nested floors give the convolution inequality at every finite X.
    for X in range(1, 1001):
        rhs = sum(b * optimum(2**j * X // 3**j, a) for j, b in enumerate(a))
        require(optimum(X, aa) >= rhs, 'finite convolution inequality')
    exact = [486, 216, 96]
    require(exact[1] == 4 * exact[0] // 9 and exact[2] == 4 * exact[1] // 9,
            'exact intermediate cutoffs')
    separate = [optimum(486, a), optimum(216, a)]
    joint = optimum(486, aa)
    require(separate == [210, 94] and joint == 306, 'genuine convolution gain')
    print(json.dumps({'finite_convolution_cutoffs_checked': [1, 1000],
                      'rounding_free_cutoffs': exact, 'separate_optima': separate,
                      'common_assignment_optimum': joint, 'strict_gain': joint - sum(separate)}))

    # Exact finite-depth shift delay, including the hostile lost-depth case.
    base, _ = full_cost_certificates(a, 12)
    shifted, _ = full_cost_certificates([0, 0, 0] + a, 12)
    require(all(shifted[h] == (base[h - 3] if h >= 3 else 0) for h in range(13)),
            'finite-depth shift delay')
    blind, _ = full_cost_certificates([0] * 10 + [1], 9)
    unshifted, _ = full_cost_certificates([1], 9)
    require(blind[-1] == 0 < unshifted[-1], 'depth-truncation hostile control')
    b = [1, 1]
    fused, _ = full_cost_certificates(convolution(a, b), 12)
    q = [Fraction(w) * Fraction(2, 3)**j / normalizer(b) for j, w in enumerate(b)]
    require(all(fused[h] >= sum(q[j] * (base[h-j] if h >= j else 0) for j in range(len(q)))
                for h in range(13)), 'depth-corrected convolution inequality')
    print(json.dumps({'shift_delay_checked_depth': 12, 'shift': 3,
                      'hostile_truncation_depth': 9, 'hostile_shift': 10,
                      'shifted_partial': str(blind[-1]), 'unshifted_partial': str(unshifted[-1]),
                      'exact_shifted_and_unshifted_limits': 'both alpha by theorem',
                      'depth_corrected_convolution_checked': True}))

    # The Følner boundary estimate is exact and independent of pairing costs.
    p = [Fraction(9, 13), Fraction(0), Fraction(4, 13)]
    for m in range(1, 21):
        u = [Fraction(1, m)] * m
        bound = sum(p[k] * min(Fraction(k, m), 1) for k in range(len(p)))
        require(tv(convolution(p, u), u) <= bound, 'uniform-window boundary estimate')
    print(json.dumps({'uniform_window_TV_checks': 20, 'test_probability_kernel': [str(x) for x in p]}))

    # Inherited finite certificate now applies also to lower logarithmic density.
    kernel = [3**j * 2**(7-j) for j in range(8)]
    cert, totals = full_cost_certificates(kernel, 20)
    require(totals[-1] == 3286041555236, 'independent inherited aggregate')
    require(cert[-1] == Fraction(821510388809, 2677850419968), 'logarithmic-density certificate')
    print(json.dumps({'inherited_eight_scale_depth': 20, 'direct_full_cost_sum': totals[-1],
                      'lower_logarithmic_density_bound': str(cert[-1]),
                      'decimal': float(cert[-1]), 'global_realization_claimed': False}))

    # All sixteen legal ten-vertex prefixes, harmonic extension to cutoff 233.
    small_harmonic = min(sum((Fraction(b[i], i) for i in range(2, 12)), Fraction(0)) for b in feasible)
    require(harmonic_optimum(11) == small_harmonic, 'exhaustive harmonic optimum')
    H233 = harmonic_optimum(233)
    largest_excess = Fraction(0)
    for bits in feasible:
        fixed_cost = sum((Fraction(bits[i], i) for i in range(2, 12)), Fraction(0))
        constrained = harmonic_optimum(233, bits)
        require(constrained <= H233 + fixed_cost + 3, 'harmonic prefix extension bound')
        largest_excess = max(largest_excess, constrained - H233)
    print(json.dumps({'harmonic_gap_nodes_checked': 232, 'harmonic_prefix_extensions': len(feasible),
                      'prefix_cutoff': 11, 'extension_cutoff': 233,
                      'unrestricted_harmonic_cost_decimal': float(H233),
                      'largest_prefix_extension_excess_decimal': float(largest_excess),
                      'all_comparisons_exact_fractions': True}))

    # The asymptotic claims about this hostile family are proved in the note;
    # these finite integer-membership controls exhibit its endpoint profile.
    phase_rows = []
    for j in [10, 20, 30]:
        X = 3**j // 2**j
        bands = phase_zero_bands(X)
        count = sum(hi-lo+1 for lo, hi in bands)
        harmonic = fsum(1 / i for lo, hi in bands for i in range(lo, hi+1))
        require(abs(count / X - (sqrt(6)-2)) <= (2*j+10) / X,
                'phase hostile endpoint count')
        phase_rows.append({'X': X, 'count': count, 'prefix_density': count/X,
                           'harmonic_density': harmonic/log(X)})
    print(json.dumps({'phase_hostile_exact_kernel_limit': 'sqrt(6)-2',
                      'phase_hostile_log_density': '1/2',
                      'membership_method': 'integer squares, not floating logarithms',
                      'finite_controls': phase_rows}))
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
