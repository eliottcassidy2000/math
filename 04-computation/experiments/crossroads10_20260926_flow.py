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


def parent(i):
    if i % 3 == 0:
        return 2*i//3
    if i % 3 == 1:
        return (2*i+1)//3
    return (2*i-1)//3


def children(i):
    return [3*i//2] if i % 2 == 0 else [(3*i-1)//2, (3*i+1)//2]


def harmonic_forest(Y, X, removed=frozenset(), prefix=None, weights=None, C=1):
    """Exact forest optimum above Y, optionally retaining a prescribed prefix."""
    zero, one = {}, {}
    for i in range(X, Y, -1):
        a, b = Fraction(0), Fraction(1, i) if weights is None else weights[i]
        if i % 2 == 0:
            j = 3*i//2
            if j <= X and (i, j) not in removed:
                a += one[j]
                b += zero[j]
        else:
            j, k = (3*i-1)//2, (3*i+1)//2
            if j <= X and (i, j) not in removed:
                a += zero[j]
                b += min(zero[j], one[j])
            if k <= X and (i, k) not in removed:
                a += min(zero[k], one[k])
                b += one[k]
        zero[i], one[i] = a, b
        require(abs(b-a) <= Fraction(3*C, i-1), 'pruned harmonic gap')
    total = Fraction(0)
    for i in range(Y+1, X+1):
        p = parent(i)
        if p > Y and (p, i) not in removed:
            continue
        allowed = [0, 1]
        if p <= Y and prefix is not None:
            if p % 2 == 0:
                allowed = [1-prefix[p]]
            elif i == (3*p-1)//2 and prefix[p] == 0:
                allowed = [0]
            elif i == (3*p+1)//2 and prefix[p] == 1:
                allowed = [1]
        total += min((zero[i], one[i])[b] for b in allowed)
    return total, {i: one[i]-zero[i] for i in zero}


def exact_phase_colors(X, seeds):
    """Boundary seeds in [1,3/2) represent exact rational geometric boundaries."""
    R = Fraction(3, 2)
    require(seeds == sorted(set(seeds)) and seeds[0] == 1 and seeds[-1] < R, 'phase seeds')
    colors, scale = {}, Fraction(1)
    for i in range(1, X+1):
        while scale*R <= i:
            scale *= R
        z = Fraction(i)/scale
        colors[i] = max(q for q, seed in enumerate(seeds) if seed <= z)
    return colors


def audit_phase_cuts():
    rows = []
    for Y, X, Q in [(11, 233, 2), (11, 233, 10), (100, 2000, 10)]:
        seeds = [Fraction(1)+Fraction(q, 2*Q) for q in range(Q)]
        colors = exact_phase_colors(X, seeds)
        cuts = {(i, j) for i in range(Y+1, X+1) for j in children(i)
                if j <= X and colors[i] != colors[j]}
        full, _ = harmonic_forest(Y, X)
        decoupled, _ = harmonic_forest(Y, X, cuts)
        gap_budget = sum((Fraction(3, j-1) for _, j in cuts), Fraction(0))
        envelope = Fraction(18*Q, Y-1)
        require(0 <= full-decoupled <= gap_budget <= envelope, 'phase-cut uniform bound')
        restored = 0
        if X == 233 and Q == 10:
            remaining = set(cuts)
            old, old_gaps = harmonic_forest(Y, X, remaining)
            for edge in sorted(cuts, key=lambda edge: edge[1], reverse=True):
                child_gap = abs(old_gaps[edge[1]])
                remaining.remove(edge)
                new, new_gaps = harmonic_forest(Y, X, remaining)
                require(0 <= new-old <= child_gap, 'whole-component reconnection')
                old, old_gaps = new, new_gaps
                restored += 1
            require(old == full, 'restoration recovered full optimum')
        rows.append({'Y': Y, 'X': X, 'phase_cells': Q, 'cut_edges': len(cuts),
                     'cost_increase_decimal': float(full-decoupled),
                     'exact_root_gap_budget_decimal': float(gap_budget),
                     'uniform_envelope': str(envelope), 'edges_restored_individually': restored})
    # Freezing the descendants after flipping just a root can be illegal.
    cut_cost, gaps = harmonic_forest(2, 8, {(3, 5)}, {2: 0})
    full_cost, _ = harmonic_forest(2, 8, prefix={2: 0})
    require(cut_cost == Fraction(1, 2) and full_cost == Fraction(33, 40), 'propagation hostile costs')
    require(gaps[5] == Fraction(13, 40), 'entire child component must be reoptimized')
    frozen = {2: 0, 3: 1, 4: 0, 5: 1, 6: 1, 7: 0, 8: 0}
    require(frozen[5] > frozen[8], 'naive root-only flip must violate a descendant clause')
    print(json.dumps({'phase_cut_controls': rows,
                      'propagation_hostile': {'cut_edge': [3, 5], 'fixed_bit2': 0,
                                              'cut_cost': str(cut_cost), 'restored_cost': str(full_cost),
                                              'whole_component_gap': str(gaps[5]),
                                              'naive_root_flip_violates': 'epsilon_5<=epsilon_8'}}))


def audit_shell_comparison():
    R = Fraction(3, 2)
    splits = []
    for Y, X in [(2, 8), (11, 233), (100, 512)]:
        lower, upper = harmonic_optimum(Y), harmonic_optimum(X)
        shell, _ = harmonic_forest(Y, X)
        defect = upper-lower-shell
        require(0 <= defect <= 3, 'harmonic interval quasi-additivity')
        splits.append({'Y': Y, 'X': X, 'join_defect_decimal': float(defect)})
    shifts = []
    m, X = 3, Fraction(233)
    a, b = int(X//R**m), int(X)
    kernel = [R**k/Fraction(m) for k in range(m)]
    for ratio in [Fraction(1), Fraction(5, 4), Fraction(3, 2)]:
        Z = ratio*X
        c, d = int(Z//R**m), int(Z)
        weights = {}
        for i in range(a+1, d+1):
            power = Fraction(1)
            while power*R <= Z/i:
                power *= R
            weights[i] = R*power/((R-1)*Z)
            require(Fraction(2, i) <= weights[i] <= Fraction(3, i), 'exact phase price range')
        fixed, _ = harmonic_forest(a, b, weights=weights, C=3)
        moved, _ = harmonic_forest(c, d, weights=weights, C=3)
        require(abs(fixed-moved) <= 12, 'same-price shifted-shell comparison')
        f = cost_from_cutoffs(d, [(int(Z//R**k), w) for k, w in enumerate(kernel)])/Z
        require(moved <= m*f+2, 'normalized uniform objective to leading shell')
        require(fixed <= m*f+14, 'closure shell upper bound')
        shifts.append({'scale_ratio': str(ratio), 'fixed_shell': [a, b], 'moved_shell': [c, d],
                       'same_price_difference_decimal': float(fixed-moved),
                       'fixed_shell_minus_m_f_decimal': float(fixed-m*f)})
    print(json.dumps({'exact_shell_split_controls': splits, 'exact_shift_price_controls': shifts,
                      'all_comparisons_exact_fractions': True}))


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
                      'decimal': float(cert[-1]), 'natural_density_realization_claimed': False}))

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
    audit_phase_cuts()
    audit_shell_comparison()
    print(json.dumps({'audited_companion_theorem': 'B_star is the attained minimum existing logarithmic density',
                      'finite_harmonic_optimum_limit': 'H(X)/log(X) -> B_star',
                      'numerical_certificate_is_only_a_lower_bound': True,
                      'natural_density_attainment_claimed': False,
                      'collatz_proof_claimed': False}))
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
