"""Exact Bellman bounds and legal finite-residue policies for global P2.

Requires NumPy for chunked integer arrays. All decisions/certificates are exact;
decimal displays alone use floating point. Run normally or with python3 -O.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
import numpy as np


def require(test, message):
    if not test:
        raise AssertionError(message)


def solve(A, b):
    n = len(b)
    a = [list(map(F, row)) + [F(value)] for row, value in zip(A, b)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        require(pivot is not None, 'exact policy matrix nonsingular')
        a[j], a[pivot] = a[pivot], a[j]
        v = a[j][j]
        a[j] = [x/v for x in a[j]]
        for i in range(n):
            if i != j and a[i][j]:
                v = a[i][j]
                a[i] = [x-v*y for x, y in zip(a[i], a[j])]
    return [row[-1] for row in a]


def policy_system(t):
    M = len(t)
    inv = pow(3, -1, M) if M > 1 else 0
    P = [[(2*r+c)*inv % M for c in [0, 1, -1]] for r in range(M)]
    A = [[0]*M for _ in range(M)]
    for r, (p0, p1, p2) in enumerate(P):
        A[r][r] += 3
        A[r][p0] += 1
        A[r][p1 if t[r] else p2] -= 1
    x = solve(A, [1+b for b in t])
    require(all(0 < z < 1 for z in x), 'strict interior state means')
    return A, P, x


def policy_bits(t, N):
    M = len(t)
    bits = [0]*(N+1)
    for n in range(2, N+1):
        choice = t[n % M]
        if n % 3 == 0:
            bits[n] = 1-bits[2*n//3]
        elif n % 3 == 1:
            bits[n] = choice*bits[(2*n+1)//3]
        else:
            bits[n] = choice+(1-choice)*bits[(2*n-1)//3]
    return bits


def map_value(n, bits):
    i = (n+1)//2
    if n % 2:
        return i-1 if bits[i] else 3*i-1
    return 3*i if bits[i] else i


def audit_policies():
    best = None
    for t in product([0, 1], repeat=4):
        _, _, x = policy_system(t)
        value = sum(x)/4
        best = value if best is None else min(best, value)
    require(best == F(35, 108), 'all sixteen modulus-four policies')
    policy = [int(r in {2, 10, 18, 26, 30}) for r in range(32)]
    A, P, x = policy_system(policy)
    y = solve(list(map(list, zip(*A))), [1]*32)
    require(all(v <= 0 if policy[r] else v >= 0 for r, v in enumerate(y)), 'exact dual signs')
    for r, (p0, p1, p2) in enumerate(P):
        require(3*x[r]+x[p0]-x[p2] >= 1, 'LP lower row')
        require(3*x[r]+x[p0]-x[p1] <= 2, 'LP upper row')
    require(sum(x) == sum(y[r]*(1+policy[r]) for r in range(32)), 'primal-dual equality')
    density = sum(x)/32
    require(density == F(8209, 25920), 'modulus-32 density')
    bits = policy_bits(policy, 30000)
    for i in range(2, 20001):
        if i % 2 == 0:
            require(bits[i]+bits[3*i//2] == 1, 'global complement clause')
        else:
            require(bits[(3*i-1)//2] <= bits[i] <= bits[(3*i+1)//2], 'global order clauses')
    for n in range(3, 40001):
        a = map_value(n, bits)
        require(a < n or map_value(a, bits) < n, 'actual two-step first descent')
    hostile = [0]*20
    for i in range(1, 20):
        q = i
        while q % 3 == 0:
            hostile[i] ^= 1
            q //= 3
    require(map_value(6, hostile) == 9 and map_value(9, hostile) == 14, 'valuation-policy hostile')
    print(json.dumps({'modulus4_policies_exhausted': 16, 'modulus4_optimum': str(best),
                      'modulus32_one_residues': [i for i, b in enumerate(policy) if b],
                      'modulus32_exact_optimum': str(density), 'decimal': float(density),
                      'exact_rational_LP_dual_verified': True,
                      'actual_sources_checked': [3, 40000], 'hostile_v3_path': [6, 9, 14]}))


def fraction_bellman(values):
    M = len(values)
    r = F(2, 3)
    answer = []
    for i in range(2*M):
        if i % 2 == 0:
            answer.append(1-r*values[(3*i//2) % M])
        else:
            answer.append(1+r*min(values[((3*i-1)//2) % M], 0)
                          +r*max(values[((3*i+1)//2) % M], 0))
    return answer


@lru_cache(None)
def direct_tree(h, i):
    """Fixed-root costs on depth h-1, with each descendant level discounted2/3."""
    if h == 1:
        return F(0), F(1)
    r = F(2, 3)
    if i % 2 == 0:
        a, b = direct_tree(h-1, 3*i//2)
        return r*b, 1+r*a
    a, b = direct_tree(h-1, (3*i-1)//2)
    c, d = direct_tree(h-1, (3*i+1)//2)
    return r*(a+min(c, d)), 1+r*(min(a, b)+d)


def audit_contraction_hostile():
    d, e = list(map(F, [2, -1])), list(map(F, [1, -2]))
    difference = [x-y for x, y in zip(fraction_bellman(d), fraction_bellman(e))]
    require(difference == [F(-2, 3), F(4, 3), F(-2, 3), F(0)], 'sup-norm hostile')
    require(max(map(abs, difference)) == F(4, 3), 'not a sup-norm contraction')
    require(sum(map(abs, difference))/4 == F(2, 3), 'Haar L1 control')
    print(json.dumps({'sup_norm_contraction_refuted': True,
                      'exact_output_difference': [str(v) for v in difference],
                      'input_L1': '1', 'output_L1': '2/3'}))


def audit_phase_sweep(N=120000):
    """Integer-only address selection; finite checks do not prove the limit."""
    levels = [1]*(N+1)
    q, next_two, next_three = 0, 2, 3
    selector_powers = {}
    for n in range(2, N+1):
        while n*next_two >= next_three:
            q += 1
            next_two *= 2
            next_three *= 3
        if q not in selector_powers:
            h = 1
            while (h+1)**3 <= q:
                h += 1
            D = (h+1)**3-h**3
            A = D*q-h**3
            selector_powers[q] = h, D-1, 2**A, 3**A
        h, power, two, three = selector_powers[q]
        levels[n] = h+int(pow(n, power)*two <= three)

    policies = {}
    values = [F(1)]
    for h in range(1, max(levels)+1):
        policies[h] = [int(v < 0) for v in values]
        values = fraction_bellman(values)
    stationary = {h: policy_bits(policy, N) for h, policy in policies.items()}
    bits = [0]*(N+1)
    patched = [0]*(N+1)
    for n in range(2, N+1):
        h = levels[n]
        policy = policies[h]
        choice = policy[n % len(policy)]
        if n % 3 == 0:
            bits[n] = 1-bits[2*n//3]
        elif n % 3 == 1:
            bits[n] = choice*bits[(2*n+1)//3]
        else:
            bits[n] = choice+(1-choice)*bits[(2*n-1)//3]
        patched[n] = stationary[h][n]

    patch_failures, first_patch_failure = 0, None
    for i in range(2, 2*N//3+1):
        if i % 2 == 0:
            require(bits[i]+bits[3*i//2] == 1, 'sweep complement clause')
            patch_ok = patched[i]+patched[3*i//2] == 1
        else:
            require(bits[(3*i-1)//2] <= bits[i] <= bits[(3*i+1)//2], 'sweep order clause')
            patch_ok = patched[(3*i-1)//2] <= patched[i] <= patched[(3*i+1)//2]
        if not patch_ok:
            patch_failures += 1
            if first_patch_failure is None:
                first_patch_failure = i
    for n in range(3, 4*N//3+1):
        a = map_value(n, bits)
        require(a < n or map_value(a, bits) < n, 'sweep actual two-step descent')
    require(first_patch_failure == 26 and [patched[26], patched[39]] == [0, 0],
            'patching stationary bits need not preserve complement clauses')
    require([bits[26], bits[39]] == [0, 1], 'recursive sweep repairs the hostile patch')
    ancestor_cases = 0
    for n in range(2, N+1):
        a, two, three = n, 1, 1
        for j in range(1, 5):
            a = (2*a+(1 if a % 3 == 1 else -1 if a % 3 == 2 else 0))//3
            two *= 2
            three *= 3
            require(abs(three*a-two*n) <= three, 'exact bounded ancestor rounding')
            ancestor_cases += 1
    checkpoints = []
    count, patch_count, differences = 0, 0, 0
    for n in range(1, N+1):
        count += bits[n]
        patch_count += patched[n]
        differences += bits[n] != patched[n]
        if n in [1000, 10000, N]:
            checkpoints.append({'X': n, 'legal_count': count, 'reference_count': patch_count,
                                'differing_bits': differences})
    print(json.dumps({'exact_phase_selector_addresses': N,
                      'selected_levels': sorted(set(levels[2:])),
                      'actual_sources_checked': [3, 4*N//3],
                      'integer_ancestor_checks': ancestor_cases,
                      'patched_reference_clause_failures': patch_failures,
                      'first_patched_reference_bad_parent': first_patch_failure,
                      'hostile_patch_bits_26_39': [patched[26], patched[39]],
                      'legal_sweep_bits_26_39': [bits[26], bits[39]],
                      'finite_counts_not_asymptotic_evidence': checkpoints}))


def certified_bellman(depth=24):
    values = np.array([1], dtype=np.int64)  # d_1=1, period1.
    denominator = 1
    old_residual = None
    direct_cases = 0
    retained = []
    chunk = 65536
    for h in range(1, depth+1):
        M = len(values)
        new_den = 3*denominator
        # Every chunk sum, including the worst absolute residual, is below2^63.
        require(6*new_den*chunk < 2**63, 'integer chunk overflow guard')
        new = np.empty(2*M, dtype=np.int64)
        residual_sum, negative_sum, new_sum = 0, 0, 0
        if h <= 8:
            for i, numerator in enumerate(values):
                a, b = direct_tree(h, i)
                require(F(int(numerator), denominator) == b-a, 'independent fixed-root tree')
                direct_cases += 1
        for start in range(0, M, chunk):
            stop = min(M, start+chunk)
            k = np.arange(start, stop, dtype=np.int64)
            even = new_den-2*values[(3*k) % M]
            odd = new_den+2*np.minimum(values[(3*k+1) % M], 0)+2*np.maximum(values[(3*k+2) % M], 0)
            new[2*start:2*stop:2] = even
            new[2*start+1:2*stop:2] = odd
            require(bool(np.all(np.abs(even) <= 3*new_den)) and bool(np.all(np.abs(odd) <= 3*new_den)),
                    'uniform bounded-potential control')
            residual_sum += int(np.abs(even-3*values[(2*k) % M]).sum())
            residual_sum += int(np.abs(odd-3*values[(2*k+1) % M]).sum())
            negative_sum += int(np.minimum(values[start:stop], 0).sum())
            new_sum += int(even.sum())+int(odd.sum())
        require(new_sum == new_den*2*M, 'exact Haar mean one')
        residual = F(residual_sum, new_den*2*M)
        require(residual <= F(2, 3)**h, 'universal contraction envelope')
        if old_residual is not None:
            require(residual <= F(2, 3)*old_residual, 'observed exact L1 contraction')
        center = (1+F(negative_sum, denominator*M))/3
        lower, upper = center-residual/2, center+residual/2
        if h in [4, 8, 12, 16, 20, 22, 24]:
            row = {'iteration': h, 'policy_period': M, 'lower_exact': str(lower),
                   'upper_exact': str(upper), 'lower_decimal': float(lower),
                   'upper_decimal': float(upper), 'interval_width_exact': str(residual)}
            retained.append(row)
            print(json.dumps(row))
        values, denominator, old_residual = new, new_den, residual
    require(lower == F(750086585295922675, 2369190669160808448), 'retained lower certificate')
    require(upper == F(750158422907862637, 2369190669160808448), 'retained upper certificate')
    print(json.dumps({'independent_fixed_root_cases': direct_cases,
                      'all_residue_states_at_final_step': 2**24,
                      'signed_64_bit_chunk_arithmetic_checked': True,
                      'finite_policy_natural_density_infimum': 'B_star',
                      'upper_natural_density_infimum': 'B_star',
                      'natural_density_attainment_proved_by_slow_phase_sweep': True,
                      'Collatz_proof_claimed': False}))


def main():
    audit_policies()
    audit_contraction_hostile()
    audit_phase_sweep()
    certified_bellman()
    print('ALL CHECKS PASSED')


if __name__ == '__main__':
    main()
