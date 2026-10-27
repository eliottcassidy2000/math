#!/usr/bin/env python3
"""Session 73: exact arithmetic carrier / gap firewall tests, stdlib only."""
import argparse
import itertools
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    answer = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            answer[p] = answer.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        answer[n] = answer.get(n, 0) + 1
    return answer


def sigma(n):
    return math.prod((p ** (e + 1) - 1) // (p - 1) for p, e in factor(n).items())


def aliquot(n):
    return sigma(n) - n


def prime(n):
    return n >= 2 and factor(n) == {n: 1}


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def squarefull(n):
    return all(e >= 2 for e in factor(n).values())


def combine(x, y):
    a, u = x
    b, v = y
    return a * b, a * v + b * u + u * v


def decode(n, mask):
    arcs = [[False] * n for _ in range(n)]
    for k, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        arcs[i][j] = bool((mask >> k) & 1)
        arcs[j][i] = not arcs[i][j]
    return arcs


def h_dp(arcs):
    n = len(arcs)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if dp[mask][last]:
                for nxt in range(n):
                    if not (mask >> nxt) & 1 and arcs[last][nxt]:
                        dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[-1])


def h_permutations(arcs):
    n = len(arcs)
    return sum(all(arcs[p[i]][p[i + 1]] for i in range(n - 1))
               for p in itertools.permutations(range(n)))


def join(a, b):
    n, m = len(a), len(b)
    c = [[False] * (n + m) for _ in range(n + m)]
    for i in range(n + m):
        for j in range(n + m):
            c[i][j] = (a[i][j] if i < n and j < n else
                       b[i-n][j-n] if i >= n and j >= n else i < n)
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--negative-control', choices=['target-only', 'drop-guard', 'omit-cube', 'count-is-weight'])
    args = ap.parse_args()
    if args.negative_control == 'target-only':
        require(aliquot(42) == aliquot(175), 'TARGET_ONLY_CONTEXT_FAILED: 54 != 73')
    if args.negative_control == 'drop-guard':
        require(aliquot(4) == 2 * aliquot(2) + sigma(2), 'PRIME_COPRIMALITY_GUARD_FAILED: 3 != 5')
    if args.negative_control == 'omit-cube':
        require(125 ** 3 <= aliquot(125) ** 4, 'PRIME_CUBE_EXCEPTION_REQUIRED: 125^3 > 31^4')
    if args.negative_control == 'count-is-weight':
        require(sum([1, 101]) <= len([1, 101]), 'EXCEPTION_COUNT_IS_NOT_WEIGHT: 102 > 2')

    # Independent sieve gives COMPLETE fibers for every target in 2..300.
    target_bound = 300
    source_bound = (target_bound - 1) ** 2
    sums = [0] * (source_bound + 1)
    for d in range(1, source_bound // 2 + 1):
        for n in range(2 * d, source_bound + 1, d):
            sums[n] += d
    fibers = defaultdict(list)
    for n in range(1, source_bound + 1):
        if 2 <= sums[n] <= target_bound:
            fibers[sums[n]].append(n)
    for n in range(1, source_bound + 1):
        if n <= 1000 or squarefull(n):
            require(sums[n] == aliquot(n), f'independent sigma disagreement at {n}')
    require(fibers[2] == fibers[5] == fibers[52] == [], 'empty-fiber control failed')
    require(fibers[6] == [6, 25], 'complete F(6) failed')
    require(fibers[7] == [8] and fibers[21] == [18, 51, 91], 'touchability firewall failed')
    require(aliquot(9) == 4 and aliquot(27) == 13, 'multiplication nonclosure witnesses failed')
    multiplication_failures = [(a, b, a*b) for a in range(2, 53) for b in range(a, 53)
                               if a*b <= 52 and fibers[a] and fibers[b] and not fibers[a*b]]
    require(multiplication_failures == [(4, 13, 52)], 'first nonclosure pair changed')

    carrier_count = 0
    for a in range(1, 151):
        for b in range(1, 151):
            if math.gcd(a, b) == 1:
                require(combine((a, aliquot(a)), (b, aliquot(b))) == (a*b, aliquot(a*b)), 'coprime carrier failed')
                carrier_count += 1
    observability_count = 0
    primes = [p for p in range(2, 102) if prime(p)]
    for a in range(1, 501):
        for p in primes:
            if a % p:
                u = aliquot(a)
                v = aliquot(a*p)
                require(v - (p+1)*u == a, 'prime-context observability failed')
                observability_count += 1
    collision_count = 0
    for m in range(2, target_bound + 1):
        for a, b in itertools.combinations(fibers[m], 2):
            p = next(p for p in primes if (a*b) % p)
            require(aliquot(a*p)-aliquot(b*p) == a-b, 'equal-target separation failed')
            collision_count += 1

    finite_bank_failures = []
    for C in [3, 5, 10]:
        failures = []
        for m in range(2, target_bound + 1):
            if not prime(m) or not fibers[m]:
                continue
            bank = False
            for a in range(2, C + 1):
                S, T = aliquot(a), sigma(a)
                if m > T and (m-T) % S == 0:
                    p = (m-T)//S
                    if prime(p) and a % p:
                        bank = True
            sf = any(squarefull(n) for n in fibers[m])
            if not bank and not sf:
                failures.append({'prime_target': m, 'complete_fiber': fibers[m]})
        require(failures, 'finite bank hostile unexpectedly absent')
        finite_bank_failures.append({'core_bank': [2,C], 'first_missed_but_touchable': failures[0]})

    # Shape dichotomy supporting the analytic O(X^(2/3)) statement.
    shape_counts = Counter()
    squarefull_prime_targets = set()
    squarefull_odd_targets = set()
    for n in range(2, source_bound + 1):
        f = factor(n)
        if not all(e >= 2 for e in f.values()):
            continue
        m = sums[n]
        omega = sum(f.values())
        if omega == 2:
            require(len(f) == 1 and m == next(iter(f)) + 1, 'prime square branch failed')
            require(not prime(m) or m == 3, 'unexpected prime-square prime image')
            shape_counts['prime_square'] += 1
        elif omega == 3:
            require(len(f) == 1 and next(iter(f)) ** 2 <= m, 'prime cube branch failed')
            shape_counts['prime_cube'] += 1
        else:
            p = min(f)
            require(p**4 <= n and n//p <= m and n**3 <= m**4, 'squarefull bulk source bound failed')
            shape_counts['omega_at_least_four'] += 1
        if m <= 1000 and prime(m):
            squarefull_prime_targets.add(m)
        if m <= 1000 and m % 2:
            squarefull_odd_targets.add(m)

    # Exact CRT density factor for fixed affine prime branches.
    # This checks LOCAL residues; no finite census is called an asymptotic proof.
    local_checks = 0
    crt_rows = []
    for a in range(2, 101):
        S, T = aliquot(a), sigma(a)
        for ell in primes:
            bad = {r for r in range(ell) if r % ell == 0 or (S*r+T) % ell == 0}
            expected = ell if S % ell == T % ell == 0 else 1 if S % ell == 0 or T % ell == 0 else 2
            require(len(bad) == expected, 'two-affine-form local exclusion count failed')
            local_checks += 1
        bank = [ell for ell in [3, 5, 7, 11] if (S*T) % ell]
        modulus = math.prod(bank)
        good_prime_residues = sum(math.gcd(r, modulus) == 1 and all((r-T) % ell for ell in bank)
                                  for r in range(modulus))
        predicted = math.prod(ell-2 for ell in bank)
        require(good_prime_residues == predicted, 'CRT prime-target survivor count failed')
        denominator = math.prod(ell-1 for ell in bank)
        if a in [3, 5, 9, 15, 21, 25]:
            crt_rows.append({'core': a, 'S': S, 'T': T, 'auxiliary_primes': bank,
                             'survivor_fraction': str(Fraction(predicted, denominator))})

    # Positive omission transport, retaining top/source n and source identity.
    omission_rows = []
    for a in [4, 6, 9, 25, 27]:
        proper = divisors(a)[:-1]
        omitted = proper[::2]
        contexts = [p for p in [2, 3, 5, 7, 11, 13] if a % p]
        for k in range(1, min(4, len(contexts)) + 1):
            P = math.prod(contexts[:k])
            lifted = [d*e for d in omitted for e in divisors(P)]
            require(len(set(lifted)) == len(omitted)*2**k, 'omission multiplicity/guard failed')
            require(sum(lifted) == sum(omitted)*sigma(P), 'weighted omission transport failed')
            exact = combine((a, aliquot(a)), (P, aliquot(P)))[1]
            truncated = combine((a, aliquot(a)-sum(omitted)), (P, aliquot(P)))[1]
            require(exact-truncated == sum(lifted), 'positive carrier omission error failed')
            omission_rows.append({'source': a, 'context': P, 'omitted_count': len(lifted),
                                  'omitted_mass': sum(lifted)})

    # Independent dynamic and permutation enumeration of all 64 labelled T4.
    h4 = Counter()
    witnesses = {}
    for mask in range(64):
        A = decode(4, mask)
        h = h_dp(A)
        require(h == h_permutations(A), 'Hamilton path independent count disagreement')
        h4[h] += 1
        witnesses.setdefault(h, mask)
    require(h4[5] == 24, 'H=5 firewall failed')
    selected = [decode(4, witnesses[h]) for h in sorted(witnesses)]
    for a in selected:
        for b in selected:
            require(h_dp(join(a,b)) == h_dp(a)*h_dp(b), 'ordinal join H product failed')
    # Find a fixed exact H=63 witness independently (not spectrum completeness).
    import random
    rng = random.Random(730021)
    h63_mask = None
    for _ in range(10000):
        mask = rng.randrange(1 << 28)
        if h_dp(decode(8, mask)) == 63:
            h63_mask = mask
            break
    require(h63_mask is not None, 'H=63 witness search exhausted its explicit finite universe')
    require(h_permutations(decode(8,h63_mask)) == 63, 'H63 permutation certificate failed')

    print(json.dumps({
        'status': 'FINITE-EXACT controls; asymptotic claims proved in REPORT.md with cited classical input',
        'complete_fibers': {'targets': [2,target_bound], 'source_bound': source_bound,
                            'selected': {str(m): fibers[m] for m in [2,3,4,5,6,7,13,21,52]}},
        'least_multiplicative_touchability_failure': multiplication_failures,
        'coprime_carrier_checks': carrier_count,
        'prime_context_observability_checks': observability_count,
        'all_complete_fiber_collision_checks': collision_count,
        'finite_bank_plus_all_squarefull_hostiles': finite_bank_failures,
        'squarefull_source_shapes_through': source_bound,
        'squarefull_shape_counts': dict(shape_counts),
        'squarefull_prime_targets_le_1000_in_source_scan': sorted(squarefull_prime_targets),
        'squarefull_odd_target_count_le_1000_in_source_scan': len(squarefull_odd_targets),
        'affine_local_residue_checks': local_checks,
        'CRT_prime_target_survivor_fractions': crt_rows,
        'omission_transport_cases': len(omission_rows),
        'omission_selected': omission_rows[-4:],
        'H4_histogram': dict(sorted(h4.items())),
        'H5_mask': witnesses[5],
        'H63_order8_mask': h63_mask,
        'mask_convention': 'unordered pairs itertools.combinations(range(n),2), low bit first; 1 means i->j',
        'aliquot_prime_uniqueness': 'NOT PROVED',
        'H_full_spectrum_completeness': 'NOT PROVED',
        'native_ROOT_positivity_or_budget': 'NOT CLAIMED',
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
