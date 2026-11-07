#!/usr/bin/env python3
"""Independent exact controls; Python3 standard library, no assertions."""
import hashlib
import itertools
import json
from formula import comb, coefficients, path_formula, require, spectrum

def ending_table(out):
    n = len(out)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    # Compute by final-edge deletion, rather than the constructor's formula.
    incoming = [sum(1 << j for j in range(n) if out[j] >> i & 1) for i in range(n)]
    for mask in range(1, 1 << n):
        if mask & (mask - 1) == 0:
            continue
        for last in range(n):
            if mask >> last & 1:
                prior = mask ^ (1 << last)
                choices = incoming[last] & prior
                total = 0
                while choices:
                    bit = choices & -choices
                    choices ^= bit
                    total += dp[prior][bit.bit_length() - 1]
                dp[mask][last] = total
    return dp

def path_count(out):
    return sum(ending_table(out)[-1])

def reverse(out):
    return [sum(1 << j for j in range(len(out)) if out[j] >> i & 1) for i in range(len(out))]

def cover_data(out):
    n = len(out)
    full = (1 << n) - 1
    end, start = ending_table(out), ending_table(reverse(out))
    K = [[0] * n for _ in range(n)]
    for mask in range(1, full):
        other = full ^ mask
        for i in range(n):
            if end[mask][i]:
                for j in range(n):
                    K[i][j] += end[mask][i] * start[other][j]
    return sum(end[-1]), start[-1], end[-1], K

def attach(out, x):
    n = len(out)
    return [row | (((x >> i) & 1) << n) for i, row in enumerate(out)] + [((1 << n) - 1) ^ x]

def attachment_formula(data, x):
    h, s, e, K = data
    n = len(s)
    return sum(s[i] * (1 - (x >> i & 1)) + e[i] * (x >> i & 1) for i in range(n)) + sum(K[i][j] * (x >> i & 1) * (1 - (x >> j & 1)) for i in range(n) for j in range(n))

def brute_paths(out):
    return sum(all(out[p[i]] >> p[i + 1] & 1 for i in range(len(out) - 1)) for p in itertools.permutations(range(len(out))))

def cycle_count(out):
    n = len(out)
    # Root at vertex zero, enumerate partial paths by a memoized suffix count.
    from functools import lru_cache
    @lru_cache(None)
    def f(last, remaining):
        if not remaining:
            return int(bool(out[last] & 1))
        choices = out[last] & remaining
        total = 0
        while choices:
            bit = choices & -choices
            choices ^= bit
            total += f(bit.bit_length() - 1, remaining ^ bit)
        return total
    return f(0, ((1 << n) - 1) ^ 1) if n >= 3 else 0

def tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        out = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
        yield out

def main():
    extension_cases, core_cases, brute_cases = 0, 0, 0
    for n in range(1, 5):
        for out in tournaments(n):
            data = cover_data(out)
            core_cases += 1
            for x in range(1 << n):
                actual = attach(out, x)
                predicted = attachment_formula(data, x)
                require(path_count(actual) == predicted, ('attachment DP', out, x))
                require(brute_paths(actual) == predicted, ('attachment permutation', out, x))
                require(predicted >= data[0], 'induced monotonicity')
                extension_cases += 1
                brute_cases += 1
    family_counts = []
    coefficient_entries = 0
    for m in range(1, 9):
        U = comb(m, 0, 0, True)[:-1]
        U = [row & ((1 << len(U)) - 1) for row in U]
        h, starts, ends, K = cover_data(U)
        require(h == (1 << m) + 1, ('base U', m))
        A = [ends[i] - starts[i] + sum(K[i]) for i in range(len(U))]
        a = m + 2
        require(A[a] == 4 * (1 << m) - 2, ('A_a', m))
        require(A[1] == 2 * (1 << m) + 2, ('A_1', m))
        require(K[a][1] == 1 << m and K[1][a] == 0, ('K_a1', m))
        for i in range(m):
            p, first, last = i + 2, int(i == 0), int(i == m - 1)
            d = (1 << (m + 1 - i)) * 3 ** i + (1 << (i + 2)) - 2 * last
            require(A[p] == d, ('A_p', m, i))
            expected = [(p, a, (1 << i) + last), (a, p, (1 << (m - 1 - i)) - first), (p, 1, (1 << i) - last), (1, p, (1 << (m - 1 - i)) + first)]
            for u, v, count in expected:
                require(K[u][v] == count, ('special K', m, u, v))
                coefficient_entries += 1
            for j in range(i + 1, m):
                r, distance = j + 2, j - i
                require(K[p][r] == (1 << (m - distance)) + int(distance == 1), ('forward K', m, i, j))
                require(K[r][p] == (1 << (m - distance)) * 3 ** (distance - 1) + (1 << (distance - 1)), ('backward K', m, i, j))
                coefficient_entries += 2
        checked = 0
        for q in range(1 << m):
            for tau in (0, 1):
                for rev in (False, True):
                    out = comb(m, q, tau, rev)
                    predicted = path_formula(m, q, tau, rev)
                    require(path_count(out) == predicted, ('comb DP', m, q, tau, rev))
                    require(cycle_count(out) == q + tau + int(not rev), ('cycle provenance control', m, q, tau, rev))
                    if m <= 3:
                        require(brute_paths(out) == predicted, ('comb permutation', m, q, tau, rev))
                        brute_cases += 1
                    checked += 1
        family_counts.append({'m': m, 'order': m + 4, 'cases': checked})
    limit = 65535
    found, scanned, max_m = spectrum(limit)
    for h, params in found.items():
        require(path_formula(*params) == h, ('Gray-code replay', h, params))
    support = sorted(found)
    low = [h for h in support if h <= 35]
    require(low == [3, 5, 9, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33], 'small exact support')
    require(11 not in found and 35 not in found, 'hostile missing family target')
    witnesses = {}
    for h in (23, 25, 27, 29, 31, 33, 37, 65535):
        if h in found:
            params = found[h]
            out = comb(*params)
            require(path_count(out) == h, ('displayed witness', h))
            witnesses[str(h)] = {'parameters_m_q_tau_reversed': list(params), 'outrows': out, 'H': h, 'C': cycle_count(out)}
    no_go = []
    for n in range(3, 13):
        out = [sum(1 << j for j in range(i + 1, n)) for i in range(n)]
        out[0] ^= 1 << (n - 1)
        out[-1] |= 1
        actual_h, actual_c = path_count(out), cycle_count(out)
        require(actual_h == (1 << (n - 2)) + 1 and actual_c == 1, ('no-go family', n))
        no_go.append({'n': n, 'H': actual_h, 'C': actual_c})
    report = {
      'status': 'FINITE-EXACT controls of all-order proofs; Chan-Pak completeness remains OPEN',
      'general_attachment': {'labeled_cores_orders_1_to_4': core_cases, 'all_attachments': extension_cases, 'all_attachments_permutation_checked': extension_cases},
      'comb_formula': {'cases': sum(x['cases'] for x in family_counts), 'by_m': family_counts, 'cover_coefficient_entries': coefficient_entries, 'total_permutation_controls': brute_cases},
      'all_length_comb_support': {'limit': limit, 'max_m_proved_sufficient': max_m, 'parameter_cases': scanned, 'attained_values': len(support), 'support_through_35': low, 'first_100_missing_odds_at_least_23': [h for h in range(23, limit + 1, 2) if h not in found][:100], 'sorted_support': support, 'sorted_support_sha256': hashlib.sha256(json.dumps(support, separators=(',', ':')).encode()).hexdigest()},
      'witnesses': witnesses,
      'nonmonotone_control': [{'parameters': [3, 6, 0, False], 'C': 7, 'H': path_formula(3, 6, 0, False)}, {'parameters': [3, 7, 0, False], 'C': 8, 'H': path_formula(3, 7, 0, False)}],
      'constant_cycle_unbounded_path_controls': no_go,
      'scope': 'The comb is one restricted construction. Its missing values are not missing ambient tournament path counts.'
    }
    print(json.dumps(report, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
