"""Small independent controls for the inherited odd-cycle counting identity.

No imports from the repository research programs. All checks use explicit
exceptions, so python -O runs the same verification. This is not a proof of
global Hamiltonian spectrum completeness or an explicit path bijection.
"""
from itertools import combinations, permutations
import json
import random


def check(ok, label):
    if not ok:
        raise RuntimeError(label)


def tournament(n, mask):
    rows = [0] * n
    for b, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if (mask >> b) & 1 else (j, i)
        rows[u] |= 1 << v
    return rows


def hamilton(rows):
    n = len(rows)
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            todo = rows[v] & ~mask
            while todo:
                b = todo & -todo
                todo -= b
                w = b.bit_length() - 1
                dp[mask | b][w] += dp[mask][v]
    return sum(dp[-1])


def cycles(rows):
    n = len(rows)
    result = []
    for length in range(3, n + 1, 2):
        for support in combinations(range(n), length):
            a, rest = support[0], support[1:]
            for tail in permutations(rest):
                word = (a,) + tail
                if all((rows[word[i]] >> word[(i + 1) % length]) & 1
                       for i in range(length)):
                    result.append((sum(1 << x for x in support), word))
    return result


def packing_profile(odd_cycles):
    counts = [1]

    def visit(start, used, size):
        for j in range(start, len(odd_cycles)):
            support = odd_cycles[j][0]
            if support & used:
                continue
            if len(counts) <= size + 1:
                counts.append(0)
            counts[size + 1] += 1
            visit(j + 1, used | support, size + 1)
    visit(0, 0, 0)
    return counts


def verify(rows):
    h = hamilton(rows)
    profile = packing_profile(cycles(rows))
    orbit_count = sum(a * (1 << (k - 1))
                      for k, a in enumerate(profile) if k)
    check(h == 1 + 2 * orbit_count, (rows, h, profile))
    return h, profile


def one_apex(bits):
    m = len(bits)
    rows = [sum(1 << j for j in range(i + 1, m)) for i in range(m)] + [0]
    for i, b in enumerate(bits):
        if b:
            rows[m] |= 1 << i
        else:
            rows[i] |= 1 << m
    return rows


def signature_index(bits):
    return sum(1 if j - i <= 2 else 1 << (j - i - 2)
               for i in range(len(bits)) for j in range(i + 1, len(bits))
               if bits[i] and not bits[j])


def main():
    spectra = {}
    tested = 0
    for n in range(1, 6):
        values = set()
        for mask in range(1 << (n * (n - 1) // 2)):
            h, _ = verify(tournament(n, mask))
            values.add(h)
            tested += 1
        spectra[n] = sorted(values)
    rng = random.Random(210)
    sampled = {}
    for n in (6, 7, 8):
        for _ in range(32):
            verify(tournament(n, rng.getrandbits(n * (n - 1) // 2)))
        sampled[n] = 32

    c3 = [2, 4, 1]
    joined = [x | (7 << 3) for x in c3] + [x << 3 for x in c3]
    h, profile = verify(joined)
    check((h, profile) == (9, [1, 2, 1]), 'C3 order-join control')

    words = 0
    for m in range(8):
        for code in range(1 << m):
            bits = [(code >> (m - 1 - i)) & 1 for i in range(m)]
            h = hamilton(one_apex(bits))
            check(h == 1 + 2 * signature_index(bits), ('word', bits, h))
            words += 1
    cap = 128
    support = {0}
    # The inherited endpoint-pair bound gives L <= 3 + floor(log2 cap).
    for m in range(2, 4 + cap.bit_length() - 1):
        for middle in range(1 << (m - 2)):
            bits = [1] + [(middle >> (m - 3 - i)) & 1
                          for i in range(m - 2)] + [0]
            index = signature_index(bits)
            if index <= cap:
                support.add(index)
    check(3 not in support and 10 not in support, 'known single-core gaps')
    check(31 in support and 6 not in support, 'known positive and hostile')
    print(json.dumps({
        'status': 'PASS',
        'all_labelled_orders': [1, 2, 3, 4, 5],
        'all_labelled_count': tested,
        'exact_small_H_spectra': spectra,
        'deterministic_sample_counts': sampled,
        'join_control': {'rows': joined, 'H': hamilton(joined),
                         'packing_profile': profile, 'index': 4},
        'one_apex_words_directly_checked': words,
        'single_core_index_cap': cap,
        'single_core_missing_indices': sorted(set(range(cap + 1)) - support),
        'scope': 'Inherited identity controls and elementary orbit-count corollary; not completeness.'
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
