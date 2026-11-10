"""Independent audit: no imports from producer, no inherited state recurrence."""
import itertools, json, hashlib, math

def check(condition):
    if not condition:
        raise RuntimeError('Independent coverage audit check failed')
from pathlib import Path

def adjacency(word):
    n = len(word) + 1
    g = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            forward = i != 0 or word[j - 1] == 'O'
            g[i if forward else j][j if forward else i] = True
    return g

def count(g):
    n = len(g)
    d = {(1 << v, v): 1 for v in range(n)}
    for m in range(1, 1 << n):
        for v in range(n):
            x = d.get((m, v), 0)
            if x:
                for w in range(n):
                    if not m >> w & 1 and g[v][w]:
                        key = (m | 1 << w, w)
                        d[key] = d.get(key, 0) + x
    return sum((d.get(((1 << n) - 1, v), 0) for v in range(n)))

def permutations(g):
    return sum((all((g[p[i]][p[i + 1]] for i in range(len(g) - 1))) for p in itertools.permutations(range(len(g)))))

def formula(r, s):
    return ((2 ** r - 1) * (2 ** s + 1) + 3) // 2

def independent_inverse(N):
    M = 2 * N - 3
    answers = []
    for r in range(3, max(3, M.bit_length() + 1)):
        q, rem = divmod(M, 2 ** r - 1)
        if rem == 0 and q >= 5 and (q - 1 & q - 2 == 0):
            answers.append((r, (q - 1).bit_length() - 1))
    return answers

def expected_bits(r, s):
    if r < s:
        return '1' * r + '0' * (s - r) + '1' * r
    if r == s:
        return '1' * (2 * r)
    return '1' + '0' * s + '1' * (r - s - 1) + '0' + '1' * s

def main():
    rows = []
    for a in range(1, 5):
        for b in range(1, 5):
            g = adjacency('O' * (a + 1) + 'Z' + 'O' * b + 'Z')
            H = count(g)
            check(H == formula(a + 2, b + 1))
            if len(g) <= 8:
                check(H == permutations(g))
            check(count([list(row) for row in zip(*g)]) == H)
            rows.append([a, b, len(g), H])
    check(count(adjacency('OOOO')) == 1)
    check(count(adjacency('OZ')) == 3)
    check(count(adjacency('OZOZ')) == 9 != formula(3, 2))
    for r in range(3, 123):
        for s in range(2, 122):
            check(bin(2 * formula(r, s) - 3)[2:] == expected_bits(r, s))
    gcd_checks = 0
    for r in range(3, 103):
        for s in range(2, 102):
            g = math.gcd(r, s)
            check(math.gcd(2 ** r - 1, 2 ** s + 1) == (2 ** g + 1 if r // g % 2 == 0 else 1))
            gcd_checks += 1
    bound = 100000
    expected = {formula(r, s): (r, s) for r in range(3, 20) for s in range(2, 20) if formula(r, s) <= bound}
    recovered = {}
    for N in range(1, bound + 1):
        inv = independent_inverse(N)
        check(len(inv) <= 1)
        if inv:
            recovered[N] = inv[0]
    check(recovered == expected)
    check(independent_inverse(7) == [] and independent_inverse(21) == [])
    counts = []
    for k in [10, 20, 40, 80, 160, 320]:
        c = sum((formula(r, s) <= 2 ** k for r in range(3, k + 2) for s in range(2, k + 2)))
        check(c == (k - 4) * (k - 3) // 2 + k // 2 - 2)
        counts.append([k, c])
    result = {'status': 'PASS', 'adjacency_DP_rows_a_b_n_H': rows, 'binary_patterns_checked': 14400, 'gcd_identity_checks': gcd_checks, 'inverse_targets_all_integers': bound, 'image_count': len(expected), 'exact_power_two_counts_k_count': counts, 'controls': {'transitive': 1, 'directed_triangle': 3, 'missing_seed_hostile': 9, 'correct_seed_minimum': 19, 'excluded_7_21': True}, 'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'scope': 'Two-burst image only; no all-tournament exclusion except inherited facts; not Chan-Pak coverage.'}
    Path(__file__).with_name('verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
