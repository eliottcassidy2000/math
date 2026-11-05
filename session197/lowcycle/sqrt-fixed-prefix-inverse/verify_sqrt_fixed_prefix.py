"""Exact fixed-prefix final-gadget inverse; no third-party packages."""
from functools import lru_cache
from pathlib import Path
import hashlib
import json

P = Path(__file__).resolve().parent

@lru_cache(None)
def g(n):
    if n < 1:
        raise ValueError('g requires a positive integer')
    if n == 1:
        return 1
    if n % 2 == 0:
        return g(n // 2)
    ell = (n + 1) // 2
    return g(ell) + ell - ell % 2


def scalar(prefix, L):
    h, e, R = prefix
    return h * g(L) + R * L - e * (L % 2)


def reverse_complement(b, k):
    q = 0
    for _ in range(k):
        q = 2 * q + (1 - (b & 1))
        b >>= 1
    return q


def inverse(prefix, target):
    h, e, R = prefix
    if not (h >= 1 and 0 <= e <= h and R >= h):
        raise ValueError('prefix inequalities are required')
    if target < 1:
        return [], 0
    K = target // R
    if K == 0:
        return [], 0
    if K == 1:
        return ([1] if scalar(prefix, 1) == target else []), 1
    k = ((K - 1).bit_length() + 1) // 2
    B = 1 << k
    table = [0] * (B + 1)
    table[1] = 1
    for n in range(2, B + 1):
        if n % 2 == 0:
            table[n] = table[n // 2]
        else:
            ell = (n + 1) // 2
            table[n] = table[ell] + ell - ell % 2
    answers = []
    evaluations = 0
    for b in range(B):
        q = reverse_complement(b, k)
        c = q % 2
        coefficient = h * q + R * B
        constant = h * (table[b + 1] - 1) + R * (b + 1) - e * ((b + 1) % 2)
        def value(x):
            return h * table[x + 1] + coefficient * x + h * c * (x % 2) + constant
        lo, hi = 0, (K - 1 - b) // B
        if hi < 0:
            continue
        while lo < hi:
            mid = (lo + hi) // 2
            evaluations += 1
            if value(mid) < target:
                lo = mid + 1
            else:
                hi = mid
        evaluations += 1
        if value(lo) == target:
            answers.append(B * lo + b + 1)
    return sorted(answers), evaluations


def gadget_graph(L):
    """Root zero; deletion transitive; literal integer gadget compiler."""
    m = (L - 1).bit_length()
    rows = [0] * (m + 1)
    for i in range(m):
        x = i + 1
        if (L - 1) >> i & 1:
            rows[x] |= 1
        else:
            rows[0] |= 1 << x
        for j in range(i + 1, m):
            rows[x] |= 1 << (j + 1)
    return rows


def seam(rows, root, L):
    a = list(rows)
    z = gadget_graph(L)
    n, m = len(a), len(z)
    for i in range(n):
        a[i] |= ((1 << m) - 1) << n
    a.extend(row << n for row in z)
    a[root] ^= 1 << n
    a[n] |= 1 << root
    return a, n


def profile(rows, root):
    n = len(rows)
    N = 1 << n
    start = [[0] * n for _ in range(N)]
    end = [[0] * n for _ in range(N)]
    count = [1] + [0] * (N - 1)
    for S in range(1, N):
        for v in range(n):
            if not S >> v & 1:
                continue
            U = S ^ (1 << v)
            if not U:
                start[S][v] = end[S][v] = 1
                continue
            for w in range(n):
                if U >> w & 1:
                    if rows[v] >> w & 1:
                        start[S][v] += start[U][w]
                    else:
                        end[S][v] += end[U][w]
        count[S] = sum(start[S])
    full = N - 1
    R = sum(start[S][root] * count[full ^ S] for S in range(1, N) if S >> root & 1)
    return count[full], end[full][root], R


def main():
    dependencies = {
        'ROOTED.md': 'ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e',
        'ARITHMETIC_INVERSE.md': '4ac1727eb27c85f8a39f5598eb9f0e1b020bc4020cbca0a94e6288cb0676c3b7',
    }
    for filename, digest in dependencies.items():
        if hashlib.sha256((P / filename).read_bytes()).hexdigest() != digest:
            raise RuntimeError(('dependency', filename))
    block_controls = monotone_controls = inverse_controls = physical_controls = 0
    for k in range(1, 9):
        B = 1 << k
        for b in range(B):
            q = reverse_complement(b, k)
            # A differently expressed bit-word check prevents reversed orientation drift.
            if q != int(''.join('1' if v == '0' else '0' for v in format(b, '0%db' % k))[::-1], 2):
                raise RuntimeError('reverse-complement orientation')
            for x in range(64):
                want = g(x + 1) + q * x + (q % 2) * (x % 2) + g(b + 1) - 1
                if g(B * x + b + 1) != want:
                    raise RuntimeError(('block', k, b, x))
                block_controls += 1
    for h in range(1, 10):
        for e in range(h + 1):
            for R in range(h, h + 5):
                prefix = h, e, R
                for target in range(0, 251):
                    got, _ = inverse(prefix, target)
                    want = [L for L in range(1, target // R + 1) if scalar(prefix, L) == target]
                    if got != want:
                        raise RuntimeError(('inverse', prefix, target, got, want))
                    inverse_controls += 1
    for h in range(1, 5):
        for e in range(h + 1):
            for R in range(h, h + 3):
                for K in range(2, 65):
                    B = 1 << (((K - 1).bit_length() + 1) // 2)
                    for b in range(B):
                        values = [scalar((h, e, R), L) for L in range(b + 1, K + 1, B)]
                        for a, z in zip(values, values[1:]):
                            if not a < z:
                                raise RuntimeError(('monotonicity', h, e, R, K, b))
                            monotone_controls += 1
    for L in range(1, 129):
        if profile(gadget_graph(L), 0)[0] != g(L):
            raise RuntimeError(('physical g', L))
        physical_controls += 1
    prefixes = [([2, 4, 1], 0, (3, 1, 3)), ([12, 5, 8, 2], 1, (5, 1, 7)),
                ([18, 20, 25, 3, 8], 0, (13, 2, 18))]
    for rows, root, expected in prefixes:
        if profile(rows, root) != expected:
            raise RuntimeError(('physical prefix', rows, root, profile(rows, root), expected))
        for L in range(1, 17):
            output, newroot = seam(rows, root, L)
            if profile(output, newroot)[0] != scalar(expected, L):
                raise RuntimeError(('physical seam', expected, L))
            physical_controls += 1
    large = []
    for (_, _, prefix), L in zip(prefixes, (100000007, 123456789, 987654321)):
        target = scalar(prefix, L)
        solutions, evaluations = inverse(prefix, target)
        if solutions != [L]:
            raise RuntimeError(('large exact inverse', prefix, target, solutions))
        large.append(dict(prefix=prefix, target=target, solutions=solutions, evaluations=evaluations))
    for bad in ((0, 0, 0), (3, -1, 3), (3, 4, 3), (3, 1, 2)):
        try:
            inverse(bad, 1)
        except ValueError:
            pass
        else:
            raise RuntimeError(('invalid prefix accepted', bad))
    result = dict(status='PASS', block_identities=block_controls, inverse_controls=inverse_controls,
                  strict_branch_comparisons=monotone_controls, physical_graph_controls=physical_controls,
                  prefix_graphs=[dict(rows=a, root=r, profile=t) for a, r, t in prefixes], large=large,
                  dependencies=dependencies)
    encoded = json.dumps(result, indent=2, sort_keys=True) + '\n'
    (P / 'sqrt_fixed_prefix_certificate.json').write_text(encoded)
    print(encoded, end='')

if __name__ == '__main__':
    main()
