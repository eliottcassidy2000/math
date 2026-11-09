# Standalone independent audit; all checks remain active under python -O.
# No producer imports.
from pathlib import Path
from itertools import permutations, product, combinations
from functools import lru_cache
from math import factorial
from collections import Counter

def decode(c, n):
    A = [[False] * n for _ in range(n)]
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j] = bool(c & 1 << b)
            A[j][i] = not A[i][j]
            b += 1
    return A

def sub(A, V):
    return [[A[i][j] for j in V] for i in V]

def count(A):
    n = len(A)

    @lru_cache(None)
    def end(S, v):
        if S == 1 << v:
            return 1
        R = S ^ 1 << v
        return sum((end(R, u) for u in range(n) if R >> u & 1 and A[u][v]))
    return sum((end((1 << n) - 1, v) for v in range(n)))

def mods(A):
    n = len(A)
    ans = []
    for r in range(2, n):
        for M in combinations(range(n), r):
            if all((len({A[v][u] for u in M}) == 1 for v in range(n) if v not in M)):
                ans.append(M)
    return ans

def prime(A):
    return len(A) >= 3 and (not mods(A))

def canon(A):
    cells = {}
    for v, row in enumerate(A):
        cells.setdefault(sum(row), []).append(v)
    candidates = product(*(permutations(cells[d]) for d in sorted(cells)))
    keys = Counter()
    for bs in candidates:
        p = sum(bs, ())
        key = tuple((A[p[i]][p[j]] for i in range(len(A)) for j in range(i + 1, len(A))))
        keys[key] += 1
    return (min(keys), next(iter(keys.values())))

def W(k):
    n = 2 * k + 1
    return [[i < j if i < n - 1 and j < n - 1 else j % 2 == 0 if i == n - 1 and j != i else i % 2 == 1 if j == n - 1 and i != j else False for j in range(n)] for i in range(n)]
import argparse
import os

def source_directory():
    parser = argparse.ArgumentParser(description='Independently certify tournament census and prime fan kernels.')
    parser.add_argument('source', nargs='?', default=os.environ.get('PRIME_AUDIT_SRC'), help='Directory containing classes_1.txt through classes_8.txt (default: parent of this script directory; PRIME_AUDIT_SRC also supported)')
    args = parser.parse_args()
    return Path(args.source).expanduser().resolve() if args.source else Path(__file__).resolve().parents[1]

def require(condition, message):
    if not condition:
        raise ValueError(message)

def main():
    SRC = source_directory()
    old = set()
    for n in range(1, 9):
        rows = [tuple(map(int, l.split()[:2])) for l in (SRC / f'classes_{n}.txt').read_text().splitlines() if l and (not l.startswith('#'))]
        keys = set()
        mass = 0
        ps = nc = crit = 0
        hs = []
        for c, h in rows:
            A = decode(c, n)
            require(count(A) == h, 'Check failed: count(A) == h')
            key, aut = canon(A)
            require(key not in keys, 'Check failed: key not in keys')
            keys.add(key)
            require(factorial(n) % aut == 0, 'Check failed: factorial(n) % aut == 0')
            mass += factorial(n) // aut
            if prime(A):
                ps += 1
                if h in old:
                    nc += 1
                    hs.append(h)
                    if not any((prime(sub(A, [u for u in range(n) if u != v])) for v in range(n))):
                        crit += 1
            if n == 8 and h == 63:
                mm = mods(A)
                require(len(mm) == 3 and all((len(m) == 2 for m in mm)), 'Check failed: len(mm) == 3 and all((len(m) == 2 for m in mm))')
                blocks = mm + [(v,) for v in range(n) if not any((v in m for m in mm))]
                Q = sub(A, [m[0] for m in blocks])
                require(canon(Q)[0] == canon(W(2))[0], 'Check failed: canon(Q)[0] == canon(W(2))[0]')
                ds = sorted((count(sub(A, [u for u in range(n) if u != v])) for v in range(n)))
                require(ds == [1, 11, 25, 25, 33, 33, 37, 37], 'Check failed: ds == [1, 11, 25, 25, 33, 33, 37, 37]')
                print('ORPHAN', c, 'modules', mm, 'deletions', ds, flush=True)
        require(mass == 2 ** (n * (n - 1) // 2), 'Check failed: mass == 2 ** (n * (n - 1) // 2)')
        print('CENSUS', n, len(rows), 'orbit_mass', mass, 'primes', ps, 'noncreator', nc, 'critical_noncreator', crit, 'noncreator_H', sorted(hs) if n <= 7 else Counter(hs), flush=True)
        old.update((h for c, h in rows))
    require(count(decode(6241, 6)) == 31, 'Check failed: count(decode(6241, 6)) == 31')
    require([sum((int(v) << j for j, v in enumerate(row))) for row in decode(6241, 6)] == [2, 12, 33, 21, 7, 27], 'Check failed: [sum((int(v) << j for j, v in enumerate(row))) for row in decode(6241, 6)] == [2, 12, 33, 21, 7, 27]')
    for k in range(1, 7):
        A = W(k)
        h = count(A)
        require(h == (4 ** (k + 1) + 6 * k + 5) // 9 and prime(A), 'Check failed: h == (4 ** (k + 1) + 6 * k + 5) // 9 and prime(A)')
        require(not any((prime(sub(A, [u for u in range(len(A)) if u != v])) for v in range(len(A)))), 'Check failed: not any((prime(sub(A, [u for u in range(len(A)) if u != v])) for v in range(len(A))))')
        print('W', 2 * k + 1, h, flush=True)
    for k in range(4, 15):
        r = k // 2
        n = k + 2
        F = [[False] * n for _ in range(n)]
        for i in range(n - 1):
            for j in range(i + 1, n - 1):
                F[i][j] = True
            F[i][-1] = i + 1 in (r, k + 1 - r, k + 1)
            F[-1][i] = not F[i][-1]
        sizes = [r - 1, 2, r - 1, 1, 1] if k % 2 == 0 else [r - 1, 1, 1, 1, r - 1, 1, 1]
        Q = W(2 if k % 2 == 0 else 3)
        labels = [b for b, s in enumerate(sizes) for _ in range(s)]
        B = [[i < j if labels[i] == labels[j] else Q[labels[i]][labels[j]] for j in range(n)] for i in range(n)]
        require(B == F, 'Check failed: B == F')
        require(count(F) == 2 ** k - 1, 'Check failed: count(F) == 2 ** k - 1')
        print('FAN', k, 'sizes', sizes, 'H', 2 ** k - 1, flush=True)
    print('ALL CHECKS PASS', flush=True)
if __name__ == '__main__':
    main()
