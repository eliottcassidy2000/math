"""Stdlib exact native two-run observer and rational kernel; no insertion assumptions."""
def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


from fractions import Fraction
from itertools import permutations, combinations

def hp_table(A):
    n = len(A)
    hp = {}
    for v in range(n):
        h = [[0] * n for _ in range(n)]
        h[v][v] = 1
        hp[1 << v] = h
    for size in range(2, n + 1):
        for c in combinations(range(n), size):
            mask = sum((1 << v for v in c))
            h = [[0] * n for _ in range(n)]
            for t in c:
                prev = hp[mask ^ 1 << t]
                for s in c:
                    h[s][t] = sum((prev[s][v] * A[v][t] for v in c if v != t))
            hp[mask] = h
    return hp

def observer(A):
    n = len(A)
    hp = hp_table(A)
    full = (1 << n) - 1
    ports = list(permutations(range(n), 2))
    idx = {x: i for i, x in enumerate(ports)}
    M = [[0] * len(ports) for _ in ports]
    for mask in range(1, full):
        a = hp[mask]
        b = hp[full ^ mask]
        aa = [(s, t, a[s][t]) for s in range(n) for t in range(n) if a[s][t]]
        bb = [(s, t, b[s][t]) for s in range(n) for t in range(n) if b[s][t]]
        for s, t, v in aa:
            for u, w, z in bb:
                M[idx[s, u]][idx[t, w]] += v * z
    return (ports, M)

def direct_observer(A):
    n = len(A)
    ports = list(permutations(range(n), 2))
    idx = {x: i for i, x in enumerate(ports)}
    M = [[0] * len(ports) for _ in ports]
    for p in permutations(range(n)):
        bad = [j for j in range(n - 1) if not A[p[j]][p[j + 1]]]
        if len(bad) > 1:
            continue
        for k in range(1, n):
            if bad and bad[0] != k - 1:
                continue
            M[idx[p[0], p[k]]][idx[p[k - 1], p[-1]]] += 1
    return (ports, M)

def rref_kernel(M):
    a = [[Fraction(x) for x in row] for row in M]
    r = 0
    piv = []
    for c in range(len(a[0])):
        q = next((i for i in range(r, len(a)) if a[i][c]), None)
        if q is None:
            continue
        a[r], a[q] = (a[q], a[r])
        z = a[r][c]
        a[r] = [v / z for v in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [x - z * y for x, y in zip(a[i], a[r])]
        piv.append(c)
        r += 1
        if r == len(a):
            break
    out = []
    for c in range(len(a[0])):
        if c not in piv:
            v = [Fraction(0)] * len(a[0])
            v[c] = 1
            for i, p in enumerate(piv):
                v[p] = -a[i][c]
            need(all((sum((x * y for x, y in zip(row, v))) == 0 for row in M)), 'exact check failed')
            out.append(v)
    return out

def det(M):
    a = [row[:] for row in M]
    n = len(a)
    prev = 1
    sign = 1
    for k in range(n - 1):
        p = next((i for i in range(k, n) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[p], a[k] = (a[k], a[p])
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = a[i][j] * pivot - a[i][k] * a[k][j]
                need(num % prev == 0, 'exact check failed')
                a[i][j] = num // prev
            a[i][k] = 0
        prev = pivot
    return sign * a[-1][-1]

def sectors(ports, M):
    pairs = [p for p in ports if p[0] < p[1]]
    idx = {p: i for i, p in enumerate(ports)}
    return tuple(([[M[idx[a]][idx[b]] + sgn * M[idx[a]][idx[b[::-1]]] for b in pairs] for a in pairs] for sgn in (1, -1)))

def inverse(M):
    n = len(M)
    a = [[Fraction(x) for x in row] + [Fraction(i == j) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        q = next((i for i in range(c, n) if a[i][c]), None)
        if q is None:
            raise ValueError('singular matrix')
        a[c], a[q] = (a[q], a[c])
        z = a[c][c]
        a[c] = [v / z for v in a[c]]
        for i in range(n):
            if i != c:
                z = a[i][c]
                a[i] = [x - z * y for x, y in zip(a[i], a[c])]
    return [row[n:] for row in a]

def matmul(A, B):
    return [[sum((x * y for x, y in zip(row, col))) for col in zip(*B)] for row in A]

def insertion_schur(A, inc):
    n = len(A)
    B = [row + [inc[u]] for u, row in enumerate(A)] + [[1 - v for v in inc] + [0]]
    ports, N = observer(B)
    idx = {p: i for i, p in enumerate(ports)}
    old = list(permutations(range(n), 2))
    new = [(n, u) for u in range(n)] + [(u, n) for u in range(n)]
    block = lambda rows, cols: [[N[idx[r]][idx[c]] for c in cols] for r in rows]
    C, D, E, F = (block(old, old), block(old, new), block(new, old), block(new, new))
    H = hp_table(A)[(1 << n) - 1]
    need(F == [[H[i % n][j % n] if i // n == j // n else 0 for j in range(2 * n)] for i in range(2 * n)], 'exact check failed')
    Finv = inverse(F)
    Z = matmul(Finv, E)
    DZ = matmul(D, Z)
    R = [[x - y for x, y in zip(r, s)] for r, s in zip(C, DZ)]
    return (old, new, N, R, Z)
