#!/usr/bin/env python3
# Checks intentionally survive python -O.
import sys
if "--fail" in sys.argv:
    raise RuntimeError("deliberate checker failure")
import itertools, math, json

def rows(n, m):
    a = [0] * n
    for k, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        u, v = (i, j) if m >> k & 1 else (j, i)
        a[u] |= 1 << v
    return a

def state(a, u):
    n = len(a)
    N = 1 << n
    dp = [[0] * n for _ in range(N)]
    H = [0] * N
    H[0] = 1
    for i in range(n):
        dp[1 << i][i] = 1
    for s in range(1, N):
        for v in range(n):
            if s >> v & 1 and s != 1 << v:
                dp[s][v] = sum((dp[s ^ 1 << v][w] for w in range(n) if a[w] >> v & 1))
        H[s] = sum(dp[s])
    rev = [sum((1 << j for j in range(n) if a[j] >> i & 1)) for i in range(n)]
    ds = [[0] * n for _ in range(N)]
    ds[1 << u][u] = 1
    for s in range(1, N):
        for v in range(n):
            if s >> v & 1 and s != 1 << v:
                ds[s][v] = sum((ds[s ^ 1 << v][w] for w in range(n) if a[w] >> v & 1))
    R = sum((sum(ds[s]) * H[N - 1 ^ s] for s in range(N)))
    return (H[-1], R, dp[-1][u], H[N - 1 ^ 1 << u])

def S(z):
    h, R, e, d = z
    return (h, 2 * R, 0, d)

def B(z):
    h, R, e, d = z
    return (h + R - e, 2 * R - d, d, d)

def aut(a, u):
    n = len(a)
    out = []
    for p in itertools.permutations(range(n)):
        if p[u] != u:
            continue
        if all((a[i] >> j & 1 == a[p[i]] >> p[j] & 1 for i in range(n) for j in range(n))):
            out.append(p)
    return len(out)

def append(a, u, b):
    n = len(a)
    z = [r | 1 << n for r in a] + [0]
    if b:
        z[u] ^= 1 << n
        z[n] |= 1 << u
    return z
if __name__ == '__main__':
    found = {}
    tested = 0
    checks = 0
    for n in range(1, 7):
        for m in range(1 << n * (n - 1) // 2):
            a = rows(n, m)
            for u in range(n):
                z = state(a, u)
                g = math.gcd(*z)
                tested += 1
                if not math.gcd(z[0], B(z)[0], B(S(z))[0], B(B(z))[0]) == g:
                    raise RuntimeError('live mathematical check failed')
                for key, zz in [('S', S(z)), ('B', B(z))]:
                    gg = math.gcd(*zz)
                    if gg > g and key not in found:
                        found[key] = {'n': n, 'mask': m, 'root': u, 'state': z, 'child': zz, 'content': g, 'child_content': gg, 'aut_root': aut(a, u), 'rows': a}
                        print(key, found[key], flush=True)
                if g > 1 and 'rigid' not in found and (aut(a, u) == 1):
                    found['rigid'] = {'n': n, 'mask': m, 'root': u, 'state': z, 'content': g, 'rows': a}
                    print('rigid', found['rigid'], flush=True)
                if n <= 4:
                    for l in range(4):
                        for w in itertools.product(range(2), repeat=l):
                            az = a
                            zz = z
                            for b in w:
                                az = append(az, u, b)
                                zz = (B if b else S)(zz)
                            if not state(az, u) == zz:
                                raise RuntimeError((n, m, u, w, state(az, u), zz))
                            checks += 1
        print('done', n, tested, flush=True)
    print(json.dumps({'tested': tested, 'physical_checks': checks, 'witnesses': found}, indent=2))
