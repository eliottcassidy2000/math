"""Independent literal ordered-cover enumeration; no imports from source probes."""
import itertools as it, math, json

def tournament(n, code):
    a = [[False] * n for _ in range(n)]
    for bit, (i, j) in enumerate(it.combinations(range(n), 2)):
        a[i][j] = bool(code & 1 << bit)
        a[j][i] = not a[i][j]
    return a

def states(a):
    n = len(a)
    v = tuple(range(n))
    H = [0] * (1 << n)
    P = []
    ends = [0] * n
    for k in range(n + 1):
        for p in it.permutations(v, k):
            if all((a[x][y] for x, y in zip(p, p[1:]))):
                mask = sum((1 << x for x in p))
                H[mask] += 1
                if p:
                    P.append((p[0], mask))
                if k == n:
                    ends[p[-1]] += 1
    full = (1 << n) - 1
    R = [0] * n
    for u, m in P:
        R[u] += H[full ^ m]
    return [(H[full], R[u], ends[u], H[full ^ 1 << u]) for u in v]

def group(a):
    n = len(a)
    return [p for p in it.permutations(range(n)) if all((a[i][j] == a[p[i]][p[j]] for i in range(n) for j in range(n)))]

def S(z):
    h, r, e, d = z
    return (h, 2 * r, 0, d)

def B(z):
    h, r, e, d = z
    return (h + r - e, 2 * r - d, d, d)

def fromrows(rows):
    return [[bool(r & 1 << j) for j in range(len(rows))] for r in rows]

def add(a, u, b):
    n = len(a)
    out = [r + [True] for r in a] + [[False] * (n + 1)]
    if b:
        out[u][n] = False
        out[n][u] = True
    return out
import sys
if '--fail' in sys.argv:
    raise RuntimeError('deliberate checker failure')
count = 0
first = {}
autchecks = 0
for n in range(1, 6):
    for code in range(1 << n * (n - 1) // 2):
        a = tournament(n, code)
        zs = states(a)
        G = group(a)
        for u, z in enumerate(zs):
            count += 1
            c = math.gcd(*z)
            g = sum((p[u] == u for p in G))
            autchecks += 1
            if not all((x % g == 0 for x in z)):
                raise RuntimeError('audit check failed')
            if not (z[0] // g not in (7, 21) and z[3] // g not in (7, 21)):
                raise RuntimeError('audit check failed')
            if not math.gcd(z[0], B(z)[0], B(S(z))[0], B(B(z))[0]) == c:
                raise RuntimeError('audit check failed')
            if c > 1 and g == 1:
                first.setdefault('rooted_rigid', [n, code, u, z, c])
            for name, f in [('S', S), ('B', B)]:
                zz = f(z)
                cc = math.gcd(*zz)
                if not cc % c == 0:
                    raise RuntimeError('audit check failed')
                if cc > c:
                    first.setdefault(name, [n, code, u, z, zz])
                for p in (3, 5, 7, 11):
                    if c % p:
                        h, r, e, d = z
                        rhs = h % p == r % p == d % p == 0 and e % p != 0 if name == 'S' else r % p == d % p == (h - e) % p == 0 and h % p != 0
                        if not (cc % p == 0) == rhs:
                            raise RuntimeError('audit check failed')
T = fromrows([24, 1, 19, 22, 2])
U = fromrows([58, 60, 57, 16, 32, 8])
tc = states(T)
uc = states(U)
if not tc[1] == uc[3] == (9, 9, 3, 3):
    raise RuntimeError('audit check failed')
if not (len(group(T)) == 1 and len(group(U)) == 9):
    raise RuntimeError('audit check failed')
if not (math.gcd(*tc[0]) == 1 and all((math.gcd(*z) == 3 for z in uc))):
    raise RuntimeError('audit check failed')
rev = lambda a: [list(row) for row in zip(*a)]
if not states(rev(T))[1] == (9, 21, 1, 3):
    raise RuntimeError('audit check failed')
if not math.gcd(*states(rev(U))[3]) == 3:
    raise RuntimeError('audit check failed')
if not states(add(T, 0, True))[0][0] == 19:
    raise RuntimeError('audit check failed')
if not states(add(rev(T), 1, True))[0][0] == 29:
    raise RuntimeError('audit check failed')
X = fromrows([62, 48, 2, 38, 44, 4])
xc = states(X)
xop = states(rev(X))
uop = states(rev(U))
if not (len(X) == len(U) == 6 and xc[2] == uc[3] == (9, 9, 3, 3)):
    raise RuntimeError('equal-order state check failed')
if not (len(group(X)) == 1 and xc[1] == (9, 11, 1, 1) and math.gcd(*xc[1]) == 1):
    raise RuntimeError('equal-order primitive reroot check failed')
if not (xop[2] == (9, 42, 0, 3) and uop[3] == (9, 36, 0, 3)):
    raise RuntimeError('equal-order converse mismatch check failed')
print(json.dumps({'literal_root_states': count, 'subgroup_checks': autchecks, 'first_occurrences': first, 'T_states': tc, 'U_states': uc, 'T_converse_root1': states(rev(T))[1], 'U_converse_root3': states(rev(U))[3], 'escape_H': [19, 29], 'X_states': xc, 'X_converse_root2': xop[2], 'equal_order': 6}, indent=2))
