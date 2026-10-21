def require(ok):
 if not ok: raise AssertionError("independent audit invariant")

"""Independent bounded-prefix auditor. No producer modules or eventual routes."""
import json, hashlib, itertools, random, sys
from pathlib import Path
from collections import defaultdict, Counter
OUT = Path(__file__).resolve().parent
ROOT = next(p for p in [OUT.parent / 'math', *OUT.parents] if (p / '05-knowledge/results/session18_native_annulus/frontier.json').is_file())
FRONTIER = ROOT / '05-knowledge/results/session18_native_annulus/frontier.json'
require(hashlib.sha256(FRONTIER.read_bytes()).hexdigest() == '391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f')
FD = json.loads(FRONTIER.read_text())
V = set(FD['frontier'])
require(len(V) == 627 and max(V) == 63665)
A = ((1, 0), (2, 1), (1, 1), (1, 2), (0, 1), (-1, 1), (-1, 0), (-2, -1), (-1, -1), (-1, -2), (0, -1), (1, -1))

def q(v):
    x, y = v
    return x * x - x * y + y * y

def b(v, w):
    return q((v[0] + w[0], v[1] + w[1])) - q(v) - q(w)

def primitive(s):
    return type(s) is int and s > 0 and (s % 2 == 1) and (s % 8 != 5)

def step(s):
    require(primitive(s))
    k = 2 if s % 8 == 1 else 1
    u = 3 * s + 1 >> k
    while u % 8 == 5:
        u = u - 1 >> 2
        k += 2
    require(primitive(u))
    return (u, k)
RANK = dict(FD['frontier_ranks'])
require(set(RANK) == V and RANK[1] == 0)
for s in V:
    if s != 1:
        u, _ = step(s)
        require(u in V and RANK[u] < RANK[s])

def charge_map(N):
    z = {s: [0, 0] for s in range(1, 4 * N, 2) if primitive(s)}
    lower = ((1, 2), (-2, -1), (1, -1))
    upper = ((0, 1), (-1, -1), (1, 0))
    for m in range(1, 4 * N, 2):
        s = m
        while s % 8 == 5:
            s = (s - 1) // 4
        a, c = (lower if m < 2 * N else upper)[m % 3]
        z[s][0] += a
        z[s][1] += c
    return {s: tuple(v) for s, v in z.items()}

class Prefix:

    def __init__(self, d):
        self.d = d
        self.cache = {s: (s, 0, 0) for s in V}
        self.edges = 0

    def route(self, s):
        path = []
        u = s
        while u not in self.cache and len(path) < self.d:
            v, k = step(u)
            path.append((u, k))
            u = v
            self.edges += 1
        if u not in self.cache:
            return None
        v, t, h = self.cache[u]
        if h + len(path) > self.d:
            return None
        for x, k in reversed(path):
            t += k
            h += 1
            self.cache[x] = (v, t, h)
        return self.cache[s]

def polygon(counts):
    total = [sum((counts[v] * v[k] for v in A)) for k in (0, 1)]
    return max((q(tuple((2 * sum((counts[A[j % 12]] * A[j % 12][k] for j in range(i, i + 6))) - total[k] for k in (0, 1)))) for i in range(12)))

def validate_rows(rows, N, d):
    seen = set()
    for s, r15, r35, u, t, h, x, y in rows:
        require(s not in seen)
        seen.add(s)
        require(primitive(s) and s < 256 * N)
        require(r15 == s % 15 and r35 == s % 35)
        require(r15 % 5 == r35 % 5 and (r15 + 15 * ((r35 - r15) % 7)) % 105 == s % 105)
        require(type(t) is int and type(h) is int and (0 <= h <= d))
        v = s
        clock = 0
        for _ in range(h):
            require(v not in V)
            v, k = step(v)
            clock += k
        require(v == u and v in V and (clock == t))

def rejected(fn):
    try:
        fn()
    except AssertionError:
        return True
    return False

def run(N, d):
    z = charge_map(64 * N)
    old = charge_map(N)
    p = Prefix(d)
    W = defaultdict(lambda: [0, 0])
    H = defaultdict(lambda: [0, 0])
    pending = Counter()
    oldpending = []
    rows = []
    digest = hashlib.sha256()
    zero = 0
    newnorm = 0
    pending_ids = []
    maxold = 0
    for s, v in z.items():
        if s < 4 * N:
            require(v == old[s])
        else:
            newnorm += q(v)
        route = p.route(s)
        if route is None:
            if s < 4 * N:
                oldpending.append(s)
            elif v != (0, 0):
                pending[v] += 1
                pending_ids.append(s)
            else:
                zero += 1
            continue
        u, t, h = route
        if s < 4 * N:
            maxold = max(maxold, h)
        row = [s, s % 15, s % 35, u, t, h, *v]
        digest.update((json.dumps(row, separators=(',', ':')) + '\n').encode())
        if len(rows) < 100 or h == d:
            rows.append(row)
        target = W if s < 4 * N else H
        sg = 1 if t % 2 == 0 else -1
        target[u][0] += sg * v[0]
        target[u][1] += sg * v[1]
    require(newnorm == 153 * N)
    require(not oldpending)
    require(all((v in A for v in pending)))
    T = {u: tuple((W[u][k] + H[u][k] for k in (0, 1))) for u in V}
    delta = sum((q(T[u]) - q(W[u]) for u in V))
    D = polygon(pending)
    L = sum((n * max((abs(b(t, v)) for t in T.values())) for v, n in pending.items()))
    validate_rows(rows, N, d)
    good = next((r for r in rows if r[5] > 0))
    controls = {}
    for name, col, amount in [('clock', 4, 2), ('mod15', 1, 5), ('mod35', 2, 5)]:
        bad = good.copy()
        bad[col] += amount
        controls[name] = rejected(lambda: validate_rows([bad], N, d))
    controls['source_reuse'] = rejected(lambda: validate_rows([good, good], N, d))
    require(all(controls.values()))
    result = dict(N=N, d=d, max_old_depth=maxold, nonzero_pending=sum(pending.values()), zero_pending=zero, delta=delta, D=D, L=L, bound=delta + D + L, budget=153 * N, pass_suffix=delta + D + L <= 153 * N, source_count=len(z), evaluated_native_edges=p.edges, arrived_rows_sha256=digest.hexdigest(), pending_source_ids=pending_ids, pending_counts={str(v): n for v, n in pending.items()}, hostile_rejections=controls)
    print({k: v for k, v in result.items() if k not in ('pending_source_ids', 'pending_counts')}, flush=True)
    return result

def algebra_controls():
    rng = random.Random(22021)
    cases = 0
    for n in range(0, 11):
        for _ in range(20):
            vs = [rng.choice(A + ((0, 0),)) for j in range(n)]
            D = max((q((sum((e * v[0] for e, v in zip(es, vs))), sum((e * v[1] for e, v in zip(es, vs))))) for es in itertools.product((-1, 1), repeat=n)))
            require(D == polygon(Counter(vs)))
            cases += 1
            T = [(2, -1), (0, 0), (-3, 2)]
            groups = [[0, 0] for _ in T]
            for v in vs:
                u = rng.randrange(4)
                sg = rng.choice((-1, 1))
                if u < 3:
                    groups[u][0] += sg * v[0]
                    groups[u][1] += sg * v[1]
            require(sum((q(g) for g in groups)) <= D)
            L = sum((max((abs(b(t, v)) for t in T)) for v in vs))
            require(sum((q(g) + b(t, g) for t, g in zip(T, groups))) <= D + L)
    require(q((2, 0)) > polygon(Counter({(1, 0): 1})))
    for _ in range(1000):
        delta, budget, arrived, J = [rng.randrange(-100, 100) for i in range(4)]
        R = delta - arrived - J
        M = budget - arrived
        require((R <= M - J) == (delta <= budget))
    return dict(brute_sign_polygon_cases=cases, random_partial_group_cases=cases, extraction_equivalences=1000, zero_and_absent_terminal=True, reuse_counterexample=True)

def general_controls():
    from math import isqrt
    vsold = [(1, 0), (1, 2)]
    vsnew = [(-1, 1), (0, 1)]
    W = [(1, -1), (-2, -1), (0, 0)]
    H = [(2, 1), (0, -1), (0, 0)]
    T = [tuple((w[k] + h[k] for k in (0, 1))) for w, h in zip(W, H)]
    base = sum((q(t) - q(w) for t, w in zip(T, W)))
    Do = polygon(Counter(vsold))
    Dn = polygon(Counter(vsnew))
    Ln = sum((max((abs(b(t, v)) for t in T)) for v in vsnew))
    Lo = sum((max((abs(b(h, v)) for h in H)) for v in vsold))
    x = 4 * Do * Dn
    r = isqrt(x)
    mixed = r + (r * r < x)
    upper = base + Dn + Ln + Lo + mixed
    high = -10 ** 10
    count = 0
    for choices in itertools.product(range(7), repeat=4):
        zo = [[0, 0] for _ in W]
        zn = [[0, 0] for _ in W]
        for i, (choice, v) in enumerate(zip(choices, vsold + vsnew)):
            if choice:
                u = (choice - 1) // 2
                sign = 1 if choice % 2 else -1
                target = zo if i < 2 else zn
                for k in (0, 1):
                    target[u][k] += sign * v[k]
        actual = sum((q(tuple((T[u][k] + zo[u][k] + zn[u][k] for k in (0, 1)))) - q(tuple((W[u][k] + zo[u][k] for k in (0, 1)))) for u in range(3)))
        expansion = base + sum((q(zn[u]) + b(T[u], zn[u]) + b(H[u], zo[u]) + b(zo[u], zn[u]) for u in range(3)))
        require(actual == expansion and actual <= upper)
        high = max(high, actual)
        count += 1
    require(q((2, 0)) - q((1, 0)) == 3 and q((1, 0)) == 1)
    return dict(assignments=count, empty_terminal_included=True, D_old=Do, D_new=Dn, L_old=Lo, L_new=Ln, mixed_ceiling=mixed, upper=upper, actual_maximum=high, ignored_old_hostile=[1, 3])
if __name__ == '__main__':
    pairs = [(4096, 104), (4096, 105), (8192, 111), (8192, 112), (7, 8), (7, 22)] if len(sys.argv) == 1 else [tuple(map(int, s.split(':'))) for s in sys.argv[1:]]
    result = {'method': 'Independent odd-row ancestor charges and strictly bounded native prefixes; no eventual routes consumed', 'algebra_controls': algebra_controls(), 'general_pending_old_controls': general_controls(), 'cases': [run(N, d) for N, d in pairs]}
    (OUT / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
