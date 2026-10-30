# Independent Session113 audit. All checks remain active under python -O.
from fractions import Fraction as F
from itertools import product, permutations
from math import gcd, factorial
import importlib.util, random, argparse
from pathlib import Path
parser = argparse.ArgumentParser(description='Independent exact LRC and composite-symmetry audit; no third-party packages required.')
parser.add_argument('packet_root', type=Path, help='Directory containing lrc/check_boundary_germ.py')
args = parser.parse_args()
s = importlib.util.spec_from_file_location('lrc', args.packet_root / 'lrc' / 'check_boundary_germ.py')
l = importlib.util.module_from_spec(s)
s.loader.exec_module(l)

def direct(V, r):
    e = sorted({F(0), F(1)} | {(k + r) / v for v in V for k in range(v)} | {(k + 1 - r) / v for v in V for k in range(v)})
    good = lambda t: all((r <= v * t % 1 <= 1 - r for v in V))
    return (sum((b - a for a, b in zip(e, e[1:]) if good((a + b) / 2)), F()), any(map(good, e)))
rng = random.Random(113)
count = 0
for _ in range(500):
    V = tuple(sorted(rng.sample(range(1, 31), rng.randrange(1, 9))))
    d = rng.randrange(3, 31)
    p = rng.randrange(1, (d + 1) // 2)
    r = F(p, d)
    A, B, marks = l.germ(V, r)
    eps = F(99, 100 * r.denominator * l.event_height(V))
    mass, nonempty = direct(V, r - eps)
    if not mass == A + eps * B:
        raise RuntimeError((V, r, eps))
    if not direct(V, r)[1] == bool(marks) == nonempty:
        raise RuntimeError('audit failed')
    if not (mass >= 2 * eps * gcd(*V) / max(V)) == nonempty:
        raise RuntimeError('audit failed')
    count += 1
print('LRC independent endpoint-midpoint checks:', count)

def compose(p, q):
    return tuple((p[q[x]] for x in range(len(p))))

def inverse(p):
    return tuple((p.index(i) for i in range(len(p))))

def closure(gs):
    n = len(gs[0])
    I = tuple(range(n))
    seen = {I}
    todo = [I]
    for p in todo:
        for a in gs:
            q = compose(a, p)
            if q not in seen:
                seen.add(q)
                todo.append(q)
    return seen
for r, m in [(3, 3), (3, 5), (5, 3), (3, 7), (3, 9), (5, 5)]:
    n = r * m
    a = tuple(((x + 1) % n for x in range(n)))
    v = tuple(((x + r) % n if x % r == 0 else x for x in range(n)))
    b = compose(v, compose(a, inverse(v)))
    G = closure([a, b])
    W = set()
    K = set()
    for t in range(r):
        for z in product(range(m), repeat=r):
            p = tuple(((x % r + t) % r + r * ((x // r + z[x % r]) % m) for x in range(n)))
            W.add(p)
            if (sum(z) - t) % gcd(r, m) == 0:
                K.add(p)
    if not G == K:
        raise RuntimeError('audit failed')
    if not b not in closure([a]):
        raise RuntimeError('audit failed')
    S = [[0] * (r + 3) for _ in range(r + 3)]
    S[0][0] = 1
    for t in range(1, r + 3):
        for j in range(1, t + 1):
            S[t][j] = S[t - 1][j - 1] + j * S[t - 1][j]
        aw = sum((sum((p[i] == i for i in range(n))) ** t for p in W)) // len(W)
        ag = sum((sum((p[i] == i for i in range(n))) ** t for p in G)) // len(G)
        formula = sum((factorial(r) // factorial(r - j) * S[t][j] * m ** (t - j) for j in range(1, min(r, t) + 1))) // r
        delta = (gcd(r, m) - 1) * factorial(r - 1) * S[t][r] * m ** (t - r) if t >= r else 0
        if not (aw, ag) == (formula, formula + delta):
            raise RuntimeError('audit failed')
    if n == 9:
        colors = {}
        c = 0
        for x in product(range(n), repeat=2):
            if x in colors:
                continue
            for p in G:
                colors[p[x[0]], p[x[1]]] = c
            c += 1
        twoclosed = []
        for p in permutations(range(n)):
            if all((colors[i, j] == colors[p[i], p[j]] for i in range(n) for j in range(n))):
                twoclosed.append(p)
        if not set(twoclosed) == W:
            raise RuntimeError('audit failed')
    print('Composite exact closure/Burnside PASS', r, m, len(G), len(W))
