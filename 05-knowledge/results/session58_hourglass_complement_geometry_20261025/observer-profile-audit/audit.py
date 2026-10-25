def _publication_require(value, message):
    if not value:
        raise RuntimeError(message)
if '--negative-check' in __import__('sys').argv:
    _publication_require(False, 'injected publication guard failure')
import json, math
from pathlib import Path
from collections import Counter, deque
from fractions import Fraction
P = Path(__file__).resolve().parent.parent / 'observer-profile-algebra'
O = Path(__file__).resolve().parent.parent / 'observer-profile-audit'
n = 6
edges = [(u, v) for u in range(n) for v in range(u + 1, n)]
indices = {e: i for i, e in enumerate(edges)}
N = 1 << len(edges)
full = N - 1
generators = []
for a in range(n - 1):
    p = list(range(n))
    p[a], p[a + 1] = (p[a + 1], p[a])
    destinations = [indices[tuple(sorted((p[u], p[v])))] for u, v in edges]
    lookup = [0] * N
    for mask in range(1, N):
        bit = mask & -mask
        lookup[mask] = lookup[mask ^ bit] ^ 1 << destinations[bit.bit_length() - 1]
    toggle = 1 << indices[a, a + 1]
    generators.append((lookup, toggle))

def all_orbits(twisted):
    owner = [-1] * N
    orbits = {}
    for start in range(N):
        if owner[start] >= 0:
            continue
        orbit = []
        queue = deque([start])
        owner[start] = start
        while queue:
            x = queue.popleft()
            orbit.append(x)
            for action, toggle in generators:
                y = action[x] ^ (toggle if twisted else 0)
                if owner[y] < 0:
                    owner[y] = start
                    queue.append(y)
        orbits[start] = orbit
    return (owner, orbits)
gc, graphs = all_orbits(False)
tc, tournaments = all_orbits(True)
gpair = [min(gc[x], gc[x ^ full]) for x in range(N)]
converse_reps = sorted({min(tc[x], tc[x ^ full]) for x in tournaments})
cols = sorted(set(gpair))
A = []
sizes = {}
for t in converse_reps:
    hist = Counter((gpair[x] for x in tournaments[t]))
    mult = math.factorial(n) // len(tournaments[t])
    _publication_require(mult * len(tournaments[t]) == math.factorial(n), 'independent audit assertion failed')
    A.append([hist[c] * mult for c in cols])
    sizes[t] = len(set(tournaments[t]) | set(tournaments[tc[t ^ full]]))

def fraction_rank(A):
    B = [list(map(Fraction, row)) for row in A]
    k = 0
    for j in range(len(B[0])):
        pivot = next((i for i in range(k, len(B)) if B[i][j]), None)
        if pivot is None:
            continue
        B[k], B[pivot] = (B[pivot], B[k])
        a = B[k][j]
        B[k] = [x / a for x in B[k]]
        for i in range(k + 1, len(B)):
            a = B[i][j]
            if a:
                B[i] = [x - a * y for x, y in zip(B[i], B[k])]
        k += 1
    return k
certificate = json.loads((P / 'direct_certificate.json').read_text())
_publication_require(cols == certificate['columns'], 'independent audit assertion failed')
_publication_require(converse_reps == certificate['representatives'], 'independent audit assertion failed')
_publication_require(A == certificate['matrix'], 'independent audit assertion failed')
_publication_require(all((sum(row) == 720 for row in A)), 'independent audit assertion failed')
r = fraction_rank(A)
_publication_require(r == 32, 'independent audit assertion failed')
w = {int(k): v for k, v in certificate['signed_kernel'].items()}
_publication_require(sum(w.values()) == 0, 'independent audit assertion failed')
_publication_require(all((sum((w[t] * A[i][j] for i, t in enumerate(converse_reps))) == 0 for j in range(len(cols)))), 'independent audit assertion failed')
pos = {int(k): Fraction(v) for k, v in certificate['mixture_positive'].items()}
neg = {int(k): Fraction(v) for k, v in certificate['mixture_negative'].items()}
_publication_require(sum(pos.values()) == sum(neg.values()) == 1, 'independent audit assertion failed')
_publication_require(not set(pos) & set(neg), 'independent audit assertion failed')
_publication_require(all((sum((pos.get(t, 0) * A[i][j] for i, t in enumerate(converse_reps))) == sum((neg.get(t, 0) * A[i][j] for i, t in enumerate(converse_reps))) for j in range(len(cols)))), 'independent audit assertion failed')
result = {'independent_method': 'Adjacent-transposition BFS, no imports from source scripts', 'graph_isomorphism_classes': len(graphs), 'tournament_isomorphism_classes': len(tournaments), 'converse_classes': len(converse_reps), 'graph_complement_pair_columns': len(cols), 'distinct_profile_rows': len(set(map(tuple, A))), 'exact_rational_rank': r, 'direct_certificate_matrix_exact_match': True, 'mixture_equality_checked_exactly': True, 'positive_support': len(pos), 'negative_support': len(neg), 'positive_mass': sum((v for v in w.values() if v > 0))}
(O / 'independent_certificate.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
