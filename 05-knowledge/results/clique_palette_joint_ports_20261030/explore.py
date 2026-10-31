"""Exploratory weighted palettes for the actual Paley seven tournament."""
from itertools import combinations, product
from functools import lru_cache
import json

N = 7
ROWS = [sum(1 << j for j in range(N) if (j-i) % N in (1, 2, 4)) for i in range(N)]
def transitive(s):
    while s:
        v = next((i for i in range(N) if (s >> i) & 1 and (ROWS[i] & s).bit_count() == s.bit_count()-1), None)
        if v is None:
            return False
        s ^= 1 << v
    return True
GOOD = [s for s in range(1, 1 << N) if transitive(s)]
CYC = [sum(1 << i for i in t) for t in combinations(range(N), 3)
       if not transitive(sum(1 << i for i in t))]
@lru_cache(None)
def palette(d):
    active = sum(1 << i for i, x in enumerate(d) if x)
    if not active:
        return 0
    pivot = (active & -active).bit_length()-1
    eligible = [s for s in GOOD if (s & active) == s and (s >> pivot) & 1]
    maximal = [s for s in eligible if not any(s != t and (s & t) == s for t in eligible)]
    return 1+min(palette(tuple(x-((s >> i) & 1) for i, x in enumerate(d))) for s in maximal)

def bound(d):
    return max(max(d), (sum(d)+2)//3,
               max((sum(d[i] for i in range(N) if (s >> i) & 1)+1)//2 for s in CYC))

if __name__ == '__main__':
    hist, failures = {}, []
    for d in product(range(4), repeat=N):
        v, b = palette(d), bound(d)
        hist[v-b] = hist.get(v-b, 0)+1
        if v != b and len(failures) < 10:
            failures.append((d, v, b))
    result = {'good_sets': len(GOOD), 'cyclic_triples': len(CYC), 'gap_histogram': hist, 'failures': failures}
    try:
        import numpy as np
        from scipy.spatial import ConvexHull
        from fractions import Fraction as F
        points = np.array([[int(s >> i & 1) for i in range(N)] for s in [0]+GOOD])
        hull = ConvexHull(points)
        facets = set()
        for row in hull.equations:
            den = max(abs(v) for v in row)
            facets.add(tuple(str(F(float(v/den)).limit_denominator(1000)) for v in row))
        result['facets'] = sorted(facets)
    except ImportError:
        pass
    print(json.dumps(result, indent=2))
