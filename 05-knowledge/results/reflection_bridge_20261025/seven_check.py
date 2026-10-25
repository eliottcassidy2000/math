"""Independent finite controls for seven.md; no producer imports."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
from collections import Counter
from pathlib import Path
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def reflect(v, mask):
    return mask if v == mask else v ^ mask


for n in range(1, 9):
    mask = (1 << n) - 1
    supports = range(1, mask + 1)
    check(all(reflect(reflect(v, mask), mask) == v for v in supports), (n, "involution"))
    check(sum(reflect(v, mask) == v for v in supports) == 1, (n, "fixed"))
    check((len(supports)-1)//2 == 2**(n-1)-1, (n, "pairs"))

points = set(range(1, 8))
lines = {frozenset((x, y, x ^ y)) for x, y in combinations(points, 2)}
check(len(lines) == 7, "Fano line count")
preserved = {line for line in lines if frozenset(reflect(v, 7) for v in line) in lines}
check(preserved == {line for line in lines if 7 in line}, "exact preserved Fano lines")
check(len(preserved) == 3, "three preserved lines")
for perm in permutations(range(3)):
    act = lambda v: sum(((v >> i) & 1) << perm[i] for i in range(3))
    check(all(act(reflect(v, 7)) == reflect(act(v), 7) for v in points), "coordinate equivariance")
    check(all(frozenset(act(v) for v in line) in lines for line in lines), "coordinate collineation")

primes = (2, 3, 5)
value = lambda v: __import__("math").prod(primes[i] for i in range(3) if (v >> i) & 1)
squarefree = {d for d in range(2, 60) if 60 % d == 0 and all(d % (p*p) != 0 for p in primes)}
check({value(v) for v in points} == squarefree, "squarefree support map")
check(60//3 == 20 and 20 not in squarefree, "ten-point complement hostile")
check(Counter(v.bit_count() for v in points) == {1: 3, 2: 3, 3: 1}, "support rank counts")

projective_rows = []
for m in range(1, 11):
    frames = set()
    for vector in product((None, *range(m)), repeat=3):
        present = [i for i, v in enumerate(vector) if v is not None]
        if not present:
            continue
        first = vector[present[0]]
        frames.add(tuple(None if v is None else (v-first) % m for v in vector))
    ranks = Counter(sum(v is not None for v in frame) for frame in frames)
    check(ranks == {1: 3, 2: 3*m, 3: m*m}, (m, "projective supports"))
    check(len(frames) == 1+(m+1)+(m+1)**2, (m, "projective count"))
    projective_rows.append({"m": m, "total": len(frames), "support_counts": dict(sorted(ranks.items()))})


def mapping(point):
    x, y, z = point
    u = 1+x*y
    return (u**3*z+y*y*u*(4+3*x*y),
            y+3*x*u*u*z+3*x*y*y*(4+3*x*y),
            2*x-3*x*x*y-x**3*z)


sigma = lambda v: (-v[0], -v[1], v[2])
tau = lambda v: (v[0], -v[1], -v[2])
w = (F(-1, 8), F(0), F(0))
p0 = (F(0), F(0), F(-1, 4))
check(mapping(w) == p0, "central leaf")
check(mapping(sigma(w)) == (F(0), F(0), F(1, 4)), "failed iterate reflection")
check(mapping(sigma(w)) != sigma(mapping(w)), "same involution fails")
for point in product((F(-1), F(-1, 2), F(0), F(1, 2), F(1)), repeat=3):
    check(mapping(sigma(point)) == tau(mapping(point)), "two-involution covariance")

result = {"status": "FINITE-EXACT", "boolean_sizes": [1, 8],
          "reflection_cycles_on_seven": {"fixed": 1, "pairs": 3},
          "fano_lines": 7, "fano_lines_preserved_by_modified_complement": 3,
          "fano_lines_destroyed": 4, "coordinate_permutations": 6,
          "projective_frames": projective_rows,
          "historical_covariance_samples": 125,
          "central_leaf_image": [str(x) for x in mapping(w)],
          "reflected_leaf_image": [str(x) for x in mapping(sigma(w))]}
text = json.dumps(result, indent=2, sort_keys=True) + "\n"
Path(__file__).with_suffix(".json").write_text(text, encoding="utf-8")
print(text, end="")
