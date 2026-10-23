"""Independent cyclic-coset Helly and subgroup completeness audit (stdlib)."""
from itertools import combinations, product
from math import gcd
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('audited_lattice', ROOT / 'lattice' / 'check.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
checks = 0


def require(value, why):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(why)


def generate(g, h, m, n):
    seen = {(0, 0)}
    stack = [(0, 0)]
    while stack:
        x, y = stack.pop()
        for a, b in (g, h):
            z = ((x + a) % m, (y + b) % n)
            if z not in seen:
                seen.add(z)
                stack.append(z)
    return frozenset(seen)


subgroups = 0
for m, n in product(range(1, 10), repeat=2):
    domain = list(product(range(m), range(n)))
    generated = {generate(g, h, m, n) for g, h in combinations(domain, 2)}
    generated.add(frozenset({(0, 0)}))
    normal = {module.subgroup(m, n, f) for f in module.forms(m, n)}
    require(generated == normal, ('two-generator subgroup census', m, n))
    subgroups += len(normal)

# Enumerate every triple of cyclic subgroup cosets through modulus 18,
# including repeated cosets, with a direct set-intersection oracle.
triples = 0
for b in range(1, 19):
    cosets = [frozenset(range(r, b, d))
              for d in range(1, b + 1) if b % d == 0 for r in range(d)]
    for x, y, z in product(cosets, repeat=3):
        intersect = x & y & z
        pairwise = bool(x & y) and bool(y & z) and bool(x & z)
        require(bool(intersect) == pairwise, ('cyclic Helly2', b, x, y, z))
        if intersect:
            require(len(intersect) == gcd(gcd(len(x), len(y)), len(z)),
                    ('constant middle multiplicity', b, x, y, z))
        triples += 1

# Minimal odd noncyclic boundary: three affine lines over F_3.
points = set(product(range(3), repeat=2))
lines = [{p for p in points if p[0] == 0},
         {p for p in points if p[1] == 0},
         {p for p in points if (p[0] + p[1]) % 3 == 1}]
require(all(lines[i] & lines[j] for i, j in combinations(range(3), 2)),
        'noncyclic pairwise positive hostile')
require(not set.intersection(*lines), 'noncyclic joint empty hostile')
print(json.dumps({'status': 'PASS', 'checks': checks,
                  'abstract_subgroups_orders_through_9': subgroups,
                  'cyclic_coset_triples_moduli_through_18': triples,
                  'noncyclic_hostile': 'three pairwise-intersecting affine lines in F_3^2'}))
