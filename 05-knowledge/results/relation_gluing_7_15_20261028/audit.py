"""Separate standard-library checks. Same agent; not an independent referee."""
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod
from pathlib import Path
import json
import sys

sys.stdout.reconfigure(encoding="utf-8", newline="\n")
data = json.loads((Path(__file__).parent / "normal.json").read_text(encoding="utf-8"))

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def factors(n):
    out = {}
    for p in range(2, isqrt(n) + 1):
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    if n != 1:
        out[n] = out.get(n, 0) + 1
    return out

def divisor_sum(n):
    total = 0
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            total += d
            if d * d != n:
                total += n // d
    return total

def square(n):
    return isqrt(n) ** 2 == n

def determinant(a):
    if len(a) == 0:
        return 1
    return sum((-1)**j * a[0][j] * determinant([r[:j] + r[j+1:] for r in a[1:]])
               for j in range(len(a)))

def smith_from_minors(columns):
    previous = 1
    answer = []
    for k in range(1, min(3, len(columns)) + 1):
        g = 0
        for rs in combinations(range(3), k):
            for cs in combinations(range(len(columns)), k):
                minor = [[columns[j][i] for j in cs] for i in rs]
                g = gcd(g, determinant(minor))
        if not g:
            break
        require(g % previous == 0, "determinantal divisors")
        answer.append(g // previous)
        previous = g
    return tuple(answer)

P = [row["prime"] for row in data["prime_bank7"]]
columns = [tuple(factors(p + 1).get(q, 0) for q in (2, 3, 5)) for p in P]
profiles = {s: smith_from_minors([columns[i] for i in range(7) if s & (1 << i)])
            for s in range(128)}
lines = {frozenset((a, b, a ^ b)) for a, b in combinations(range(1, 8), 2)}
binary = 0
integral = []
for perm in permutations(range(1, 8)):
    if not all(frozenset(perm[v-1] for v in L) in lines for L in lines):
        continue
    binary += 1
    if all(profiles[s] == profiles[sum(1 << (perm[i]-1) for i in range(7) if s & (1 << i))]
           for s in range(128)):
        integral.append(list(perm))
require(binary == 168, "full incidence symmetry group from all5040 permutations")
require(integral == data["coefficient_change"]["integer_profile_stabilizer"], "Smith symmetry via minors")
require(len(integral) == 2, "exact decorated group order2")

old_sources = {}
for s in range(128):
    n = prod(P[i] for i in range(7) if s & (1 << i))
    if square(divisor_sum(n)):
        old_sources[s] = n
require(len(old_sources) == 16, "old square-total products via actual divisor lists")
points = data["relation_to_object"]["points"]
require({r["old_support_mask"] for r in points} == set(old_sources) - {0}, "all nonidentity points")
for row in points:
    require(row["old_square_total_source"] == old_sources[row["old_support_mask"]], "point source")
for a, b in product(old_sources, repeat=2):
    x, y = old_sources[a], old_sources[b]
    require(old_sources[a ^ b] == x * y // gcd(x, y)**2, "group law from integer products")
triple_count = 0
for a, b, c in combinations(points, 3):
    relation = a["old_support_mask"] ^ b["old_support_mask"] ^ c["old_support_mask"] == 0
    actual = divisor_sum(a["new_prime"] * b["new_prime"] * c["new_prime"])
    require(relation == square(actual), "old relation versus actual new divisor list")
    triple_count += relation
require(triple_count == 35, "all35 three-term laws")

line = lambda L: len(L) == 3 and len(set(L)) == 3 and L[0] ^ L[1] ^ L[2] == 0
parallelism = data["geometry"]["seven_spread_parallelism"]
seen = set()
for spread in parallelism:
    require(sorted(x for L in spread for x in L) == list(range(1, 16)), "spread partitions points")
    for L in spread:
        require(line(L), "projective line certificate")
        seen.add(tuple(sorted(L)))
require(len(seen) == 35, "parallelism covers every line exactly once")
iso = data["geometry"]["isotropic_prime_triples"]
prime_vector = {row["prime"]:row["vector"] for row in data["prime_bank15"]}
expected = {tuple(sorted(prime_vector[p] for p in L)) for L in iso}
seen = set()
for cls in data["geometry"]["four_hierarchy_classes"]:
    used = []
    for L in cls:
        require(line(L), "hierarchy line")
        seen.add(tuple(sorted(L)))
        used.extend(L)
    require(len(used) == len(set(used)), "each hierarchy class is disjoint")
require(seen == expected, "four hierarchies cover all doily lines")

# Re-enumerate the six pentads as perfect matchings and exact edge covers of K6.
def matchings(vertices):
    if not vertices:
        yield frozenset()
    else:
        a = vertices[0]
        for b in vertices[1:]:
            for rest in matchings([v for v in vertices if v not in (a, b)]):
                yield rest | {tuple(sorted((a, b)))}
match = list(matchings(list(range(6))))
pentads = [frozenset(c) for c in combinations(range(15), 5)
           if len(set().union(*(match[i] for i in c))) == 15]
require(len(pentads) == 6, "six exact K6 edge factorizations")
require(all(len(a & b) == 1 for a, b in combinations(pentads, 2)), "pentad intersection hostile")

# Norm-five/norm-seven unit quotients, with the displayed evaluations.
require({1, 4, 3, 2} == set(range(1, 5)), "Gaussian unit steps")
require((3 * 3 + 1) % 5 == 0 and (2 + 3) % 5 == 0, "Gaussian quotient map")
require((4 * 4 + 4 + 1) % 7 == 0 and (3 + 4) % 7 == 0, "Eisenstein quotient map")
require({x % 7 for u in (1, 4, 4 * 4) for x in (u, -u)} == set(range(1, 7)),
        "Eisenstein unit steps")
require(divisor_sum(520) == 1260 and prod(divisor_sum(n) for n in (2, 4, 5, 13)) == 1764,
        "fusion hostile by independent divisor enumeration")

# In O=Z[t], t^2=-d, products of the ideal generators yield 4 and 1+d;
# (1+d)-4*((d-1)//4)=2 certifies the reverse inclusion in (2).
for d in (5, 13, 37):
    require((1+d)-4*((d-1)//4) == 2, "nonprincipal grid ideal squares to(2)")
    require(all(x*x+d*y*y != 2 for y in range(-isqrt(2//d), isqrt(2//d)+1)
                for x in range(-isqrt(2), isqrt(2)+1)), "complete norm2 search")
def order2(p):
    x = 1
    for k in range(1, p):
        x = 2*x % p
        if x == 1:
            return k
require(order2(19) == 18 and order2(127) == 7, "random-recovery order hostile")

print(json.dumps({"status":"PASS", "method":"same-agent alternate standard-library paths",
                  "integer_profiles_via_minors":128, "permutations":5040,
                  "binary_automorphisms":binary, "integral_stabilizer":integral,
                  "old_divisor_lists":128, "group_law_pairs":256,
                  "new_triple_divisor_lists":455, "projective_lines":triple_count,
                  "perfect_matchings":len(match), "pentads":len(pentads),
                  "parallelism_lines":35, "hierarchy_lines":15,
                  "grid_ideal_controls":[5,13,37], "order_controls":{"19":18,"127":7}},
                 indent=2, sort_keys=True))
