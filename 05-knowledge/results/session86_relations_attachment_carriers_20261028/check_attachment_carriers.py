#!/usr/bin/env python3
"""Small independent controls for Session86's typed attachment carrier.

Standard library only. Every mathematical test survives python -O.
"""
import itertools as it
import json


checks = 0


def require(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)


def multiply(a, b):
    return tuple(sum(a[2*i+k] * b[2*k+j] for k in range(2))
                 for i in range(2) for j in range(2))


tables = list(it.product(range(2), repeat=4))
for a, b, c in it.product(tables, repeat=3):
    require(multiply(multiply(a, b), c) == multiply(a, multiply(b, c)),
            "weighted relation composition is not associative")
for a in it.product(range(3), repeat=4):
    for i, j in it.product(range(2), repeat=2):
        # Pinning the two exposed ports extracts exactly one coefficient.
        probe = tuple(int((k, l) == (i, j))
                      for k in range(2) for l in range(2))
        require(sum(x*y for x, y in zip(a, probe)) == a[2*i+j],
                "singleton context fails to extract a coefficient")


def project(relation, scope):
    return {tuple(x[i] for i in scope) for x in relation}


parity = {}
for n in (3, 4, 5):
    even = {x for x in it.product(range(2), repeat=n) if sum(x) % 2 == 0}
    odd = set(it.product(range(2), repeat=n)) - even
    for size in range(n):
        for scope in it.combinations(range(n), size):
            require(project(even, scope) == project(odd, scope),
                    "proper marginal distinguishes parity")
    require({x[1:]+x[:1] for x in even} == even, "even loses cyclic invariance")
    require({x[1:]+x[:1] for x in odd} == odd, "odd loses cyclic invariance")
    require(even != odd, "full parity relations unexpectedly agree")
    parity[n] = {"each_relation_size": len(even),
                 "all_proper_marginals_equal": True,
                 "both_cyclically_invariant": True}


def matchings(vertices):
    vertices = tuple(vertices)
    if not vertices:
        yield frozenset()
        return
    a = vertices[0]
    for b in vertices[1:]:
        rest = tuple(v for v in vertices if v not in (a, b))
        for m in matchings(rest):
            yield m | {tuple(sorted((a, b)))}


matching_counts = {}
for r in range(1, 5):
    universe = list(matchings(range(2*r)))
    expected = 1
    for j in range(1, 2*r, 2):
        expected *= j
    require(len(universe) == expected, "wrong perfect matching universe")
    for m, n in it.combinations(universe, 2):
        edge = next(iter(n-m))
        require(edge in n and edge not in m, "matching guard hostile absent")
    for m in universe:
        partners = {a: b for edge in m for a, b in (edge, edge[::-1])}
        for a, b in it.combinations(range(2*r), 2):
            if partners[a] == b:
                continue
            removed = {tuple(sorted((a, partners[a]))),
                       tuple(sorted((b, partners[b])))}
            nxt = (m-removed) | {tuple(sorted((partners[a], partners[b])))}
            require(len(nxt) == r-1 and
                    set(it.chain.from_iterable(nxt)) == set(range(2*r))-{a,b},
                    "matching update is not closed")
    matching_counts[r] = len(universe)


def c3_action_fixed(action, element):
    return sum(action(element, x) == x for x in range(3))


trivial = [c3_action_fixed(lambda g, x: x, g) for g in range(3)]
regular = [c3_action_fixed(lambda g, x: (g+x) % 3, g) for g in range(3)]
require(trivial[0] == regular[0] == 3, "same total mass control fails")
require(sum(trivial)//3 == 3 and sum(regular)//3 == 1,
        "hidden action fails to change orbit count")

# All element-fixed counts can agree while subgroup-fixed counts differ.
g9 = tuple(it.product(range(3), repeat=2))
add9 = lambda x, y: ((x[0]+y[0]) % 3, (x[1]+y[1]) % 3)
lines9 = sorted({frozenset(((k*x) % 3, (k*y) % 3) for k in range(3))
                  for x, y in g9 if (x, y) != (0, 0)}, key=lambda s: sorted(s))
cosets9 = [sorted({frozenset(add9(g, h) for h in line) for g in g9},
                  key=lambda s: sorted(s)) for line in lines9]
x9 = tuple((0, g) for g in g9) + tuple((1, i) for i in range(3))
y9 = tuple((i, coset) for i, cosets in enumerate(cosets9) for coset in cosets)


def act_x9(g, x):
    tag, point = x
    return (0, add9(g, point)) if tag == 0 else x


def act_y9(g, y):
    tag, coset = y
    return tag, frozenset(add9(g, point) for point in coset)


char_x9 = [sum(act_x9(g, x) == x for x in x9) for g in g9]
char_y9 = [sum(act_y9(g, y) == y for y in y9) for g in g9]
full_x9 = sum(all(act_x9(g, x) == x for g in g9) for x in x9)
full_y9 = sum(all(act_y9(g, y) == y for g in g9) for y in y9)
require(len(lines9) == 4 and len(x9) == len(y9) == 12, "wrong C3^2 sets")
require(char_x9 == char_y9 == [12]+[3]*8, "C3^2 characters differ")
require((full_x9, full_y9) == (3, 0), "C3^2 full fixed counts do not differ")
require(sum(char_x9)//9 == sum(char_y9)//9 == 4, "C3^2 orbit counts differ")


def tournament(odd_block):
    # A+ = C3 -> singleton. A- is its exact converse on the same four labels.
    # The triangle is local 0,1,2 and the sink is local 3. HOSTILES.md uses
    # the distinct source-at-local-0 convention. Its named-port-8 marker
    # count is not asserted here without transporting the labels.
    a = [[0]*12 for _ in range(12)]
    for u, v in it.combinations(range(12), 2):
        bu, bv = u//4, v//4
        iu, iv = u % 4, v % 4
        if bu != bv:
            forward = (bv-bu) % 3 == 1
        else:
            forward = iv == 3 or (iu < 3 and iv < 3 and (iv-iu) % 3 == 1)
            if bu == odd_block:
                forward = not forward
        a[u][v] = int(forward)
        a[v][u] = int(not forward)
    return a


def cyclic_triples(a):
    return {triple for triple in it.combinations(range(len(a)), 3)
            if all(sum(a[u][v] for v in triple if v != u) == 1 for u in triple)}


def hamilton_count(a):
    n = len(a)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            prior = mask ^ (1 << v)
            if prior:
                dp[mask][v] = sum(dp[prior][u] for u in range(n)
                                  if (prior >> u) & 1 and a[u][v])
    return sum(dp[-1])


def automorphism_count(a):
    n = len(a)
    degree = [sum(row) for row in a]
    candidates = {u: [v for v in range(n) if degree[u] == degree[v]]
                  for u in range(n)}
    order = sorted(range(n), key=lambda u: (len(candidates[u]), u))
    image = {}
    used = set()

    def search(k):
        if k == n:
            return 1
        u = order[k]
        total = 0
        for v in candidates[u]:
            if v in used or any(a[u][w] != a[v][z] for w, z in image.items()):
                continue
            image[u] = v
            used.add(v)
            total += search(k+1)
            used.remove(v)
            del image[u]
        return total

    return search(0)


plus, minus = tournament(None), tournament(2)
cplus, cminus = cyclic_triples(plus), cyclic_triples(minus)
hplus, hminus = hamilton_count(plus), hamilton_count(minus)
aplus, aminus = automorphism_count(plus), automorphism_count(minus)
require(cplus == cminus and len(cplus) == 67, "81/27 cyclic hypergraphs differ")
require(hplus == hminus == 462915, "81/27 Hamilton count mismatch")
require((aplus, aminus) == (81, 27), "81/27 automorphism control mismatch")

print(json.dumps({"status": "PASS", "explicit_checks": checks,
                  "boolean_2_by_2_composition_triples": 16**3,
                  "parity_hostiles": parity,
                  "matching_states": matching_counts,
                  "hidden_C3_action": {"trivial_element_fixed_counts": trivial,
                                       "regular_element_fixed_counts": regular,
                                       "orbit_counts": [3, 1]},
                  "element_character_loss_C3_squared": {
                      "common_character": char_x9,
                      "common_orbit_count": 4,
                      "full_group_fixed_counts": [full_x9, full_y9]},
                  "tournament_81_27": {"same_labelled_cyclic_triples": 67,
                                       "Hamilton_counts": [hplus, hminus],
                                       "automorphism_orders": [aplus, aminus],
                                       "Hamilton_word_orbits":
                                       [hplus//aplus, hminus//aminus]}},
                 indent=2, sort_keys=True))
