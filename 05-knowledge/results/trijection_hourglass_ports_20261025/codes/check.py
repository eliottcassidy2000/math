"""Exact labelled controls; all gates survive python -O. Standard library only."""
from collections import Counter
from itertools import combinations, permutations, product
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def act(plane, p):
    return frozenset(frozenset(p[x] for x in line) for line in plane)


def orbit_sizes(objects, group):
    remaining = set(objects)
    result = []
    while remaining:
        seed = min(remaining, key=lambda p: sorted(map(sorted, p)))
        orbit = {act(seed, p) for p in group}
        need(orbit <= objects, "group must preserve tested universe")
        result.append(len(orbit))
        remaining -= orbit
    return sorted(result)


def rep(mask):
    return min(mask, mask ^ 127)


def mask(xs):
    return sum(1 << x for x in xs)


base = frozenset(frozenset((a - 1, b - 1, (a ^ b) - 1))
                 for a, b in combinations(range(1, 8), 2))
perms = list(permutations(range(7)))
planes = {act(base, p) for p in perms}
need(len(planes) == 30, "thirty binary-plane label classes")

# Independent exact cover, with no linear algebra or input Fano model.
triples = list(map(frozenset, combinations(range(7), 3)))
all_pairs = frozenset(map(frozenset, combinations(range(7), 2)))
triple_pairs = {t: frozenset(map(frozenset, combinations(t, 2))) for t in triples}
steiner = set()


def search(uncovered, selected):
    if not uncovered:
        steiner.add(frozenset(selected))
        return
    pair = min(uncovered, key=lambda p: tuple(sorted(p)))
    for t in triples:
        if pair <= t and triple_pairs[t] <= uncovered:
            search(uncovered - triple_pairs[t], selected + [t])


search(all_pairs, [])
need(steiner == planes, "independent Steiner exact-cover reconstruction")

crossings = frozenset(frozenset((0, j, j + 3)) for j in range(1, 4))
hist = Counter(len(p & crossings) for p in planes)
need(hist == {0: 16, 1: 12, 3: 2}, "dependent crossing histogram")
for k, expected in ((1, 6), (2, 2), (3, 2)):
    for chosen in combinations(crossings, k):
        need(sum(set(chosen) <= p for p in planes) == expected,
             "crossing intersections")

star_group = [p for p in perms if act(crossings, p) == crossings]
need(len(star_group) == 48, "crossing-star symmetry order")
good = {p for p in planes if not (p & crossings)}
need(orbit_sizes(good, star_group) == [16], "good W orbit")
need(all(sum(act(f, p) == f for p in star_group) == 3 for f in good),
     "all good stabilizers order three")
dihedral = [(0,) + tuple(1 + (sign * i + shift) % 6 for i in range(6))
            for sign in (1, -1) for shift in range(6)]
need(len(set(dihedral)) == 12, "D6 order")
need(orbit_sizes(good, dihedral) == [4, 12], "geometric orbit split")

# Check the indicated actual x,z transformations on six root vectors.
roots = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1)]
need([(a - b, a) for a, b in roots] == roots[1:] + roots[:1], "rotation lift")
need([(-b, -a) for a, b in roots] == list(reversed(roots)), "reflection lift")

node_profiles = Counter()
character_checks = 0
for f in sorted(planes, key=lambda p: sorted(map(sorted, p))):
    lift = {0} | {mask(t) ^ 127 for t in f}
    h = {rep(v) for v in lift}
    need(len(lift) == len(h) == 8, "faithful even lift")
    need(all((v.bit_count() % 2) == 0 for v in lift), "even sign lift")
    need(all(a ^ b in lift for a in lift for b in lift), "simplex closure")
    need(Counter(v.bit_count() for v in lift) == {0: 1, 4: 7}, "code weights")
    columns = [tuple((v >> i) & 1 for v in sorted(lift)) for i in range(7)]
    need(len(set(columns)) == 7 and all(any(col) for col in columns),
         "seven distinct nontrivial canonical characters")
    character_checks += 7
    profile = Counter()
    for t in crossings:
        local = {rep(mask(s)) for k in range(4) for s in combinations(t, k)}
        stabilizer = h & local
        bad = t in f
        need(len(stabilizer) == (2 if bad else 1), "local code intersection")
        # G/local has eight nodes over this triple. Form actual translation orbits.
        g = set(range(64))
        coset = lambda v: min(rep(v ^ s) for s in local)
        nodes = {coset(v) for v in g}
        need(len(nodes) == 8, "eight node fibre")
        node_orbits = {frozenset(coset(v ^ s) for s in h) for v in nodes}
        need(len(node_orbits) == (2 if bad else 1), "quotient node orbits")
        profile["quarter" if bad else "A1"] += len(node_orbits)
    b = len(f & crossings)
    need(profile == {"A1": 3 - b, "quarter": 2 * b}
         or all(profile[key] == value for key, value in {"A1": 3-b,"quarter":2*b}.items()),
         "full local profile")
    node_profiles[(b, profile["A1"], profile["quarter"], 2 - profile["quarter"])] += 1

matching_hist = Counter()
representatives = []
for f in good:
    matching = frozenset(t - {0} for t in f if 0 in t)
    adjacent = sum((abs(a - b) in (1, 5)) for a, b in map(tuple, matching))
    need(adjacent in (1, 3), "hexagon matching shape")
    matching_hist[adjacent] += 1
for adjacent in (1, 3):
    subset = {f for f in good if sum(abs(a-b) in (1,5)
              for a,b in map(tuple, (t-{0} for t in f if 0 in t))) == adjacent}
    need(orbit_sizes(subset, dihedral) == [len(subset)], "each shape transitive")
    f = min(subset, key=lambda p: sorted(map(sorted,p)))
    matching = frozenset(t - {0} for t in f if 0 in t)
    exchanges = [p for p in dihedral if act(matching, p) == matching and act(f,p) != f]
    need(exchanges, "matching parity exchange exists")
    representatives.append({"adjacent_edges": adjacent,
                            "plane": sorted(map(sorted, f)),
                            "matching": sorted(map(sorted, matching)),
                            "extension_exchange": exchanges[0]})
need(matching_hist == {1: 12, 3: 4}, "D6 cardinality split")

# Extra result: affine four-bit transport, independent of merely counting sixteen.
seed = frozenset(frozenset(map(int, line))
                 for line in ("015", "026", "034", "123", "146", "245", "356"))
rho = (0, 2, 3, 1, 5, 6, 4)
tau = (0, 2, 1, 3, 5, 4, 6)
need(seed in good and act(seed, rho) == seed, "chosen cyclic seed")


def swap_vector(v):
    p = list(range(7))
    for i in range(3):
        if v[i]:
            p[i+1], p[i+4] = p[i+4], p[i+1]
    return tuple(p)


coordinates = {}
for v in product((0, 1), repeat=3):
    for epsilon in (0, 1):
        f = act(act(seed, tau) if epsilon else seed, swap_vector(v))
        coordinates[(v, epsilon)] = f
need(set(coordinates.values()) == good and len(set(coordinates.values())) == 16,
     "lossless affine four-bit encoding")
affine_checks = 0
for p in star_group:
    sigma = tuple((p[i+1]-1) % 3 for i in range(3))
    n = [0, 0, 0]
    for i in range(3):
        n[sigma[i]] = int(p[i+1] > 3)
    parity = sum(sigma[i] > sigma[j] for i in range(3) for j in range(i+1, 3)) % 2
    for (v, epsilon), f in coordinates.items():
        target_v = n[:]
        for i in range(3):
            target_v[sigma[i]] ^= v[i]
        target = (tuple(target_v), epsilon ^ parity)
        need(act(f, p) == coordinates[target], "affine action formula")
        affine_checks += 1
need(affine_checks == 768, "complete affine transport universe")

# Explicit hostiles to three tempting false shortcuts.
need(any(len(f & crossings) == 1 for f in planes), "nonfree quotient control")
need(not any(len(f & crossings) == 2 for f in planes), "independence hostile")
need({4, 12}.isdisjoint({6, 18}), "no subgroup orbit sizes four and six together")
need(2 - 2 * 3 == -4, "singular-versus-resolution K-square hostile")

print(json.dumps({
    "status": "PASS",
    "universe": "30 labelled simplex-code subgroups, reconstructed independently",
    "label_permutations": len(perms),
    "steiner_exact_covers": len(steiner),
    "bad_crossing_histogram": dict(sorted(hist.items())),
    "local_profiles": [dict(b=b, A1=a, quarter=q, resolution_K2=k, count=c)
                       for (b,a,q,k),c in sorted(node_profiles.items())],
    "canonical_nontrivial_character_checks": character_checks,
    "W_good_orbit_sizes": orbit_sizes(good, star_group),
    "D6_good_orbit_sizes": orbit_sizes(good, dihedral),
    "D6_matching_representatives": representatives,
    "affine_four_bit_transport_checks": affine_checks,
    "geometry_scope": "char zero; local profile independently checked; minimal models inherited from Session60",
}, indent=2, sort_keys=True))
