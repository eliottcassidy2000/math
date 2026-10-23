#!/usr/bin/env python3
"""Independent exhaustive Paley-7 audit. No producer imports; no Git access."""
import hashlib
import json
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path

V = tuple(range(7))
QR = {1, 2, 4}
ID = V


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def arc(x, y):
    return (y - x) % 7 in QR


def compose(p, q):
    return tuple(p[q[x]] for x in V)


def inverse(p):
    return tuple(p.index(x) for x in V)


def order(p):
    q, n = p, 1
    while q != ID:
        q, n = compose(p, q), n + 1
    return n


def affine(a, b):
    return tuple((a * x + b) % 7 for x in V)


def presentation(path):
    return tuple(int(arc(path[i], path[j])) for i in V for j in range(i + 1, 7))


def is_blue(path):
    return all(arc(path[i], path[j]) == arc(path[6 - j], path[6 - i])
               for i in V for j in range(i + 1, 7))


all_perms = tuple(permutations(V))
autos = {p for p in all_perms
         if all(arc(p[x], p[y]) == arc(x, y) for x in V for y in V if x != y)}
antis = {p for p in all_perms
         if all(arc(p[x], p[y]) == arc(y, x) for x in V for y in V if x != y)}
full_affine = {affine(a, b) for a in range(1, 7) for b in V}
need(autos == {affine(a, b) for a in QR for b in V}, "automorphisms")
need(antis == {affine(a, b) for a in {3, 5, 6} for b in V}, "anti-automorphisms")
need(autos | antis == full_affine, "signed group is AGL(1,7)")
need({order(p) for p in full_affine} == {1, 2, 3, 6, 7}, "element orders")
anti_involutions = {p for p in antis if order(p) == 2}
need(anti_involutions == {affine(6, b) for b in V}, "seven anti-involutions")

# Independently establish the affine coset action, including its faithful kernel.
stabilizer = {affine(a, 0) for a in range(1, 7)}
cosets = {x: {compose(affine(1, x), h) for h in stabilizer} for x in V}
need(len(set.union(*cosets.values())) == 42, "seven disjoint six-element cosets")
for g in full_affine:
    for x in V:
        need({compose(g, h) for h in cosets[x]} == cosets[g[x]], "left coset action")
need({g for g in full_affine if all(g[x] == x for x in V)} == {ID}, "coset kernel")

# Conjugation on C7 fixes its identity and kills every translation.
translations = {affine(1, b) for b in V}
conjugation_kernel = set()
for a in range(1, 7):
    for b in V:
        g = affine(a, b)
        images = [compose(compose(g, affine(1, k)), inverse(g)) for k in V]
        need(images == [affine(1, a * k % 7) for k in V], "conjugation multiplier")
        if images == [affine(1, k) for k in V]:
            conjugation_kernel.add(g)
need(conjugation_kernel == translations, "conjugation is not affine point action")

paths = {p for p in all_perms if all(arc(p[i], p[i + 1]) for i in range(6))}
need(len(paths) == 189, "Hamiltonian path count")

# A separate subset dynamic program checks 189 without permutation filtering.
dp = {(1 << v, v): 1 for v in V}
for mask in range(1, 1 << 7):
    for v in V:
        count = dp.get((mask, v), 0)
        for w in V:
            if count and not mask & (1 << w) and arc(v, w):
                key = (mask | (1 << w), w)
                dp[key] = dp.get(key, 0) + count
endpoint_counts = [dp.get((127, v), 0) for v in V]
need(endpoint_counts == [27] * 7, "independent path dynamic program")

presentation_counts = Counter(presentation(p) for p in paths)
blue_paths = {p for p in paths if is_blue(p)}
blue_presentations = {presentation(p) for p in blue_paths}
need(len(presentation_counts) == 9, "nine presentations")
need(set(presentation_counts.values()) == {21}, "free automorphism action")
need(len(blue_paths) == 63 and len(blue_presentations) == 3, "blue fiber")
need(all(sum(tuple(g[p[6-i]] for i in V) == p for g in antis) == int(p in blue_paths)
         for p in paths), "unique anti-involution for each blue path")
fixed_by_anti = {g: sum(tuple(g[p[6-i]] for i in V) == p for p in paths) for g in antis}
need(Counter(fixed_by_anti.values()) == Counter({0: 14, 9: 7}), "signed fixed points")

unseen = set(paths)
signed_orbits = []
while unseen:
    p = min(unseen)
    orbit = {tuple(g[p[i]] for i in V) for g in autos}
    orbit |= {tuple(g[p[6-i]] for i in V) for g in antis}
    need(orbit <= paths, "signed path action preserves directed paths")
    signed_orbits.append(len(orbit))
    unseen -= orbit
need(sorted(signed_orbits) == [21, 21, 21, 42, 42, 42], "signed path orbits")

# Separate direct powerset check of the ordinary/signed regular-subset claim.
regular_subsets = {}
for label, group in (("ordinary", autos), ("signed", full_affine)):
    by_size = Counter()
    for mask in range(128):
        subset = {v for v in V if mask & (1 << v)}
        if sum({g[x] for x in subset} == subset for g in group) == 1:
            by_size[len(subset)] += 1
    regular_subsets[label] = dict(sorted(by_size.items()))
need(regular_subsets == {"ordinary": {2: 21, 3: 21, 4: 21, 5: 21}, "signed": {}},
     "regular vertex subsets")

# Optional independent odd-cycle inventory checks the path-count identity.
cycles_by_length = {}
cycle_supports = []
for length in (3, 5, 7):
    cycles = []
    for support in combinations(V, length):
        for tail in permutations(support[1:]):
            cycle = (support[0],) + tail
            if all(arc(cycle[i], cycle[(i + 1) % length]) for i in range(length)):
                cycles.append(cycle)
                cycle_supports.append(set(cycle))
    cycles_by_length[str(length)] = len(cycles)
disjoint_odd_cycle_pairs = sum(not a & b for a, b in combinations(cycle_supports, 2))
need(cycles_by_length == {"3": 14, "5": 42, "7": 24}, "directed odd cycles")
need(disjoint_odd_cycle_pairs == 7, "disjoint odd-cycle pairs")
need(1 + 2 * sum(cycles_by_length.values()) + 4 * disjoint_odd_cycle_pairs == 189,
     "odd-cycle path-count control")

root = Path(__file__).resolve().parent.parent
sources = {}
for name in ("PALEY_MONSTER_BRIDGE.md", "check.py", "exact.json", "optimized.json"):
    path = root / "paley" / name
    if path.exists():
        sources[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
producer_path = root / "paley" / "exact.json"
if producer_path.exists():
    producer = json.loads(producer_path.read_text())
    need((producer["H"], producer["Aut"], producer["h"], producer["b"])
         == (len(paths), len(autos), len(presentation_counts), len(blue_presentations)),
         "producer census cross-check")
    need({presentation(tuple(p)) for p in producer["normalized_presentations"]}
         == set(presentation_counts), "producer full presentation census")
    need({presentation(tuple(p)) for p in producer["blue_presentations"]}
         == blue_presentations, "producer blue presentation census")
receipt = {
    "status": "PASS",
    "imports_producer_code": False,
    "P7": {
        "quadratic_residues": sorted(QR), "hamiltonian_paths": len(paths),
        "per_endpoint_counts": endpoint_counts,
        "automorphisms": len(autos), "anti_automorphisms": len(antis),
        "affine_signed_group_order": len(full_affine),
        "anti_involutions": len(anti_involutions),
        "presentations": len(presentation_counts),
        "blue_paths": len(blue_paths), "blue_presentations": len(blue_presentations),
        "blue_presentation_full_edge_bitstrings": sorted(''.join(map(str, x)) for x in blue_presentations),
        "black_presentations": len(presentation_counts) - len(blue_presentations),
        "paired_black_degree": (len(presentation_counts) - len(blue_presentations)) // 2,
        "signed_path_orbit_sizes": sorted(signed_orbits),
        "anti_fixed_path_distribution": dict(sorted(Counter(fixed_by_anti.values()).items())),
        "directed_odd_cycles": cycles_by_length,
        "disjoint_odd_cycle_pairs": disjoint_odd_cycle_pairs,
        "regular_vertex_subsets_by_size": regular_subsets,
    },
    "action_distinction": {
        "coset_degree": 7, "coset_point_stabilizer_order": 6,
        "faithful_affine_quotient_order": 42,
        "conjugation_image_order": 6,
        "conjugation_kernel_in_affine_quotient_order": 7,
        "conjugation_fixes_identity": True,
    },
    "producer_files_present": sources,
}
print(json.dumps(receipt, indent=2, sort_keys=True))
