#!/usr/bin/env python3
"""Exact G21 action on all labelled Fano systems; independent two-way generation."""
import hashlib
import itertools as it
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(msg)


def mask(xs):
    return sum(1 << x for x in xs)


def vertices(s):
    return [i for i in range(7) if s >> i & 1]


def move(s, p):
    return mask(p[i] for i in vertices(s))


def act(F, p):
    return tuple(sorted(move(t, p) for t in F))


def matchings(xs):
    if not xs:
        yield ()
        return
    a = xs[0]
    for b in xs[1:]:
        for rest in matchings([x for x in xs[1:] if x != b]):
            yield tuple(sorted(((a, b),) + rest))


def main():
    triples = [mask(t) for t in it.combinations(range(7), 3)]
    pairs = list(it.combinations(range(7), 2))
    pi = {p: i for i, p in enumerate(pairs)}
    covered = {t: sum(1 << pi[p] for p in it.combinations(vertices(t), 2)) for t in triples}
    complete = (1 << 21) - 1
    exact = set()

    def extend(used, lines):
        if used == complete:
            exact.add(tuple(sorted(lines)))
            return
        missing = next(i for i in range(21) if not (used >> i & 1))
        for t in triples:
            if covered[t] >> missing & 1 and not covered[t] & used:
                extend(used | covered[t], lines + [t])

    extend(0, [])
    linear = set()
    for labels in it.permutations(range(1, 8)):
        lines = tuple(sorted(mask(t) for t in it.combinations(range(7), 3)
                             if labels[t[0]] ^ labels[t[1]] ^ labels[t[2]] == 0))
        linear.add(lines)
    need(exact == linear and len(exact) == 30, "independent exact-cover and vector-space census")
    for F in exact:
        need(len(F) == 7 and all(sum((covered[t] >> i) & 1 for t in F) == 1 for i in range(21)), "literal Fano pairs")
    group = [(a, b, tuple((a * x + b) % 7 for x in range(7))) for a in (1, 2, 4) for b in range(7)]
    gp = {p for _, _, p in group}
    need(len(gp) == 21, "faithful affine group")
    for p, q in it.product(gp, repeat=2):
        need(tuple(p[q[x]] for x in range(7)) in gp, "group composition")
    for F, (_, _, p) in it.product(exact, group):
        need(act(F, p) in exact, "preserved Fano predicate")
    D = {1, 2, 4}
    Dminus = {3, 5, 6}
    fixed = [tuple(sorted(mask((x + b) % 7 for x in d) for b in range(7))) for d in (D, Dminus)]
    need(all(F in exact and all(act(F, p) == F for p in gp) for F in fixed), "two fixed difference systems")
    need(set(fixed) == {F for F in exact if all(act(F, p) == F for p in gp)}, "exact fixed systems")
    todo = set(exact)
    orbits = []
    while todo:
        F = min(todo)
        orb = {act(F, p) for p in gp}
        stabilizer = [[a, b] for a, b, p in group if act(F, p) == F]
        need(len(orb) * len(stabilizer) == 21, "orbit stabilizer")
        profile = sorted({(len(set(fixed[0]) & set(K)), len(set(fixed[1]) & set(K))) for K in orb})
        need(len(profile) == 1, "fixed-system intersections invariant")
        orbits.append({"representative_lines": [vertices(t) for t in F], "size": len(orb),
                       "stabilizer_affine_ab": stabilizer, "intersections_with_fixed": profile[0]})
        todo -= orb
    orbits.sort(key=lambda r: (r["size"], r["representative_lines"]))
    need([r["size"] for r in orbits] == [1, 1, 7, 7, 7, 7], "correct Fano orbit profile")
    addresses = {}
    for F in exact - set(fixed):
        stabilizer = [p for p in gp if act(F, p) == F]
        centres = [c for c in range(7) if all(p[c] == c for p in stabilizer)]
        need(len(stabilizer) == 3 and len(centres) == 1, "unique stabilizer centre")
        profile = (len(set(fixed[0]) & set(F)), len(set(fixed[1]) & set(F)))
        address = (*profile, centres[0])
        need(address not in addresses, "four-copy address injective")
        addresses[address] = F
        for _, _, p in group:
            target = act(F, p)
            st = [q for q in gp if act(target, q) == target]
            need(all(q[p[centres[0]]] == p[centres[0]] for q in st), "centre address equivariance")
    types = {(0, 1), (1, 0), (1, 3), (3, 1)}
    need(set(addresses) == {(a, b, c) for a, b in types for c in range(7)}, "four copies of natural seven-set")
    element_fix = []
    for a, b, p in group:
        count = sum(act(F, p) == F for F in exact)
        expected = 30 if (a, b) == (1, 0) else 2 if a == 1 else 6
        need(count == expected, "element fixed-point formula")
        element_fix.append({"a": a, "b": b, "fixed_systems": count})
    need(sum(r["fixed_systems"] for r in element_fix) == 21 * len(orbits), "Burnside exact count")
    # C3 fixed systems: three invariant pairings through the fixed point, two completions each.
    sigma = tuple(2 * x % 7 for x in range(7))
    invariant = []
    for M in matchings(list(range(1, 7))):
        image = tuple(sorted(tuple(sorted((sigma[a], sigma[b]))) for a, b in M))
        if image != M:
            continue
        pencil = {mask((0, a, b)) for a, b in M}
        completions = [F for F in exact if pencil <= set(F)]
        need(len(completions) == 2 and all(act(F, sigma) == F for F in completions), "two invariant completions")
        invariant.append({"pairs": M, "completions": [[vertices(t) for t in F] for F in sorted(completions)]})
    need(len(invariant) == 3, "three C3-invariant pairings")
    # All five cyclic orbits on triples; difference-set condition singles out two.
    todo = set(triples)
    cyclic = []
    triple_affine_sizes = []
    affine_todo = set(triples)
    while affine_todo:
        t = min(affine_todo)
        orb = {move(t, p) for p in gp}
        triple_affine_sizes.append(len(orb))
        affine_todo -= orb
    need(sorted(triple_affine_sizes) == [7, 7, 21], "triple orbit profile differs")
    while todo:
        t = min(todo)
        orb = {move(t, p) for a, _, p in group if a == 1}
        B = vertices(t)
        differences = Counter((x - y) % 7 for x in B for y in B if x != y)
        good = differences == Counter(range(1, 7))
        need((tuple(sorted(orb)) in exact) == good, "cyclic Fano iff difference set")
        cyclic.append({"base_triple": B, "orbit_size": len(orb), "differences": dict(sorted(differences.items())), "Fano": good})
        todo -= orb
    need(len(cyclic) == 5 and sum(r["Fano"] for r in cyclic) == 2, "exact two cyclic Fano systems")
    # The remaining28 are not invariant under the full168 automorphisms of one fixed Fano system.
    aut = [p for p in it.permutations(range(7)) if act(fixed[0], p) == fixed[0]]
    need(len(aut) == 168, "full Fano automorphism size")
    p = next(p for p in aut if act(fixed[1], p) != fixed[1])
    image = act(fixed[1], p)
    need(image in exact and image not in fixed, "two-deletion not invariant under full Fano group")
    result = {"status": "FINITE-EXACT; structural proof RESERVED pending audit", "universe": "all30 labelled Fano systems on fixed point setF7; all21 affine maps; all35 triples",
              "orbits": orbits, "element_fixed_counts": element_fix,
              "fixed_systems": [[vertices(t) for t in F] for F in fixed],
              "four_copy_addresses": [{"type": [a, b], "centre": c, "lines": [vertices(t) for t in addresses[a, b, c]]}
                                      for a, b, c in sorted(addresses)],
              "C3_invariant_pencils": invariant, "cyclic_triple_orbits": cyclic,
              "triple_affine_orbit_sizes": sorted(triple_affine_sizes),
              "larger_group_hostile": {"permutation_images_0_to_6": p, "fixed_first_system": True,
                                        "second_system_image": [vertices(t) for t in image]},
              "checks": CHECKS, "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
