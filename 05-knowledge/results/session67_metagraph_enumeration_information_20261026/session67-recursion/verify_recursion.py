#!/usr/bin/env python3
"""Standalone exact checks for rooted augmentation and SCC interval insertion.

Only Python's standard library is used; no repository code or executable is
imported. All checks remain enabled with python -O. Default max order is 7.
Representation: rows[v] has bit u iff v -> u. Printed codes use pairs i<j
in lexicographic order, bit one iff i -> j. Canonical labels sort outdegrees
and then minimize this code within equal-degree blocks.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import factorial, gcd, prod
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


@lru_cache(None)
def pairs(n):
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


@lru_cache(None)
def decode(n, code):
    rows = [0] * n
    for k, (i, j) in enumerate(pairs(n)):
        if code >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def encode(rows, order=None):
    if order is None:
        order = tuple(range(len(rows)))
    return sum(((rows[order[i]] >> order[j]) & 1) << k
               for k, (i, j) in enumerate(pairs(len(rows))))


def degree_orders(rows, root=None):
    buckets = defaultdict(list)
    for v, row in enumerate(rows):
        if v != root:
            key = (row.bit_count(),) if root is None else (
                (rows[root] >> v) & 1, row.bit_count())
            buckets[key].append(v)
    blocks = [buckets[k] for k in sorted(buckets)]
    for choices in product(*(permutations(b) for b in blocks)):
        yield (() if root is None else (root,)) + tuple(
            v for choice in choices for v in choice)


@lru_cache(None)
def canonical(rows):
    return min((encode(rows, order), order) for order in degree_orders(rows))


@lru_cache(None)
def root_key(rows, root):
    return min(encode(rows, order) for order in degree_orders(rows, root))


@lru_cache(None)
def automorphisms(rows):
    n = len(rows)
    baseline = tuple(sorted(range(n), key=lambda v: (rows[v].bit_count(), v)))
    base_code = encode(rows, baseline)
    result = []
    for target in degree_orders(rows):
        if encode(rows, target) == base_code:
            p = [0] * n
            for u, v in zip(baseline, target):
                p[u] = v
            result.append(tuple(p))
    return tuple(result)


def move_mask(mask, p):
    return sum(((mask >> v) & 1) << p[v] for v in range(len(p)))


@lru_cache(None)
def subset_orbits(rows):
    group = automorphisms(rows)
    return tuple(mask for mask in range(1 << len(rows))
                 if mask == min(move_mask(mask, p) for p in group))


def extension(rows, mask):
    # mask is the outgoing neighborhood of the new, final vertex.
    n = len(rows)
    return tuple(row | ((1 - ((mask >> v) & 1)) << n)
                 for v, row in enumerate(rows)) + (mask,)


def induced(rows, vertices):
    return tuple(sum(((rows[v] >> u) & 1) << j
                     for j, u in enumerate(vertices)) for v in vertices)


@lru_cache(None)
def sccs(rows):
    n = len(rows)
    reach = [rows[v] | (1 << v) for v in range(n)]
    for k in range(n):
        for v in range(n):
            if reach[v] >> k & 1:
                reach[v] |= reach[k]
    unseen = set(range(n))
    components = []
    while unseen:
        v = min(unseen)
        component = tuple(u for u in sorted(unseen)
                          if (reach[v] >> u & 1) and (reach[u] >> v & 1))
        unseen.difference_update(component)
        components.append(component)
    return tuple(sorted(components, key=lambda c: -sum(
        (rows[c[0]] >> d[0]) & 1 for d in components if c != d)))


def predicted_sccs(rows, mask):
    cs = sccs(rows)
    n, k = len(rows), len(cs)
    left = next((i for i, c in enumerate(cs)
                 if any(mask >> v & 1 for v in c)), k)
    right = next((i for i in range(k - 1, -1, -1)
                  if any(not (mask >> v & 1) for v in cs[i])), -1)
    if left > right:
        need(left == right + 1, "acyclic insertion is a single SCC cut")
        return cs[:left] + ((n,),) + cs[left:]
    merged = tuple(sorted((n,) + tuple(v for c in cs[left:right + 1] for v in c)))
    return cs[:left] + (merged,) + cs[right + 1:]


def cycle_type(p):
    unseen = set(range(len(p)))
    lengths = []
    while unseen:
        v = min(unseen)
        length = 0
        while v in unseen:
            unseen.remove(v)
            length += 1
            v = p[v]
        lengths.append(length)
    return tuple(sorted(lengths))


@lru_cache(None)
def core_order(rows):
    n = len(rows)
    homogeneous_pairs = sum(
        all(((rows[u] >> w) & 1) == ((rows[v] >> w) & 1)
            for w in range(n) if w not in (u, v))
        for u, v in pairs(n))
    return n - homogeneous_pairs


def odd_partitions(n, least=1):
    if not n:
        yield ()
    for k in range(least, n + 1, 2):
        for tail in odd_partitions(n - k, k):
            yield (k,) + tail


def burnside(n, rooted=False):
    total = Fraction()
    for part in odd_partitions(n):
        multiplicities = Counter(part)
        denominator = prod(k ** m * factorial(m) for k, m in multiplicities.items())
        exponent = sum((k - 1) // 2 for k in part) + sum(
            gcd(part[i], part[j]) for i in range(len(part)) for j in range(i + 1, len(part)))
        total += Fraction((1 << exponent) * (multiplicities[1] if rooted else 1), denominator)
    need(total.denominator == 1, "Burnside integrality")
    return total.numerator


def strong_orbit_prediction(rows):
    cs = sccs(rows)
    if not cs:
        return 1
    es = [len(subset_orbits(induced(rows, c))) for c in cs]
    if len(cs) == 1:
        return es[0] - 2
    return (es[0] - 1) * (es[-1] - 1) * prod(es[1:-1])


def verify(max_order, exhaustive_order):
    levels = [{0: ()}]
    summaries = []
    extension_stats = {}
    all_scc_trials = 0
    for n in range(1, max_order + 1):
        accepted = Counter()
        root_presentations = Counter()
        reciprocal_root_mass = defaultdict(Fraction)
        orbit_size_mass = defaultdict(Fraction)
        raw_symmetry_mass = defaultdict(Fraction)
        all_children = set()
        strong_rooted_extensions = 0
        for parent_code, parent in levels[-1].items():
            gp = automorphisms(parent)
            reps = subset_orbits(parent)
            cycle_burnside = Fraction(sum(1 << len(cycle_type(p)) for p in gp), len(gp))
            need(cycle_burnside == len(reps), "subset orbit Burnside")
            child_set = set()
            child_aut_profile = Counter()
            core_mass = defaultdict(Fraction)
            strong_here = 0
            for mask in reps:
                child = extension(parent, mask)
                code, order = canonical(child)
                canon = decode(n, code)
                root = order.index(n - 1)
                group = automorphisms(canon)
                orbit = {p[root] for p in group}
                root_keys = tuple(root_key(canon, v) for v in range(n))
                r = len(set(root_keys))
                stabilizer = sum(move_mask(mask, p) == mask for p in gp)
                need(len(group) == len(orbit) * stabilizer, "root stabilizer compatibility")
                root_presentations[code] += 1
                reciprocal_root_mass[code] += Fraction(1, r)
                orbit_size_mass[code] += Fraction(len(orbit), n)
                raw_symmetry_mass[code] += Fraction(len(gp), stabilizer) * Fraction(len(group), n * len(gp))
                all_children.add(code)
                child_set.add(code)
                child_aut_profile[len(group)] += 1
                core_mass[core_order(canon)] += Fraction(len(group), n * stabilizer)
                if root_keys[root] == min(root_keys):
                    accepted[code] += 1
                if len(sccs(child)) == 1:
                    strong_here += 1
            need(strong_here == strong_orbit_prediction(parent), "factorized strong rooted extension count")
            strong_rooted_extensions += strong_here
            extension_stats[(n - 1, parent_code)] = {
                "distinct_children": len(child_set),
                "child_codes": sorted(child_set),
                "rooted_extension_aut_profile": dict(sorted(child_aut_profile.items())),
                "strong_rooted_extensions": strong_here,
                "symmetry_weighted_core_transition": {k: str(v) for k, v in sorted(core_mass.items())},
            }
        need(set(accepted) == all_children, "canonical-root augmentation covers every child")
        need(set(accepted.values()) == {1}, "canonical-root augmentation has no duplicates")
        bank = {code: decode(n, code) for code in sorted(accepted)}
        need(len(bank) == burnside(n), "class count independently matches odd-partition Burnside")
        rooted_count = 0
        strong_count = 0
        strong_rooted_count = 0
        for code, rows in bank.items():
            r = len({root_key(rows, v) for v in range(n)})
            need(root_presentations[code] == r, "individual root-forgetting fibre size")
            need(reciprocal_root_mass[code] == orbit_size_mass[code] == raw_symmetry_mass[code] == 1,
                 "each child has exact mass one in all three repairs")
            rooted_count += r
            if len(sccs(rows)) == 1:
                strong_count += 1
                strong_rooted_count += r
        need(rooted_count == burnside(n, rooted=True), "rooted Burnside comparison")
        need(strong_rooted_count == strong_rooted_extensions, "strong rooted count two routes")
        need(sum(factorial(n) // len(automorphisms(rows)) for rows in bank.values()) == 1 << len(pairs(n)),
             "full labelled orbit mass")
        levels.append(bank)
        if n <= exhaustive_order:
            exhaustive = {canonical(decode(n, code))[0] for code in range(1 << len(pairs(n)))}
            need(exhaustive == set(bank), "independent exhaustive labelled class set")
        summaries.append({"n": n, "classes": len(bank), "rooted": rooted_count,
                          "strong": strong_count, "strong_rooted": strong_rooted_count,
                          "naive_divide_n": str(Fraction(rooted_count, n))})
    # Exhaust all labelled parents through the specified bound and every insertion.
    for n in range(exhaustive_order + 1):
        for code in range(1 << len(pairs(n))):
            rows = decode(n, code)
            for mask in range(1 << n):
                need(sccs(extension(rows, mask)) == predicted_sccs(rows, mask), "exact SCC interval update")
                all_scc_trials += 1
    a = [1] + [s["classes"] for s in summaries]
    s = [0] + [row["strong"] for row in summaries]
    b = [0] + [row["rooted"] for row in summaries]
    q = [0] + [row["strong_rooted"] for row in summaries]
    for n in range(1, max_order + 1):
        need(a[n] == sum(s[k] * a[n-k] for k in range(1, n+1)), "SCC sequence identity A=1/(1-S)")
        need(b[n] == sum(a[i] * q[j] * a[n-i-j] for i in range(n) for j in range(1, n-i+1)),
             "rooted SCC identity B=A^2 Q")
    # Minimal naive-card hostile: both vertex orbits of T2 delete to K1.
    t2 = next(iter(levels[2].values()))
    deletion_cards = [canonical(induced(t2, tuple(u for u in range(2) if u != v)))[0] for v in range(2)]
    need(len(set(deletion_cards)) == 1 and len({root_key(t2, v) for v in range(2)}) == 2,
         "unrooted minimum-card selector is not a root-orbit selector")
    # Find equal complete vertex-action cycle indices with unequal one-step responses.
    hostile = None
    for n in range(1, max_order):
        by_index = defaultdict(list)
        for code, rows in levels[n].items():
            counts = Counter(cycle_type(p) for p in automorphisms(rows))
            cycle_index = tuple(sorted((part, str(Fraction(m, len(automorphisms(rows))))) for part, m in counts.items()))
            by_index[cycle_index].append(code)
        for cycle_index, codes in by_index.items():
            for i, left in enumerate(codes):
                for right in codes[i+1:]:
                    ls, rs = extension_stats[(n, left)], extension_stats[(n, right)]
                    if ls["distinct_children"] != rs["distinct_children"]:
                        hostile = {"n": n, "parent_codes": [left, right],
                                   "parent_rows": [levels[n][left], levels[n][right]],
                                   "cycle_index": cycle_index, "responses": [ls, rs]}
                        break
                if hostile:
                    break
            if hostile:
                break
        if hostile:
            break
    need(hostile is not None, "cycle-index growth hostile found")
    # Stronger family hostile: same core order, SCC-size sequence and cycle index.
    family_hostile = None
    for n in range(1, max_order):
        families = defaultdict(list)
        for code, rows in levels[n].items():
            group = automorphisms(rows)
            ci = tuple(sorted((part, str(Fraction(m, len(group))))
                              for part, m in Counter(cycle_type(p) for p in group).items()))
            key = (core_order(rows), tuple(map(len, sccs(rows))), ci)
            families[key].append(code)
        for key, codes in families.items():
            for i, left in enumerate(codes):
                for right in codes[i+1:]:
                    ls, rs = extension_stats[(n, left)], extension_stats[(n, right)]
                    if ls["symmetry_weighted_core_transition"] != rs["symmetry_weighted_core_transition"]:
                        family_hostile = {"n": n, "parent_codes": [left, right],
                                          "parent_rows": [levels[n][left], levels[n][right]],
                                          "core_order": key[0], "scc_sizes": key[1], "cycle_index": key[2],
                                          "responses": [ls, rs]}
                        break
                if family_hostile:
                    break
            if family_hostile:
                break
        if family_hostile:
            break
    return {"status": "PASS", "max_order": max_order,
            "exhaustive_labelled_order": exhaustive_order, "scc_insertion_trials": all_scc_trials,
            "table": summaries, "minimal_equal_cycle_index_unequal_distinct_children": hostile,
            "minimal_equal_core_scc_cycle_index_unequal_weighted_core_transition": family_hostile,
            "minimum_unrooted_card_selector_failure_order": 2,
            "minimum_naive_root_division_failure_order": 3}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=7)
    parser.add_argument("--exhaustive-order", type=int, default=5)
    args = parser.parse_args()
    need(args.max_order >= 5 and 0 <= args.exhaustive_order <= min(6, args.max_order), "supported verification bounds")
    print(json.dumps(verify(args.max_order, args.exhaustive_order), sort_keys=True, indent=2))
