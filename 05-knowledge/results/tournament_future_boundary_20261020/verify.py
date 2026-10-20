#!/usr/bin/env python3
"""Exact future-response Burnside DP and marked tribonacci defect checks.

Standard library only. All checks remain active with python -O.
Run: python3 verify.py --max-order 100 [--large-order 200]
The optional large run extends the computed sequence, not the independent
external-data control (which ends at 76).
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
from pathlib import Path
import argparse
import json

CHECKS = Counter()


def require(condition, label):
    CHECKS[label] += 1
    if not condition:
        raise RuntimeError(label)


def profile_dp(bound, *, core=False, no_fixed=False, ascending=False):
    """Return all coefficients through bound, aggregating future responses.

    For core=True use p_l=x^l/(1+x^l) in Z_T; otherwise p_l=x^l.
    no_fixed=True omits l=1. State degree includes substitution excess;
    the response profile counts only actual permutation cycles.
    """
    lengths = list(range(3 if no_fixed else 1, bound + 1, 2))
    if not ascending:
        lengths.reverse()
    states = {(0, (0,) * len(lengths)): Q(1)}
    rows = []
    for index, ell in enumerate(lengths):
        future = lengths[index + 1:]
        new = defaultdict(Q)
        transitions = 0
        for (degree, response), weight in states.items():
            admissible = [j for j in lengths[index:] if j <= bound - degree]
            old = dict(zip(admissible, response))
            charge = old.get(ell, 0)
            for k in range((bound - degree) // ell + 1):
                exponent = k * (ell - 1) // 2 + ell * k * (k - 1) // 2 + k * charge
                factor = weight * Q(1 << exponent, ell ** k * factorial(k))
                excesses = range((bound - degree) // ell - k + 1) if core and k else (0,)
                for excess in excesses:
                    new_degree = degree + ell * (k + excess)
                    profile = tuple(old[j] + k * gcd(ell, j)
                                    for j in future if j <= bound - new_degree)
                    sign_weight = (-1) ** excess * comb(k + excess - 1, excess) if k else 1
                    new[new_degree, profile] += factor * sign_weight
                    transitions += 1
        states = {key: value for key, value in new.items() if value}
        rows.append({"cycle_length": ell, "states": len(states), "transitions": transitions})
    answer = [Q(0)] * (bound + 1)
    for (degree, profile), weight in states.items():
        require(not profile, "terminal profile empty")
        answer[degree] += weight
    return answer, {"peak_states": max([1] + [row["states"] for row in rows]),
                    "transitions": sum(row["transitions"] for row in rows),
                    "rows": rows}


def odd_partitions(n, minimum=1):
    if not n:
        yield ()
    for part in range(minimum, n + 1, 2):
        for rest in odd_partitions(n - part, part):
            yield (part,) + rest


def cycle_weight(parts):
    exponent = sum((a - 1) // 2 for a in parts)
    exponent += sum(gcd(a, b) for a, b in combinations(parts, 2))
    centralizer = 1
    for a, multiplicity in Counter(parts).items():
        centralizer *= a ** multiplicity * factorial(multiplicity)
    return Q(1 << exponent, centralizer)


def direct_series(bound, core=False, no_fixed=False):
    """Independent partition sum, with direct polynomial substitution."""
    result = [Q(0)] * (bound + 1)
    for mass in range(bound + 1):
        for parts in odd_partitions(mass, 3 if no_fixed else 1):
            polynomial = [Q(1)] + [Q(0)] * (bound - mass)
            if core:
                for ell, k in Counter(parts).items():
                    new = [Q(0)] * len(polynomial)
                    for offset, value in enumerate(polynomial):
                        for excess in range((bound - mass - offset) // ell + 1):
                            new[offset + ell * excess] += (
                                value * (-1) ** excess * comb(k + excess - 1, excess))
                    polynomial = new
            weight = cycle_weight(parts)
            for offset, coefficient in enumerate(polynomial):
                result[mass + offset] += weight * coefficient
    return result


def permutation_for(parts):
    perm = []
    cursor = 0
    for ell in parts:
        perm.extend(list(range(cursor + 1, cursor + ell)) + [cursor])
        cursor += ell
    return perm


def fixed_core_count(parts):
    """Direct orientation-orbit enumeration, not the gcd exponent formula."""
    perm = permutation_for(parts)
    n = len(perm)
    unseen = set(combinations(range(n), 2))
    orbits = []
    while unseen:
        a, b = min(unseen)
        seed = (a, b)
        orbit = []
        while True:
            orbit.append((a, b))
            unseen.remove(tuple(sorted((a, b))))
            a, b = perm[a], perm[b]
            if (a, b) == seed:
                break
            require((a, b) != seed[::-1], "odd permutation never reverses pair")
        orbits.append(orbit)
    count = 0
    for mask in range(1 << len(orbits)):
        adjacency = [0] * n
        for index, orbit in enumerate(orbits):
            for a, b in orbit:
                if (mask >> index) & 1:
                    a, b = b, a
                adjacency[a] |= 1 << b
        if core_order(adjacency) == n:
            count += 1
    centralizer = 1
    for ell, k in Counter(parts).items():
        centralizer *= ell ** k * factorial(k)
    return Q(count, centralizer)


def core_order(adjacency):
    n = len(adjacency)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            a = parent[a]
        return a

    for a, b in combinations(range(n), 2):
        if not (adjacency[a] ^ adjacency[b]) & ~((1 << a) | (1 << b)):
            parent[find(a)] = find(b)
    return len({find(a) for a in range(n)})


def compositions(n, previous_one=False):
    if n == 0:
        yield ()
    for first in range(1, n + 1):
        if first == 1 and previous_one:
            continue
        for rest in compositions(n - first, first == 1):
            yield (first,) + rest


def decode(parts):
    cursor = sum(parts)
    path = []
    for length in parts:
        path.extend(range(cursor - length, cursor))
        cursor -= length
    return tuple(path)


def tribonacci_arc(a, b):
    return b == a + 1 or b <= a - 2


def line_child(parts):
    path = decode(parts)
    n = len(path)
    edges = {tuple(sorted(edge)) for edge in zip(path, path[1:])}
    adjacency = [0] * n
    for a, b in combinations(range(n), 2):
        # Converse T*, then toggle the chosen Hamiltonian-path arcs.
        if (b == a + 1) ^ ((a, b) in edges):
            adjacency[b] |= 1 << a
        else:
            adjacency[a] |= 1 << b
    return adjacency


def cap_expansion(parts):
    capped = tuple(min(3, length) for length in parts)
    prototype = line_child(capped)
    # Vertices have increasing numeric order, opposite to the run order.
    fibers = []
    pointer = 0
    for length in reversed(parts):
        if length <= 2:
            fibers.extend([pointer + i] for i in range(length))
        else:
            fibers.extend(([pointer], list(range(pointer + 1, pointer + length - 1)),
                           [pointer + length - 1]))
        pointer += length
    expanded = [0] * sum(parts)
    for i, fiber in enumerate(fibers):
        for a, b in combinations(fiber, 2):
            expanded[a] |= 1 << b
        for j in range(i + 1, len(fibers)):
            for a in fiber:
                for b in fibers[j]:
                    if (prototype[i] >> j) & 1:
                        expanded[a] |= 1 << b
                    else:
                        expanded[b] |= 1 << a
    return capped, prototype, expanded


def minimum_feedback(adjacency):
    """Independent subset DP: choose the last vertex of each total order."""
    n = len(adjacency)
    best = [0] * (1 << n)
    for mask in range(1, 1 << n):
        candidates = []
        bits = mask
        while bits:
            bit = bits & -bits
            bits -= bit
            vertex = bit.bit_length() - 1
            rest = mask ^ bit
            candidates.append(best[rest] + (adjacency[vertex] & rest).bit_count())
        best[mask] = min(candidates)
    return best[-1]


def path_checks():
    totals = [1]
    census = []
    for n in range(1, 15):
        codes = list(compositions(n))
        paths = {decode(code) for code in codes}
        require(len(paths) == len(codes), "run code injective")
        for code in codes:
            path = decode(code)
            require(all(tribonacci_arc(a, b) for a, b in zip(path, path[1:])),
                    "decoded path admissible")
            require((path[0], path[-1]) == (n - code[0], code[-1] - 1),
                    "source and target ports")
            require(tuple(n - 1 - x for x in reversed(path)) == decode(code[::-1]),
                    "reflection reverses run code")
            capped, prototype, expanded = cap_expansion(code)
            actual = line_child(code)
            require(actual == expanded, "exact labeled line-child reconstruction")
            require(core_order(actual) == core_order(prototype), "canonical core survives stretching")
            require(core_order(actual) <= 3 * len(code) - 2, "sharp run core bound")
            if len(code) == 2 and min(code) >= 2:
                require(core_order(actual) == 3, "two-run triangle core")
                require(all((actual[0] >> j) & 1 for j in range(1, n-1))
                        and all((actual[j] >> (n-1)) & 1 for j in range(1, n-1))
                        and actual[n-1] & 1, "two-run marked triangle orientation")
            if min(code) >= 2:
                blocks = []
                cursor = 0
                for length in reversed(code):
                    blocks.append((cursor, cursor+length-1))
                    cursor += length
                used_edges = set()
                order = list(range(n))
                for i in range(len(blocks)-1):
                    left, right = blocks[i]
                    next_left, next_right = blocks[i+1]
                    triangle = ((left, right), (right, next_right), (next_right, left))
                    require(all((actual[a] >> b) & 1 for a, b in triangle),
                            "feedback lower-bound triangle")
                    edge_set = {tuple(sorted(e)) for e in triangle}
                    require(not edge_set & used_edges, "feedback triangles edge disjoint")
                    used_edges |= edge_set
                    order[right], order[next_left] = order[next_left], order[right]
                cost = sum((actual[b] >> a) & 1 for a, b in combinations(order, 2))
                require(cost == len(code)-1, "feedback upper-bound explicit order")
                if n <= 10:
                    require(minimum_feedback(actual) == len(code)-1, "independent feedback subset DP")
        if n <= 8:
            brute = {path for path in permutations(range(n))
                     if all(tribonacci_arc(a, b) for a, b in zip(path, path[1:]))}
            require(paths == brute, "independent full permutation path enumeration")
        totals.append(len(codes))
        if n >= 3:
            require(totals[n] == sum(totals[n-3:n]), "inherited tribonacci recurrence")
        blue = sum(code == code[::-1] for code in codes)
        require(blue == totals[n // 2 + 1], "inherited blue law")
        profile = Counter((len(code), code.count(1)) for code in codes)
        for r in range(1, n + 1):
            for s in range(r + 1):
                large = r - s
                expected = 0
                if large and n - 2*r + s >= 0 and s <= large + 1:
                    expected = comb(large + 1, s) * comb(n-r-1, large-1)
                if large == 0 and r == s == n == 1:
                    expected = 1
                require(profile[r, s] == expected, "bigraded run count")
        census.append({"n": n, "paths": len(codes), "blue": blue})
    type_rows = []
    for r in range(1, 9):
        codes = [c for c in product((1, 2, 3), repeat=r)
                 if not any(a == b == 1 for a, b in zip(c, c[1:]))]
        maximum = max(core_order(line_child(code)) for code in codes)
        require(maximum == 3*r-2, "sharp run bound attained")
        if r >= 3:
            require(core_order(line_child((3,) * r)) == 3*r-2, "all-three extremal carrier")
        type_rows.append({"runs": r, "types": len(codes), "maximum_core_order": maximum})
    return {"census": census, "type_rows": type_rows}


def graph_checks():
    all_edges = set(combinations(range(4), 2))
    classes = Counter()
    for edges in combinations(all_edges, 3):
        degree = sorted(sum(v in edge for edge in edges) for v in range(4))
        classes[tuple(degree)] += 1
    require(classes == {(1, 1, 2, 2): 12, (1, 1, 1, 3): 4, (0, 2, 2, 2): 4},
            "complete four-vertex three-edge universe")
    star = {(0, 1), (0, 2), (0, 3)}
    triangle = {(1, 2), (1, 3), (2, 3)}
    require(all_edges - star == triangle, "star triangle complement")
    # Historical clique-quotient hostile at fugacity 2.
    require(Q(1 + 3*2 + 2**2) != Q((1 + 2*2)**2, 1+2),
            "clique separator needs extension vector")
    # A claw inside the odd-interval conflict graph of T*_9.
    intervals = [set(range(9)), set(range(3)), set(range(3, 6)), set(range(6, 9))]
    conflict = {(a, b) for a, b in combinations(range(4), 2)
                if intervals[a] & intervals[b]}
    require(conflict == star, "tribonacci odd-cycle conflict claw")
    for edges, kernel in [(star, (1, -1, -1, -1)),
                          (triangle, (1, 0, 0, 0)),
                          ({(0, 1), (1, 2), (2, 3)}, (1, -1, 1, -1))]:
        require(all(kernel[a] + kernel[b] == 0 for a, b in edges), "located signless kernel")
    return {str(key): value for key, value in sorted(classes.items())}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=100)
    parser.add_argument("--large-order", type=int, default=0)
    args = parser.parse_args()
    require(args.max_order >= 24, "audit order at least twenty-four")
    bound = args.max_order
    tournaments, tstats = profile_dp(bound)
    cores, cstats = profile_dp(bound, core=True)
    boundary, bstats = profile_dp(bound, core=True, no_fixed=True)
    for row in (tournaments, cores):
        require(all(x.denominator == 1 and x >= 0 for x in row), "class coefficients nonnegative integral")
    require(all(x >= 0 for x in boundary), "boundary coefficients nonnegative")
    for core, no_fixed, values in [(False, False, tournaments),
                                  (True, False, cores), (True, True, boundary)]:
        independent = direct_series(18, core, no_fixed)
        require(independent == values[:19], "direct partition/substitution control")
        alternate, _ = profile_dp(24, core=core, no_fixed=no_fixed, ascending=True)
        require(alternate == values[:25], "opposite elimination order")
    plain = direct_series(min(bound, 60))
    require(plain == tournaments[:len(plain)], "uncompressed Burnside through sixty")
    for n in range(7):
        fixed = sum((fixed_core_count(parts) for parts in odd_partitions(n)), Q())
        free = sum((fixed_core_count(parts) for parts in odd_partitions(n, 3)), Q())
        require(fixed == cores[n], "pair-orbit direct core Burnside")
        require(free == boundary[n], "pair-orbit direct fixed-free boundary")
    require(cores[:9] == list(map(Q, (1, 1, 0, 1, 2, 5, 26, 255, 4602))),
            "inherited core bank through eight")
    external = Path(__file__).with_name("a000568-reference.txt")
    external_count = 0
    for line in external.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        n, value = map(int, line.split())
        if n <= bound:
            require(tournaments[n] == value, "OEIS independent table")
            external_count += 1
    # Histories with the same mass and cycle count need different responses at 3.
    alpha, beta = (11, 5), (9, 7)
    require(sum(alpha) == sum(beta) and len(alpha) == len(beta), "hostile equal coarse data")
    require(sum(gcd(a, 3) for a in alpha) != sum(gcd(a, 3) for a in beta),
            "hostile distinguishes future three-cycle")
    require(sum(gcd(a, 1) for a in alpha) == sum(gcd(a, 1) for a in beta),
            "histories safely merge after three eliminated")
    result = {"scope": {"maximum_order": bound, "direct_Burnside": min(bound, 60),
                        "direct_core_substitution": 18, "orientation_orbits": 6,
                        "external_table_terms": external_count},
              "tournaments": list(map(str, tournaments)), "core_classes": list(map(str, cores)),
              "fixed_free_core_mass": list(map(str, boundary)),
              "statistics": {"tournaments": tstats, "core_classes": cstats, "boundary": bstats},
              "paths": path_checks(), "four_vertex_graphs": graph_checks(),
              "reference_sha256": sha256(external.read_bytes()).hexdigest()}
    require(sha256(str(tournaments[100]).encode()).hexdigest()
            == "1d290d31b74fb08a9e73a6de3e842fc4a1d5bca8c1657107e7f14c0a0471e9c4"
            if bound >= 100 else True, "incoming independently recorded A100")
    if args.large_order:
        large, stats = profile_dp(args.large_order)
        require(large[:bound+1] == tournaments, "bound extension preserves earlier coefficients")
        require(all(x.denominator == 1 for x in large), "large coefficients integral")
        result["large"] = {"order": args.large_order, "coefficients": list(map(str, large)), "statistics": stats}
    result["checks"] = dict(sorted(CHECKS.items()))
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
