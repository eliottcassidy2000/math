#!/usr/bin/env python3
"""Exact independent paths for signed deletion-card recurrence, stdlib only.

No imports of repository producers. All checks survive python -O.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial, gcd
import json

CHECKS = Counter()


def check(ok, tag):
    CHECKS[tag] += 1
    if not ok:
        raise RuntimeError(tag)


@lru_cache(None)
def pairs(n):
    return tuple(combinations(range(n), 2))


@lru_cache(None)
def perms(n):
    return tuple(permutations(range(n)))


@lru_cache(None)
def actions(n):
    index = {edge: i for i, edge in enumerate(pairs(n))}
    out = []
    for p in perms(n):
        edgebits = tuple(1 << index[tuple(sorted((p[a], p[b])))]
                         for a, b in pairs(n))
        reverse = sum(edgebits[i] for i, (a, b) in enumerate(pairs(n))
                      if p[a] > p[b])
        out.append((p, edgebits, reverse))
    return tuple(out)


def image(mask, edgebits):
    result = 0
    while mask:
        bit = mask & -mask
        result |= edgebits[bit.bit_length() - 1]
        mask -= bit
    return result


@lru_cache(None)
def catalogue(n, tournament=False):
    unseen = set(range(1 << len(pairs(n))))
    representatives = []
    canonical = {}
    while unseen:
        rep = min(unseen)
        representatives.append(rep)
        orbit = {image(rep, bits) ^ (reverse if tournament else 0)
                 for _, bits, reverse in actions(n)}
        for mask in orbit:
            canonical[mask] = rep
        unseen.difference_update(orbit)
    return tuple(representatives), canonical


def rows(n, mask):
    out = [0] * n
    for k, (a, b) in enumerate(pairs(n)):
        if mask >> k & 1:
            out[a] |= 1 << b
            out[b] |= 1 << a
    return out


def induced_mask(adjacency, vertices):
    return sum(1 << k for k, (a, b) in enumerate(pairs(len(vertices)))
               if adjacency[vertices[a]] >> vertices[b] & 1)


def sign(n, mask, p):
    return (-1) ** sum(1 for k, (a, b) in enumerate(pairs(n))
                       if mask >> k & 1 and p[a] > p[b])


def edge_cycle_sign(n, mask, p):
    """Independent edge-orbit route: odd iff an orbit returns reversed."""
    remaining = {edge for k, edge in enumerate(pairs(n)) if mask >> k & 1}
    parity = 0
    while remaining:
        a, b = first = min(remaining)
        while True:
            remaining.discard(tuple(sorted((a, b))))
            a, b = p[a], p[b]
            if tuple(sorted((a, b))) == first:
                parity ^= (a > b)
                break
    return (-1) ** parity


@lru_cache(None)
def automorphisms(n, mask):
    return tuple(p for p, bits, _ in actions(n) if image(mask, bits) == mask)


def even(n, mask):
    return all(sign(n, mask, p) == 1 for p in automorphisms(n, mask))


def vertex_orbits(n, group):
    unseen = set(range(n))
    out = []
    while unseen:
        v = min(unseen)
        orbit = {p[v] for p in group}
        out.append(tuple(sorted(orbit)))
        unseen.difference_update(orbit)
    return out


def subset_image(subset, p):
    return sum(1 << p[v] for v in range(len(p)) if subset >> v & 1)


def subset_orbits(n, group):
    unseen = set(range(1 << n))
    out = []
    while unseen:
        subset = min(unseen)
        out.append(subset)
        unseen.difference_update(subset_image(subset, p) for p in group)
    return out


def extension(n, mask, subset):
    adjacency = rows(n, mask)
    adjacency.append(subset)
    for v in range(n):
        if subset >> v & 1:
            adjacency[v] |= 1 << n
    return induced_mask(adjacency, tuple(range(n + 1)))


def card_maps(n, parent, subset, child, w):
    """Order-n marked-card isomorphisms, without calling child Aut/canon."""
    source = rows(n, parent)
    target = rows(n + 1, child)
    available = tuple(u for u in range(n + 1) if u != w)
    out = []
    for destinations in permutations(available):
        if any(bool(subset >> a & 1) != bool(target[w] >> destinations[a] & 1)
               for a in range(n)):
            continue
        if any(bool(source[a] >> b & 1) !=
               bool(target[destinations[a]] >> destinations[b] & 1)
               for a, b in pairs(n)):
            continue
        out.append(destinations + (w,))
    return out


def components(n, mask):
    adjacency = rows(n, mask)
    unseen = set(range(n))
    out = []
    while unseen:
        reached = {min(unseen)}
        todo = list(reached)
        while todo:
            u = todo.pop()
            for v in tuple(unseen - reached):
                if adjacency[u] >> v & 1:
                    reached.add(v)
                    todo.append(v)
        unseen.difference_update(reached)
        out.append(tuple(sorted(reached)))
    return out


def odd_partitions(total, minimum=1):
    if total == 0:
        yield ()
    for first in range(minimum, total + 1, 2):
        for rest in odd_partitions(total - first, first):
            yield (first,) + rest


def burnside(n):
    result = Fraction(0)
    for partition in odd_partitions(n):
        counts = Counter(partition)
        denominator = 1
        for length, multiplicity in counts.items():
            denominator *= length ** multiplicity * factorial(multiplicity)
        exponent = sum((a - 1) // 2 for a in partition)
        exponent += sum(gcd(a, b) for a, b in combinations(partition, 2))
        result += Fraction(1 << exponent, denominator)
    check(result.denominator == 1, "Burnside_integrality")
    return result.numerator


def main():
    limit = 6
    graph_counts, labelled_even, even_counts, tournament_counts = [], [], [], []
    odd_connected_root_orbits = [0] * (limit + 1)
    odd_root_orbits = [0] * (limit + 1)
    even_root_orbits = [0] * (limit + 1)
    point_tournaments = [0] * (limit + 1)
    for n in range(limit + 1):
        representatives, _ = catalogue(n)
        graph_counts.append(len(representatives))
        even_counts.append(sum(even(n, mask) for mask in representatives))
        labelled_even.append(sum(factorial(n) // len(automorphisms(n, mask))
                                 for mask in representatives if even(n, mask)))
        tournament_reps, _ = catalogue(n, True)
        tournament_counts.append(len(tournament_reps))
        check(even_counts[-1] == tournament_counts[-1] == burnside(n),
              "three_count_paths")
        for tournament in tournament_reps:
            group = tuple(p for p, bits, reverse in actions(n)
                          if image(tournament, bits) ^ reverse == tournament)
            point_tournaments[n] += len(vertex_orbits(n, group))
        for mask in representatives:
            group = automorphisms(n, mask)
            global_even = even(n, mask)
            for p in group:
                check(sign(n, mask, p) == edge_cycle_sign(n, mask, p),
                      "independent_edge_sign")
            adjacency = rows(n, mask)
            pieces = components(n, mask)
            odd_pieces = [piece for piece in pieces
                          if not even(len(piece), induced_mask(adjacency, piece))]
            check(global_even == (not odd_pieces), "component_evenness")
            for orbit in vertex_orbits(n, group):
                v = orbit[0]
                stabilizer = [p for p in group if p[v] == v]
                relative = all(sign(n, mask, p) == 1 for p in stabilizer)
                if relative:
                    transport = {(u, w): sign(n, mask, next(p for p in group if p[u] == w))
                                 for u in orbit for w in orbit}
                    for u in orbit:
                        for w in orbit:
                            for t in orbit:
                                check(transport[u, w] * transport[w, t] == transport[u, t],
                                      "signed_transport_composition")
                        check(sum(transport[u, w] for w in orbit) ==
                              (len(orbit) if global_even else 0), "signed_matrix_row_sum")
                if global_even:
                    even_root_orbits[n] += 1
                elif relative:
                    odd_root_orbits[n] += 1
                    odd_connected_root_orbits[n] += (len(pieces) == 1)
                    check(len(odd_pieces) == 1 and v in odd_pieces[0],
                          "single_odd_component_localization")
                    piece = odd_pieces[0]
                    localmask = induced_mask(adjacency, piece)
                    localgroup = automorphisms(len(piece), localmask)
                    localv = piece.index(v)
                    localorbit = {p[localv] for p in localgroup}
                    check(len(orbit) == len(localorbit), "local_orbit_size")
                    check(all(sign(len(piece), localmask, p) == 1
                              for p in localgroup if p[localv] == localv),
                          "local_root_evenness")
                if relative and not global_even:
                    signs = [sign(n, mask, next(p for p in group if p[v] == w))
                             for w in orbit]
                    check(sum(signs) == 0, "balanced_odd_orbit")
        check(point_tournaments[n] == even_root_orbits[n] + odd_root_orbits[n],
              "relative_even_rooted_identity")
        check(odd_root_orbits[n] == sum(even_counts[j] *
              odd_connected_root_orbits[n-j] for j in range(n + 1)),
              "defect_series_factorization")

    recurrence = []
    for n in range(limit):
        N = n + 1
        total_uniform = Fraction(0)
        total_filtered = 0
        parent_even_only = 0
        charges = Counter()
        child_canon = catalogue(N)[1]
        rooted = 0
        eligible = 0
        negative_transports = 0
        for parent in catalogue(n)[0]:
            group = automorphisms(n, parent)
            for subset in subset_orbits(n, group):
                rooted += 1
                K = [p for p in group if subset_image(subset, p) == subset]
                relative = all(sign(n, parent, p) == 1 for p in K)
                eligible += relative
                child = extension(n, parent, subset)
                childgroup = automorphisms(N, child)
                check(len(K) == sum(p[n] == n for p in childgroup),
                      "root_stabilizer_identification")
                check(relative == all(sign(N, child, p) == 1
                                      for p in childgroup if p[n] == n),
                      "root_sign_restriction")
                signed_sum = 0
                chosen_sum = 0
                root_orbit = set()
                for w in range(N):
                    maps = card_maps(n, parent, subset, child, w)
                    expected = {p for p in childgroup if p[n] == w}
                    check(set(maps) == expected, "card_maps_vs_child_automorphisms")
                    if not maps:
                        continue
                    root_orbit.add(w)
                    check(len(maps) == len(K), "card_transport_torsor")
                    signed_values = [sign(N, child, p) for p in maps]
                    signed_sum += sum(signed_values)
                    if relative:
                        check(len(set(signed_values)) == 1,
                              "transporter_choice_independence")
                        chosen_sum += signed_values[0]
                        negative_transports += (signed_values[0] == -1)
                    else:
                        check(sum(signed_values) == 0, "nonrelative_fibre_cancellation")
                charge = Fraction(signed_sum, len(K))
                expected_charge = len(root_orbit) if even(N, child) else 0
                check(charge == expected_charge, "root_charge_dichotomy")
                check(charge == (chosen_sum if relative else 0), "filtered_equals_uniform")
                total_uniform += charge
                total_filtered += chosen_sum if relative else 0
                parent_even_only += chosen_sum if relative and even(n, parent) else 0
                charges[child_canon[child]] += charge
        check(total_uniform / N == even_counts[N] == Fraction(total_filtered, N),
              "all_parent_recurrence")
        for child, charge in charges.items():
            check(charge == (N if even(N, child) else 0), "unrooted_child_total")
        recurrence.append({"child_order": N, "rooted_graph_types": rooted,
                           "relative_even_rooted_types": eligible,
                           "negative_transport_vertices": negative_transports,
                           "signed_charge": int(total_uniform),
                           "A000568": int(total_uniform / N),
                           "even_parent_only_weight": str(Fraction(parent_even_only, N))})

    # Smallest distinctions and obligatory odd-parent branch.
    check(tournament_counts[2] == even_counts[2] == 1 and labelled_even[2] == 1
          and 2 ** len(pairs(2)) == 2, "labelled_vs_unlabelled_hostile")
    check(even(3, 3) and not even(3, 7), "path_even_triangle_odd")
    check(all(row.bit_count() % 2 == 0 for row in rows(3, 7)) and
          any(row.bit_count() % 2 for row in rows(3, 3)), "Eulerian_not_sign_even")
    check(recurrence[2]["even_parent_only_weight"] == "4/3", "odd_parent_is_necessary")
    check(extension(2, 1, 1) == 3 and even(3, 3), "odd_K2_parent_to_even_P3")
    check(len(card_maps(2, 1, 1, 3, 2)) == 1 and
          len(automorphisms(2, 1)) == 2, "neighbourhood_mark_is_necessary")

    values = [burnside(n) for n in range(13)]
    check(values[:7] == [1, 1, 1, 2, 4, 12, 56], "classical_initial_values")
    print(json.dumps({"status": "PASS", "universe": "all graph/tournament types n=0..6; all attachments to graph parents n=0..5",
                      "checks": sum(CHECKS.values()), "checks_by_kind": dict(sorted(CHECKS.items())),
                      "graph_types": graph_counts, "sign_even_types": even_counts,
                      "labelled_sign_even_graphs": labelled_even,
                      "tournament_types": tournament_counts,
                      "pointed_tournament_types": point_tournaments,
                      "odd_graph_relative_even_root_orbits": odd_root_orbits,
                      "odd_connected_relative_even_root_orbits": odd_connected_root_orbits,
                      "card_recurrence": recurrence,
                      "classical_Burnside_A000568_through_12": values},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
