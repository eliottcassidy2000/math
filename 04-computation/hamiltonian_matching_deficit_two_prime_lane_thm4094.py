#!/usr/bin/env python3
"""Exact referee for THM-4094.

The proof of the infinite statements is in the theorem.  This companion
independently checks the minimal insertion hostiles, every deletion pair
through order five, an actual support/multiplicity collision, the first
finite-prefix multiplicative boundary, the two-prime-lane decomposition
through one million, and finite witnesses used in the compactness argument.

Only the Python standard library is used.  Every executable gate goes through
``require`` so ``python -O`` cannot disable the audit.
"""

from hashlib import sha256
from itertools import permutations
import json
from math import prod


def require(condition, label):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + label)


def pair_bit(n, i, j):
    require(0 <= i < j < n, "pair_bit domain")
    return i * (2 * n - i - 1) // 2 + (j - i - 1)


def edge(mask, n, i, j):
    """Return whether the tournament encoded by ``mask`` has arc i -> j."""
    require(0 <= i < n and 0 <= j < n and i != j, "edge domain")
    if i > j:
        return not edge(mask, n, j, i)
    return bool((mask >> pair_bit(n, i, j)) & 1)


def mask_from_arc_set(n, arcs):
    arcs = set(arcs)
    mask = 0
    for i in range(n):
        for j in range(i + 1, n):
            forward = (i, j) in arcs
            backward = (j, i) in arcs
            require(forward != backward, "exactly one orientation per pair")
            if forward:
                mask |= 1 << pair_bit(n, i, j)
    return mask


def hamilton_paths(mask, n):
    return tuple(
        path
        for path in permutations(range(n))
        if all(edge(mask, n, path[i], path[i + 1]) for i in range(n - 1))
    )


def is_strong(mask, n):
    for start in range(n):
        seen = {start}
        stack = [start]
        while stack:
            i = stack.pop()
            for j in range(n):
                if j not in seen and edge(mask, n, i, j):
                    seen.add(j)
                    stack.append(j)
        if len(seen) != n:
            return False
    return True


def deletion_incidence(mask, n, vertex):
    """Return the full old/new deletion incidence and its exact deficit."""
    old_vertices = tuple(i for i in range(n) if i != vertex)
    old_paths = tuple(
        path
        for path in permutations(old_vertices)
        if all(edge(mask, n, path[i], path[i + 1]) for i in range(n - 2))
    )
    new_paths = hamilton_paths(mask, n)
    fibers = {path: [] for path in old_paths}
    orphans = []
    for new_path in new_paths:
        reduced = tuple(x for x in new_path if x != vertex)
        if reduced in fibers:
            fibers[reduced].append(new_path)
        else:
            orphans.append(new_path)
    require(all(fibers[path] for path in old_paths), "left-total insertion incidence")
    type_ii_total = sum(len(fibers[path]) - 1 for path in old_paths)
    deficit = type_ii_total + len(orphans)
    require(len(new_paths) - len(old_paths) == deficit, "deletion deficit identity")
    return old_paths, new_paths, fibers, tuple(orphans), type_ii_total, deficit


def sieve_smallest_prime(limit):
    smallest = list(range(limit + 1))
    if limit >= 1:
        smallest[1] = 1
    p = 2
    while p * p <= limit:
        if smallest[p] == p:
            for multiple in range(p * p, limit + 1, p):
                if smallest[multiple] == multiple:
                    smallest[multiple] = p
        p += 1
    return smallest


def is_prime(n, smallest):
    return n >= 2 and smallest[n] == n


def prime_factors(n, smallest):
    factors = []
    while n > 1:
        p = smallest[n]
        factors.append(p)
        n //= p
    return factors


def product_closure(generators, cap):
    reached = {1}
    changed = True
    while changed:
        changed = False
        for a in tuple(reached):
            for b in generators:
                value = a * b
                if value <= cap and value not in reached:
                    reached.add(value)
                    changed = True
    return reached


def two_lane_decomposition(n, smallest):
    """Construct the proof decomposition into p, 7p, 49, 63, and 343 atoms."""
    require(n % 2 == 1 and n not in {7, 21}, "two-lane input is allowed")
    valuation = 0
    unit = n
    while unit % 7 == 0:
        valuation += 1
        unit //= 7
    factors = prime_factors(unit, smallest)
    if valuation == 0:
        return factors
    if valuation % 2 == 0:
        return [49] * (valuation // 2) + factors
    if valuation >= 3:
        return [343] + [49] * ((valuation - 3) // 2) + factors
    require(valuation == 1 and unit not in {1, 3}, "excluded 7 and 21 boundary")
    if all(p == 3 for p in factors):
        require(len(factors) >= 2, "63 carry exponent")
        return [63] + [3] * (len(factors) - 2)
    index = next(i for i, p in enumerate(factors) if p >= 5)
    prime = factors.pop(index)
    return [7 * prime] + factors


def main():
    ledger = {}

    transitive_three = mask_from_arc_set(3, {(0, 1), (0, 2), (1, 2)})
    cyclic_three = mask_from_arc_set(3, {(0, 1), (1, 2), (2, 0)})
    minimal_rows = []
    print("MINIMAL INSERTION HOSTILES (delete v=2; common base arc 0->1)")
    for name, mask in (("transitive", transitive_three), ("C3", cyclic_three)):
        old, new, fibers, orphans, type_ii, deficit = deletion_incidence(mask, 3, 2)
        degrees = tuple(len(fibers[path]) for path in old)
        row = {
            "name": name,
            "mask": mask,
            "H_old": len(old),
            "H_new": len(new),
            "degrees": degrees,
            "orphans": len(orphans),
            "type_ii": type_ii,
            "deficit": deficit,
        }
        minimal_rows.append(row)
        print(
            name,
            "mask=", mask,
            "H_old=", len(old),
            "H_new=", len(new),
            "left_degrees=", degrees,
            "type_ii=", type_ii,
            "orphans=", len(orphans),
            "deficit=", deficit,
        )
        print("  old_paths=", old)
        print("  new_paths=", new)
        print("  orphan_paths=", orphans)
    require(minimal_rows[0]["deficit"] == 0, "transitive hostile deficit")
    require(minimal_rows[1]["deficit"] == 2, "C3 hostile deficit")
    ledger["minimal_rows"] = minimal_rows

    print("\nFULL INSERTION-DEFICIT AUDIT")
    insertion_rows = []
    for n in range(2, 6):
        pairs = 0
        failures = 0
        positive_orphan_pairs = 0
        paper_identity_failures = 0
        orphan_but_balanced = 0
        for mask in range(1 << (n * (n - 1) // 2)):
            for vertex in range(n):
                pairs += 1
                old, _, fibers, orphans, type_ii, deficit = deletion_incidence(
                    mask, n, vertex
                )
                positive_orphan_pairs += bool(orphans)
                paper_identity_failures += type_ii != len(orphans)
                orphan_but_balanced += bool(orphans) and type_ii == len(orphans)
                failures += deficit < 0 or not all(fibers[path] for path in old)
        row = {
            "n": n,
            "pairs": pairs,
            "failures": failures,
            "positive_orphan_pairs": positive_orphan_pairs,
            "paper_identity_failures": paper_identity_failures,
            "orphan_but_balanced": orphan_but_balanced,
        }
        insertion_rows.append(row)
        print(
            "n=", n,
            "pairs=", pairs,
            "failures=", failures,
            "positive_orphan_pairs=", positive_orphan_pairs,
            "paper_identity_failures=", paper_identity_failures,
            "orphan_but_balanced=", orphan_but_balanced,
        )
        require(failures == 0, "full insertion audit n=" + str(n))
    require(
        [row["pairs"] for row in insertion_rows] == [4, 24, 256, 5120],
        "insertion universe sizes",
    )
    ledger["insertion_rows"] = insertion_rows

    strong_h9_mask = None
    for mask in range(1 << 10):
        if is_strong(mask, 5) and len(hamilton_paths(mask, 5)) == 9:
            strong_h9_mask = mask
            break
    require(strong_h9_mask is not None, "strong order-five H=9 witness")
    join_arcs = set()
    for shift in (0, 3):
        join_arcs |= {
            (shift, shift + 1),
            (shift + 1, shift + 2),
            (shift + 2, shift),
        }
    join_arcs |= {(i, j) for i in range(3) for j in range(3, 6)}
    c3_join_c3_mask = mask_from_arc_set(6, join_arcs)
    join_h = len(hamilton_paths(c3_join_c3_mask, 6))
    require(join_h == 9, "C3 order-join C3 has H=9")
    require(not is_strong(c3_join_c3_mask, 6), "C3 join C3 is reducible")
    support_row = {
        "strong_n5_mask": strong_h9_mask,
        "strong_H": 9,
        "join_n6_mask": c3_join_c3_mask,
        "join_H": join_h,
    }
    ledger["support_collision"] = support_row
    print("\nSUPPORT VERSUS OBJECT-MULTIPLICITY HOSTILE")
    print("strong_n5_mask=", strong_h9_mask, "H=", 9)
    print(
        "C3_join_C3_mask=", c3_join_c3_mask,
        "H=", join_h,
        "strong=", is_strong(c3_join_c3_mask, 6),
    )
    print("two distinct objects; support={9}; object-term multiset=[9,9]")

    smallest = sieve_smallest_prime(1_000_000)
    prefix = {m for m in range(1, 610, 2) if m not in {7, 21}}
    reached = product_closure(prefix, 5000)
    allowed = [m for m in range(1, 5001, 2) if m not in {7, 21}]
    missing = [m for m in allowed if m not in reached]
    require(611 in reached, "611 is forced as 13*47")
    require(is_prime(613, smallest) and 613 not in reached, "613 prefix boundary")
    require(missing[0] == 613, "first prefix-closure miss")
    prefix_row = {
        "base_count": len(prefix),
        "base_max": max(prefix),
        "first_twenty_missing": missing[:20],
    }
    ledger["prefix_closure"] = prefix_row
    print("\nFINITE-PREFIX / MULTIPLICATIVE-CLOSURE HOSTILE")
    print("base_count=", len(prefix), "base_max=", max(prefix))
    print("611_factorization=13*47; reached=", 611 in reached)
    print("613_prime=", is_prime(613, smallest), "reached=", 613 in reached)
    print("first_20_allowed_missing_through_5000=", missing[:20])

    lane_failures = []
    atom_kinds = set()
    checked = 0
    for n in range(1, 1_000_001, 2):
        if n in {7, 21}:
            continue
        checked += 1
        decomposition = two_lane_decomposition(n, smallest)
        if prod(decomposition) != n:
            lane_failures.append((n, decomposition))
        for atom in decomposition:
            if is_prime(atom, smallest) and atom != 7:
                atom_kinds.add("ordinary-prime")
            elif atom in {49, 63, 343}:
                atom_kinds.add("finite-7-carry")
            elif atom % 7 == 0 and is_prime(atom // 7, smallest) and atom // 7 != 3:
                atom_kinds.add("7-prime")
            else:
                lane_failures.append((n, atom))
    require(not lane_failures, "two-prime-lane decomposition")
    boundary_inputs = (49, 63, 343, 611, 613, 623, 9261)
    boundary_decompositions = {
        str(n): two_lane_decomposition(n, smallest) for n in boundary_inputs
    }
    lane_row = {
        "checked": checked,
        "cap": 1_000_000,
        "failures": len(lane_failures),
        "atom_kinds": sorted(atom_kinds),
        "boundary_decompositions": boundary_decompositions,
    }
    ledger["two_lane"] = lane_row
    print(
        "two_lane_decomposition_checked_allowed_odds_through=", 1_000_000,
        "checked=", checked,
        "failures=", len(lane_failures),
        "atom_kinds=", sorted(atom_kinds),
    )
    print("boundary_decompositions=", boundary_decompositions)

    compactness_rows = []
    print("\nCOMPACTNESS FINITE-SUBSET WITNESSES (C3 order joins)")
    for threshold in (1, 2, 10, 100, 1000, 10000, 1_000_000):
        copies = 0
        h_value = 1
        while h_value < threshold:
            copies += 1
            h_value *= 3
        require(h_value >= threshold, "compactness finite witness")
        vertices = 1 if copies == 0 else 3 * copies
        row = {
            "threshold": threshold,
            "copies": copies,
            "vertices": vertices,
            "H": h_value,
        }
        compactness_rows.append(row)
        print(
            "threshold=", threshold,
            "copies=", copies,
            "vertices=", vertices,
            "H=", h_value,
        )
    ledger["compactness_rows"] = compactness_rows

    semantic_bytes = json.dumps(
        ledger, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    semantic_digest = sha256(semantic_bytes).hexdigest()
    print("\nsemantic_sha256=", semantic_digest)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
