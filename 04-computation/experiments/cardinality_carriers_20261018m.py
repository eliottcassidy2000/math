#!/usr/bin/env python3
"""Exact labelled carriers behind the fifteen atlas; standard library only.

Universe: all even subsets of 2g+2 labels modulo complement, g=1,2,3.
No arithmetic, geometric, or rationality filter is silently imposed.
Two independent constructions: isotropic-span enumeration and affine-plane
permutation orbits. All checks survive python -O. No novelty is claimed for
the classical symplectic, Fano, or Hamming-code correspondences.
"""
import argparse
from collections import Counter
from itertools import combinations, permutations
import json
from math import factorial, prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def matchings(labels):
    if not labels:
        yield ()
        return
    first = labels[0]
    for j in labels[1:]:
        rest = tuple(x for x in labels if x not in (first, j))
        for tail in matchings(rest):
            yield ((first, j),) + tail


def relabel(mask, perm):
    return sum(1 << perm[i] for i in range(len(perm)) if mask >> i & 1)


def carrier(g):
    n = 2 * g + 2
    full = (1 << n) - 1
    canonical = lambda a: min(a, a ^ full)
    vectors = sorted({canonical(a) for a in range(full + 1)
                      if a.bit_count() % 2 == 0})
    pairing = lambda a, b: (a & b).bit_count() % 2
    require(len(vectors) == 2 ** (2 * g), "quotient dimension")
    require(all(pairing(a, a) == 0 for a in vectors), "alternating")
    require(all(a == 0 or any(pairing(a, b) for b in vectors)
                for a in vectors), "nondegenerate")
    spaces = {frozenset([0])}
    level_counts = [1]
    for dimension in range(1, g + 1):
        next_spaces = set()
        for space in spaces:
            for v in vectors:
                if v in space or any(pairing(v, x) for x in space):
                    continue
                next_spaces.add(space | frozenset(canonical(v ^ x) for x in space))
        spaces = next_spaces
        level_counts.append(len(spaces))
        require(all(len(s) == 2 ** dimension for s in spaces), "span size")
    require(len(spaces) == prod(2 ** i + 1 for i in range(1, g + 1)),
            "independent Lagrangian incidence formula")
    matching_images = {}
    for matching in matchings(tuple(range(n))):
        space = frozenset([0])
        for i, j in matching:
            v = canonical((1 << i) | (1 << j))
            space |= frozenset(canonical(x ^ v) for x in space)
        require(space in spaces, "matching isotropy")
        require(space not in matching_images, "matching injection")
        matching_images[space] = matching
    require(len(matching_images) == prod(range(1, n, 2)), "double factorial")
    pair_vectors = {canonical((1 << i) | (1 << j)) for i, j in combinations(range(n), 2)}
    pair_hist = Counter(len(s & pair_vectors) for s in spaces)
    return canonical, spaces, matching_images, {
        "g": g, "labels": n, "nonzero_vectors": len(vectors) - 1,
        "isotropic_subspaces_by_dimension": level_counts,
        "lagrangians": len(spaces), "matching_images": len(matching_images),
        "missing": len(spaces) - len(matching_images),
        "pair_class_histogram": dict(sorted(pair_hist.items())),
    }


def set_partitions(labels):
    if not labels:
        yield ()
        return
    first, *rest = labels
    for partition in set_partitions(rest):
        yield ((first,),) + partition
        for i in range(len(partition)):
            yield partition[:i] + ((first,) + partition[i],) + partition[i + 1:]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = [carrier(g) for g in (1, 2, 3)]
    canonical, spaces, images, _ = data[1]
    duads = {sum(1 << i for i in pair) for pair in combinations(range(6), 2)}
    fixed = []
    for name, perm in [("transposition", (1, 0, 2, 3, 4, 5)),
                       ("triple_transposition", (1, 0, 3, 2, 5, 4))]:
        fixed.append({"cycle_type": name,
                      "duads": sum(relabel(v, perm) == v for v in duads),
                      "matchings": sum(frozenset(canonical(relabel(v, perm)) for v in s) == s
                                       for s in spaces)})
    require([(r["duads"], r["matchings"]) for r in fixed] == [(7, 3), (3, 7)],
            "hostile: untwisted S6 point-line bijection")
    for perm in permutations(range(6)):
        for space, matching in images.items():
            moved_pairs = [canonical((1 << perm[i]) | (1 << perm[j])) for i, j in matching]
            moved_space = frozenset(canonical(relabel(v, perm)) for v in space)
            require(moved_space == frozenset([0] + moved_pairs), "S6 equivariance")
    require(all(sum(canonical(v) in s for s in spaces) == 3 for v in duads),
            "doily incidence in complement-quotient coordinates")

    # Mark label 5: remove its pair, leaving an A5 double transposition
    # with the partner of 5 fixed. Conversely, adjoin that pair again.
    involutions = set()
    for matching in images.values():
        perm = list(range(5))
        for i, j in matching:
            if j != 5 and i != 5:
                perm[i], perm[j] = j, i
        involutions.add(tuple(perm))
    require(len(involutions) == 15, "pointed A5 bijection")
    require(all(sum(i == p[i] for i in range(5)) == 1 for p in involutions),
            "one marked fixed point")

    canonical8, spaces8, images8, _ = data[2]
    full8 = 255
    missing = spaces8 - images8.keys()
    designs = set()
    for space in missing:
        blocks = frozenset(v for x in space - {0} for v in (x, x ^ full8))
        require(len(blocks) == 14 and all(v.bit_count() == 4 for v in blocks),
                "four-subset lift")
        for triple in combinations(range(8), 3):
            mask = sum(1 << i for i in triple)
            require(sum(mask & b == mask for b in blocks) == 1, "S(3,4,8)")
        code = blocks | {0, full8}
        require(all((a ^ b) in code for a in code for b in code), "Hamming closure")
        require(all((a & b).bit_count() % 2 == 0 for a in code for b in code),
                "self-dual by size and orthogonality")
        require(Counter(v.bit_count() for v in code) == {0: 1, 4: 14, 8: 1},
                "extended Hamming weight distribution")
        designs.add(blocks)
    # Independent: all affine planes a.x=b in the fixed affine space F2^3.
    affine = frozenset(sum(1 << x for x in range(8) if (a & x).bit_count() % 2 == b)
                       for a in range(1, 8) for b in (0, 1))
    orbit = {frozenset(relabel(block, perm) for block in affine)
             for perm in permutations(range(8))}
    require(orbit == designs and len(orbit) == 30, "independent affine orbit")
    require(factorial(8) // len(orbit) == 8 * 168, "AGL stabilizer")
    punctured = {v & 127 for v in affine | {0, full8}}
    require(len(punctured) == 16, "Hamming puncture is injective")
    balls = [v ^ e for v in punctured for e in [0] + [1 << i for i in range(7)]]
    require(len(balls) == len(set(balls)) == 128, "perfect packing")

    # Bell(4)=15 is labelled; quotienting labels gives p(4)=5.
    partitions4 = list(set_partitions(list(range(4))))
    shapes = Counter(tuple(sorted(map(len, p), reverse=True)) for p in partitions4)
    require(len(partitions4) == 15 and len(shapes) == 5, "partition typing")
    require(len(list(set_partitions(list(range(5))))) == 52, "Bell neighbour hostile")
    # A split symplectic four-space has 9 product lines and 6 graph lines.
    # Fix nonorthogonal u,v. Their plane is nondegenerate.
    u, v = 3, 5
    plane = {0, u, v, u ^ v}
    orthogonal = {x for x in range(32) if x.bit_count() % 2 == 0
                  and (x & u).bit_count() % 2 == 0 and (x & v).bit_count() % 2 == 0}
    require(len(plane) == len(orthogonal) == 4 and plane & orthogonal == {0},
            "symplectic split")
    intersections = Counter((len(s & plane), len(s & orthogonal)) for s in spaces)
    require(intersections == {(2, 2): 9, (1, 1): 6}, "product versus graph kernels")
    # Independent signed permutation matrices for two-qubit X^x Z^z.
    symplectic_basis = (3, 5, 15, 23)
    pauli = {}
    for label in range(16):
        vector = 0
        for i, basis_vector in enumerate(symplectic_basis):
            if label >> i & 1:
                vector ^= basis_vector
        x = (label & 1) | ((label >> 1) & 2)
        z = ((label >> 1) & 1) | ((label >> 2) & 2)
        pauli[vector] = (x, z)
    require(len(pauli) == 16, "Pauli symplectic basis")
    for va, (xa, za) in pauli.items():
        for vb, (xb, zb) in pauli.items():
            sign = -1 if (va & vb).bit_count() % 2 else 1
            for t in range(4):
                ab = (-1) ** ((zb & t).bit_count() + (za & (t ^ xb)).bit_count())
                ba = (-1) ** ((za & t).bit_count() + (zb & (t ^ xa)).bit_count())
                require(ab == sign * ba, "Pauli commutator from actual matrix action")
    # Same F21 action on Fano flags and Paley arcs, with chosen base flag/arc.
    difference_set = {0, 1, 3}
    fano_lines = {frozenset((x + b) % 7 for x in difference_set) for b in range(7)}
    flag_to_arc = {}
    for a in (1, 2, 4):
        for b in range(7):
            line = frozenset((a * x + b) % 7 for x in difference_set)
            require(line in fano_lines and b in line, "F21 preserves Fano incidence")
            flag_to_arc[(b, line)] = (b, (a + b) % 7)
    require(len(flag_to_arc) == len(set(flag_to_arc.values())) == 21,
            "regular flag and arc actions")
    for a in (1, 2, 4):
        for b in range(7):
            for (point, line), (u, v) in flag_to_arc.items():
                moved = ((a * point + b) % 7,
                         frozenset((a * x + b) % 7 for x in line))
                require(flag_to_arc[moved] == ((a * u + b) % 7, (a * v + b) % 7),
                        "F21 equivariance after base choice")
    receipt = {
        "status": "FINITE-EXACT; classical carriers independently enumerated",
        "universe": "all even subsets modulo complement for 4,6,8 labelled points",
        "genera": [d[3] for d in data], "fixed_points_hostile": fixed,
        "s6_equivariance_cases": 720 * 15, "doily_point_degrees": [3] * 15,
        "marked_A5_involutions": len(involutions),
        "F21_flag_arc_equivariance_cases": 21 * 21,
        "Pauli_matrix_commutator_cases": 16 * 16 * 4,
        "genus2_split_intersections": {str(k): v for k, v in sorted(intersections.items())},
        "missing_genus3_affine_designs": len(designs),
        "affine_permutations_checked": factorial(8), "affine_stabilizer": 1344,
        "extended_Hamming_weights": {"0": 1, "4": 14, "8": 1},
        "punctured_Hamming_perfect_packing": {"codewords": 16, "ball_size": 8, "total": 128},
        "set_partition_4_shapes": {str(k): v for k, v in sorted(shapes.items())},
        "Bell_5_hostile": 52, "nonzero_F2_5": 31,
    }
    encoded = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end="")


if __name__ == "__main__":
    main()
