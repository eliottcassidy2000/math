#!/usr/bin/env python3
"""Exact verifier for the universal weighted r=6 finite-horn atlas.

The old finite-horn loop was indexed by a seven-speed core P and then by
six killers.  This verifier changes the vertex set.  Its vertices are 35
rational-time obligations (q,a), all of which are safe for *every* speed
1,...,12.  A candidate killer k is incident to the obligations on which it
is unsafe.  Integer weights give a dual set-cover certificate:

    total obligation weight = 505,
    every k in [92,332] kills weight at most 84.

Consequently six killers kill weight at most 504, so at least one rational
time survives.  All arithmetic below is integer arithmetic; no optimizer is
used by the verifier.  The weights were discovered from the exact rational
dual of the finite set-cover LP and then frozen here.

Tournament Analysis is deliberately only a diagnostic quotient here.  On
killer vertices use load(k)-load(l) as the pairwise observable and break a
load tie by k.  The resulting tournament is transitive.  It retains the
single-killer capacity used by the proof but destroys obligation ownership
and higher intersections; the weighted incidence hypergraph is the actual
proof-facing object.
"""

from collections import Counter


def least_abs_residue(x: int, q: int) -> int:
    r = x % q
    return min(r, q - r)


def threshold(q: int) -> int:
    return (q + 13) // 14


# (q, a, integer weight).  The common LP-dual denominator was 84.
ATLAS = (
    (26, 2, 81),
    (27, 2, 7),
    (28, 2, 3),
    (28, 13, 7),
    (40, 3, 10),
    (41, 3, 3),
    (41, 19, 5),
    (42, 13, 1),
    (53, 4, 16),
    (55, 4, 2),
    (56, 13, 21),
    (68, 5, 3),
    (68, 21, 21),
    (69, 5, 1),
    (69, 16, 12),
    (70, 27, 20),
    (79, 6, 3),
    (81, 25, 5),
    (82, 19, 6),
    (83, 6, 1),
    (83, 32, 5),
    (84, 13, 8),
    (93, 7, 9),
    (94, 7, 13),
    (97, 7, 1),
    (97, 30, 2),
    (98, 15, 6),
    (105, 8, 4),
    (106, 49, 46),
    (107, 33, 10),
    (109, 42, 49),
    (110, 17, 25),
    (111, 8, 1),
    (111, 17, 64),
    (112, 43, 34),
)

KILLER_MIN = 92
KILLER_MAX = 332
CAPACITY = 84


def kill_load(k: int) -> int:
    return sum(
        weight
        for q, a, weight in ATLAS
        if least_abs_residue(k * a, q) < threshold(q)
    )


def main() -> None:
    total_weight = sum(weight for _, _, weight in ATLAS)
    assert len(ATLAS) == 35
    assert total_weight == 505
    assert max(q for q, _, _ in ATLAS) == 112

    # Every atlas time is safe for the entire possible core {1,...,12};
    # hence it is safe for every seven-speed subcore.
    core_slacks = []
    for q, a, _ in ATLAS:
        h = threshold(q)
        slack = min(least_abs_residue(p * a, q) - h for p in range(1, 13))
        assert slack >= 0
        core_slacks.append(slack)

    loads = {k: kill_load(k) for k in range(KILLER_MIN, KILLER_MAX + 1)}
    max_load = max(loads.values())
    maximizers = [k for k, load in loads.items() if load == max_load]
    assert max_load == CAPACITY
    assert 6 * max_load == 504 < total_weight

    histogram = Counter(loads.values())
    raw_tie_pairs = sum(n * (n - 1) // 2 for n in histogram.values())

    print("### universal weighted r=6 finite-horn atlas ###")
    print(f"atlas obligations: {len(ATLAS)}")
    print(f"modulus range: {min(q for q, _, _ in ATLAS)}..{max(q for q, _, _ in ATLAS)}")
    print(f"total integer weight: {total_weight}")
    print(f"minimum all-core safety slack: {min(core_slacks)}")
    print(f"candidate killers checked exactly: {KILLER_MIN}..{KILLER_MAX} ({len(loads)})")
    print(f"maximum killed weight: {max_load}")
    print("maximizing killers: " + " ".join(map(str, maximizers)))
    print(f"six-killer capacity: 6*{max_load} = {6 * max_load} < {total_weight}")
    print("CONCLUSION: every six candidate killers leave an atlas obligation safe")
    print()
    print("### frozen obligation table (q a weight core_slack) ###")
    for (q, a, weight), slack in zip(ATLAS, core_slacks):
        print(f"{q:3d} {a:3d} {weight:3d} {slack:3d}")
    print()
    print("### killer-load histogram (load: count) ###")
    for load in sorted(histogram):
        print(f"{load:3d}: {histogram[load]:3d}")
    print()
    print("### Tournament Analysis (lossy diagnostic) ###")
    print("vertices: the 241 candidate killers")
    print("observable: load(k)-load(l); tie gauge: increasing integer k")
    print("tournament: transitive; directed 3-cycles: 0; SCCs: 241")
    print("score histogram after tie gauge: every score 0..240 occurs once")
    print("Hamiltonian paths: 1")
    print(f"edges changed by reversing only the tie gauge: {raw_tie_pairs}")
    print("proof-facing object: weighted obligation--killer incidence hypergraph")
    print("destroyed by tournament quotient: obligation owners and higher overlaps")


if __name__ == "__main__":
    main()
