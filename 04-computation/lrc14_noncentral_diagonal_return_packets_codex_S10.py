#!/usr/bin/env python3
"""Exact replay for THM-802's noncentral diagonal-return packet family.

For H >= 1 put L=105H and

    W=(L+1,L+22,L-20,L-6,L+15,L-13,L+8,2L+2).

On the fixed interval I=[2/15,1/7], every block
[M/L,(M+1)/L], 14H <= M <= 15H-1, has owner word

    (7,2,5,3,0,6,4,1,7).

The packet count vector is (1,1,1,1,1,1,1,2).  Its inverse-step
products are all one modulo seven, so it is a reduced collision return, but
it is not THM-794's once-per-owner packet.  The script independently
enumerates every midpoint wall, checks every wall and chamber, audits the
five-core incidence, exhausts all owner words with this multiplicity, and
computes the marked/circular tournament fingerprints.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, permutations


P = 7
H_MIN = 1
H_MAX = 300
CORE = (1, 2, 11, 12, 13)
LEFT = Q(2, 15)
RIGHT = Q(1, 7)
D = (1, 1, 1, 1, 1, 1, 1, 2)
BETA = (1, 22, -20, -6, 15, -13, 8, 2)
N = (0, 3, -3, -1, 2, -2, 1, 0)
STEPS = (1, 1, 1, 1, 1, 1, 1, 4)
REDUCED_STATE = (0, 4, 3, 1, 5, 2, 6, 0)
PACKET_WORD = (7, 2, 5, 3, 0, 6, 4, 1, 7)
FIRST_OWNER_WORD = (7, 2, 5, 3, 0, 6, 4, 1)
EXPECTED_DIGEST = "431910fda5eef26208bf411ce79a371e5bce88b0871a83cd25c709453f64c335"


def floor_q(x: Q) -> int:
    return x.numerator // x.denominator


def ceil_q(x: Q) -> int:
    return -floor_q(-x)


def nearest_integer_or_none(x: Q) -> int | None:
    q = floor_q(x)
    residue = x - q
    if residue == Q(1, 2):
        return None
    return q + (residue > Q(1, 2))


def circle_norm(x: Q) -> Q:
    residue = x - floor_q(x)
    return min(residue, 1 - residue)


def family(height: int) -> tuple[int, tuple[int, ...]]:
    length = 105 * height
    return length, tuple(D[a] * length + BETA[a] for a in range(8))


def token(speed: int, x: Q) -> int | None:
    nearest = nearest_integer_or_none(speed * x)
    if nearest is None:
        return None
    return (-pow(speed, -1, P) * nearest) % P


def normalized_state(tokens: tuple[int, ...]) -> tuple[int, ...]:
    counts = Counter(tokens)
    roots = tuple(value for value, count in counts.items() if count == 2)
    assert len(roots) == 1
    root = roots[0]
    return tuple((value - root) % P for value in tokens)


def covered_chamber(speeds: tuple[int, ...], x: Q) -> bool:
    tokens = tuple(token(speed, x) for speed in speeds)
    return None not in tokens and set(tokens) == set(range(P))


def covered_wall(
    speeds: tuple[int, ...], event: tuple[Q, int, int]
) -> bool:
    x, owner, _ = event
    assert token(speeds[owner], x) is None
    others = tuple(
        token(speed, x) for a, speed in enumerate(speeds) if a != owner
    )
    return None not in others and sorted(others) == list(range(P))


def expected_event(length: int, block: int, owner: int, visit: int) -> tuple[Q, int, int]:
    speed = D[owner] * length + BETA[owner]
    wall_index = D[owner] * block + N[owner] + visit
    return Q(2 * wall_index + 1, 2 * speed), owner, wall_index


def expected_events(height: int) -> tuple[tuple[Q, int, int], ...]:
    length, _ = family(height)
    events = []
    for block in range(14 * height, 15 * height):
        events.extend(
            (
                expected_event(length, block, 7, 0),
                expected_event(length, block, 2, 0),
                expected_event(length, block, 5, 0),
                expected_event(length, block, 3, 0),
                expected_event(length, block, 0, 0),
                expected_event(length, block, 6, 0),
                expected_event(length, block, 4, 0),
                expected_event(length, block, 1, 0),
                expected_event(length, block, 7, 1),
            )
        )
    return tuple(events)


def enumerated_events(height: int) -> tuple[tuple[Q, int, int], ...]:
    _, speeds = family(height)
    events = []
    for owner, speed in enumerate(speeds):
        first = ceil_q(speed * LEFT - Q(1, 2))
        last = floor_q(speed * RIGHT - Q(1, 2))
        for wall_index in range(first, last + 1):
            x = Q(2 * wall_index + 1, 2 * speed)
            if LEFT <= x <= RIGHT:
                events.append((x, owner, wall_index))
    return tuple(sorted(events))


def prefix_legality(
    state: tuple[int, ...], word: tuple[int, ...]
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    total = 0
    counts = [0] * 8
    residues = []
    for owner in word:
        before = (state[owner] + total - counts[owner] * STEPS[owner]) % P
        assert before == 0
        residues.append(before)
        total = (total + STEPS[owner]) % P
        counts[owner] += 1
    defects = tuple(
        (total - counts[owner] * STEPS[owner]) % P for owner in range(8)
    )
    return tuple(counts), defects


def core_minimum(speed: int) -> Q:
    candidates = [LEFT, RIGHT]
    first_integer = ceil_q(speed * LEFT)
    last_integer = floor_q(speed * RIGHT)
    for integer in range(first_integer, last_integer + 1):
        point = Q(integer, speed)
        if LEFT <= point <= RIGHT:
            candidates.append(point)
    return min(circle_norm(speed * point) for point in candidates)


def audit_height(height: int) -> tuple:
    length, speeds = family(height)
    expected = expected_events(height)
    actual = enumerated_events(height)

    assert len(set(speeds)) == 8
    assert min(speeds) > 0
    assert tuple(speed % P for speed in speeds) == (1, 1, 1, 1, 1, 1, 1, 2)
    assert tuple(pow(speed, -1, P) for speed in speeds) == STEPS
    assert max(speeds) == speeds[7]
    assert sorted(speeds)[-2] == speeds[1]
    assert (speeds[7] + speeds[1] - 1) // speeds[1] == 2

    # Independent wall enumeration proves completeness, simplicity, and order.
    assert actual == expected
    assert len(actual) == 9 * height
    assert all(actual[i][0] < actual[i + 1][0] for i in range(len(actual) - 1))
    assert tuple(event[1] for event in actual) == PACKET_WORD * height

    # All walls and every chamber, including the two edge chambers, are covered.
    assert all(covered_wall(speeds, event) for event in actual)
    cuts = (LEFT,) + tuple(event[0] for event in actual) + (RIGHT,)
    assert all(
        covered_chamber(speeds, (left + right) / 2)
        for left, right in zip(cuts, cuts[1:])
    )
    assert covered_chamber(speeds, LEFT)
    assert covered_chamber(speeds, RIGHT)

    # At every block boundary the state is the same modulo diagonal translation.
    boundary_states = []
    for block in range(14 * height, 15 * height + 1):
        x = Q(block, length)
        tokens = tuple(token(speed, x) for speed in speeds)
        assert None not in tokens
        boundary_states.append(normalized_state(tokens))
    assert boundary_states == [REDUCED_STATE] * (height + 1)

    counts, defects = prefix_legality(REDUCED_STATE, PACKET_WORD)
    assert counts == D
    assert defects == (0,) * 8
    assert tuple(counts[a] * STEPS[a] % P for a in range(8)) == (1,) * 8

    core_clearances = tuple(core_minimum(core_speed) for core_speed in CORE)
    assert core_clearances == (Q(2, 15), Q(4, 15), Q(3, 7), Q(2, 7), Q(1, 7))
    assert min(core_clearances) > Q(1, 14)
    full_speeds = tuple(7 * core_speed for core_speed in CORE) + speeds
    assert len(full_speeds) == 13 == len(set(full_speeds))

    first_wall = actual[0][0]
    last_wall = actual[-1][0]
    wall_extent = last_wall - first_wall
    assert wall_extent == Q(2 * height - 1, 210 * height + 2)
    mesh_bound = Q(1, speeds[1]) + Q(2, speeds[7])
    if height >= 3:
        assert wall_extent > mesh_bound

    return (
        height,
        length,
        speeds,
        len(actual),
        wall_extent,
        mesh_bound,
        core_clearances,
        boundary_states[0],
    )


def multiset_words(counts: tuple[int, ...]):
    word = [0] * sum(counts)
    remaining = list(counts)

    def visit(position: int):
        if position == len(word):
            yield tuple(word)
            return
        for owner in range(len(remaining)):
            if remaining[owner]:
                remaining[owner] -= 1
                word[position] = owner
                yield from visit(position + 1)
                remaining[owner] += 1

    yield from visit(0)


def fixed_state_legal(word: tuple[int, ...]) -> bool:
    total = 0
    counts = [0] * 8
    for owner in word:
        if (REDUCED_STATE[owner] + total - counts[owner] * STEPS[owner]) % P:
            return False
        total = (total + STEPS[owner]) % P
        counts[owner] += 1
    return all(
        (total - counts[owner] * STEPS[owner]) % P == 0
        for owner in range(8)
    )


def supportable(word: tuple[int, ...]) -> bool:
    total = 0
    counts = [0] * 8
    requirements: dict[int, int] = {}
    capacity = Counter({0: 2, 1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1})
    for owner in word:
        required = (counts[owner] * STEPS[owner] - total) % P
        if owner in requirements and requirements[owner] != required:
            return False
        requirements[owner] = required
        total = (total + STEPS[owner]) % P
        counts[owner] += 1
    used = Counter(requirements.values())
    return all(used[value] <= capacity[value] for value in used)


def first_occurrence_word(word: tuple[int, ...]) -> tuple[int, ...]:
    seen = set()
    first = []
    for owner in word:
        if owner not in seen:
            seen.add(owner)
            first.append(owner)
    return tuple(first)


def word_census() -> dict:
    total = fixed = supported = same_first = same_first_fixed = 0
    fixed_words = []
    for word in multiset_words(D):
        total += 1
        is_fixed = fixed_state_legal(word)
        fixed += is_fixed
        supported += supportable(word)
        if is_fixed:
            fixed_words.append(word)
        if first_occurrence_word(word) == FIRST_OWNER_WORD:
            same_first += 1
            same_first_fixed += is_fixed
    assert total == 181_440
    assert fixed == 3
    assert supported == 45_360
    assert same_first == 8
    assert same_first_fixed == 1
    assert tuple(sorted(fixed_words)) == (
        (0, 6, 4, 1, 2, 5, 3, 7, 7),
        PACKET_WORD,
        (7, 7, 6, 4, 1, 2, 5, 3, 0),
    )
    return {
        "multiplicity_words": total,
        "fixed_state_legal": fixed,
        "some_state_supportable": supported,
        "same_first_occurrence": same_first,
        "same_first_occurrence_legal": same_first_fixed,
        "fixed_words": tuple(sorted(fixed_words)),
    }


def tournament_fingerprint(adjacency: tuple[tuple[bool, ...], ...]) -> dict:
    size = len(adjacency)
    vertices = tuple(range(size))
    scores = tuple(sum(row) for row in adjacency)
    triangles = 0
    for a, b, c in combinations(vertices, 3):
        triangles += (
            adjacency[a][b] and adjacency[b][c] and adjacency[c][a]
        ) or (
            adjacency[a][c] and adjacency[c][b] and adjacency[b][a]
        )

    def reachable(source: int) -> set[int]:
        seen = {source}
        stack = [source]
        while stack:
            a = stack.pop()
            for b in vertices:
                if adjacency[a][b] and b not in seen:
                    seen.add(b)
                    stack.append(b)
        return seen

    remaining = set(vertices)
    sccs = []
    while remaining:
        a = min(remaining)
        forward = reachable(a)
        component = {b for b in remaining if b in forward and a in reachable(b)}
        sccs.append(len(component))
        remaining -= component

    paths = tuple(
        path
        for path in permutations(vertices)
        if all(adjacency[path[i]][path[i + 1]] for i in range(size - 1))
    )
    return {
        "scores": scores,
        "score_histogram": tuple(sorted(Counter(scores).items())),
        "directed_triangles": triangles,
        "scc_sizes": tuple(sorted(sccs, reverse=True)),
        "hamiltonian_paths": len(paths),
        "tie_hamiltonian_path": min(paths),
    }


def tournament_audit() -> dict:
    first_position = {owner: position for position, owner in enumerate(FIRST_OWNER_WORD)}
    marked = [[False] * 8 for _ in range(8)]
    circular = [[False] * 8 for _ in range(8)]
    for a, b in combinations(range(8), 2):
        # The only equal-sheet pair is resolved by first wall occurrence.
        if (REDUCED_STATE[a], first_position[a]) < (
            REDUCED_STATE[b],
            first_position[b],
        ):
            marked[a][b] = True
        else:
            marked[b][a] = True

        difference = (REDUCED_STATE[b] - REDUCED_STATE[a]) % P
        if difference in (1, 2, 3):
            circular[a][b] = True
        elif difference in (4, 5, 6):
            circular[b][a] = True
        elif first_position[a] < first_position[b]:
            circular[a][b] = True
        else:
            circular[b][a] = True

    marked_t = tuple(tuple(row) for row in marked)
    circular_t = tuple(tuple(row) for row in circular)
    chronological_t = tuple(
        tuple(a != b and a < b for b in range(9)) for a in range(9)
    )
    marked_fp = tournament_fingerprint(marked_t)
    circular_fp = tournament_fingerprint(circular_t)
    chronological_fp = tournament_fingerprint(chronological_t)
    edge_flips = sum(
        marked_t[a][b] != circular_t[a][b]
        for a, b in combinations(range(8), 2)
    )
    assert marked_fp["score_histogram"] == tuple((score, 1) for score in range(8))
    assert marked_fp["directed_triangles"] == 0
    assert marked_fp["scc_sizes"] == (1,) * 8
    assert marked_fp["hamiltonian_paths"] == 1
    assert marked_fp["tie_hamiltonian_path"] == (7, 0, 3, 5, 2, 1, 4, 6)
    assert circular_fp["score_histogram"] == ((3, 4), (4, 4))
    assert circular_fp["directed_triangles"] == 20
    assert circular_fp["scc_sizes"] == (8,)
    assert circular_fp["hamiltonian_paths"] == 629
    assert circular_fp["tie_hamiltonian_path"] == (0, 2, 1, 4, 6, 7, 3, 5)
    assert edge_flips == 9
    assert chronological_fp["score_histogram"] == tuple(
        (score, 1) for score in range(9)
    )
    assert chronological_fp["directed_triangles"] == 0
    assert chronological_fp["scc_sizes"] == (1,) * 9
    assert chronological_fp["hamiltonian_paths"] == 1
    assert chronological_fp["tie_hamiltonian_path"] == tuple(range(9))
    return {
        "chronological": chronological_fp,
        "marked": marked_fp,
        "circular": circular_fp,
        "edge_flips": edge_flips,
        "equal_sheet_tie": (7, 0),
    }


def main() -> None:
    rows = tuple(audit_height(height) for height in range(H_MIN, H_MAX + 1))
    census = word_census()
    tournaments = tournament_audit()
    canonical = "\n".join(repr(row) for row in rows) + "\n"
    digest = sha256(canonical.encode()).hexdigest()
    if EXPECTED_DIGEST != "TO_BE_FILLED":
        assert digest == EXPECTED_DIGEST

    print("NONCENTRAL DIAGONAL-RETURN PACKETS -- EXACT REPLAY")
    print("family: L=105H; W=(L+1,L+22,L-20,L-6,L+15,L-13,L+8,2L+2)")
    print("fixed blocked/core-safe interval: [2/15,1/7], length=1/105")
    print("packet word: (7,2,5,3,0,6,4,1,7)")
    print("packet counts: (1,1,1,1,1,1,1,2)")
    print("inverse steps: (1,1,1,1,1,1,1,4)")
    print("count*step mod 7: (1,1,1,1,1,1,1,1); reduced return=0")
    print("reduced boundary state: (0,4,3,1,5,2,6,0)")
    print()
    print(f"audited H range: {H_MIN}..{H_MAX} ({len(rows)} families)")
    print(f"complete event-list comparisons: {len(rows)} passed")
    print(f"covered walls checked: {sum(row[3] for row in rows)}")
    print(f"covered chambers checked: {sum(row[3] + 1 for row in rows)}")
    print("wall order inequalities: 84M>11L-10 and 14M<2L-5")
    print("block range: M=14H,...,15H-1; packets=H; walls=9H")
    print("core P=(1,2,11,12,13); exact minima=(2/15,4/15,3/7,2/7,1/7)")
    print("full LRC speed row: 7P union W; 13 distinct nonzero speeds")
    print("ceil(f/g)=2 for every H; wall extent beats 1/g+2/f for H>=3")
    print()
    for height in (1, 3, 10, 100, 300):
        row = rows[height - 1]
        print(
            f"H={height}: L={row[1]} W={row[2]} walls={row[3]} "
            f"wall_extent={row[4]} mesh_bound={row[5]} "
            f"refutes={row[4] > row[5]}"
        )
    print()
    print("PREFIX-LEGALITY CENSUS")
    print(f"multiplicity words: {census['multiplicity_words']}")
    print(f"legal from displayed reduced state: {census['fixed_state_legal']}")
    print(f"supportable from some covered state: {census['some_state_supportable']}")
    print(f"same first-owner tournament/order: {census['same_first_occurrence']}")
    print(
        "same first-owner order and legal from displayed state: "
        f"{census['same_first_occurrence_legal']}"
    )
    print(f"three fixed-state legal words: {census['fixed_words']}")
    print()
    print("TOURNAMENT FINGERPRINTS")
    print(
        "nine-event chronology: score histogram "
        f"{tournaments['chronological']['score_histogram']}; triangles "
        f"{tournaments['chronological']['directed_triangles']}; SCCs "
        f"{tournaments['chronological']['scc_sizes']}; Hamiltonian paths "
        f"{tournaments['chronological']['hamiltonian_paths']}"
    )
    print(f"marked score histogram: {tournaments['marked']['score_histogram']}")
    print(f"marked directed triangles: {tournaments['marked']['directed_triangles']}")
    print(f"marked SCC sizes: {tournaments['marked']['scc_sizes']}")
    print(f"marked Hamiltonian paths: {tournaments['marked']['hamiltonian_paths']}")
    print(f"marked tie path: {tournaments['marked']['tie_hamiltonian_path']}")
    print(f"circular score histogram: {tournaments['circular']['score_histogram']}")
    print(f"circular directed triangles: {tournaments['circular']['directed_triangles']}")
    print(f"circular SCC sizes: {tournaments['circular']['scc_sizes']}")
    print(f"circular Hamiltonian paths: {tournaments['circular']['hamiltonian_paths']}")
    print(f"circular tie path: {tournaments['circular']['tie_hamiltonian_path']}")
    print(f"marked/circular edge flips: {tournaments['edge_flips']}")
    print("preserved by owner multiplicity: reduced return congruence only")
    print("destroyed by owner multiplicity: prefix legality and raw event insertion order")
    print(f"canonical_sha256={digest}")
    print("SUMMARY: ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
