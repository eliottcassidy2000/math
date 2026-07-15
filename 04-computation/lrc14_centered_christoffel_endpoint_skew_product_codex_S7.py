#!/usr/bin/env python3
"""Exact audit for THM-778's centered-Christoffel endpoint skew product.

The script treats endpoint events, rather than runners, as the primary base
objects.  For two owner clocks it verifies the centered mechanical-word
formula, its Euclidean parity cocycle, and the exact odd/odd tie count.  For a
whole speed family it reconstructs the global simultaneous-wall order from
pairwise centered Beatty ranks, then drives the prime-seven token fibre by the
THM-773 translations.

Tournament Analysis uses owner clocks as vertices only after this richer base
has been built.  The pairwise observable is which owner has the earlier next
endpoint.  Every instantaneous tournament is transitive; its labelled
Hamiltonian path changes through the centered-Christoffel word.  Thus the
ordinary isomorphism class deliberately forgets the theorem's payload.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import floor, gcd
from pathlib import Path

from lrc14_prime7_sheet_monodromy_metagraph_codex_S6 import (
    ALL_PERMS,
    circular_tournament,
    encode_fixed_path,
    least_hamiltonian_path,
)


P = 7


def ceil_fraction(x: Fraction) -> int:
    return -floor(-x)


def continued_fraction(num: int, den: int) -> tuple[int, ...]:
    """Simple continued fraction of num/den, including a possible a_0."""
    assert num >= 0 and den > 0
    digits = []
    while den:
        a, num, den = num // den, den, num % den
        digits.append(a)
    return tuple(digits)


def from_continued_fraction(digits: tuple[int, ...]) -> tuple[int, int]:
    num, den = 1, 0
    for a in reversed(digits):
        num, den = a * num + den, num
    return num, den


def centered_count(u: int, v: int, i: int) -> int:
    """Number of v-midpoints strictly before u's i-th midpoint."""
    assert u > 0 and v > 0 and 0 <= i < u
    return ceil_fraction(Fraction(v * (2 * i + 1) - u, 2 * u))


def centered_height(num: int, den: int, phase: int, i: int) -> int:
    """F^phase_(num,den)(i), the two-coset centered Beatty height."""
    assert num >= 0 and den > 0 and phase in (0, 1)
    return ceil_fraction(Fraction(num * i, den) + Fraction(num, 2 * den) - Fraction(phase, 2))


def centered_increment_word(num: int, den: int, phase: int = 1) -> tuple[int, ...]:
    return tuple(
        centered_height(num, den, phase, i + 1)
        - centered_height(num, den, phase, i)
        for i in range(den)
    )


def pair_blocks_direct(u: int, v: int) -> tuple[str, ...]:
    events: dict[Fraction, list[str]] = defaultdict(list)
    for i in range(u):
        events[Fraction(2 * i + 1, 2 * u)].append("u")
    for j in range(v):
        events[Fraction(2 * j + 1, 2 * v)].append("v")
    return tuple("".join(events[x]) for x in sorted(events))


def pair_blocks_formula(u: int, v: int) -> tuple[str, ...]:
    blocks: list[str] = []
    used_v = 0
    for i in range(u):
        count = centered_count(u, v, i)
        while used_v < count:
            blocks.append("v")
            used_v += 1
        tied = used_v < v and v * (2 * i + 1) == u * (2 * used_v + 1)
        if tied:
            blocks.append("uv")
            used_v += 1
        else:
            blocks.append("u")
    blocks.extend("v" for _ in range(v - used_v))
    return tuple(blocks)


def run_length_word(blocks: tuple[str, ...]) -> tuple[tuple[str, int], ...]:
    if not blocks:
        return ()
    runs: list[tuple[str, int]] = []
    symbol = blocks[0]
    length = 1
    for block in blocks[1:]:
        if block == symbol:
            length += 1
        else:
            runs.append((symbol, length))
            symbol, length = block, 1
    runs.append((symbol, length))
    return tuple(runs)


def pair_packet(u: int, v: int) -> dict:
    g = gcd(u, v)
    p, q = u // g, v // g
    blocks = pair_blocks_formula(u, v)
    reduced = pair_blocks_formula(p, q)
    assert blocks == reduced * g
    digits = continued_fraction(q, p)
    assert from_continued_fraction(digits) == (q, p)
    ties = sum(block == "uv" for block in blocks)
    expected_ties = g if p % 2 == q % 2 == 1 else 0
    assert ties == expected_ties
    a, r = divmod(q, p)
    increments = centered_increment_word(q, p)
    assert set(increments) <= {a, a + bool(r)}
    assert sum(value - a for value in increments) == r
    return {
        "u": u,
        "v": v,
        "gcd": g,
        "reduced_ratio_v_over_u": [q, p],
        "continued_fraction": digits,
        "euclid_depth": len(digits),
        "tie_blocks": ties,
        "block_count": len(blocks),
        "run_count": len(run_length_word(blocks)),
        "runs": run_length_word(blocks),
    }


def pair_audit(limit: int) -> dict:
    pairs = formula_failures = repetition_failures = tie_failures = 0
    balanced_failures = phase_failures = 0
    cf_depths: Counter[int] = Counter()
    run_counts: Counter[int] = Counter()
    for u in range(1, limit + 1):
        for v in range(1, limit + 1):
            pairs += 1
            direct = pair_blocks_direct(u, v)
            formula = pair_blocks_formula(u, v)
            formula_failures += direct != formula
            g = gcd(u, v)
            p, q = u // g, v // g
            repetition_failures += formula != pair_blocks_formula(p, q) * g
            predicted_ties = g if p % 2 == q % 2 == 1 else 0
            tie_failures += sum(block == "uv" for block in formula) != predicted_ties

            digits = continued_fraction(q, p)
            assert from_continued_fraction(digits) == (q, p)
            cf_depths[len(digits)] += 1
            run_counts[len(run_length_word(formula))] += 1

            increments = centered_increment_word(q, p)
            for length in range(1, p + 1):
                doubled = increments + increments
                sums = [sum(doubled[start : start + length]) for start in range(p)]
                if max(sums) - min(sums) > 1:
                    balanced_failures += 1
                    break

            for phase in (0, 1):
                a, r = divmod(q, p)
                next_phase = (a - phase) % 2
                offset = (a - phase + next_phase) // 2
                for i in range(-p, 2 * p + 1):
                    lhs = centered_height(q, p, phase, i)
                    rhs = a * i + offset + centered_height(r, p, next_phase, i)
                    if lhs != rhs:
                        phase_failures += 1
                        break

    assert not any(
        (formula_failures, repetition_failures, tie_failures, balanced_failures, phase_failures)
    )
    return {
        "ordered_pairs": pairs,
        "limit": limit,
        "formula_failures": formula_failures,
        "gcd_repetition_failures": repetition_failures,
        "tie_formula_failures": tie_failures,
        "cyclic_balance_failures": balanced_failures,
        "euclidean_phase_cocycle_failures": phase_failures,
        "cf_depth_histogram": dict(sorted(cf_depths.items())),
        "run_count_range": [min(run_counts), max(run_counts)],
    }


def token(w: int, x: Fraction) -> int:
    return (-pow(w, -1, P) * floor(w * x + Fraction(1, 2))) % P


def root_polynomial(position: tuple[int, ...]) -> tuple[int, ...]:
    coeff = [1]
    for k in position:
        nxt = [0] * (len(coeff) + 1)
        for i, value in enumerate(coeff):
            nxt[i] = (nxt[i] - k * value) % P
            nxt[i + 1] = (nxt[i + 1] + value) % P
        coeff = nxt
    return tuple(coeff)


def polynomial_remainder(coefficients: tuple[int, ...]) -> tuple[int, ...]:
    work = list(coefficients)
    for degree in range(len(work) - 1, P - 1, -1):
        lead = work[degree] % P
        work[degree] = 0
        work[degree - P + 1] = (work[degree - P + 1] + lead) % P
    return tuple((work[i] if i < len(work) else 0) % P for i in range(P))


def covered(position: tuple[int, ...]) -> bool:
    return not any(polynomial_remainder(root_polynomial(position)))


def event_rank(speeds: tuple[int, ...], owner: int, local_index: int) -> int:
    w = speeds[owner]
    return local_index + sum(
        centered_count(w, other, local_index)
        for b, other in enumerate(speeds)
        if b != owner
    )


def direct_event_blocks(speeds: tuple[int, ...]) -> tuple[tuple[Fraction, tuple[tuple[int, int], ...]], ...]:
    events: dict[Fraction, list[tuple[int, int]]] = defaultdict(list)
    for owner, w in enumerate(speeds):
        for i in range(w):
            events[Fraction(2 * i + 1, 2 * w)].append((owner, i))
    return tuple((x, tuple(events[x])) for x in sorted(events))


def ranked_event_blocks(speeds: tuple[int, ...]) -> tuple[tuple[int, tuple[tuple[int, int], ...]], ...]:
    ranks: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for owner, w in enumerate(speeds):
        for i in range(w):
            ranks[event_rank(speeds, owner, i)].append((owner, i))
    return tuple((rank, tuple(ranks[rank])) for rank in sorted(ranks))


def next_event_path(speeds: tuple[int, ...], next_indices: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sorted(
            range(len(speeds)),
            key=lambda a: (Fraction(2 * next_indices[a] + 1, 2 * speeds[a]), a),
        )
    )


def path_edge_flips(first: tuple[int, ...], second: tuple[int, ...]) -> int:
    rank_first = {owner: i for i, owner in enumerate(first)}
    rank_second = {owner: i for i, owner in enumerate(second)}
    return sum(
        (rank_first[a] < rank_first[b]) != (rank_second[a] < rank_second[b])
        for a in range(len(first))
        for b in range(a + 1, len(first))
    )


def family_movie(speeds: tuple[int, ...]) -> dict:
    assert all(w % P for w in speeds)
    direct = direct_event_blocks(speeds)
    ranked = ranked_event_blocks(speeds)
    rank_times: dict[int, Fraction] = {}
    for cumulative, (x, events) in zip((rank for rank, _ in ranked), direct):
        rank_times[cumulative] = x
        assert tuple(sorted(events)) == tuple(sorted(dict(ranked)[cumulative]))
        assert cumulative == sum(len(old_events) for old_x, old_events in direct if old_x < x)
    assert len(rank_times) == len(direct)
    total_events = sum(speeds)
    mirror_failures = []
    rank_by_events = {frozenset(events): rank for rank, events in ranked}
    for block_index, (x, events) in enumerate(direct):
        mirror_events = frozenset((a, speeds[a] - 1 - i) for a, i in events)
        mirror_x, observed_mirror = direct[-1 - block_index]
        rank = rank_by_events[frozenset(events)]
        mirror_rank = rank_by_events[mirror_events]
        if (
            mirror_x != 1 - x
            or frozenset(observed_mirror) != mirror_events
            or rank + mirror_rank != total_events - len(events)
        ):
            mirror_failures.append(block_index)
    assert not mirror_failures

    position = [0] * len(speeds)
    next_indices = [0] * len(speeds)
    paths = []
    exact_chambers = []
    covered_walls = []
    transition_failures = []
    for block_index, (x, events) in enumerate(direct):
        event_owners = {e[0] for e in events}
        wall_position = tuple(k for a, k in enumerate(position) if a not in event_owners)
        if covered(wall_position):
            record = {
                "x": x,
                "block_index": block_index,
                "rank": rank_by_events[frozenset(events)],
                "owner_indices": tuple(a for a, _ in events),
                "owner_speeds": tuple(speeds[a] for a, _ in events),
                "remaining_tokens": wall_position,
            }
            if len(events) == 1 and len(wall_position) == P:
                owner = events[0][0]
                adjacency = circular_tournament(wall_position)
                path = least_hamiltonian_path(adjacency)
                record.update(
                    {
                        "a267_mask": encode_fixed_path(adjacency, path),
                        "least_hamiltonian_path": path,
                        "duplicate_sheet_before": position[owner],
                        "duplicate_sheet_after": (
                            position[owner] - pow(speeds[owner], -1, P)
                        )
                        % P,
                    }
                )
            covered_walls.append(record)

        for owner, local_index in events:
            assert next_indices[owner] == local_index
            position[owner] = (position[owner] - pow(speeds[owner], -1, P)) % P
            next_indices[owner] += 1

        right = direct[block_index + 1][0] if block_index + 1 < len(direct) else Fraction(1)
        midpoint = (x + right) / 2
        expected = tuple(token(w, midpoint) for w in speeds)
        if tuple(position) != expected:
            transition_failures.append(block_index)
        if covered(tuple(position)):
            exact_chambers.append(
                {
                    "left": x,
                    "right": right,
                    "position": tuple(position),
                }
            )
        paths.append(next_event_path(speeds, tuple(next_indices)))

    assert tuple(position) == (P - 1,) * len(speeds)
    assert not transition_failures
    assert len(paths) == len(direct)
    edge_flips = [
        path_edge_flips(paths[i - 1], paths[i])
        for i in range(len(paths))
    ]
    covered_owner_word = tuple(
        tuple(wall["owner_speeds"])
        for wall in covered_walls
    )
    assert covered_owner_word == tuple(reversed(covered_owner_word))
    covered_mask_word = tuple(
        wall.get("a267_mask")
        for wall in covered_walls
    )
    covered_duplicate_word = tuple(
        (wall.get("duplicate_sheet_before"), wall.get("duplicate_sheet_after"))
        for wall in covered_walls
    )
    covered_wall_gaps = []
    for left_wall, right_wall in zip(covered_walls, covered_walls[1:]):
        left_index = left_wall["block_index"]
        right_index = right_wall["block_index"]
        between = direct[left_index + 1 : right_index]
        owner_counts = [0] * len(speeds)
        for _, block_events in between:
            for owner, _ in block_events:
                owner_counts[owner] += 1
        covered_wall_gaps.append(
            {
                "left": left_wall["x"],
                "right": right_wall["x"],
                "wall_blocks": len(between),
                "individual_events": sum(owner_counts),
                "owner_counts": tuple(owner_counts),
                "owner_counts_mod7": tuple(count % P for count in owner_counts),
            }
        )

    pair_packets = [
        pair_packet(speeds[a], speeds[b])
        for a in range(len(speeds))
        for b in range(a + 1, len(speeds))
    ]
    rank_digest = hashlib.sha256(
        json.dumps(
            [
                [rank, [[a, i] for a, i in events]]
                for rank, events in ranked
            ],
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    return {
        "speeds": speeds,
        "individual_events": total_events,
        "wall_blocks": len(direct),
        "simultaneous_blocks": sum(len(events) > 1 for _, events in direct),
        "max_block_size": max(len(events) for _, events in direct),
        "rank_reconstruction_failures": 0,
        "mirror_rank_failures": mirror_failures,
        "rank_digest": rank_digest,
        "exact_chamber_count": len(exact_chambers),
        "exact_chambers": exact_chambers,
        "covered_wall_count": len(covered_walls),
        "covered_walls": covered_walls,
        "covered_wall_owner_word": covered_owner_word,
        "covered_wall_mask_word": covered_mask_word,
        "covered_wall_duplicate_word": covered_duplicate_word,
        "covered_wall_gaps": covered_wall_gaps,
        "covered_wall_palindrome": True,
        "global_carry": tuple(position),
        "transition_failures": transition_failures,
        "next_event_tournament": {
            "iso_class": "transitive",
            "score_histogram": {i: 1 for i in range(len(speeds))},
            "directed_3cycles": 0,
            "scc_sizes": (1,) * len(speeds),
            "hamiltonian_paths_per_chamber": 1,
            "distinct_labelled_paths": len(set(paths)),
            "total_edge_flips": sum(edge_flips),
            "edge_flip_histogram": dict(sorted(Counter(edge_flips).items())),
            "tie_hamiltonian_path": tuple(range(len(speeds))),
        },
        "pair_cf_depth_histogram": dict(
            sorted(Counter(packet["euclid_depth"] for packet in pair_packets).items())
        ),
        "pair_tie_blocks": sum(packet["tie_blocks"] for packet in pair_packets),
        "pair_run_count_range": [
            min(packet["run_count"] for packet in pair_packets),
            max(packet["run_count"] for packet in pair_packets),
        ],
    }


def reflection_mask_audit() -> dict:
    """Does owner-labelled reflection descend to the 25 mask indices?"""
    images: dict[int, set[int]] = defaultdict(set)
    for position in ALL_PERMS:
        adjacency = circular_tournament(position)
        path = least_hamiltonian_path(adjacency)
        mask = encode_fixed_path(adjacency, path)
        reflected = tuple(P - 1 - k for k in position)
        reflected_adjacency = circular_tournament(reflected)
        reflected_path = least_hamiltonian_path(reflected_adjacency)
        reflected_mask = encode_fixed_path(reflected_adjacency, reflected_path)
        images[mask].add(reflected_mask)
    image_sizes = Counter(len(targets) for targets in images.values())
    assert len(images) == 25
    assert sum(size == 1 for size in map(len, images.values())) == 9
    assert max(map(len, images.values())) == 7
    return {
        "assignments": len(ALL_PERMS),
        "source_masks": len(images),
        "single_valued_source_masks": sum(len(targets) == 1 for targets in images.values()),
        "maximum_reflected_images": max(map(len, images.values())),
        "image_size_histogram": dict(sorted(image_sizes.items())),
        "images": {mask: tuple(sorted(targets)) for mask, targets in sorted(images.items())},
        "descends_to_mask_function": False,
    }


def render(
    pair_result: dict,
    movies: list[dict],
    sample_pairs: list[dict],
    mask_reflection: dict,
) -> str:
    lines = [
        "THM-778 CENTERED-CHRISTOFFEL ENDPOINT SKEW PRODUCT",
        "=" * 78,
        "PAIRWISE EXACT AUDIT",
        f"  ordered pairs u,v <= {pair_result['limit']}: {pair_result['ordered_pairs']}",
        f"  direct/formula failures: {pair_result['formula_failures']}",
        f"  gcd repetition failures: {pair_result['gcd_repetition_failures']}",
        f"  odd/odd tie-count failures: {pair_result['tie_formula_failures']}",
        f"  cyclic-balance failures: {pair_result['cyclic_balance_failures']}",
        f"  Euclidean parity-cocycle failures: {pair_result['euclidean_phase_cocycle_failures']}",
        f"  CF-depth histogram: {pair_result['cf_depth_histogram']}",
        "",
        "SAMPLE CENTERED WORD PACKETS",
    ]
    for packet in sample_pairs:
        lines.append(
            "  "
            f"({packet['u']},{packet['v']}): gcd={packet['gcd']} "
            f"v/u={packet['reduced_ratio_v_over_u'][0]}/{packet['reduced_ratio_v_over_u'][1]} "
            f"CF={tuple(packet['continued_fraction'])} ties={packet['tie_blocks']} "
            f"blocks={packet['block_count']} runs={packet['run_count']}"
        )
        runs = tuple(tuple(x) for x in packet["runs"])
        if len(runs) > 20:
            run_display = runs[:8] + (("...", len(runs) - 16),) + runs[-8:]
        else:
            run_display = runs
        lines.append(f"    run word={run_display}")

    lines.extend(["", "GLOBAL RANK RECONSTRUCTION + PRIME-SEVEN FIBRE"])
    for movie in movies:
        lines.append(
            "  "
            f"W={tuple(movie['speeds'])}: individual={movie['individual_events']} "
            f"walls={movie['wall_blocks']} simultaneous={movie['simultaneous_blocks']} "
            f"max_block={movie['max_block_size']} rank_failures={movie['rank_reconstruction_failures']} "
            f"mirror_failures={movie['mirror_rank_failures']}"
        )
        lines.append(
            "    "
            f"exact_chambers={movie['exact_chamber_count']} covered_walls={movie['covered_wall_count']} "
            f"carry={tuple(movie['global_carry'])} transition_failures={movie['transition_failures']}"
        )
        if movie["covered_walls"]:
            lines.append(
                "    "
                f"covered-wall owner word={tuple(movie['covered_wall_owner_word'])} "
                f"palindrome={movie['covered_wall_palindrome']}"
            )
            lines.append(
                "    "
                f"covered-wall a267 mask word={tuple(movie['covered_wall_mask_word'])}"
            )
            lines.append(
                "    "
                f"covered-wall duplicate-root word={tuple(tuple(x) for x in movie['covered_wall_duplicate_word'])}"
            )
            lines.append(
                "    "
                f"between-wall block word={tuple(gap['wall_blocks'] for gap in movie['covered_wall_gaps'])}"
            )
            lines.append(
                "    "
                f"between-wall event word={tuple(gap['individual_events'] for gap in movie['covered_wall_gaps'])}"
            )
        lines.append(
            "    "
            f"rank_sha256={movie['rank_digest']} CF-depths={movie['pair_cf_depth_histogram']} "
            f"pair_ties={movie['pair_tie_blocks']} pair_runs={movie['pair_run_count_range']}"
        )
        tournament = movie["next_event_tournament"]
        lines.append(
            "    next-event tournament: "
            f"iso={tournament['iso_class']} labelled_paths={tournament['distinct_labelled_paths']} "
            f"edge_flips={tournament['total_edge_flips']} "
            f"flip_hist={tournament['edge_flip_histogram']} c3=0 SCC=singletons Hpaths=1"
        )
        for wall in movie["covered_walls"][:8]:
            lines.append(
                "    covered wall "
                f"x={wall['x']} owners={tuple(wall['owner_speeds'])} "
                f"tokens={tuple(wall['remaining_tokens'])} "
                f"mask={wall.get('a267_mask')} "
                f"duplicate={wall.get('duplicate_sheet_before')}->{wall.get('duplicate_sheet_after')}"
            )

    lines.extend(
        [
            "",
            "REFLECTION VERSUS THE 25-MASK QUOTIENT",
            f"  assignments audited: {mask_reflection['assignments']}",
            f"  masks with unique reflected image: {mask_reflection['single_valued_source_masks']}/25",
            f"  maximum reflected images from one mask: {mask_reflection['maximum_reflected_images']}",
            f"  image-size histogram: {mask_reflection['image_size_histogram']}",
            f"  descends to a mask function: {mask_reflection['descends_to_mask_function']}",
            "",
            "PRESERVATION AUDIT",
            "  exact base carrier: owner-local event + centered Beatty rank + simultaneous block",
            "  exact fibre carrier: owner token + inverse step mod 7 + cover polynomial",
            "  CF role: recursive address of each pairwise midpoint cutting sequence",
            "  iso-node loss: every next-event tournament is transitive, so the node forgets the entire labelled word",
            "  challenged vertices: endpoint events and Euclidean cells carry more proof data than runner vertices",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=80)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    pair_result = pair_audit(args.limit)
    families = [
        (1, 4, 5, 6, 8, 9, 10),
        (12, 38, 72, 96, 151, 169, 188),
        (108, 169, 143, 213, 206, 197, 30, 162),
        (1, 2, 3, 4, 5, 8, 10),
        (1, 2, 3, 4, 5, 8, 11),
    ]
    movies = [family_movie(speeds) for speeds in families]
    seed_wall = next(
        wall
        for wall in movies[2]["covered_walls"]
        if wall["x"] == Fraction(19, 216)
    )
    assert seed_wall["owner_speeds"] == (108,)
    assert movies[0]["exact_chamber_count"] == 2
    assert movies[1]["exact_chamber_count"] == 6
    duplicate_word = movies[2]["covered_wall_duplicate_word"]
    assert duplicate_word[:5] == duplicate_word[5:]
    gap_blocks = tuple(gap["wall_blocks"] for gap in movies[2]["covered_wall_gaps"])
    gap_events = tuple(gap["individual_events"] for gap in movies[2]["covered_wall_gaps"])
    gap_counts = tuple(gap["owner_counts"] for gap in movies[2]["covered_wall_gaps"])
    assert gap_blocks == tuple(reversed(gap_blocks))
    assert gap_events == tuple(reversed(gap_events))
    assert gap_counts == tuple(reversed(gap_counts))

    sample_pairs = [
        pair_packet(1, 10),
        pair_packet(4, 9),
        pair_packet(12, 188),
        pair_packet(108, 162),
    ]
    mask_reflection = reflection_mask_audit()
    output = render(pair_result, movies, sample_pairs, mask_reflection)
    print(output, end="")
    if args.out:
        args.out.write_text(output)
    if args.json:
        payload = {
            "schema_version": 1,
            "theorem": "THM-778",
            "pair_audit": pair_result,
            "sample_pairs": sample_pairs,
            "movies": movies,
            "reflection_mask_audit": mask_reflection,
        }
        args.json.write_text(json.dumps(payload, separators=(",", ":"), default=str) + "\n")


if __name__ == "__main__":
    main()
