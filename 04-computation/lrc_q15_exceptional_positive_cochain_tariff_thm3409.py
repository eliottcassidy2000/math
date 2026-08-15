#!/usr/bin/env python3
"""Exact companion for THM-3409: the exceptional q=15 cochain tariff.

The owner edge U=(1,7,8,9,11,13) is the unique rank-six q=15 physical
edge outside the capped common-mode-centre locus in the pinned precursor
audits.  This companion independently resolves its entire physical-time
locus into strict event cells.  At every representative time it enumerates
every THM-3398 selected mode whose open interval contains that time, rather
than accepting only one witness returned by a search.

All arithmetic is integral or Fraction-exact.  Runtime gates survive
``python -O``.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import itertools
from collections import Counter
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODE_PATH = ROOT / "04-computation/lrc_general_finite_mode_sheet_cover_cochain_thm3398.py"
PINS = (
    (
        "THM-3398",
        ROOT / "01-canon/theorems/THM-3398-general-finite-mode-sheet-cover-cochain.md",
        "01901da2bb382184cfe4466550afe79255598f580f00a761fc32731a52ec9378",
    ),
    (
        "THM-3398-script",
        MODE_PATH,
        "82929cbf6903701533c1b1f6ebed143e5c8f9edc570dfe2895cf8db70e478da9",
    ),
    (
        "THM-3398-output",
        ROOT / "05-knowledge/results/lrc_general_finite_mode_sheet_cover_cochain_thm3398.out",
        "ab25331039813f8c83626a66d0d0d8157e8b3826a76fccc690452a2cdad3169b",
    ),
    (
        "THM-3405",
        ROOT / "01-canon/theorems/THM-3405-common-centre-gcd-gauge-and-boolean-half-twist.md",
        "d3e7dbeeb85c6f897bd9e31270bd0b6602ae4feac3b46a45eb5ce23ae5d24fe0",
    ),
    (
        "mobile-boundary-script",
        ROOT / "04-computation/lrc_mobile_common_centre_crt_rank_probe_20260815.py",
        "dd157efbd3bd7da34b75cba7f30fb55cfa0381bb6a74ece7abade7f4a2c439fa",
    ),
    (
        "mobile-boundary-output",
        ROOT / "05-knowledge/results/lrc_mobile_common_centre_crt_rank_probe_20260815.out",
        "b874f11a8f81604e8bfefe01a6437bfd10a21b909c3ee3e37d0c8a714528594f",
    ),
    (
        "full-physical-script",
        ROOT / "04-computation/lrc_q8_q15_full_physical_clutter_audit_20260815.py",
        "e54b77eeae05484cbbfacd904e850815f7e78a5e3306f21ad87a68ffbfae9e2e",
    ),
    (
        "full-physical-output",
        ROOT / "05-knowledge/results/lrc_q8_q15_full_physical_clutter_audit_20260815.out",
        "8b8cc6f45ab8b14b0ba26afb29748ccf3bc08f4f103581ba9afd7167391e8008",
    ),
)

Q = 15
SPEEDS = (1, 7, 8, 9, 11, 13)
EXPECTED_POSITIVE_BLOCKS = (
    (0, 1),
    (5, 7),
    (8, 10, 12),
    (4, 9, 14),
    (2, 6, 13),
    (3, 11),
)
EXPECTED_POSITIVE_COCHAIN = (
    2, 2, 3, 3, 4,
    -2, 3, -1, 2,
    6, 2, 6,
    -6, -3,
    5,
)
EXPECTED_NEGATIVE_BLOCKS = (
    (0, 1),
    (9, 11),
    (4, 6, 8),
    (2, 7, 12),
    (3, 10, 14),
    (5, 13),
)
EXPECTED_SEMANTIC_SHA256 = "80ecee351d7648a6a4c00f879a649d8fd4dbe87119929e5a5ded72d693c9859d"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_repr(value: object) -> str:
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def load_modes():
    spec = importlib.util.spec_from_file_location("thm3398_modes_for_3409", MODE_PATH)
    require(spec is not None and spec.loader is not None, "mode module spec")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def containing_centre(mode: tuple[object, ...], speed: int, time: Fraction):
    """Return the unique lifted mode centre containing ``time``, if any."""
    residue = int(mode[1])
    width = int(mode[2])
    shifted = speed * time - Fraction(residue, 2 * Q)
    floor = shifted.numerator // shifted.denominator
    centres = []
    for tooth in range(floor - 1, floor + 3):
        if abs(shifted - tooth) < Fraction(width, 14 * Q):
            centres.append(
                Fraction(residue, 2 * Q * speed) + Fraction(tooth, speed)
            )
    require(len(centres) <= 1, (speed, time, mode, centres))
    return None if not centres else centres[0]


def exact_block(module, speed: int, time: Fraction) -> tuple[int, ...]:
    return tuple(
        sheet
        for sheet in range(Q)
        if module.source_danger(Q, speed, time, sheet)
    )


def is_cover(module, time: Fraction) -> bool:
    return all(
        any(module.source_danger(Q, speed, time, sheet) for speed in SPEEDS)
        for sheet in range(Q)
    )


def cochain_values(centres: tuple[Fraction, ...]) -> tuple[int, ...]:
    values = []
    for left, right in itertools.combinations(range(len(SPEEDS)), 2):
        value = (
            2
            * Q
            * SPEEDS[left]
            * SPEEDS[right]
            * (centres[left] - centres[right])
        )
        require(value.denominator == 1, (left, right, value))
        values.append(value.numerator)
    return tuple(values)


def edge_map(values: tuple[int, ...]) -> dict[tuple[int, int], int]:
    return dict(zip(itertools.combinations(range(len(SPEEDS)), 2), values))


def check_cochain(module, modes, values) -> None:
    mapped = edge_map(values)
    for left, right in itertools.combinations(range(len(SPEEDS)), 2):
        require(
            mapped[(left, right)]
            in module.gap_values(
                Q, SPEEDS[left], SPEEDS[right], modes[left], modes[right]
            ),
            ("pair fibre", left, right, mapped[(left, right)]),
        )
    for left, middle, right in itertools.combinations(range(len(SPEEDS)), 3):
        circulation = (
            SPEEDS[right] * mapped[(left, middle)]
            - SPEEDS[middle] * mapped[(left, right)]
            + SPEEDS[left] * mapped[(middle, right)]
        )
        require(circulation == 0, ("triangle", left, middle, right, circulation))
    compiled = module.complete_cochain(Q, SPEEDS, modes)
    require(compiled is not None, ("complete cochain", modes, values))
    require(tuple(value for _, _, value in compiled) == values, (compiled, values))
    realized = module.realize_cochain(Q, SPEEDS, modes, compiled)
    require(is_cover(module, realized), ("realized cover", realized))


def translated_blocks(blocks: tuple[tuple[int, ...], ...], shift: int):
    return tuple(
        tuple(sorted((sheet + shift) % Q for sheet in block))
        for block in blocks
    )


def canonical_blocks(blocks: tuple[tuple[int, ...], ...]):
    return min(translated_blocks(blocks, shift) for shift in range(Q))


def reflected_blocks(blocks: tuple[tuple[int, ...], ...]):
    return tuple(
        tuple(sorted((1 - sheet) % Q for sheet in block))
        for block in blocks
    )


def negated_blocks(blocks: tuple[tuple[int, ...], ...]):
    return tuple(
        tuple(sorted((-sheet) % Q for sheet in block))
        for block in blocks
    )


def connected_tree(indices: tuple[int, ...], edges) -> bool:
    parents = list(range(len(SPEEDS)))

    def root(vertex: int) -> int:
        while parents[vertex] != vertex:
            parents[vertex] = parents[parents[vertex]]
            vertex = parents[vertex]
        return vertex

    for index in indices:
        left, right = edges[index]
        root_left = root(left)
        root_right = root(right)
        if root_left == root_right:
            return False
        parents[root_left] = root_right
    return len({root(vertex) for vertex in range(len(SPEEDS))}) == 1


def tree_tariff(values: tuple[int, ...]):
    edges = tuple(itertools.combinations(range(len(SPEEDS)), 2))
    rows = []
    for indices in itertools.combinations(range(len(edges)), len(SPEEDS) - 1):
        if not connected_tree(indices, edges):
            continue
        absolute = tuple(abs(values[index]) for index in indices)
        rows.append((sum(absolute), max(absolute), indices))
    require(len(rows) == 6 ** 4, ("Cayley tree count", len(rows)))
    minimum_l1 = min(row[0] for row in rows)
    minimum_linf = min(row[1] for row in rows)
    minimum_l1_count = sum(row[0] == minimum_l1 for row in rows)
    minimum_linf_count = sum(row[1] == minimum_linf for row in rows)
    owner_seven_star = tuple(
        index for index, edge in enumerate(edges) if 1 in edge
    )
    star_absolute = tuple(abs(values[index]) for index in owner_seven_star)
    require((sum(star_absolute), max(star_absolute)) == (10, 3), star_absolute)
    return (
        len(rows),
        minimum_l1,
        minimum_l1_count,
        minimum_linf,
        minimum_linf_count,
        tuple(edges[index] for index in owner_seven_star),
        tuple(values[index] for index in owner_seven_star),
    )


def packet_census(module):
    scale, samples = module.event_samples(Q, SPEEDS)
    require(len(samples) % 2 == 0, len(samples))
    boundary_count = len(samples) // 2
    events = tuple(sample // 2 for sample in samples[:boundary_count])
    require(tuple(sorted(events)) == events, "event ordering")

    universe = frozenset(range(Q))
    packets = {}
    cover_samples = 0
    boundary_covers = 0
    interval_checks = 0
    packet_instances = 0

    for sample_index, sample in enumerate(samples):
        time = Fraction(sample, 2 * scale)
        physical_cover = is_cover(module, time)
        if physical_cover:
            cover_samples += 1
            boundary_covers += int(sample_index < boundary_count)

        active_banks = []
        actual_blocks = []
        for speed in SPEEDS:
            actual = exact_block(module, speed, time)
            actual_blocks.append(actual)
            active = []
            for mode in module.owner_modes(Q, speed):
                centre = containing_centre(mode, speed, time)
                direct = all(sheet in actual for sheet in mode[0])
                require(bool(centre is not None) == direct, (
                    "interval/direct mismatch", speed, time, mode, actual, centre,
                ))
                interval_checks += 1
                if centre is not None:
                    active.append((mode, centre))
            active_banks.append(tuple(active))

        packet_count_at_sample = 0
        for chosen in itertools.product(*active_banks):
            modes = tuple(item[0] for item in chosen)
            blocks = tuple(tuple(sorted(mode[0])) for mode in modes)
            if frozenset().union(*(mode[0] for mode in modes)) != universe:
                continue
            require(physical_cover, ("selected cover without physical cover", time))
            centres = tuple(item[1] for item in chosen)
            values = cochain_values(centres)
            check_cochain(module, modes, values)
            require(blocks == tuple(actual_blocks), (
                "nonmaximal covering packet", time, blocks, actual_blocks,
            ))

            require(sample_index >= boundary_count, ("boundary packet", time))
            cell_index = sample_index - boundary_count
            left = events[cell_index]
            right = events[(cell_index + 1) % boundary_count]
            if cell_index + 1 == boundary_count:
                right += scale
            cell = (Fraction(left, scale), Fraction(right, scale))
            require(cell[0] < time < cell[1], (cell, time))

            mode_data = tuple(
                (int(mode[1]), int(mode[2]), int(mode[3]), int(mode[4]))
                for mode in modes
            )
            key = (blocks, mode_data, centres, values)
            packets.setdefault(key, []).append((time, cell))
            packet_instances += 1
            packet_count_at_sample += 1
        require(packet_count_at_sample == int(physical_cover), (
            "packet multiplicity at sample", time, physical_cover, packet_count_at_sample,
        ))

    require((scale, boundary_count, len(samples)) == (15_135_120, 1_260, 2_520), (
        scale, boundary_count, len(samples),
    ))
    require((cover_samples, boundary_covers) == (30, 0), (cover_samples, boundary_covers))
    require(packet_instances == len(packets) == 30, (packet_instances, len(packets)))
    require(all(len(occurrences) == 1 for occurrences in packets.values()), "duplicate packet")

    widths = Counter(
        cell[1] - cell[0]
        for occurrences in packets.values()
        for _, cell in occurrences
    )
    total_measure = sum(
        (width * count for width, count in widths.items()), Fraction(0)
    )
    require(widths == Counter({Fraction(1, 3_696): 30}), widths)
    require(total_measure == Fraction(5, 616), total_measure)

    orbit_bank = {}
    for key, occurrences in packets.items():
        blocks, _, _, values = key
        orbit_key = (canonical_blocks(blocks), values)
        orbit_bank.setdefault(orbit_key, []).extend(occurrences)
    orbit_records = tuple(
        sorted(
            (
                blocks,
                values,
                len(occurrences),
                sum((cell[1] - cell[0] for _, cell in occurrences), Fraction(0)),
            )
            for (blocks, values), occurrences in orbit_bank.items()
        )
    )
    expected_orbits = (
        (
            EXPECTED_POSITIVE_BLOCKS,
            EXPECTED_POSITIVE_COCHAIN,
            15,
            Fraction(5, 1_232),
        ),
        (
            EXPECTED_NEGATIVE_BLOCKS,
            tuple(-value for value in EXPECTED_POSITIVE_COCHAIN),
            15,
            Fraction(5, 1_232),
        ),
    )
    require(orbit_records == expected_orbits, (orbit_records, expected_orbits))
    require(
        reflected_blocks(EXPECTED_POSITIVE_BLOCKS) == EXPECTED_NEGATIVE_BLOCKS,
        "canonical time-reversal block action",
    )

    packet_by_time = {
        occurrences[0][0]: key for key, occurrences in packets.items()
    }
    require(len(packet_by_time) == len(packets), "packet time collision")
    for time, key in packet_by_time.items():
        reversed_key = packet_by_time.get((-time) % 1)
        require(reversed_key is not None, ("missing reversed packet", time))
        require(reversed_key[0] == negated_blocks(key[0]), (
            "physical time-reversal block action", time, key[0], reversed_key[0],
        ))
        require(reversed_key[3] == tuple(-value for value in key[3]), (
            "physical time-reversal cochain action", time, key[3], reversed_key[3],
        ))

    for blocks, values, _, _ in orbit_records:
        orbit_blocks = {
            translated_blocks(blocks, shift) for shift in range(Q)
        }
        actual_orbit_blocks = {
            key[0]
            for key in packets
            if key[3] == values
        }
        require(len(orbit_blocks) == Q, ("nonfree translation orbit", blocks))
        require(actual_orbit_blocks == orbit_blocks, (
            "translation orbit mismatch", blocks, actual_orbit_blocks, orbit_blocks,
        ))

    norm_profile = Counter(
        (
            sum(abs(value) for value in key[3]),
            max(abs(value) for value in key[3]),
            sum(value * value for value in key[3]),
        )
        for key in packets
    )
    require(norm_profile == Counter({(50, 6, 206): 30}), norm_profile)

    packet_rows = tuple(
        sorted(
            (
                key[0],
                key[1],
                key[2],
                key[3],
                occurrences[0][0],
                occurrences[0][1],
            )
            for key, occurrences in packets.items()
        )
    )
    return (
        scale,
        boundary_count,
        len(samples),
        interval_checks,
        cover_samples,
        boundary_covers,
        len(packets),
        tuple(sorted(widths.items())),
        total_measure,
        orbit_records,
        tuple(sorted(norm_profile.items())),
        digest_repr(packet_rows),
        digest_repr(tuple((row[4], row[5]) for row in packet_rows)),
    )


def main() -> None:
    for name, path, expected in PINS:
        actual = sha256_lf(path)
        require(actual == expected, ("dependency changed", name, path, actual, expected))

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    module = load_modes()
    census = packet_census(module)
    tariff = tree_tariff(EXPECTED_POSITIVE_COCHAIN)
    semantic = digest_repr((Q, SPEEDS, census, tariff))
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    (
        scale,
        boundaries,
        samples,
        interval_checks,
        cover_samples,
        boundary_covers,
        packet_count,
        widths,
        total_measure,
        orbits,
        norms,
        packet_digest,
        cell_digest,
    ) = census
    print("THM-3409 Q15 EXCEPTIONAL POSITIVE-COCHAIN RIGIDITY AND TARIFF")
    print(f"source_sha256_lf={sha256_lf(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINS)}")
    print("status=PROVED exact event-cell classification plus VERIFIED-EXACT companion;independently_audited")
    print(f"edge=q{Q}:{SPEEDS};inherited_boundary=unique_rank6_physical_edge_outside_capped_zero_cochain_slice")
    print(f"event_universe=(scale,boundaries,samples,interval_checks)=({scale},{boundaries},{samples},{interval_checks})")
    print(f"cover_locus=(open_cells,boundary_covers,width_histogram,total_measure)=({cover_samples},{boundary_covers},{widths},{total_measure})")
    print(f"selected_packets={packet_count};all_maximal_exact_partitions=YES;packet_sha256={packet_digest};cell_sha256={cell_digest}")
    print(f"translation_orbits=(canonical_blocks,lex_pair_cochain,count,measure)={orbits}")
    print("time_reversal=t_to_minus_t_and_ell_to_minus_ell;canonical_representative_map=ell_to_1_minus_ell;cochain_sign_flip=YES")
    print(f"complete_tariff=(L1,Linf,L2_squared)_histogram={norms}")
    print(f"spanning_tree_tariff=(tree_count,minL1,minL1_count,minLinf,minLinf_count,owner7_star_edges,owner7_star_values)={tariff}")
    print("scope=fixed_edge_and_owner_order;one_selected_mode_per_distinct_owner;strict_open_physical_cover;no_LRC14_decrement;not_a_uniform_all_edges_tariff")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
