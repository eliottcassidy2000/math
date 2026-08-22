#!/usr/bin/env python3
"""Exact mobile common-centre rank probe and general zero-centre CRT cover.

A mobile common-centre certificate chooses one THM-3398 selected-block mode
per owner, allows their mode-centre lattices to share an arbitrary rational
point c, and asks their sheet blocks to cover Z/qZ.  This differs from
THM-3401's fixed physical source time zero.  The literal finite probe uses
owners 1,...,14 for q=15,...,28.  Separately, an elementary construction gives
an abstract zero-centre cover for every q>=15 from prime-kernel modes and one
unit trimode per sign class.

This is an unnumbered PROVED-ELEMENTARY / FINITE-EXACT complement to THM-3401,
not an alternate reading of its fixed-zero rank.  MISTAKE-384 records the
former prose-level identification of fixed zero with the whole zero-cochain
locus; the theorem's precise ranks and proof are unchanged.
Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODE_PATH = ROOT / "04-computation/lrc_general_finite_mode_sheet_cover_cochain_thm3398.py"
PINNED = (
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
        "THM-3401",
        ROOT / "01-canon/theorems/THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight.md",
        "dae088cbda12fb64d24f84ab26a6879e94939e04cb03601d8fb996a48c077716",
    ),
    (
        "THM-3401-script",
        ROOT / "04-computation/lrc_centered_transverse_sheet_cover_rank_thm3401.py",
        "fff146868de4b1ec5993ed404fca4200e8e3eac47a7cb902ff51d559eef228e0",
    ),
    (
        "THM-3401-output",
        ROOT / "05-knowledge/results/lrc_centered_transverse_sheet_cover_rank_thm3401.out",
        "12f2cf337c982068c8cfa0cb351a2772f05e55c1be18df4d0414f9a7251327dd",
    ),
    (
        "full-q8-q15-script",
        ROOT / "04-computation/lrc_q8_q15_full_physical_clutter_audit_20260815.py",
        "e54b77eeae05484cbbfacd904e850815f7e78a5e3306f21ad87a68ffbfae9e2e",
    ),
    (
        "full-q8-q15-output",
        ROOT / "05-knowledge/results/lrc_q8_q15_full_physical_clutter_audit_20260815.out",
        "8b8cc6f45ab8b14b0ba26afb29748ccf3bc08f4f103581ba9afd7167391e8008",
    ),
)

EXPECTED_MOBILE_COMMON_CENTRE_RANKS = (
    (15, 6),
    (16, 4),
    (17, 8),
    (18, 4),
    (19, 9),
    (20, 6),
    (21, 8),
    (22, 6),
    (23, 6),
    (24, 6),
    (25, 7),
    (26, 8),
    (27, 9),
    (28, 8),
)
EXPECTED_FIXED_ZERO_RANKS = (
    (15, 6),
    (16, 5),
    (17, 8),
    (18, 5),
    (19, 9),
    (20, 6),
    (21, 8),
    (22, 7),
    (23, 11),
    (24, 6),
    (25, 11),
    (26, 8),
    (27, 10),
    (28, 8),
)
EXPECTED_SEMANTIC_DIGEST = "a9a06cdc130ec2eb1d8292166fd1e4ca27f16b4b0a3210821b7c976d2464a467"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


def load_modes():
    spec = spec_from_file_location("thm3398_modes", MODE_PATH)
    require(spec is not None and spec.loader is not None, "mode module spec")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prime_divisors(value):
    result = []
    divisor = 2
    remaining = value
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            result.append(divisor)
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        result.append(remaining)
    return tuple(result)


def euler_phi(value):
    result = value
    for prime in prime_divisors(value):
        result = result // prime * (prime - 1)
    return result


def sign_representatives(q):
    return tuple(
        residue
        for residue in range(1, (q + 1) // 2)
        if gcd(residue, q) == 1
    )


def selected_mode(module, q, speed, start, size):
    matches = tuple(
        mode
        for mode in module.owner_modes(q, speed)
        if int(mode[3]) == start % int(mode[6]) and int(mode[4]) == size
    )
    require(len(matches) == 1, ("selected mode", q, speed, start, size, matches))
    return matches[0]


def canonical_crt_cover(module, q):
    primes = prime_divisors(q)
    kernel_speeds = () if len(primes) == 1 and primes[0] == q else tuple(q // p for p in primes)
    unit_speeds = sign_representatives(q)
    speeds = tuple(sorted(kernel_speeds + unit_speeds))
    require(len(speeds) == len(set(speeds)), (q, "owner collision", speeds))

    modes_by_speed = {}
    for speed in kernel_speeds:
        prime = q // speed
        mode = selected_mode(module, q, speed, 0, 1)
        require(set(mode[0]) == set(range(0, q, prime)), (q, speed, mode[0]))
        modes_by_speed[speed] = mode
    for speed in unit_speeds:
        mode = selected_mode(module, q, speed, q - 1, 3)
        inverse = pow(speed, -1, q)
        expected = {0, inverse, (-inverse) % q}
        require(set(mode[0]) == expected, (q, speed, mode[0], expected))
        modes_by_speed[speed] = mode

    modes = tuple(modes_by_speed[speed] for speed in speeds)
    require(all(int(mode[1]) == 0 for mode in modes), (q, "nonzero centre"))
    require(set().union(*(set(mode[0]) for mode in modes)) == set(range(q)), (q, "cover"))
    for left, right in combinations(range(len(speeds)), 2):
        require(
            0 in module.gap_values(q, speeds[left], speeds[right], modes[left], modes[right]),
            (q, "zero gap", speeds[left], speeds[right]),
        )

    nonunits = {residue for residue in range(q) if gcd(residue, q) != 1}
    kernel_union = set().union(
        *(set(modes_by_speed[speed][0]) for speed in kernel_speeds)
    ) if kernel_speeds else {0}
    require(kernel_union == nonunits, (q, "kernel union", kernel_union, nonunits))
    unit_blocks = tuple(set(modes_by_speed[speed][0]) - {0} for speed in unit_speeds)
    require(set().union(*unit_blocks) == set(range(q)) - nonunits, (q, "unit union"))
    require(sum(map(len, unit_blocks)) == euler_phi(q), (q, "unit overlap"))

    expected_rank = euler_phi(q) // 2 + (0 if not kernel_speeds else len(primes))
    require(len(speeds) == expected_rank, (q, "rank formula", speeds, expected_rank))
    return (
        q,
        speeds,
        tuple(tuple(sorted(mode[0])) for mode in modes),
        expected_rank,
        len(primes),
        euler_phi(q) // 2,
    )


def mode_centres(q, speed, mode):
    h = int(mode[1])
    return tuple(
        Fraction(2 * q * tooth + h, 2 * q * speed) % 1
        for tooth in range(speed)
    )


def mobile_common_centre_atlas(module, q, owners):
    atlas = {}
    for speed in owners:
        for mode in module.owner_modes(q, speed):
            for centre in mode_centres(q, speed, mode):
                atlas.setdefault(centre, {}).setdefault(speed, []).append(mode)

    reduced = {}
    for centre, owner_banks in atlas.items():
        reduced[centre] = {}
        for speed, bank in owner_banks.items():
            unique = {}
            for mode in bank:
                unique[frozenset(mode[0])] = mode
            blocks = tuple(unique)
            require(
                all(left <= right or right <= left for left, right in combinations(blocks, 2)),
                (q, speed, centre, "nonchain common-centre modes", blocks),
            )
            largest = max(unique.values(), key=lambda mode: (len(mode[0]), -int(mode[4]), -int(mode[3])))
            reduced[centre][speed] = largest
    return reduced


def mobile_common_centre_rank(module, q):
    owners = tuple(range(1, 15))
    full = (1 << q) - 1
    atlas = mobile_common_centre_atlas(module, q, owners)
    best = len(owners) + 1
    maximalized_certificates = []
    cover_capable_centres = []

    for centre in sorted(atlas):
        owner_modes = atlas[centre]
        masks = {
            speed: sum(1 << sheet for sheet in owner_modes[speed][0])
            for speed in owner_modes
        }
        all_available = 0
        for mask in masks.values():
            all_available |= mask
        if all_available != full:
            continue
        cover_capable_centres.append(centre)
        available = tuple(sorted(masks))
        for size in range(1, min(best, len(available)) + 1):
            found_at_size = []
            for speeds in combinations(available, size):
                covered = 0
                for speed in speeds:
                    covered |= masks[speed]
                if covered == full:
                    found_at_size.append(speeds)
            if found_at_size:
                if size < best:
                    best = size
                    maximalized_certificates = []
                if size == best:
                    for speeds in found_at_size:
                        maximalized_certificates.append(
                            (
                                centre,
                                speeds,
                                tuple(tuple(sorted(owner_modes[speed][0])) for speed in speeds),
                                tuple(
                                    (
                                        int(owner_modes[speed][1]),
                                        int(owner_modes[speed][2]),
                                        int(owner_modes[speed][3]),
                                        int(owner_modes[speed][4]),
                                    )
                                    for speed in speeds
                                ),
                            )
                        )
                break

    require(best <= len(owners), (q, "no mobile common-centre cover"))
    maximalized_certificates = tuple(sorted(set(maximalized_certificates), key=repr))
    require(
        all(len(item[1]) == best for item in maximalized_certificates),
        (q, "certificate rank"),
    )
    canonical = maximalized_certificates[0]
    cover_capable_twists = tuple(
        sorted({Fraction(q) * centre % 1 for centre in cover_capable_centres})
    )
    minimum_twists = tuple(
        sorted({Fraction(q) * item[0] % 1 for item in maximalized_certificates})
    )
    return (
        q,
        best,
        len(maximalized_certificates),
        canonical,
        tuple(sorted({item[1] for item in maximalized_certificates})),
        cover_capable_twists,
        minimum_twists,
    )


def common_centre_witness(module, q, centre, speeds):
    atlas = mobile_common_centre_atlas(module, q, speeds)
    require(centre in atlas, (q, centre, "missing centre"))
    require(all(speed in atlas[centre] for speed in speeds), (q, centre, speeds))
    modes = tuple(atlas[centre][speed] for speed in speeds)
    require(set().union(*(set(mode[0]) for mode in modes)) == set(range(q)), (q, "cover"))
    require(
        all(
            0 in module.gap_values(q, speeds[left], speeds[right], modes[left], modes[right])
            for left, right in combinations(range(len(speeds)), 2)
        ),
        (q, "nonzero cochain"),
    )
    return (
        q,
        centre,
        speeds,
        tuple(tuple(sorted(mode[0])) for mode in modes),
        tuple((int(mode[1]), int(mode[2]), int(mode[3]), int(mode[4])) for mode in modes),
    )


def physical_edges_at_rank(module, q, owners, rank):
    scale, samples = module.event_samples(q, owners)
    bank = {}
    for speed in owners:
        sheets = []
        for sheet in range(q):
            bits = 0
            for index, sample in enumerate(samples):
                if module.integer_sample_danger(q, speed, sample, scale, sheet):
                    bits |= 1 << index
            sheets.append(bits)
        bank[speed] = tuple(sheets)

    sample_universe = (1 << len(samples)) - 1
    edges = []
    for speeds in combinations(owners, rank):
        simultaneous = sample_universe
        for sheet in range(q):
            blocked = 0
            for speed in speeds:
                blocked |= bank[speed][sheet]
            simultaneous &= blocked
            if not simultaneous:
                break
        if simultaneous:
            edges.append(speeds)
    return tuple(edges)


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))
    module = load_modes()

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

    crt_records = tuple(canonical_crt_cover(module, q) for q in range(15, 101))
    require(all(max(record[1]) <= 14 for record in crt_records[:14]), "literal owner cap q15..28")
    require(crt_records[0][1] == (1, 2, 3, 4, 5, 7), crt_records[0])
    require(crt_records[1][1] == (1, 3, 5, 7, 8), crt_records[1])
    require(crt_records[3][1] == (1, 5, 6, 7, 9), crt_records[3])

    common_centre_records = tuple(mobile_common_centre_rank(module, q) for q in range(15, 29))
    common_centre_ranks = tuple((q, rank) for q, rank, _, _, _, _, _ in common_centre_records)
    require(
        common_centre_ranks == EXPECTED_MOBILE_COMMON_CENTRE_RANKS,
        ("mobile common-centre ranks", common_centre_ranks),
    )

    speed_sets = {q: sets for q, _, _, _, sets, _, _ in common_centre_records}
    require((1, 2, 3, 4, 5, 7) in speed_sets[15], speed_sets[15])
    require((2, 6, 10, 14) in speed_sets[16], speed_sets[16])
    require((2, 10, 12, 14) in speed_sets[18], speed_sets[18])

    crt_q15_q28 = tuple((record[0], record[3]) for record in crt_records[:14])
    require(
        crt_q15_q28 == EXPECTED_FIXED_ZERO_RANKS,
        ("THM-3401 fixed-zero ranks", crt_q15_q28),
    )
    savings = tuple(
        (q, crt_rank - common_centre_rank_value)
        for (q, crt_rank), (_, common_centre_rank_value) in zip(crt_q15_q28, common_centre_ranks)
    )
    compression_support = tuple(q for q, saving in savings if saving)
    twist_profiles = tuple(
        (q, cover_twists, minimum_twists)
        for q, _, _, _, _, cover_twists, minimum_twists in common_centre_records
    )
    require(
        all(cover_twists == (Fraction(0), Fraction(1, 2)) for _, cover_twists, _ in twist_profiles),
        ("cover-capable twists", twist_profiles),
    )
    require(
        all(minimum_twists == (Fraction(1, 2),) for q, _, minimum_twists in twist_profiles if q in compression_support),
        ("strict-gap minimum twists", twist_profiles),
    )
    harmonic_support_mass = sum((Fraction(1, q) for q in compression_support), Fraction(0))
    harmonic_savings_mass = sum(
        (Fraction(saving, q) for q, saving in savings), Fraction(0)
    )
    canonical_excesses = tuple(
        (q, sum(map(len, canonical[2])) - q)
        for q, _, _, canonical, _, _, _ in common_centre_records
    )
    perfect_partition_support = tuple(q for q, excess in canonical_excesses if excess == 0)
    rank_four_support = tuple(q for q, rank in common_centre_ranks if rank == 4)
    rank_six_support = tuple(q for q, rank in common_centre_ranks if rank == 6)

    outside_pool_controls = (
        common_centre_witness(module, 25, Fraction(1, 50), (1, 9, 10, 11, 19, 21)),
        common_centre_witness(module, 27, Fraction(1, 54), (3, 15, 18, 21)),
    )
    require(all(any(speed > 14 for speed in row[2]) for row in outside_pool_controls), outside_pool_controls)

    q15_physical_rank6 = physical_edges_at_rank(module, 15, tuple(range(1, 15)), 6)
    q15_physical_digest = sha256(repr(q15_physical_rank6).encode("ascii")).hexdigest()
    require(len(q15_physical_rank6) == 157, len(q15_physical_rank6))
    require(
        q15_physical_digest == "cc2eba81d46e6ae4e2f3705fc74b9705527200a5d89b393c8aabcff6f9d1fc6d",
        q15_physical_digest,
    )
    q15_common_centre_sets = set(speed_sets[15])
    q15_exceptional = tuple(edge for edge in q15_physical_rank6 if edge not in q15_common_centre_sets)
    require(len(q15_exceptional) == 1, q15_exceptional)
    require(
        not q15_common_centre_sets - set(q15_physical_rank6),
        "common-centre nonphysical q15 edge",
    )
    q15_common_centre_gcd_profile = tuple(
        sorted(
            Counter(
                tuple(sorted((gcd(speed, 15) for speed in edge), reverse=True))
                for edge in q15_common_centre_sets
            ).items()
        )
    )
    require(
        q15_common_centre_gcd_profile == (((5, 3, 1, 1, 1, 1), 156),),
        q15_common_centre_gcd_profile,
    )
    exceptional_gcd_profile = tuple(
        sorted((gcd(speed, 15) for speed in q15_exceptional[0]), reverse=True)
    )
    require(exceptional_gcd_profile == (3, 1, 1, 1, 1, 1), exceptional_gcd_profile)
    exceptional_witness = module.mode_cover_all_firing(15, q15_exceptional[0])
    require(exceptional_witness is not None, q15_exceptional)
    exceptional_order, exceptional_modes, exceptional_cochain, exceptional_time = exceptional_witness
    exceptional_profile = (
        q15_exceptional[0],
        exceptional_order,
        tuple(tuple(sorted(mode[0])) for mode in exceptional_modes),
        tuple((int(mode[1]), int(mode[2]), int(mode[3]), int(mode[4])) for mode in exceptional_modes),
        exceptional_cochain,
        exceptional_time,
        sum(abs(edge[2]) for edge in exceptional_cochain),
        max(abs(edge[2]) for edge in exceptional_cochain),
    )
    require(any(edge[2] for edge in exceptional_cochain), exceptional_profile)

    semantic = ExactDigest()
    semantic.add(("crt", crt_records))
    semantic.add(("fixed_zero", EXPECTED_FIXED_ZERO_RANKS))
    semantic.add(("mobile_common_centre", common_centre_records))
    semantic.add(("twist_profiles", twist_profiles))
    semantic.add(("outside_pool_controls", outside_pool_controls))
    semantic.add(("savings", savings))
    semantic.add(("harmonic", compression_support, harmonic_support_mass, harmonic_savings_mass))
    semantic.add(("partition", canonical_excesses, perfect_partition_support, rank_four_support, rank_six_support))
    semantic.add(
        (
            "q15_mobile_common_centre_boundary",
            q15_physical_digest,
            q15_common_centre_gcd_profile,
            exceptional_profile,
        )
    )
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("LRC ZERO-CENTRE CRT AND MOBILE COMMON-CENTRE RANK EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=PROVED-ELEMENTARY zero-centre CRT upper cover for every q>=15;FINITE-EXACT mobile_common_centre_rank_owner_pool_1..14_q15..28;MISTAKE384_fixed_zero_scope_repair;unnumbered_complement_to_THM3401")
    print("mobile_common_centre_definition=one_selected_mode_per_owner;mode_centre_lattices_share_arbitrary_c;blocks_cover_Zmodq;not_fixed_source_time_zero")
    print("crt=one_kernel_owner_q/p_per_prime_divisor_plus_one_unit_triphase_per_unit_sign_class;all_h=0;all_pij=0")
    print("crt_rank=omega(q)+phi(q)/2_for_composite_q;phi(q)/2_for_prime_q;verified_q15..100")
    print(f"fixed_zero_ranks_q15_q28={crt_q15_q28}")
    print(f"mobile_common_centre_ranks_owner_pool_1..14_q15_q28={common_centre_ranks}")
    print(f"fixed_zero_minus_mobile_rank_gap={savings}")
    print(f"twist_profiles_qc_mod1=(q,cover_capable_twists,minimum_twists)={twist_profiles}")
    print(f"outside_owner_pool_hostiles={outside_pool_controls}")
    print(f"compression_support={compression_support};harmonic_support_mass={harmonic_support_mass};harmonic_savings_mass={harmonic_savings_mass}")
    print(f"canonical_block_excesses={canonical_excesses};perfect_partition_support={perfect_partition_support};rank4_support={rank_four_support};rank6_support={rank_six_support}")
    print(f"q15_rank6_physical_vs_mobile_common_centre=157_vs_{len(q15_common_centre_sets)};common_centre_gcd_profile={q15_common_centre_gcd_profile};unique_positive_cochain_edge={exceptional_profile}")
    for q, rank, count, canonical, speed_sets_at_q, cover_twists, minimum_twists in common_centre_records:
        speed_set_digest = sha256(repr(speed_sets_at_q).encode("ascii")).hexdigest()
        print(f"q={q};mobile_common_centre_rank_owner_pool_1..14={rank};maximalized_certificate_count={count};speed_set_count={len(speed_sets_at_q)};cover_capable_twists={cover_twists};minimum_twists={minimum_twists};speed_set_sha256={speed_set_digest};canonical={canonical}")
    print("dilation_controls=q8_to_q16_(1,3,5,7)->(2,6,10,14);q9_to_q18_(1,5,6,7)->(2,10,12,14)")
    print("scope=literal_owner_pool_1..14_for_mobile_common_centre_minima;abstract_positive_owners_for_general_CRT;THM3401_is_distinct_fixed_zero_subslice_of_zero_cochains;no_core_rescue_or_LRC14_ledger_consequence")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
