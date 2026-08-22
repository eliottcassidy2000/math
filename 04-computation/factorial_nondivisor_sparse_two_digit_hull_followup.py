#!/usr/bin/env python3
"""Exact audit for the nondivisor residue-hull follow-up to THM-3483.

This companion does not extend the d-boundary.  It mines the three stored
THM-3483 census blocks, checks the p=3 and p=5 automatic-admissibility
tables, verifies the sparse two-unit-digit hull formula on a declared grid,
and replays the positive row d=6590 and the hostile row d=6518.
"""

import ast
from collections import Counter
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RHO_NAME = "factorial_nondivisor_residue_digit_pair_compiler_thm3483.py"
RHO_SHA256 = "9e37ead620f141617a9c6d51c182e09c034945793092e56e39fb061254662723"
ADAPTIVE_NAME = "factorial_adaptive_rho_block_6000.py"
ADAPTIVE_SHA256 = "b65edcf2870714ca57456b8297afdd05284a09b302ec4b84d2e57829520c94d1"
OUTPUTS = (
    ("factorial_adaptive_rho_block_4000_thm3483.out",
     "48733decf5874c197b989d7731f0864082a65fe7b9056af67e0149d1e8a94896"),
    ("factorial_adaptive_rho_block_6000.out",
     "ab629edc04e31d1889741688897bfe60f5249df60b64116391217393962b1ddf"),
    ("factorial_adaptive_rho_block_10000.out",
     "18b131aed2f380b1c1bace8beeb8488ced0e24599f4b7484d66a14e5869c0d22"),
)
GRID = ((3, 8), (5, 6), (7, 5), (11, 4))


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def file_digest(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def semantic_digest(value):
    data = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(data).hexdigest()


def load_pinned(name, expected, module_name):
    path = HERE / name
    actual = file_digest(path)
    require(actual == expected, (name, actual, expected))
    spec = importlib.util.spec_from_file_location(module_name, path)
    require(spec is not None and spec.loader is not None, ("bad import", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_assignment(text, label, terminal=None):
    line = next(row for row in text.splitlines() if row.startswith(label + "="))
    value = line[len(label) + 1:]
    if terminal is not None:
        value = value.split(terminal, 1)[0]
    return ast.literal_eval(value)


def load_census():
    killers, packets = [], {}
    for name, expected in OUTPUTS:
        path = ROOT / "05-knowledge" / "results" / name
        actual = file_digest(path)
        require(actual == expected, (name, actual, expected))
        text = path.read_text(encoding="utf-8")
        killers.extend(parse_assignment(text, "rho_killers", " histogram="))
        packets.update(parse_assignment(text, "rho_needed_divisor_packets"))
    require(len(killers) == 281, len(killers))
    require(len(packets) == 281, len(packets))
    return tuple(killers), packets


def heights(prime, exponent):
    power = prime ** exponent
    return 2 * (power - 1) // (prime - 1)


def expected_sparse(prime, large_exponent, small_exponent):
    large = prime ** large_exponent
    small = prime ** small_exponent
    small_height = heights(prime, small_exponent)
    large_height = heights(prime, large_exponent)
    f_hull = (
        (0, 0),
        (small, small_height),
        (large + small, large_height + small_height),
    )
    g_hull = (
        (0, 0),
        (1, 0),
        (small + 1, small_height),
        (large + small + 1, large_height + small_height),
    )
    blocks = (
        (str(Fraction(small_height, small)), small, small),
        (str(Fraction(large_height, large)), large, large),
    )
    degrees = (0, small, large, large + small)
    return f_hull, g_hull, blocks, degrees


def sparse_representation(d, prime):
    n, nonzero, exponent = d - 2, [], 0
    while n:
        digit = n % prime
        if digit:
            nonzero.append((exponent, digit))
        n //= prime
        exponent += 1
    if len(nonzero) != 2 or any(digit != 1 for _, digit in nonzero):
        return None
    small, large = nonzero
    if small[0] < 1:
        return None
    return large[0], small[0]


def main():
    rho_module = load_pinned(RHO_NAME, RHO_SHA256, "sparse_rho")
    adaptive = load_pinned(ADAPTIVE_NAME, ADAPTIVE_SHA256, "sparse_adaptive")
    killers, packets = load_census()

    # For G, n == 1 and d == 2 modulo p.  These are the complete residue
    # tables j_0=0,...,p-1.  F needs only the THM-3161 vertex residues 0,h.
    g_tables, f_vertex_tables = {}, {}
    for prime, expected_g in ((3, (2, 1, 2)), (5, (3, 1, 1, 2, 2))):
        d, n_g, n_f = prime + 2, prime + 1, prime
        g_table = tuple(rho_module.rho(n_g, j, d, prime) for j in range(prime))
        h = (prime - 1) // 2
        f_table = tuple(rho_module.rho(n_f, j, d, prime) for j in (0, h))
        require(g_table == expected_g, (prime, g_table, expected_g))
        require(f_table == (1, 1), (prime, f_table))
        g_tables[prime], f_vertex_tables[prime] = g_table, f_table

    killer_histogram = tuple(sorted(Counter(prime for _, prime in killers).items()))
    require(killer_histogram == (
        (3, 31), (5, 23), (7, 81), (11, 83), (13, 31),
        (17, 18), (19, 9), (23, 4), (29, 1),
    ), killer_histogram)
    lane_summary = []
    for prime in (3, 5):
        lane = tuple((d, killer) for d, killer in killers if d % prime == 2)
        at_prime = tuple(d for d, killer in lane if killer == prime)
        require(all(d % prime == 2 for d, killer in killers if killer == prime), prime)
        lane_summary.append((prime, len(lane), len(at_prime), len(lane) - len(at_prime)))
    require(tuple(lane_summary) == ((3, 66, 31, 35), (5, 38, 23, 15)),
            lane_summary)

    # Direct finite controls for the all-height sparse theorem.
    cells = []
    for prime, max_large_exponent in GRID:
        for large_exponent in range(2, max_large_exponent + 1):
            for small_exponent in range(1, large_exponent):
                large = prime ** large_exponent
                small = prime ** small_exponent
                n, d = large + small, large + small + 2
                f_hull = rho_module.raw_hull(n, prime)
                g_hull = rho_module.raw_hull(n + 1, prime)
                degrees, blocks = rho_module.pair_degrees(f_hull, g_hull)
                expected = expected_sparse(prime, large_exponent, small_exponent)
                require((f_hull, g_hull, blocks, degrees) == expected,
                        (prime, large_exponent, small_exponent,
                         f_hull, g_hull, blocks, degrees, expected))
                f_rhos = tuple(rho_module.rho(n, j, d, prime) for j, _ in f_hull)
                g_rhos = tuple(rho_module.rho(n + 1, j, d, prime)
                                for j, _ in g_hull)
                require(f_rhos == (1, 1, 1),
                        (prime, large_exponent, small_exponent, f_rhos))
                require(g_rhos == (pow(2, -1, prime), 1, 1, 1),
                        (prime, large_exponent, small_exponent, g_rhos))
                cells.append((prime, large_exponent, small_exponent,
                              f_hull, g_hull, blocks, degrees, f_rhos, g_rhos))
    require(len(cells) == 59, len(cells))

    sparse_census = tuple(
        (d, prime, sparse_representation(d, prime))
        for d, killer in killers
        for prime in (3, 5)
        if sparse_representation(d, prime) is not None
    )
    require(sparse_census == ((6590, 3, (8, 3)),), sparse_census)

    positive_d, positive_prime = 6590, 3
    positive_expected = expected_sparse(positive_prime, 8, 3)
    positive_packet = packets[positive_d]
    positive_intersection = tuple(sorted(set(positive_packet) & set(positive_expected[3])))
    require(positive_packet == tuple(599 * t for t in range(1, 11)), positive_packet)
    require(positive_intersection == (), positive_intersection)
    require(all(3 ** 3 < degree < 3 ** 8 for degree in positive_packet),
            positive_packet)

    # Congruence-only hostile: p=3 is admissible at d=6518, but its integer
    # slope block preserves the complete six-degree divisor packet.  p=29 is
    # the unique final killer in the stored census.
    divisor, rho_check, independent, coefficient, digital = adaptive.load_engines()
    hostile = adaptive.scan_row(divisor, rho_check, coefficient, digital, 6518)
    require(hostile[1:3] == ("rho", 29), hostile[1:3])
    require(hostile[-1] == (3087, 3430, 4802, 5145, 5488, 5831), hostile[-1])
    hostile_p3 = next(record for record in hostile[4] if record[0] == 3)
    hostile_p29 = next(record for record in hostile[4] if record[0] == 29)
    require(hostile_p3[1] == "admissible" and hostile_p3[-1] == hostile[-1],
            hostile_p3)
    require(hostile_p29[-2:] == ((5831,), ()), hostile_p29)

    # The strict K>s hypothesis avoids digit collision.  At p=3,K=s=2 the
    # integer-slope block fills every degree, not merely {0,9,18}.
    boundary_n, boundary_p = 18, 3
    boundary_f = rho_module.raw_hull(boundary_n, boundary_p)
    boundary_g = rho_module.raw_hull(boundary_n + 1, boundary_p)
    boundary_degrees, boundary_blocks = rho_module.pair_degrees(boundary_f, boundary_g)
    require(boundary_degrees == tuple(range(19)), boundary_degrees)

    # A hostile to the tempting predecessor lane d-1=5q, q==2 (mod 3):
    # d=56 has q=11, is automatically p=3 admissible, and retains all four
    # degrees forced by the q-divisor barcode.
    five_q_d, five_q_q = 56, 11
    five_q_n = five_q_d - 2
    five_q_f = rho_module.raw_hull(five_q_n, 3)
    five_q_g = rho_module.raw_hull(five_q_n + 1, 3)
    five_q_degrees, five_q_blocks = rho_module.pair_degrees(five_q_f, five_q_g)
    five_q_candidates = tuple(five_q_q * t for t in range(1, 5))
    five_q_intersection = tuple(
        degree for degree in five_q_candidates if degree in set(five_q_degrees)
    )
    require(five_q_intersection == five_q_candidates, five_q_intersection)

    # The p=5 analogue d-1=3q, q==2 (mod 5), also needs a digit sidecar.
    three_q_d, three_q_q = 22, 7
    three_q_n = three_q_d - 2
    three_q_f = rho_module.raw_hull(three_q_n, 5)
    three_q_g = rho_module.raw_hull(three_q_n + 1, 5)
    three_q_degrees, three_q_blocks = rho_module.pair_degrees(three_q_f, three_q_g)
    three_q_candidates = (three_q_q, 2 * three_q_q)
    three_q_intersection = tuple(
        degree for degree in three_q_candidates if degree in set(three_q_degrees)
    )
    require(three_q_intersection == three_q_candidates, three_q_intersection)

    semantic = (
        tuple(sorted(g_tables.items())), tuple(sorted(f_vertex_tables.items())),
        killer_histogram, tuple(lane_summary), tuple(cells), sparse_census,
        (positive_d, positive_packet, positive_expected, positive_intersection),
        (6518, hostile[-1], hostile_p3, hostile_p29),
        (boundary_p, boundary_n, boundary_f, boundary_g,
         boundary_blocks, boundary_degrees),
        (five_q_d, five_q_q, five_q_f, five_q_g, five_q_blocks,
         five_q_candidates, five_q_intersection),
        (three_q_d, three_q_q, three_q_f, three_q_g, three_q_blocks,
         three_q_candidates, three_q_intersection),
    )

    print("FACTORIAL NONDIVISOR SPARSE TWO-DIGIT HULL FOLLOW-UP")
    print("scope=exact-support quadratic resonance pair; no d-boundary extension")
    print(f"rho_G_tables={tuple(sorted(g_tables.items()))}")
    print(f"rho_F_vertex_tables={tuple(sorted(f_vertex_tables.items()))}")
    print(f"rho_killer_histogram={killer_histogram} total={len(killers)}")
    print(f"automatic_lane_summary={tuple(lane_summary)} schema=(p,lane_rows,killed_at_p,not_killed_at_p)")
    print(f"sparse_grid={GRID} cells={len(cells)} status=exact_control_for_proved_formula")
    print(f"sparse_census_rows={sparse_census}")
    print(f"d6590_packet={positive_packet} local_degrees={positive_expected[3]} intersection={positive_intersection}")
    print(f"d6518_packet={hostile[-1]} p3_blocks={hostile_p3[4]} p3_post={hostile_p3[-1]}")
    print(f"d6518_p29_pre_post={hostile_p29[-2:]} final_killer={hostile[2]}")
    print(f"collision_boundary=(p={boundary_p},n={boundary_n},blocks={boundary_blocks},degree_count={len(boundary_degrees)})")
    print(f"five_q_hostile=(d={five_q_d},q={five_q_q},blocks={five_q_blocks},candidates={five_q_candidates},intersection={five_q_intersection})")
    print(f"three_q_hostile=(d={three_q_d},q={three_q_q},blocks={three_q_blocks},candidates={three_q_candidates},intersection={three_q_intersection})")
    print(f"semantic_sha256={semantic_digest(semantic)}")


if __name__ == "__main__":
    main()
