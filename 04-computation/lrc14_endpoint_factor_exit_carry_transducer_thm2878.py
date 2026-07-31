#!/usr/bin/env python3
"""Exact companion for THM-2878 endpoint-factor carry transducer.

The companion answers four questions left by the THM-2874
full-pattern-complement addendum.

1. What are all lawful representatives of one of the 169 endpoint
   characters, and what does that force on coordinates ell_0 and ell_5?
2. Can any address, representative gauge, source/target origin, residue
   shift, or orientation realize the literal (guard safe, q5 dangerous)
   corner on the selected atom orbit?
3. Does adjoining one binary tag create a faithful four-corner factor
   square?
4. Can a V4 filler carry the order-thirteen ancestry holonomy, or is the
   nonsplit C169 carry fibre still independently necessary?

It proves no row exclusion and no LRC(14) conclusion.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
RESULTS = ROOT / "05-knowledge/results"
sys.path.insert(0, str(COMP))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


PINNED = {
    COMP / "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
    COMP / "lrc14_literal_fixed_sheet_allocation_thm2806.py":
        "311d0d85500f0c65945ebe5913f09d34a16293119c942b42eeaa854fbf85f71e",
    COMP / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.py":
        "2227f59c717095da0f2042096ada145de4e3661530c9aa2cc9020f42c8237a8b",
    RESULTS / "lrc14_literal_fixed_sheet_allocation_thm2806.out":
        "a970776ed95128b5745c1fd370af768778b409d931a57b68006e04271e00813f",
    RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out":
        "424fd2e83a618f862a5ee1b5f073a93fe236d10cdc5412eab1b54dec5e537eac",
}
for path, expected in PINNED.items():
    require(lf_hash(path) == expected, f"pinned dependency changed: {path.name}")


import lrc14_extended_carrier_endpoint_lib as endpoint_base
import lrc14_literal_fixed_sheet_allocation_thm2806 as allocation


P = 13
T = endpoint_base.T
UNIT = T // P
W = tuple(endpoint_base.WMOD)
SPEEDS = tuple(endpoint_base.endpoint.W)
SOURCE_ATOM = tuple(allocation.ATOM_INTERVAL)
TARGET_ATOM = tuple(
    endpoint + allocation.physical.SHIFT for endpoint in SOURCE_ATOM
)
BASE_ATOMS = (SOURCE_ATOM, TARGET_ATOM)
CORNER_NAMES = {
    ("S", "S"): "SS",
    ("D", "S"): "DS",
    ("D", "D"): "DD",
    ("S", "D"): "SD",
}
FACTOR_NAMES = (
    "guard", "u1", "u2", "u3", "u4", "u5", "target_a", "target_b",
    "deep_C3",
)


def rank_mod(matrix, prime):
    work = [[entry % prime for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(pivot_row, rows)
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            entry * inverse % prime for entry in work[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            scalar = work[row][column]
            work[row] = [
                (left - scalar * right) % prime
                for left, right in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def shifted_atom(atom, residue):
    left = (atom[0] + residue * UNIT) % T
    right = (atom[1] + residue * UNIT) % T
    require(left < right, "selected atom crossed the circle seam")
    return left, right


def factor_state(atom, index, ell_i):
    speed = SPEEDS[index]
    if index == 0:
        denominator = 91
        low = -13 - 7 * ell_i
        high = 13 - 7 * ell_i
    else:
        denominator = 182
        low = -13 - 14 * ell_i
        high = 13 - 14 * ell_i
    require(T % (denominator * speed) == 0, "comb grid no longer resolves")
    unit = T // (denominator * speed)
    step = denominator * unit
    length = (high - low) * unit
    base = (low % denominator) * unit
    left, right = atom
    first = (left - base) // step
    hits = []
    for translate in range(first - 2, first + 4):
        danger_left = base + translate * step
        danger_right = danger_left + length
        if max(left, danger_left) < min(right, danger_right):
            hits.append((danger_left, danger_right))
    if not hits:
        return "S"
    if (
        len(hits) == 1
        and hits[0][0] <= left
        and right <= hits[0][1]
    ):
        return "D"
    return "P"


def literal_pattern(atom, residue, ell):
    moved = shifted_atom(atom, residue)
    return "".join(
        factor_state(moved, index, ell[index])
        for index in range(len(SPEEDS))
    )


def normalized_value(table, key):
    if 0 in key:
        return 0
    return table[key]


def v4_cohomology():
    """Return normalized H^2(V4,F13) ranks for the trivial action."""
    nonzero = (1, 2, 3)
    one_keys = nonzero
    two_keys = tuple((g, h) for g in nonzero for h in nonzero)
    one_index = {key: index for index, key in enumerate(one_keys)}
    two_index = {key: index for index, key in enumerate(two_keys)}

    delta_one = []
    for g, h in two_keys:
        row = [0] * len(one_keys)
        row[one_index[g]] += 1
        row[one_index[h]] += 1
        if g ^ h:
            row[one_index[g ^ h]] -= 1
        delta_one.append(row)

    delta_two = []
    for g in range(4):
        for h in range(4):
            for k in range(4):
                row = [0] * len(two_keys)
                terms = (
                    (+1, (h, k)),
                    (-1, (g ^ h, k)),
                    (+1, (g, h ^ k)),
                    (-1, (g, h)),
                )
                for sign, key in terms:
                    if 0 not in key:
                        row[two_index[key]] += sign
                delta_two.append(row)
    b2_rank = rank_mod(delta_one, P)
    z2_dimension = len(two_keys) - rank_mod(delta_two, P)
    require(
        b2_rank == z2_dimension == 3,
        "normalized V4 cohomology ceased to vanish",
    )
    return len(one_keys), len(two_keys), b2_rank, z2_dimension


def cyclic_carry_obstruction():
    """Check that the C169 carry cocycle is not a split C13 coboundary."""
    one_keys = tuple(range(1, P))
    one_index = {key: index for index, key in enumerate(one_keys)}
    delta_one = []
    carry = []
    for h in range(P):
        for k in range(P):
            row = [0] * len(one_keys)
            if h:
                row[one_index[h]] += 1
            if k:
                row[one_index[k]] += 1
            total = (h + k) % P
            if total:
                row[one_index[total]] -= 1
            delta_one.append(row)
            carry.append((h + k) // P)

    # The standard carry is a normalized two-cocycle.
    cocycle_checks = 0
    for h in range(P):
        for k in range(P):
            for ell in range(P):
                left = (
                    (h + k) // P
                    + (((h + k) % P) + ell) // P
                )
                right = (
                    (k + ell) // P
                    + (h + ((k + ell) % P)) // P
                )
                require(left == right, "cyclic carry cocycle changed")
                cocycle_checks += 1

    rank = rank_mod(delta_one, P)
    augmented = rank_mod(
        [
            row + [carry[index]]
            for index, row in enumerate(delta_one)
        ],
        P,
    )
    require(
        rank == 11 and augmented == 12,
        "C169 carry unexpectedly became a split C13 coboundary",
    )
    return rank, augmented, cocycle_checks


def forward_cycle_steps(vertices):
    return tuple(
        (vertices[(index + 1) % len(vertices)] - vertices[index]) % P
        for index in range(len(vertices))
    )


def lifted_cycle_carry(start, steps):
    ancestry = 0
    residue = start
    for step in steps:
        ancestry += (residue + step) // P
        residue = (residue + step) % P
    return ancestry, residue


def main() -> None:
    require(
        P == 13
        and W == (1, 1, 1, 1, 1, 1, 0, 0, 0)
        and endpoint_base.endpoint.gf13_rank(
            endpoint_base.RELATION_ROWS
        ) == 6,
        "endpoint relation packet changed",
    )
    relation_rows = tuple(tuple(row) for row in endpoint_base.RELATION_ROWS)
    for vector in (
        W,
        *endpoint_base.REPS.values(),
    ):
        require(
            all(
                sum(left * right for left, right in zip(row, vector)) % P
                == 0
                for row in relation_rows
            ),
            "purported endpoint character left L-perp",
        )

    representatives = {}
    pair_projection = {}
    for address, canonical in endpoint_base.REPS.items():
        require(
            canonical[0] == canonical[5] == 0,
            "canonical guard/q5 coordinates changed",
        )
        for gauge in range(P):
            ell = tuple(
                (canonical[index] + gauge * W[index]) % P
                for index in range(len(W))
            )
            representatives[address, gauge] = ell
            pair_projection[address, gauge] = (ell[0], ell[5])
    require(
        len(set(representatives.values())) == P**3
        and set(pair_projection.values()) == {(s, s) for s in range(P)}
        and all(
            sum(pair_projection[key] == (s, s) for key in pair_projection)
            == P**2
            for s in range(P)
        ),
        "lawful representative census or diagonal projection changed",
    )

    expected_zero_orbit = (
        "SSSSSSSSD",
        "SSDSSSSSD",
        "SSDSSSSSD",
        "SSSSSSSSD",
        "SDSSSSSSD",
        "DDSSSSSSD",
        "DSSSSDSSD",
        "DSSSSDSSD",
        "DSSSSSSSD",
        "SSSSDSSSD",
        "SSSSDSSSD",
        "SSSSSSSSD",
        "SSSDSSSSD",
    )
    zero = (0,) * len(W)
    zero_orbits = tuple(
        tuple(literal_pattern(atom, q, zero) for q in range(P))
        for atom in BASE_ATOMS
    )
    require(
        all(orbit == expected_zero_orbit for orbit in zero_orbits)
        and expected_zero_orbit[3] == expected_zero_orbit[11]
        == "SSSSSSSSD"
        and expected_zero_orbit[7] == "DSSSSDSSD",
        "source/target zero-chart factor orbit changed",
    )

    # Exhaust all 36 two-factor projections over the full representative
    # gauge and residue plane at the distinguished address.  Source and
    # target atoms give the same table.
    zero_gauge_patterns = []
    for atom in BASE_ATOMS:
        table = {}
        for gauge in range(P):
            ell = tuple(gauge * value % P for value in W)
            for residue in range(P):
                table[gauge, residue] = literal_pattern(atom, residue, ell)
        zero_gauge_patterns.append(table)
    require(
        zero_gauge_patterns[0] == zero_gauge_patterns[1],
        "source/target pair-projection atlases differ",
    )
    pattern_table = zero_gauge_patterns[0]
    pair_images = {}
    for left in range(len(W)):
        for right in range(left + 1, len(W)):
            pair_images[left, right] = {
                (pattern[left], pattern[right])
                for pattern in pattern_table.values()
            }
    pair_size_histogram = {
        size: sum(len(image) == size for image in pair_images.values())
        for size in range(1, 5)
    }
    full_pairs = tuple(
        pair for pair, image in pair_images.items() if len(image) == 4
    )
    three_missing_sd = tuple(
        pair
        for pair, image in pair_images.items()
        if len(image) == 3 and ("S", "D") not in image
    )
    three_missing_dd = tuple(
        pair
        for pair, image in pair_images.items()
        if len(image) == 3 and ("D", "D") not in image
    )
    axis_pairs = tuple(
        pair for pair, image in pair_images.items() if len(image) == 2
    )
    constant_pairs = tuple(
        pair for pair, image in pair_images.items() if len(image) == 1
    )
    require(
        pair_size_histogram == {1: 3, 2: 18, 3: 14, 4: 1}
        and full_pairs == ((0, 1),)
        and three_missing_sd == ((0, 5),)
        and len(three_missing_dd) == 13
        and len(axis_pairs) == 18
        and constant_pairs == ((6, 7), (6, 8), (7, 8)),
        "36-pair classification changed",
    )
    danger_arcs = tuple(
        tuple(
            residue
            for residue, pattern in enumerate(expected_zero_orbit)
            if pattern[index] == "D"
        )
        for index in range(len(W))
    )
    require(
        danger_arcs
        == (
            (5, 6, 7, 8),
            (4, 5),
            (1, 2),
            (12,),
            (9, 10),
            (6, 7),
            (),
            (),
            tuple(range(P)),
        ),
        "factor danger arcs changed",
    )

    # The guard/u1 full square is unique at the distinguished zero address,
    # not globally over all 169 canonical addresses.  Audit every address
    # before stating the sharp scope.  Source and target give the same
    # q-labelled word orbit at every address.
    address_full_pairs = {}
    for address, canonical in endpoint_base.REPS.items():
        alpha, beta = address
        source_orbit = tuple(
            literal_pattern(SOURCE_ATOM, residue, canonical)
            for residue in range(P)
        )
        target_orbit = tuple(
            literal_pattern(TARGET_ATOM, residue, canonical)
            for residue in range(P)
        )
        require(
            source_orbit == target_orbit,
            "source/target all-address word orbits differ",
        )
        address_danger_sets = tuple(
            frozenset(
                residue
                for residue, pattern in enumerate(source_orbit)
                if pattern[index] == "D"
            )
            for index in range(len(W))
        )
        require(
            address_danger_sets[:6]
            == (
                frozenset((5, 6, 7, 8)),
                frozenset(((4 + alpha) % P, (5 + alpha) % P)),
                frozenset(((1 + beta) % P, (2 + beta) % P)),
                frozenset((12,)),
                frozenset((9, 10)),
                frozenset((6, 7)),
            )
            and all(
                danger in (frozenset(), frozenset(range(P)))
                for danger in address_danger_sets[6:8]
            )
            and address_danger_sets[8] == frozenset(range(P)),
            "all-address moving-arc formula changed",
        )
        source_pair_images = {
            (left, right): {
                (pattern[left], pattern[right])
                for pattern in source_orbit
            }
            for left in range(len(W))
            for right in range(left + 1, len(W))
        }
        target_pair_images = {
            (left, right): {
                (pattern[left], pattern[right])
                for pattern in target_orbit
            }
            for left in range(len(W))
            for right in range(left + 1, len(W))
        }
        require(
            source_pair_images == target_pair_images,
            "source/target all-address pair atlases differ",
        )
        address_full_pairs[address] = tuple(
            pair
            for pair, image in source_pair_images.items()
            if len(image) == 4
        )
    address_full_set_histogram = Counter(address_full_pairs.values())
    address_full_count_histogram = dict(sorted(Counter(
        len(pairs) for pairs in address_full_pairs.values()
    ).items()))
    address_full_pair_occurrences = dict(sorted(Counter(
        pair
        for pairs in address_full_pairs.values()
        for pair in pairs
    ).items()))
    require(
        len(address_full_set_histogram) == 38
        and address_full_count_histogram
        == {0: 54, 1: 64, 2: 38, 3: 10, 4: 3}
        and address_full_pairs[(0, 0)] == ((0, 1),)
        and address_full_pairs[(1, 0)] == ((1, 5),)
        and address_full_pair_occurrences
        == {
            (0, 1): 26,
            (0, 2): 26,
            (1, 2): 26,
            (1, 4): 26,
            (1, 5): 26,
            (2, 4): 26,
            (2, 5): 26,
        },
        "all-address full-pair census changed",
    )

    # At the zero address the unique honest two-bit square is guard versus
    # u1.  There are three clean copies in that q-orbit, according to which
    # all-safe/deep-danger representative is chosen for the SS corner.
    square_corner_order = (
        ("S", "S"), ("S", "D"), ("D", "D"), ("D", "S")
    )
    clean_residue_squares = tuple(
        (safe_residue, 4, 5, 8) for safe_residue in (0, 3, 11)
    )
    for square in clean_residue_squares:
        patterns = tuple(expected_zero_orbit[q] for q in square)
        require(
            tuple((pattern[0], pattern[1]) for pattern in patterns)
            == square_corner_order
            and all(
                pattern[2:8] == "SSSSSS" and pattern[8] == "D"
                for pattern in patterns
            ),
            "clean guard/u1 residue square changed",
        )
        steps = forward_cycle_steps(square)
        carry, returned = lifted_cycle_carry(square[0], steps)
        require(
            sum(steps) == P and carry == 1 and returned == square[0],
            "clean residue square stopped winding once",
        )

    # The same clean square occurs at fixed physical residue q7 after
    # adjoining the representative-gauge coordinate s.  It is a square on
    # the marked representative pullback, not a descended old-Ghat object.
    clean_q7_gauge_squares = tuple(
        (safe_gauge, 10, 11, 1) for safe_gauge in (4, 6, 9)
    )
    for square in clean_q7_gauge_squares:
        patterns = tuple(pattern_table[gauge, 7] for gauge in square)
        require(
            tuple((pattern[0], pattern[1]) for pattern in patterns)
            == square_corner_order
            and all(
                pattern[2:8] == "SSSSSS" and pattern[8] == "D"
                for pattern in patterns
            ),
            "fixed-q7 clean guard/u1 gauge square changed",
        )
        steps = forward_cycle_steps(square)
        carry, returned = lifted_cycle_carry(square[0], steps)
        require(
            sum(steps) == P and carry == 1 and returned == square[0],
            "clean representative square stopped winding once",
        )

    # Parent-proposed faithful parity rechart:
    #
    #     j = q5_bit XOR u1_bit,       q5_bit = j XOR u1_bit.
    #
    # It completes the guard/q5 punctured projection without forgetting the
    # original q5 truth.  At canonical q7, u1 is safe, hence j=q5.
    def parity_corner(pattern):
        guard = int(pattern[0] == "D")
        u1 = int(pattern[1] == "D")
        q5 = int(pattern[5] == "D")
        j = q5 ^ u1
        return guard, j, u1, q5

    parity_image = {
        parity_corner(pattern)[:2] for pattern in pattern_table.values()
    }
    require(
        parity_image == {(0, 0), (0, 1), (1, 0), (1, 1)}
        and all(
            q5 == (j ^ u1)
            for guard, j, u1, q5 in map(
                parity_corner, pattern_table.values()
            )
        ),
        "faithful q5/u1 parity rechart changed",
    )
    canonical_q7_parity = parity_corner(pattern_table[0, 7])
    require(
        canonical_q7_parity == (1, 1, 0, 1),
        "canonical q7 no longer has j=q5=danger",
    )

    # This projected square contains q7 and keeps the deep-C3 factor
    # dangerous on the same source/target atom chart.
    parity_residue_square = (3, 4, 7, 8)
    parity_gauge_square = (9, 10, 0, 1)
    for coordinate_square, fetch in (
        (
            parity_residue_square,
            lambda coordinate: expected_zero_orbit[coordinate],
        ),
        (
            parity_gauge_square,
            lambda coordinate: pattern_table[coordinate, 7],
        ),
    ):
        patterns = tuple(fetch(coordinate) for coordinate in coordinate_square)
        require(
            tuple(parity_corner(pattern)[:2] for pattern in patterns)
            == ((0, 0), (0, 1), (1, 1), (1, 0))
            and all(pattern[8] == "D" for pattern in patterns)
            and any(
                parity_corner(pattern) == canonical_q7_parity
                for pattern in patterns
            ),
            "q7 parity square changed",
        )
        steps = forward_cycle_steps(coordinate_square)
        carry, returned = lifted_cycle_carry(
            coordinate_square[0], steps
        )
        require(
            steps == (1, 3, 1, 8)
            and sum(steps) == P
            and carry == 1
            and returned == coordinate_square[0],
            "q7 parity square winding changed",
        )

    fixed_q7_parity_patterns = tuple(
        pattern_table[gauge, 7] for gauge in parity_gauge_square
    )
    fixed_q7_e3_vertices = tuple(
        gauge
        for gauge, pattern in zip(
            parity_gauge_square, fixed_q7_parity_patterns
        )
        if pattern[:8] == "SSSSSSSS" and pattern[8] == "D"
    )
    lawfully_cotransported_patterns = tuple(
        pattern_table[gauge, (7 - gauge) % P]
        for gauge in parity_gauge_square
    )
    require(
        fixed_q7_e3_vertices == (9,)
        and len(set(lawfully_cotransported_patterns)) == 1
        and lawfully_cotransported_patterns[0] == expected_zero_orbit[7],
        "fixed-fibre versus descended representative typing changed",
    )

    # The full word is not blind to the basepoint if its *directed unit
    # transitions* are retained.  Factor u3 (index 3, speed 40) is dangerous
    # only at q12, so its unique D->S exit is exactly q12->q0.
    marker_index = 3

    def marker_exits(gauge, start, length):
        total = 0
        for offset in range(length):
            source = pattern_table[gauge, (start + offset) % P]
            target = pattern_table[gauge, (start + offset + 1) % P]
            total += source[marker_index] == "D" and target[marker_index] == "S"
        return total

    canonical_marker_checks = 0
    marked_marker_checks = 0
    for residue in range(P):
        for step in range(P):
            require(
                marker_exits(0, residue, step)
                == (residue + step) // P,
                "canonical u3 marker stopped computing carry",
            )
            canonical_marker_checks += 1
            for gauge in range(P):
                # Under lawful marked transport ell->ell+sW, the physical
                # residue is q-s.  The shifted singleton marker then computes
                # the original q-carry, not a chart-dependent carry.
                require(
                    marker_exits(
                        gauge, (residue - gauge) % P, step
                    )
                    == (residue + step) // P,
                    "marked-gauge u3 carry marker lost covariance",
                )
                marked_marker_checks += 1
    require(
        canonical_marker_checks == P**2
        and marked_marker_checks == P**3
        and marker_exits(0, 3, 4) == 0
        and marker_exits(0, 3, 8) == 0
        and marker_exits(0, 11, 9) == 1,
        "direct/via carry-marker triangle changed",
    )

    # The positive event count is additive as a path-semigroup observable.
    # Comparing concatenation with the reduced step gives exactly the central
    # carry in the THM-2851 composition law.
    semigroup_checks = 0
    reduced_composition_checks = 0
    for residue in range(P):
        for first_step in range(P):
            for second_step in range(P):
                concatenated = marker_exits(
                    0, residue, first_step + second_step
                )
                split = (
                    marker_exits(0, residue, first_step)
                    + marker_exits(
                        0,
                        (residue + first_step) % P,
                        second_step,
                    )
                )
                require(
                    concatenated == split,
                    "positive marker count lost concatenation additivity",
                )
                semigroup_checks += 1
                reduced_step = (first_step + second_step) % P
                central = (first_step + second_step) // P
                require(
                    split
                    == marker_exits(0, residue, reduced_step) + central,
                    "marker count stopped deriving the reduced lift law",
                )
                reduced_composition_checks += 1

    # Exhaust all nine bits and both oriented toggle events at the zero
    # address.  Then repeat at every canonical address.  Shifted u1/u2
    # events can coincide with carry at some addresses, but u3 D->S is the
    # unique candidate present uniformly at all 169 addresses.
    def carry_event_candidates(orbit):
        candidates = []
        for index in range(len(W)):
            for orientation in ("S_to_D", "D_to_S"):
                matches = True
                for residue in range(P):
                    for step in range(P):
                        count = 0
                        for offset in range(step):
                            source = orbit[(residue + offset) % P][index]
                            target = orbit[
                                (residue + offset + 1) % P
                            ][index]
                            count += (
                                source == orientation[0]
                                and target == orientation[5]
                            )
                        if count != (residue + step) // P:
                            matches = False
                if matches:
                    candidates.append((index, orientation))
        return tuple(candidates)

    full_cycle_toggle_counts = []
    for index in range(len(W)):
        enter = exit_count = 0
        for residue in range(P):
            source = expected_zero_orbit[residue][index]
            target = expected_zero_orbit[(residue + 1) % P][index]
            enter += source == "S" and target == "D"
            exit_count += source == "D" and target == "S"
        full_cycle_toggle_counts.append((enter, exit_count))
    oriented_event_candidates = carry_event_candidates(expected_zero_orbit)
    require(
        tuple(full_cycle_toggle_counts)
        == ((1, 1), (1, 1), (1, 1), (1, 1), (1, 1), (1, 1),
            (0, 0), (0, 0), (0, 0))
        and oriented_event_candidates == ((3, "D_to_S"),),
        "zero-address oriented-factor carry detector uniqueness changed",
    )

    address_marker_candidates = {
        address: carry_event_candidates(tuple(
            literal_pattern(SOURCE_ATOM, residue, canonical)
            for residue in range(P)
        ))
        for address, canonical in endpoint_base.REPS.items()
    }
    address_marker_candidate_histogram = Counter(
        address_marker_candidates.values()
    )
    expected_marker_candidate_histogram = Counter({
        ((3, "D_to_S"),): 121,
        ((1, "D_to_S"), (3, "D_to_S")): 11,
        ((1, "S_to_D"), (3, "D_to_S")): 11,
        ((2, "D_to_S"), (3, "D_to_S")): 11,
        ((2, "S_to_D"), (3, "D_to_S")): 11,
        (
            (1, "D_to_S"), (2, "D_to_S"), (3, "D_to_S")
        ): 1,
        (
            (1, "D_to_S"), (2, "S_to_D"), (3, "D_to_S")
        ): 1,
        (
            (1, "S_to_D"), (2, "D_to_S"), (3, "D_to_S")
        ): 1,
        (
            (1, "S_to_D"), (2, "S_to_D"), (3, "D_to_S")
        ): 1,
    })
    address_marker_occurrences = dict(sorted(Counter(
        candidate
        for candidates in address_marker_candidates.values()
        for candidate in candidates
    ).items()))
    extra_marker_addresses = sum(
        candidates != ((3, "D_to_S"),)
        for candidates in address_marker_candidates.values()
    )
    require(
        address_marker_candidate_histogram
        == expected_marker_candidate_histogram
        and address_marker_candidates[(0, 0)] == ((3, "D_to_S"),)
        and address_marker_candidates[(7, 0)]
        == ((1, "D_to_S"), (3, "D_to_S"))
        and all(
            (3, "D_to_S") in candidates
            for candidates in address_marker_candidates.values()
        )
        and extra_marker_addresses == 48
        and address_marker_occurrences
        == {
            (1, "D_to_S"): 13,
            (1, "S_to_D"): 13,
            (2, "D_to_S"): 13,
            (2, "S_to_D"): 13,
            (3, "D_to_S"): 169,
        },
        "all-address oriented-factor carry census changed",
    )
    for (alpha, beta), candidates in address_marker_candidates.items():
        expected_candidates = []
        if alpha == 9:
            expected_candidates.append((1, "S_to_D"))
        if alpha == 7:
            expected_candidates.append((1, "D_to_S"))
        if beta == 12:
            expected_candidates.append((2, "S_to_D"))
        if beta == 10:
            expected_candidates.append((2, "D_to_S"))
        expected_candidates.append((3, "D_to_S"))
        require(
            candidates == tuple(expected_candidates),
            "wrap-boundary formula stopped explaining marker candidates",
        )

    def reversed_marker_exits(start, length):
        total = 0
        for offset in range(length):
            source = expected_zero_orbit[(start - offset) % P]
            target = expected_zero_orbit[(start - offset - 1) % P]
            total += (
                source[marker_index] == "D"
                and target[marker_index] == "S"
            )
        return total

    require(
        reversed_marker_exits(7, 9) == 1
        and marker_exits(0, 11, 9) == 1,
        "truth-word reversal hostile changed",
    )
    require(
        all(
            representative[marker_index] == 0
            for representative in endpoint_base.REPS.values()
        ),
        "u3 marker acquired address dependence",
    )

    direct_path = tuple(range(3, 8))
    via_path = tuple(range(3, 12)) + (12,) + tuple(range(0, 8))
    direct_words = tuple(expected_zero_orbit[q] for q in direct_path)
    via_words = tuple(expected_zero_orbit[q] for q in via_path)

    def guard_q5_level(pattern):
        pair = pattern[0] + pattern[5]
        return {"SS": 0, "DS": 1, "DD": 2}[pair]

    direct_levels = tuple(map(guard_q5_level, direct_words))
    via_levels = tuple(map(guard_q5_level, via_words))

    def up_down(levels):
        up = down = 0
        for left, right in zip(levels, levels[1:]):
            up += max(0, right - left)
            down += max(0, left - right)
        return up, down

    require(
        direct_levels == (0, 0, 1, 2, 2)
        and via_levels
        == (0, 0, 1, 2, 2, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 2, 2)
        and up_down(direct_levels) == (2, 0)
        and up_down(via_levels) == (4, 2)
        and direct_levels[-1] - direct_levels[0]
        == via_levels[-1] - via_levels[0] == 2,
        "guard/q5 filtration path census changed",
    )
    consecutive_word_pairs = tuple(
        (
            expected_zero_orbit[q],
            expected_zero_orbit[(q + 1) % P],
        )
        for q in range(P)
    )
    marker_pair = (
        expected_zero_orbit[12],
        expected_zero_orbit[0],
    )
    require(
        consecutive_word_pairs.count(marker_pair) == 1
        and consecutive_word_pairs[12] == marker_pair,
        "full-word q12/q0 basepoint marker ceased to be unique",
    )

    # All four fixed-q7 truth-polarized endpoint restrictions are the same
    # physical atom.  Therefore their raw endpoint coefficient is identical:
    # the Boolean square is coefficient-flat before any external ancestry
    # lift is attached.
    q7_atom = shifted_atom(TARGET_ATOM, 7)
    q7_endpoint_values = tuple(
        allocation.endpoint_sum(
            ((*q7_atom, 1),), -endpoint_base.Y0, embedding
        )
        for embedding in endpoint_base.endpoint.MODS
    )
    require(
        all(q7_endpoint_values),
        "fixed-q7 endpoint atom coefficient vanished",
    )

    fixed_counts = {"SS": 0, "DS": 0, "DD": 0, "SD": 0}
    transported_counts = {"SS": 0, "DS": 0, "DD": 0, "SD": 0}
    covariance_checks = 0
    for atom in BASE_ATOMS:
        for address, canonical in endpoint_base.REPS.items():
            for gauge in range(P):
                ell = representatives[address, gauge]
                for residue in range(P):
                    fixed = literal_pattern(atom, residue, ell)
                    fixed_corner = fixed[0] + fixed[5]
                    fixed_counts[fixed_corner] += 1

                    # THM-2763/2806's lawful marked representative transport:
                    # ell -> ell+sW is accompanied by atom q -> q-s.
                    transported = literal_pattern(
                        atom, (residue - gauge) % P, ell
                    )
                    canonical_pattern = literal_pattern(
                        atom, residue, canonical
                    )
                    require(
                        transported == canonical_pattern,
                        "factorwise marked representative covariance failed",
                    )
                    transported_corner = transported[0] + transported[5]
                    transported_counts[transported_corner] += 1
                    covariance_checks += 1
    expected_counts = {
        "SS": 39546,
        "DS": 8788,
        "DD": 8788,
        "SD": 0,
    }
    require(
        fixed_counts == transported_counts == expected_counts
        and covariance_checks == 2 * P**4,
        "guard/q5 corner census changed",
    )

    # Since every residue is already exhausted, reversing the residue
    # orientation, adding an arbitrary residue offset, and swapping either
    # source/target origin only permute the rows above.  Record the direct
    # nesting law on the zero chart.
    zero_pair_orbit = tuple(
        pattern[0] + pattern[5] for pattern in expected_zero_orbit
    )
    guard_danger = tuple(
        q for q, corner in enumerate(zero_pair_orbit) if corner[0] == "D"
    )
    q5_danger = tuple(
        q for q, corner in enumerate(zero_pair_orbit) if corner[1] == "D"
    )
    require(
        zero_pair_orbit
        == ("SS", "SS", "SS", "SS", "SS", "DS", "DD",
            "DD", "DS", "SS", "SS", "SS", "SS")
        and guard_danger == (5, 6, 7, 8)
        and q5_danger == (6, 7)
        and set(q5_danger) < set(guard_danger),
        "literal danger nesting changed",
    )

    # A binary decoration does not change the literal projection.  XORing
    # the decoration into the guard label manufactures all four *formal*
    # labels, but precisely by changing the guard predicate.
    physical_corners = ((0, 0), (1, 0), (1, 1))
    decorated = tuple(
        (corner, bit) for corner in physical_corners for bit in (0, 1)
    )
    literal_projection = {corner for corner, _bit in decorated}
    xor_projection = {
        (corner[0] ^ bit, corner[1]) for corner, bit in decorated
    }
    xor_guard_mismatches = sum(
        (corner[0] ^ bit) != corner[0] for corner, bit in decorated
    )
    require(
        literal_projection == set(physical_corners)
        and (0, 1) not in literal_projection
        and xor_projection == {(0, 0), (1, 0), (1, 1), (0, 1)}
        and xor_guard_mismatches == 3,
        "binary-decoration hostile changed",
    )

    v4_dimensions = v4_cohomology()
    split_vertex_phase_checks = 0
    for phase_10 in range(P):
        for phase_11 in range(P):
            for phase_01 in range(P):
                phase = {
                    0: 0,
                    1: phase_10,
                    3: phase_11,
                    2: phase_01,
                }
                boundary = (
                    (phase[1] - phase[0])
                    + (phase[3] - phase[1])
                    - (phase[2] - phase[0])
                    - (phase[3] - phase[2])
                ) % P
                require(boundary == 0, "split V4 vertex phase curved")
                split_vertex_phase_checks += 1
    require(
        split_vertex_phase_checks == P**3,
        "split V4 phase census changed",
    )

    carry_rank, carry_augmented_rank, carry_checks = (
        cyclic_carry_obstruction()
    )
    q3 = (0, 3)

    def natural_lift(state, step):
        ancestry, residue = state
        return (
            ancestry + (residue + step) // P,
            (residue + step) % P,
        )

    q11 = natural_lift(q3, 8)
    q7_via = natural_lift(q11, 9)
    q7_direct = natural_lift(q3, 4)
    require(
        q11 == (0, 11)
        and q7_via == (1, 7)
        and q7_direct == (0, 7)
        and q7_via[0] - q7_direct[0] == 1,
        "q3/q11/q7 ancestry carry changed",
    )
    ancestry_holonomy_exponent = 3 * (
        q7_via[0] - q7_direct[0]
    ) % P
    require(
        ancestry_holonomy_exponent == 3,
        "character-three ancestry holonomy changed",
    )

    smallest_phase_fibre = next(
        order for order in range(2, P + 1) if order % P == 0
    )
    require(
        smallest_phase_fibre == P
        and all(order % P for order in range(2, P)),
        "smaller phase fibre acquired order-thirteen character",
    )

    semantic_output = (
        RESULTS / "lrc14_q3_q11_q7_bockstein_holonomy_thm2851.out"
    ).read_text(encoding="utf-8").replace("\r\n", "\n")
    require(
        "semantic_leg=THM2835_PINNED:"
        "449_QA(q11,a)_to_QAB(q7,a+1)" in semantic_output
        and "law=L_k*L_h=T^floor((h+k)/13)*L_(h+k_mod13)"
        in semantic_output
        and "faithful_first_carry_sidecar_min_states=13" in semantic_output
        and "faithful_q_and_carry_joint_min_states=169" in semantic_output,
        "pinned QA/QAB ancestry typing changed",
    )

    print("THM-2878 ENDPOINT-FACTOR EXIT CARRY TRANSDUCER")
    print(
        f"endpoint_relation_rank=6; Lperp_size={len(representatives)}; "
        "Ghat_size=169; representative_law=ell(a,b)+sW"
    )
    print(
        f"Wmod13={W}; canonical_(ell0,ell5)=(0,0); "
        "all_lawful_(ell0,ell5)={(s,s):s in F13}; "
        "each_diagonal_pair_multiplicity=169"
    )
    print(
        f"source_zero_orbit={zero_orbits[0]}; "
        f"target_zero_orbit={zero_orbits[1]}"
    )
    print(
        f"all_36_pair_sizes_zero_address={pair_size_histogram}; "
        f"unique_full_pair_zero_address={full_pairs}; "
        f"three_corner_missing_SD={three_missing_sd}; "
        f"three_corner_missing_DD={three_missing_dd}; "
        f"axis_pairs={axis_pairs}; constant_pairs={constant_pairs}"
    )
    print(
        f"all_address_full_pair_count_histogram="
        f"{address_full_count_histogram}; "
        f"distinct_full_pair_sets={len(address_full_set_histogram)}; "
        f"full_pairs_at_(0,0)={address_full_pairs[(0, 0)]}; "
        f"full_pairs_at_(1,0)={address_full_pairs[(1, 0)]}; "
        f"pair_occurrences={address_full_pair_occurrences}"
    )
    print(
        "all_address_arc_law=u1={4+alpha,5+alpha};"
        "u2={1+beta,2+beta};other_nonconstant_arcs_fixed;"
        "target_bits_constant_per_address;"
        "full_pairs_are_crossing_arcs;"
        "each_of_7_pair_types_has_2_boundary_placements*13=26"
    )
    print(
        f"factor_order={FACTOR_NAMES}; danger_arcs={danger_arcs}; "
        "unique_crossing_zero_address=guard{5,6,7,8}_with_u1{4,5}"
    )
    print(
        f"clean_guard_u1_residue_squares={clean_residue_squares}; "
        f"clean_fixed_q7_gauge_squares={clean_q7_gauge_squares}; "
        "every_boundary_winds_once"
    )
    print(
        "faithful_parity_rechart=j=q5_XOR_u1; inverse=q5=j_XOR_u1; "
        f"parity_residue_square={parity_residue_square}; "
        f"fixed_q7_parity_gauge_square={parity_gauge_square}; "
        f"fixed_q7_words={fixed_q7_parity_patterns}; "
        f"fixed_q7_PAT_E3_vertices={fixed_q7_e3_vertices}"
    )
    print(
        "fixed_q7_square_scope=representative-fibre pullback; lawful "
        "co-transport q->q-s collapses its four words to "
        f"{lawfully_cotransported_patterns[0]}; "
        f"common_atom_endpoint_values={q7_endpoint_values}; "
        "truth-polarized coefficient boundary=flat"
    )
    print(
        f"guard_danger_residues={guard_danger}; "
        f"q5_danger_residues={q5_danger}; "
        "nesting=q5_danger_strict_subset_guard_danger"
    )
    print(
        f"universe=2_atoms*169_addresses*13_representative_gauges*"
        f"13_residues={covariance_checks}; "
        f"fixed_corner_counts={fixed_counts}; "
        f"transported_corner_counts={transported_counts}"
    )
    print(
        "origin_residue_hostile=source/target swap, q->-q+c, and address "
        "changes only permute an exhausted universe; literal_SD_count=0"
    )
    print(
        f"marked_gauge_covariance_checks={covariance_checks}; "
        "law=(ell,atom_q)->(ell+sW,atom_(q-s)); "
        "all_nine_factor_patterns_preserved=1"
    )
    print(
        f"binary_fibre=U3xC2_has_{len(decorated)}_states; "
        f"literal_projection={sorted(literal_projection)}; "
        f"xor_regrading_projection={sorted(xor_projection)}; "
        f"xor_guard_mismatches={xor_guard_mismatches}/6"
    )
    print(
        f"guard_q5_filtration_direct={direct_levels}; "
        f"via={via_levels}; direct_(up,down)={up_down(direct_levels)}; "
        f"via_(up,down)={up_down(via_levels)}; net_thresholds=2_both"
    )
    print(
        f"directed_basepoint_marker=factor{marker_index}_{FACTOR_NAMES[marker_index]}"
        "_D_to_S; danger_arc={12}; edge=q12_to_q0; "
        f"canonical_checks={canonical_marker_checks}; "
        f"marked_gauge_checks={marked_marker_checks}; "
        "unique_among_18_oriented_toggles_zero_address="
        f"{oriented_event_candidates}"
    )
    print(
        "address_uniform_marker=(3,D_to_S); "
        f"candidate_set_histogram={dict(address_marker_candidate_histogram)}; "
        f"extra_candidate_addresses={extra_marker_addresses}; "
        f"candidate_occurrences={address_marker_occurrences}; "
        "minimal_extra_witness_address_(7,0)=((1,D_to_S),(3,D_to_S))"
    )
    print(
        "extra_marker_boundary_law=u1_D_to_S@alpha7,"
        "u1_S_to_D@alpha9,u2_D_to_S@beta10,u2_S_to_D@beta12;"
        "extra_addresses=26+26-4=48;u3_fixed_by_ell3=0"
    )
    print(
        f"positive_path_semigroup_additivity={semigroup_checks}; "
        f"reduced_lift_compositions={reduced_composition_checks}; "
        "law=L_kL_h=T^floor((h+k)/13)L_(h+k_mod13); "
        "reversal_hostile=positive_and_reversed_paths_both_have_one_D_to_S_"
        "exit; inherited_positive_h_orientation_is_load_bearing"
    )
    print(
        "vertex_coboundary_boundary=signed_enter_minus_exit_cancels; "
        f"full_cycle_(enter,exit)_per_factor={tuple(full_cycle_toggle_counts)}; "
        "positive_D_to_S_count_is_noninvertible_path-semigroup_observable"
    )
    print(
        f"V4_normalized_cochains=(C1={v4_dimensions[0]},"
        f"C2={v4_dimensions[1]},B2_rank={v4_dimensions[2]},"
        f"Z2_dim={v4_dimensions[3]},H2_dim=0); "
        f"split_vertex_phase_boundaries_flat={split_vertex_phase_checks}"
    )
    print(
        f"C13_carry=(delta1_rank={carry_rank},"
        f"augmented_rank={carry_augmented_rank},"
        f"cocycle_checks={carry_checks}); "
        "triangle=q3--8->q11--9->q7 versus q3--4->q7; "
        f"carry=1; chi3_holonomy_exponent={ancestry_holonomy_exponent}"
    )
    print(
        f"smallest_phase_fibre_with_order13={smallest_phase_fibre}; "
        "minimal_faithful_extra_state=initial_C13_ancestry_coordinate; "
        "the_full_word_supplies_its_carry_cocycle_but_not_its_initial_"
        "ancestry_state"
    )
    print(
        "verdict=literal_guard_q5_SD_is_globally_absent, but guard_u1_is_the_"
        "unique_honest_square_at_the_zero_address and "
        "(guard,j=q5_XOR_u1,u1) is a faithful q7-retaining zero-address "
        "rechart; its vertex coefficient is flat, while the directed u3 "
        "exit is the unique address-uniform marker and reconstructs exactly "
        "the nonsplit carry cocycle"
    )
    print(
        "scope=exact selected source/target atom orbit and endpoint-character "
        "representative gauge; inherited positive-h orientation and QA/QAB "
        "typing pinned, but no map from the marker/initial ancestry fibre to "
        "a common physical QA-to-QAB current, E3 contraction, row exclusion, "
        "or LRC14 proof"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
