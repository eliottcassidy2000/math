#!/usr/bin/env python3
"""Independent hostile audit of the r=5/U_full owner-node coupling.

This checker fixes its object from THM-2471/2594 and THM-3514 before
consulting the candidate implementation.  It imports only those canonical
dependencies, then uses a separately written route:

* rebuild the unsplit r=5 source profiles with an independent transfer-fold;
* desheet U_full endpoint intervals directly, without candidate atom buckets;
* pull Q back under x=13y and integrate it with a prefix primitive, rather
  than the candidate's per-segment x-sweep;
* compare literal guarded endpoint sets pointwise with restored unguarded
  masks; and
* invert and Fourier analyse the resulting 7 x 13 table.

The primary hostile is structural.  At w_u=(y+u)/13, the U_full owner-in
condition is

    ||13 w_u|| < 1/14  <=>  ||y|| < 1/14.

It is therefore exactly the source cell ell=0 almost everywhere.  Complete
Fourier support can only be the spectrum of delta_0(ell) tensor R(t), not
genuine two-coordinate cell/residue mixing.

The surviving object is an explicit one-common-base owner-node integrand
candidate.  It is not an exact address C(a;X,m), a temporal/source-arrival
intertwiner, U_clock chronology, a physical current, a row exclusion, or an
LRC(14) result.
"""

from __future__ import annotations

from bisect import bisect_right
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    "04-computation/"
    "lrc_r5_ufull_owner_node_common_base_current_"
    "independent_audit_20260816.py"
)
OUTPUT = (
    "05-knowledge/results/"
    "lrc_r5_ufull_owner_node_common_base_current_"
    "independent_audit_20260816.out"
)

SOURCE_PATH = ROOT / "04-computation/lrc14_stage2_theta_contraction_opus_20260728.py"
SOURCE_SHA256 = "09c43af0a0a56c7a0833bbfd13ed6a96bc5a7a3718aa1bc6b77a144bde101a06"
ENDPOINT_AUDIT_PATH = (
    ROOT
    / "04-computation/"
    "lrc_ufull_owner_boundary_k4xf13_endpoint_factorization_"
    "independent_audit_20260816.py"
)
ENDPOINT_AUDIT_SHA256 = "f89be10c65bb77270199f9399b155d5a2c82c0da121b3e8589fe3c1f7e9824fc"
ENDPOINT_PARENT_SHA256 = "ad2a620cdc238f28e3384698b2c612f38cdf2566bd56b76d1cbabcc03107ec0b"

P = 13
Q = 7
JOINT_COORDINATE = 9684279613402457983920
JOINT_ORDER = 125895634974231953790960
JOINT_PRIME = 755373809845391722745761
JOINT_GENERATOR = 23
JOINT_ROOT = 148035889
ORDER_FACTORS = (
    (2, 4), (3, 3), (5, 1), (7, 2), (11, 1), (13, 8),
    (53, 1), (61, 1), (131, 1), (313, 1),
)
PRIME_MINUS_ONE_FACTORS = (
    (2, 5), (3, 4), (5, 1), (7, 2), (11, 1), (13, 8),
    (53, 1), (61, 1), (131, 1), (313, 1),
)
Y_FREQUENCY = 57122
CHARACTER_EXPONENT = P * Y_FREQUENCY  # exponent on the C-grid in an order-13C root
SOURCE_PROFILE_SHA256 = "2de29f9be7fd16ceb4be5d15f7a71aa3d09f2907ec39f7af0b2017eadf3c18d2"
SOURCE_TOTAL_NUMERATOR = 168908463464745122312762880
EXPECTED_ROLE_VALUES = (
    125385278409587426725290,
    657486478079327229022863,
    223272610175651920448188,
)
EXPECTED_ERASURE_VALUES = (
    471060960989539924924555,
    594285905723663416170205,
    632148865111268231500111,
)
EXPECTED_FIXED = 317699132065964946247468
EXPECTED_ERASURE_FIXED = 454454155013190282848607
EXPECTED_SEMANTIC_SHA256 = "88d4be52bcb16a52ab2656ff0c0b6bf70e33a2174652ee7e1df62376426f24e6"

CONTROL_TRIPLES = ((0, 0, 0), (0, 1, 6), (1, 0, 6), (6, 6, 12), (12, 12, 0))
DANGER_EXPECTED = {
    "left": frozenset((12, 0, 1)),
    "right": frozenset((11, 12, 0)),
}


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def load_module(path: Path, name: str, expected_hash: str):
    require(lf_sha256(path) == expected_hash, (name, "dependency hash drift"))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, (name, "module loader"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SRC = load_module(SOURCE_PATH, "independent_owner_source", SOURCE_SHA256)
END_AUDIT = load_module(
    ENDPOINT_AUDIT_PATH, "independent_owner_endpoint", ENDPOINT_AUDIT_SHA256
)
EP = END_AUDIT.M
require(lf_sha256(END_AUDIT.PARENT_PATH) == ENDPOINT_PARENT_SHA256,
        "THM-3514 endpoint parent hash drift")


def factor_product(factors: tuple[tuple[int, int], ...]) -> int:
    value = 1
    for prime, power in factors:
        value *= prime**power
    return value


def is_small_prime(number: int) -> bool:
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 1
    return True


def split_field_certificate() -> tuple[int, int, int]:
    require(factor_product(PRIME_MINUS_ONE_FACTORS) == JOINT_PRIME - 1,
            "incomplete p-1 factorization")
    require(all(is_small_prime(prime) for prime, _power in PRIME_MINUS_ONE_FACTORS),
            "nonprime factor in p-1 certificate")
    require(pow(JOINT_GENERATOR, JOINT_PRIME - 1, JOINT_PRIME) == 1,
            "Lucas Fermat gate")
    for prime, _power in PRIME_MINUS_ONE_FACTORS:
        require(
            gcd(
                pow(JOINT_GENERATOR, (JOINT_PRIME - 1) // prime, JOINT_PRIME) - 1,
                JOINT_PRIME,
            ) == 1,
            ("Lucas primitive gate", prime),
        )
    require(pow(JOINT_GENERATOR, 6, JOINT_PRIME) == JOINT_ROOT,
            "declared root is not generator^6")
    require(factor_product(ORDER_FACTORS) == JOINT_ORDER, "order factorization")
    require(pow(JOINT_ROOT, JOINT_ORDER, JOINT_PRIME) == 1, "root power")
    for prime, _power in ORDER_FACTORS:
        require(pow(JOINT_ROOT, JOINT_ORDER // prime, JOINT_PRIME) != 1,
                ("root loses order", prime))
    return JOINT_PRIME, JOINT_GENERATOR, JOINT_ROOT


def merge_intervals(intervals: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    answer: list[list[int]] = []
    for left, right in sorted(intervals):
        require(left < right, ("empty interval", left, right))
        if answer and left <= answer[-1][1]:
            if right > answer[-1][1]:
                answer[-1][1] = right
        else:
            answer.append([left, right])
    return tuple((left, right) for left, right in answer)


def profile_value(profile: tuple[tuple[int, ...], tuple[int, ...]], point: int) -> int:
    starts, values = profile
    return values[bisect_right(starts, point) - 1]


def fold_weighted(
    pieces: list[tuple[int, int, int]], multiplier: int, grid: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Independently construct the numerator profile of multiplier*P_m h."""
    constant = 0
    events: dict[int, int] = {}
    input_mass = 0
    for left, right, weight in pieces:
        if not weight:
            continue
        require(0 <= left < right <= grid, ("fold piece", left, right))
        input_mass += weight * (right - left)
        quotient, remainder = divmod(multiplier * (right - left), grid)
        constant += quotient * weight
        if not remainder:
            continue
        start = multiplier * left % grid
        stop = start + remainder
        if stop <= grid:
            events[start] = events.get(start, 0) + weight
            events[stop] = events.get(stop, 0) - weight
        else:
            events[start] = events.get(start, 0) + weight
            events[grid] = events.get(grid, 0) - weight
            events[0] = events.get(0, 0) + weight
            events[stop - grid] = events.get(stop - grid, 0) - weight
    current = constant
    starts: list[int] = []
    values: list[int] = []
    for position in sorted(set(events) | {0}):
        current += events.get(position, 0)
        if starts and values[-1] == current:
            continue
        starts.append(position)
        values.append(current)
    require(starts and starts[0] == 0 and all(value >= 0 for value in values),
            "fold profile shape")
    mass = 0
    for index, left in enumerate(starts):
        right = starts[index + 1] if index + 1 < len(starts) else grid
        mass += values[index] * (right - left)
    require(mass == multiplier * input_mass, ("fold mass", mass, input_mass))
    return tuple(starts), tuple(values)


def intersect_profile_with_set(
    profile: tuple[tuple[int, ...], tuple[int, ...]],
    intervals: list[tuple[int, int]],
    grid: int,
) -> list[tuple[int, int, int]]:
    starts, values = profile
    pieces: list[tuple[int, int, int]] = []
    for left, right in intervals:
        index = bisect_right(starts, left) - 1
        cursor = left
        while cursor < right:
            stop = starts[index + 1] if index + 1 < len(starts) else grid
            stop = min(stop, right)
            if values[index]:
                pieces.append((cursor, stop, values[index]))
            cursor = stop
            index += 1
    return pieces


def pull_profile_to_root(
    profile: tuple[tuple[int, ...], tuple[int, ...]], root: int, grid: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    require(grid % P == 0, "source grid is not 13-divisible")
    starts, values = profile
    sheet = grid // P
    lower = root * sheet
    upper = lower + sheet
    index = bisect_right(starts, lower) - 1
    out_starts: list[int] = []
    out_values: list[int] = []
    cursor = lower
    while cursor < upper:
        stop = starts[index + 1] if index + 1 < len(starts) else grid
        stop = min(stop, upper)
        mapped = P * cursor - root * grid
        value = values[index]
        if not out_starts or out_values[-1] != value:
            out_starts.append(mapped)
            out_values.append(value)
        cursor = stop
        index += 1
    require(out_starts[0] == 0, ("root profile start", root))
    return tuple(out_starts), tuple(out_values)


def source_profiles():
    grid = SRC.T_DEN
    e_intervals = SRC.build_set(SRC.PAT_E, SRC.ZELL)
    q_intervals = SRC.build_set(SRC.PAT_QB, SRC.ZELL)
    require(Fraction(SRC.measure(e_intervals), grid) == Fraction(1882176, 28589561),
            "source E measure")
    require(Fraction(SRC.measure(q_intervals), grid)
            == Fraction(4839079319, 190921088358), "source Q measure")

    packet = fold_weighted(
        [(left, right, 1) for left, right in e_intervals], SRC.RPKT, grid
    )
    f_pieces = intersect_profile_with_set(packet, q_intervals, grid)
    f_mass = sum(weight * (right - left) for left, right, weight in f_pieces)
    require(Fraction(f_mass, SRC.RPKT * grid) == SRC.RHO_B, "source f anchor")
    whole_u = fold_weighted(f_pieces, SRC.DCOLL, grid)
    whole_v = fold_weighted(
        [(left, right, 1) for left, right in e_intervals], SRC.DCOLL, grid
    )

    # Independent fold agrees with the canonical inherited constructor.
    require(whole_u == tuple(tuple(row) for row in SRC.weighted_fold(
        f_pieces, SRC.DCOLL, grid)), "source U fold disagreement")
    require(whole_v == tuple(tuple(row) for row in SRC.weighted_fold(
        [(left, right, 1) for left, right in e_intervals], SRC.DCOLL, grid
    )), "source V fold disagreement")

    roots_u = tuple(pull_profile_to_root(whole_u, root, grid) for root in range(P))
    roots_v = tuple(pull_profile_to_root(whole_v, root, grid) for root in range(P))
    require(roots_u == tuple(
        tuple(tuple(row) for row in SRC.extract_window(
            whole_u[0], whole_u[1], root, P, grid
        )) for root in range(P)
    ), "source U root pull disagreement")
    require(roots_v == tuple(
        tuple(tuple(row) for row in SRC.extract_window(
            whole_v[0], whole_v[1], root, P, grid
        )) for root in range(P)
    ), "source V root pull disagreement")

    profile_digest = digest_json((roots_u, roots_v))
    require(profile_digest == SOURCE_PROFILE_SHA256,
            ("source profile digest", profile_digest))
    boundaries = sorted(
        {0, grid}
        | {position for profile in roots_u + roots_v for position in profile[0]}
    )
    total = diagonal = 0
    support_types: dict[tuple[int, int], int] = {}
    for left, right in zip(boundaries, boundaries[1:]):
        u_values = tuple(profile_value(profile, left) for profile in roots_u)
        v_values = tuple(profile_value(profile, left) for profile in roots_v)
        same = sum(u * v for u, v in zip(u_values, v_values))
        require(same == 0, ("pointwise same-root source", left, right))
        total += sum(u_values) * sum(v_values) * (right - left)
        diagonal += same * (right - left)
        support = (
            sum(value != 0 for value in u_values),
            sum(value != 0 for value in v_values),
        )
        support_types[support] = support_types.get(support, 0) + right - left
    require(total == SOURCE_TOTAL_NUMERATOR, ("source total", total))
    require(diagonal == 0, "source diagonal")
    require(Fraction(total, SRC.RPKT * SRC.DCOLL**2 * grid) == 169 * SRC.I5,
            "169 I5 source anchor")
    require(support_types == {
        (1, 12): 42548128262640,
        (2, 11): 255288769575840,
    }, ("source support types", support_types))
    return roots_u, roots_v, tuple(boundaries), profile_digest, total, support_types


def scale_profile(
    profile: tuple[tuple[int, ...], tuple[int, ...]], numerator: int, denominator: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    require(numerator % denominator == 0, ("profile scale", numerator, denominator))
    scale = numerator // denominator
    return tuple(position * scale for position in profile[0]), profile[1]


def endpoint_word_and_grid() -> tuple[tuple[int, ...], int]:
    word = EP.to_current(EP.U_FULL_REL)
    endpoint_grid = 182 * EP.lcm_tuple(word)
    require(word == (1, 183, 27, 131, 53, 313, 13, 2197, 742586),
            ("U_full word", word))
    require(endpoint_grid == 483730250419703196, "endpoint grid")
    return word, endpoint_grid


def endpoint_events(
    word: tuple[int, ...], endpoint_grid: int, alpha: int, beta: int,
    literal_tau: int | None,
) -> tuple[dict[int, int], int, int]:
    pattern = dict(EP.PATTERN_E)
    if literal_tau is None:
        require(pattern.pop(EP.GUARD) == "guard_safe", "guard deletion type")
        tau = 0
    else:
        tau = literal_tau
    shift = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    intervals = EP.fast_build_set(word, endpoint_grid, pattern, shift)
    scale = JOINT_COORDINATE // endpoint_grid
    require(scale * endpoint_grid == JOINT_COORDINATE, "endpoint coordinate scale")
    sheet_width = endpoint_grid // P
    events: dict[int, int] = {}
    mapped = 0
    for left, right in intervals:
        first = min(P - 1, left // sheet_width)
        last = min(P - 1, (right - 1) // sheet_width)
        for sheet in range(first, last + 1):
            piece_left = max(left, sheet * sheet_width)
            piece_right = min(right, (sheet + 1) * sheet_width)
            if piece_left >= piece_right:
                continue
            y_left = (P * piece_left - sheet * endpoint_grid) * scale
            y_right = (P * piece_right - sheet * endpoint_grid) * scale
            require(0 <= y_left < y_right <= JOINT_COORDINATE,
                    ("desheeted endpoint", sheet, y_left, y_right))
            bit = 1 << sheet
            events[y_left] = events.get(y_left, 0) ^ bit
            events[y_right] = events.get(y_right, 0) ^ bit
            mapped += 1
    return events, len(intervals), mapped


def chamber(left: int, right: int) -> str:
    twice_midpoint = left + right
    if Q * twice_midpoint < 2 * JOINT_COORDINATE:
        return "left"
    if Q * twice_midpoint < 12 * JOINT_COORDINATE:
        return "middle"
    return "right"


def cell_index(left: int, right: int) -> int:
    band = Q * (left + right) // JOINT_COORDINATE
    require(0 <= band < 2 * Q, ("cell band", band))
    return 0 if band in (0, 2 * Q - 1) else (band + 1) // 2


def danger_arcs() -> dict[str, frozenset[int]]:
    # Midpoints avoid every null boundary.  The literal H guard is unsafe when
    # ||(y+sheet)/13|| < 1/7.
    probes = {"left": Fraction(1, 28), "right": Fraction(27, 28)}
    derived = {}
    for name, y in probes.items():
        unsafe = set()
        for sheet in range(P):
            value = (y + sheet) / P
            fractional = value - int(value)
            distance = min(fractional, 1 - fractional)
            if distance < Fraction(1, 7):
                unsafe.add(sheet)
        derived[name] = frozenset(unsafe)
    require(derived == DANGER_EXPECTED, ("guard danger arcs", derived))
    return derived


def guard_mask(name: str, tau: int, danger: dict[str, frozenset[int]]) -> int:
    require(name in danger, ("guard chamber", name))
    mask = 0
    for sheet in range(P):
        if (sheet + tau) % P not in danger[name]:
            mask |= 1 << sheet
    return mask


def fixed_boundaries(source_boundaries: tuple[int, ...], source_grid: int) -> set[int]:
    scale = JOINT_COORDINATE // source_grid
    answer = {position * scale for position in source_boundaries}
    answer.update(index * (JOINT_COORDINATE // (2 * Q)) for index in range(2 * Q + 1))
    answer.update(index * (JOINT_COORDINATE // P) for index in range(P + 1))
    answer.update((0, JOINT_COORDINATE))
    return answer


def literal_guard_restoration(
    word: tuple[int, ...], endpoint_grid: int, alpha: int, beta: int, tau: int,
    boundaries: set[int], danger: dict[str, frozenset[int]],
) -> tuple[int, int]:
    unguarded_pattern = dict(EP.PATTERN_E)
    require(unguarded_pattern.pop(EP.GUARD) == "guard_safe",
            "reference guard deletion type")
    unguarded_shift = (0, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    direct_shift = (tau, -alpha % P, -beta % P, 0, 0, 0, 0, alpha, beta)
    fast_unguarded = EP.fast_build_set(
        word, endpoint_grid, unguarded_pattern, unguarded_shift
    )
    reference_unguarded = EP.reference_build_set(
        word, endpoint_grid, unguarded_pattern, unguarded_shift
    )
    fast_direct = EP.fast_build_set(
        word, endpoint_grid, EP.PATTERN_E, direct_shift
    )
    reference_direct = EP.reference_build_set(
        word, endpoint_grid, EP.PATTERN_E, direct_shift
    )
    require(merge_intervals(fast_unguarded) == tuple(reference_unguarded),
            ("independent unguarded interval engine", alpha, beta))
    require(merge_intervals(fast_direct) == tuple(reference_direct),
            ("independent guarded interval engine", alpha, beta, tau))

    unguarded, unguarded_count, _mapped = endpoint_events(
        word, endpoint_grid, alpha, beta, None
    )
    direct, direct_count, _direct_mapped = endpoint_events(
        word, endpoint_grid, alpha, beta, tau
    )
    positions = sorted(set(unguarded) | set(direct) | boundaries)
    unguarded_mask = direct_mask = 0
    tested = 0
    for left, right in zip(positions, positions[1:]):
        unguarded_mask ^= unguarded.get(left, 0)
        direct_mask ^= direct.get(left, 0)
        if left == right:
            continue
        name = chamber(left, right)
        expected = 0 if name == "middle" else (
            unguarded_mask & guard_mask(name, tau, danger)
        )
        require(direct_mask == expected,
                ("literal guard restoration", alpha, beta, tau, left, right,
                 direct_mask, expected))
        tested += 1
    unguarded_mask ^= unguarded.get(positions[-1], 0)
    direct_mask ^= direct.get(positions[-1], 0)
    require(unguarded_mask == direct_mask == 0, "guard masks did not close")
    return unguarded_count + direct_count, tested


class HarmonicPrimitive:
    def __init__(self, word: tuple[int, ...], endpoint_grid: int):
        zero = (0,) * 9
        q_intervals = EP.fast_build_set(word, endpoint_grid, EP.PATTERN_QA, zero)
        reference_q = EP.reference_build_set(
            word, endpoint_grid, EP.PATTERN_QA, zero
        )
        require(merge_intervals(q_intervals) == tuple(reference_q),
                "independent Q interval engine")
        require(Fraction(EP.validate_intervals(
            q_intervals, endpoint_grid, "independent U_full Q"
        ), endpoint_grid) == Fraction(197889477091847, 6201669877175682),
                "endpoint Q measure")
        scale = JOINT_COORDINATE // (P * endpoint_grid)
        require(scale * P * endpoint_grid == JOINT_COORDINATE,
                "Q(13y) pullback scale")
        pulled = []
        for branch in range(P):
            offset = branch * endpoint_grid
            pulled.extend(
                ((left + offset) * scale, (right + offset) * scale)
                for left, right in q_intervals
            )
        self.intervals = merge_intervals(pulled)
        self.starts = tuple(left for left, _right in self.intervals)
        self.ends = tuple(right for _left, right in self.intervals)
        prefix = [0]
        for left, right in self.intervals:
            contribution = (
                pow(JOINT_ROOT, CHARACTER_EXPONENT * left % JOINT_ORDER, JOINT_PRIME)
                - pow(JOINT_ROOT, CHARACTER_EXPONENT * right % JOINT_ORDER, JOINT_PRIME)
            ) % JOINT_PRIME
            prefix.append((prefix[-1] + contribution) % JOINT_PRIME)
        self.prefix = tuple(prefix)

    def value(self, point: int) -> int:
        require(0 <= point <= JOINT_COORDINATE, ("primitive point", point))
        if point == JOINT_COORDINATE:
            return self.prefix[-1]
        index = bisect_right(self.starts, point) - 1
        if index < 0:
            return 0
        if point < self.ends[index]:
            return (
                self.prefix[index]
                + pow(
                    JOINT_ROOT,
                    CHARACTER_EXPONENT * self.starts[index] % JOINT_ORDER,
                    JOINT_PRIME,
                )
                - pow(
                    JOINT_ROOT,
                    CHARACTER_EXPONENT * point % JOINT_ORDER,
                    JOINT_PRIME,
                )
            ) % JOINT_PRIME
        return self.prefix[index + 1]


def integrate_character_bank(
    word: tuple[int, ...], endpoint_grid: int,
    source_u: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
    source_v: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
    source_grid: int, source_boundaries: tuple[int, ...],
    harmonic: HarmonicPrimitive, danger: dict[str, frozenset[int]],
):
    scaled_u = tuple(scale_profile(profile, JOINT_COORDINATE, source_grid)
                     for profile in source_u)
    scaled_v = tuple(scale_profile(profile, JOINT_COORDINATE, source_grid)
                     for profile in source_v)
    base_boundaries = fixed_boundaries(source_boundaries, source_grid)
    zeta = pow(JOINT_ROOT, JOINT_ORDER // P, JOINT_PRIME)
    gamma = []
    erased = []
    diagonal = []
    totals = [0, 0, 0, 0]
    support_by_cell = [0] * Q

    for alpha in range(P):
        for beta in range(P):
            events, interval_count, mapped_count = endpoint_events(
                word, endpoint_grid, alpha, beta, None
            )
            positions = sorted(set(events) | base_boundaries)
            mask = 0
            coupled_rows = [[0] * Q for _tau in range(P)]
            erased_rows = [[0] * Q for _tau in range(P)]
            diagonal_rows = [[0] * Q for _tau in range(P)]
            active = q_active = 0
            primitive_left = harmonic.value(positions[0])
            for position_index, (left, right) in enumerate(zip(positions, positions[1:])):
                mask ^= events.get(left, 0)
                primitive_right = harmonic.value(right)
                jump = (primitive_right - primitive_left) % JOINT_PRIME
                primitive_left = primitive_right
                if left == right or not mask:
                    continue
                active += 1
                name = chamber(left, right)
                require(name != "middle",
                        ("owner support entered middle chamber", alpha, beta, left, right))
                ell = cell_index(left, right)
                require(ell == 0,
                        ("owner pullback escaped cell zero", alpha, beta, left, right, ell))
                if not jump:
                    continue
                q_active += 1
                u_values = tuple(profile_value(profile, left) for profile in scaled_u)
                v_values = tuple(profile_value(profile, left) for profile in scaled_v)
                for tau in range(P):
                    selected = mask & guard_mask(name, tau, danger)
                    left_sum = right_sum = same_sum = count = 0
                    for sheet in range(P):
                        if not ((selected >> sheet) & 1):
                            continue
                        left_sum += u_values[sheet]
                        right_sum += v_values[sheet]
                        same_sum += u_values[sheet] * v_values[sheet]
                        count += 1
                    require(same_sum == 0,
                            ("same-root current", alpha, beta, tau, left, right))
                    coupled_rows[tau][ell] = (
                        coupled_rows[tau][ell] + left_sum * right_sum * jump
                    ) % JOINT_PRIME
                    erased_rows[tau][ell] = (
                        erased_rows[tau][ell] + count * count * jump
                    ) % JOINT_PRIME
                    diagonal_rows[tau][ell] = (
                        diagonal_rows[tau][ell] + same_sum * jump
                    ) % JOINT_PRIME
            mask ^= events.get(positions[-1], 0)
            require(mask == 0, ("endpoint mask closure", alpha, beta, mask))
            phase = pow(zeta, beta, JOINT_PRIME)
            for tau in range(P):
                row = tuple(phase * value % JOINT_PRIME for value in coupled_rows[tau])
                erow = tuple(phase * value % JOINT_PRIME for value in erased_rows[tau])
                drow = tuple(phase * value % JOINT_PRIME for value in diagonal_rows[tau])
                require(all(value == 0 for value in drow), "same-root gamma row")
                require(all(value == 0 for value in row[1:]), "coupled row left cell zero")
                require(all(value == 0 for value in erow[1:]), "erased row left cell zero")
                gamma.append(row)
                erased.append(erow)
                diagonal.append(drow)
                for ell, value in enumerate(row):
                    support_by_cell[ell] += int(value != 0)
            totals[0] += interval_count
            totals[1] += mapped_count
            totals[2] += active
            totals[3] += q_active
    require(len(gamma) == len(erased) == len(diagonal) == P**3, "character count")
    require(all(all(value == 0 for value in row) for row in diagonal),
            "same-root bank")
    return tuple(gamma), tuple(erased), tuple(diagonal), tuple(totals), tuple(support_by_cell)


def inverse_table(gamma, zeta: int):
    normalizer = pow(P**3, -1, JOINT_PRIME)
    table = [[0] * P for _ell in range(Q)]
    index = 0
    for alpha in range(P):
        for _beta in range(P):
            for tau in range(P):
                row = gamma[index]
                index += 1
                for relation_t in range(P):
                    phase = pow(zeta, -(alpha + tau * relation_t) % P, JOINT_PRIME)
                    for ell in range(Q):
                        table[ell][relation_t] = (
                            table[ell][relation_t] + row[ell] * phase
                        ) % JOINT_PRIME
    return tuple(tuple(value * normalizer % JOINT_PRIME for value in row) for row in table)


def fourier_2d(matrix, eta: int, zeta: int):
    return tuple(tuple(sum(
        matrix[ell][relation_t]
        * pow(eta, -h * ell % Q, JOINT_PRIME)
        * pow(zeta, -k * relation_t % P, JOINT_PRIME)
        for ell in range(Q) for relation_t in range(P)
    ) % JOINT_PRIME for k in range(P)) for h in range(Q))


def support_shape(spectrum) -> tuple[int, int, int, int, int]:
    dc = int(spectrum[0][0] != 0)
    cell_axis = sum(spectrum[h][0] != 0 for h in range(1, Q))
    residue_axis = sum(spectrum[0][k] != 0 for k in range(1, P))
    mixed = sum(spectrum[h][k] != 0 for h in range(1, Q) for k in range(1, P))
    return dc + cell_axis + residue_axis + mixed, dc, cell_axis, residue_axis, mixed


def rank_mod(matrix) -> int:
    rows = [list(row) for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, JOINT_PRIME)
        rows[rank] = [value * inverse % JOINT_PRIME for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (value - factor * pivot_value) % JOINT_PRIME
                for value, pivot_value in zip(rows[row], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def interaction(matrix):
    inv_q = pow(Q, -1, JOINT_PRIME)
    inv_p = pow(P, -1, JOINT_PRIME)
    inv_qp = pow(Q * P, -1, JOINT_PRIME)
    row_sums = tuple(sum(row) % JOINT_PRIME for row in matrix)
    column_sums = tuple(sum(matrix[ell][relation_t] for ell in range(Q)) % JOINT_PRIME
                        for relation_t in range(P))
    grand = sum(row_sums) % JOINT_PRIME
    return tuple(tuple((
        matrix[ell][relation_t]
        - row_sums[ell] * inv_p
        - column_sums[relation_t] * inv_q
        + grand * inv_qp
    ) % JOINT_PRIME for relation_t in range(P)) for ell in range(Q))


def role_value(table, relation_t: int) -> int:
    return sum(table[ell][relation_t] for ell in range(Q)) % JOINT_PRIME


def main() -> None:
    field = split_field_certificate()
    source_u, source_v, source_boundaries, profile_digest, source_total, support_types = (
        source_profiles()
    )
    word, endpoint_grid = endpoint_word_and_grid()
    source_grid = SRC.T_DEN
    require(JOINT_COORDINATE == lcm(source_grid, P * endpoint_grid),
            "joint coordinate")
    require(JOINT_ORDER == lcm(source_grid, P * P * endpoint_grid),
            "joint order")
    require(JOINT_ORDER == P * JOINT_COORDINATE, "joint order/coordinate ratio")

    danger = danger_arcs()
    boundaries = fixed_boundaries(source_boundaries, source_grid)
    guard_records = tuple(
        (triple, *literal_guard_restoration(
            word, endpoint_grid, *triple, boundaries, danger
        ))
        for triple in CONTROL_TRIPLES
    )
    harmonic = HarmonicPrimitive(word, endpoint_grid)
    gamma, erased, diagonal, work_counts, support_by_cell = integrate_character_bank(
        word, endpoint_grid, source_u, source_v, source_grid, source_boundaries,
        harmonic, danger,
    )

    zeta = pow(JOINT_ROOT, JOINT_ORDER // P, JOINT_PRIME)
    eta = pow(JOINT_ROOT, JOINT_ORDER // Q, JOINT_PRIME)
    coupled_table = inverse_table(gamma, zeta)
    erased_table = inverse_table(erased, zeta)
    coupled_interaction = interaction(coupled_table)
    coupled_spectrum = fourier_2d(coupled_table, eta, zeta)
    erased_spectrum = fourier_2d(erased_table, eta, zeta)
    interaction_spectrum = fourier_2d(coupled_interaction, eta, zeta)

    ranks = (
        rank_mod(coupled_table), rank_mod(erased_table), rank_mod(coupled_interaction)
    )
    require(ranks == (1, 1, 1), ("coordinate ranks", ranks))
    require(all(all(value == 0 for value in row) for row in coupled_table[1:]),
            "coupled inverse escaped cell zero")
    require(all(all(value == 0 for value in row) for row in erased_table[1:]),
            "erased inverse escaped cell zero")
    residue_profile = coupled_table[0]
    residue_spectrum = tuple(sum(
        residue_profile[relation_t] * pow(zeta, -k * relation_t % P, JOINT_PRIME)
        for relation_t in range(P)
    ) % JOINT_PRIME for k in range(P))
    require(all(value != 0 for value in residue_spectrum),
            "one-dimensional residue mode vanished")
    require(all(
        coupled_spectrum[h][k] == residue_spectrum[k]
        for h in range(Q) for k in range(P)
    ), "delta-cell Fourier factorization")

    inv_q = pow(Q, -1, JOINT_PRIME)
    mean = sum(residue_profile) * pow(P, -1, JOINT_PRIME) % JOINT_PRIME
    require(all(
        coupled_interaction[ell][relation_t] == (
            ((1 if ell == 0 else 0) - inv_q)
            * (residue_profile[relation_t] - mean)
        ) % JOINT_PRIME
        for ell in range(Q) for relation_t in range(P)
    ), "rank-one ANOVA factorization")

    shapes = (
        support_shape(coupled_spectrum),
        support_shape(erased_spectrum),
        support_shape(interaction_spectrum),
    )
    require(shapes == (
        (91, 1, 6, 12, 72),
        (91, 1, 6, 12, 72),
        (72, 0, 0, 0, 72),
    ), ("spectral shapes", shapes))

    role_values = (role_value(coupled_table, 1), role_value(coupled_table, 0))
    bridge = (role_values[0] - role_values[1]) % JOINT_PRIME
    erased_values = (role_value(erased_table, 1), role_value(erased_table, 0))
    erased_bridge = (erased_values[0] - erased_values[1]) % JOINT_PRIME
    require((*role_values, bridge) == EXPECTED_ROLE_VALUES,
            ("candidate role comparison", role_values, bridge))
    require((*erased_values, erased_bridge) == EXPECTED_ERASURE_VALUES,
            ("candidate erasure comparison", erased_values, erased_bridge))

    relation_t = 6
    fixed_profile = tuple(coupled_table[ell][relation_t] for ell in range(Q))
    erased_fixed_profile = tuple(erased_table[ell][relation_t] for ell in range(Q))
    require(fixed_profile == (EXPECTED_FIXED,) + (0,) * 6,
            ("fixed relation profile", fixed_profile))
    require(erased_fixed_profile == (EXPECTED_ERASURE_FIXED,) + (0,) * 6,
            ("erased fixed relation", erased_fixed_profile))
    fixed_f7 = tuple(sum(
        fixed_profile[ell] * pow(eta, -h * ell % Q, JOINT_PRIME)
        for ell in range(Q)
    ) % JOINT_PRIME for h in range(Q))
    require(fixed_f7 == (EXPECTED_FIXED,) * Q, "fixed F7 delta spectrum")

    comparison = {
        "gamma": digest_json(gamma),
        "erased_gamma": digest_json(erased),
        "coupled_table": digest_json(coupled_table),
        "erased_table": digest_json(erased_table),
        "coupled_spectrum": digest_json(coupled_spectrum),
        "interaction_spectrum": digest_json(interaction_spectrum),
    }
    support_record = tuple(sorted(support_types.items()))
    semantic = digest_json((
        field, JOINT_COORDINATE, JOINT_ORDER, profile_digest, source_total,
        support_record, guard_records, work_counts, support_by_cell, ranks, shapes,
        role_values, bridge, erased_values, erased_bridge, fixed_profile,
        erased_fixed_profile, fixed_f7, comparison,
    ))
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic digest", semantic))

    print("== independent hostile audit: r=5/U_full owner-node common base ==")
    print(f"dependencies=((THM2594_source,{SOURCE_SHA256}),(THM3514_audit,{ENDPOINT_AUDIT_SHA256}),(THM3514_parent,{ENDPOINT_PARENT_SHA256}))")
    print(f"common_scales=(source={source_grid},endpoint={endpoint_grid},coordinate={JOINT_COORDINATE},order={JOINT_ORDER})")
    print(f"split_field=(prime={field[0]},generator={field[1]},root={field[2]},lucas=PASS,exact_order=PASS)")
    print(f"frequency_descent=(-13*w_u+742599*w_q)=57122*y+integer;character_exponent_on_C={CHARACTER_EXPONENT};phase=zeta13^beta")
    print(f"source_rebuild=(profile_sha256={profile_digest},service_numerator={source_total},169I5=PASS,same_root_pointwise_zero=PASS,support_types={support_types})")
    print(f"literal_guard_restoration={guard_records}: PASS")
    print("primary_hostile=U_full_owner_in_at_w_u_is_||13w_u||<1/14=||y||<1/14=cell_0_a.e.: PASS")
    print(f"cell_support=(character_rows={support_by_cell},coordinate_ranks={ranks})")
    print(f"work_counts_(endpoint_intervals,mapped_intervals,active_segments,q_active_segments)={work_counts}")
    print(f"inverse_roles=(q_H={role_values[0]},q_q5={role_values[1]},bridge={bridge})")
    print(f"source_erasure=(q_H={erased_values[0]},q_q5={erased_values[1]},bridge={erased_bridge})")
    print(f"spectral_shapes_(total,dc,F7axis,F13axis,mixed)=(coupled={shapes[0]},erased={shapes[1]},ANOVA={shapes[2]})")
    print("separable_mechanism=table=delta_0(cell)*R(residue);ANOVA=(delta_0-1/7)*(R-mean(R));full_support_is_not_two_coordinate_mixing")
    print(f"fixed_(1,0,6)=(profile={fixed_profile},F7={fixed_f7})")
    print(f"comparison_sha256={comparison}")
    print(f"semantic_sha256={semantic}")
    print("verdict=ACCEPT narrow one-base owner-node integrand and nonzero 13-residue role bridge; REJECT genuine 7x13 spectral closure/current interpretation")
    print("typing=all root sums are inside one y-integral and no Fubini mismatch occurs; common root gauge and linked-node observable are imposed candidates, not a proved source/arrival or U_clock transport")
    print("scope=no exact C(a;X,m),no physical current,no arrival/source-time identification,no U_clock chronology,no row exclusion,no LRC(14)")
    print("all exact hostile checks passed")


if __name__ == "__main__":
    main()
