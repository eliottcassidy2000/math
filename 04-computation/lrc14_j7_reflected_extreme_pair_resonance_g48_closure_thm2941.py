#!/usr/bin/env python3
r"""Close the reflected extreme-pair rays 4, 5, 6 and the g>=48 tail.

This verifier works inside THM-2941's sufficient reflected ``k=1`` family,
after the arbitrary-level bank has left 561 bodies.  The two bodies whose
same-level graph is not K6 are disjoint from this residual, so an uncertified
packet on a residual body has six distinct levels.

Write the physical extreme levels as

    m=gP,  Q=gR,  gcd(P,R)=1,  R/P>3.

The low-phase condition is ``P+R<=7``.  In this ratio range its only primitive
channels are (1,4), (1,5), and (1,6).  For one of these rays, let ``a`` be the
label at level g and ``b`` the label at level Rg.  On a body-safe cell
``[j,j+1]``, the primitive fibre has transverse coordinate

    z(s)=(b-Ra)(j+s)/L,  0<=s<=1.

For R<=6 this fibre is a single periodized trapezoid.  If

    M=max(|h*j-kL|, |h*(j+1)-kL|),  h=|b-Ra|,

then its cellwise floor is

    1/(7R)                         if M <= (R-1)L/14,
    ((R+1)L/14-M)/(R L)           if (R-1)L/14 < M < (R+1)L/14,
    0                              otherwise.

The exact optimizer checks safe-range endpoints and the integers adjacent to
``kL/h-1/2``.  Across all 561*30 physical orientations its floors are at
least 3/112, 2/105, and 1/42 for R=4,5,6.

For actual scale g, subdividing by ``u=(r+x)/g`` leaves the primitive fibre
after dropping the two ``x/(gL)`` perturbations.  Their two-clause symmetric
difference is at most ``4(a+b)/(gL)``.  Exact distinct-level debt then proves
analytic tails starting at g=17,22,19.  The finite prefixes g=2..16,
g=2..21, and g=2..18 are checked directly at the same selected cells.  The
direct engine is an independent integer two-pointer formula for the full
teeth on a body-safe cell; it is compared with the promoted Fraction interval
engine on all 50,490 scale-two controls.

Every other reduced channel with R/P>3 has P+R>=8, so the inherited primitive
fibre floor is 1/105.  The same affine transport and a termwise worst debt
comparison prove all such channels when g>=48.  The surrogate invoice is
negative in exactly one row at g=47, making 48 sharp for this uniform invoice,
not for the physical pair.

Together with the Q>=23m cone, the remaining reflected certificate locus is

    D>=6, m>=2, 3<Q/m<23, Q/m not in {4,5,6}, gcd(m,Q)<=47.

This is a sufficient-family reduction, not a physical-survivor census or a
proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
Q23_SOURCE = ROOT / "04-computation/lrc14_j7_reflected_extreme_pair_q23m_cone_closure_thm2941.py"
Q23_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_q23m_cone_closure_thm2941.out"
LOW_SOURCE = ROOT / "04-computation/lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
LOW_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_resonance_g48_closure_thm2941.out"

EXPECTED_Q23_SOURCE_SHA256 = "135e94154c400cb612c3f93d03557b121ddafe51f5f6cf064c057848acb72a1d"
EXPECTED_Q23_OUTPUT_SHA256 = "e97426ea5016b23beeee3e00ee24439bee37caddbc0528fc3ec53a112fc93131"
EXPECTED_Q23_SEMANTIC = "88644e8527610dcfd9c1108889a533ceee5883aeb56a84ed0da5edb6891c7e5e"
EXPECTED_LOW_SOURCE_SHA256 = "416c36f16f7c821feb8d260882711d2717069147b8604a93ba60432785cf1d1c"
EXPECTED_LOW_OUTPUT_SHA256 = "1a2863c10f32c0580ab7a5ac866df29f0ba38def1b9b713551cc27ae76676601"
EXPECTED_LOW_SEMANTIC = "0063d2beb7c627f1836576081689720d41f45a078b4b7b3e06aeb81574039248"

BODY_COUNT = 561
ORIENTATION_COUNT = 30
RAYS = (4, 5, 6)
TAIL_START = {4: 17, 5: 22, 6: 19}
FINITE_ROW_COUNT = BODY_COUNT * ORIENTATION_COUNT * sum(TAIL_START[r] - 2 for r in RAYS)
PROFILE_COUNT = BODY_COUNT * ORIENTATION_COUNT * len(RAYS)
ENGINE_CONTROL_COUNT = PROFILE_COUNT
HIGH_FLOOR = F(1, 105)
GENERIC_TAIL_START = 48
H = (1, 2, 3, 4, 6, 12)

EXPECTED_PROFILE_MINIMA = {
    4: (F(3, 112), H, (5, 4), 155, 39, 42, -42, 168),
    5: (F(2, 105), H, (5, 3), 155, 52, 56, -56, 168),
    6: (F(1, 42), (1, 2, 3, 4, 5, 6), (0, 1), 779, 4, 244, -4, 840),
}
EXPECTED_TAIL_MINIMA = {
    4: (
        F(1830813628341083957, 4081224172991372132880), H, (5, 4),
        17, 155, F(3, 112),
    ),
    5: (
        F(949961811207548632, 1178577790606909135905), H, (5, 3),
        22, 155, F(2, 105),
    ),
    6: (
        F(249759800474998027, 949285825547992885050), H, (5, 4),
        19, 155, F(1, 42),
    ),
}
EXPECTED_TAIL_HOSTILES = {
    4: (-F(12848328910247272, 10754449400337062805), H, (5, 4), 16, 155),
    5: (-F(2045008180052, 33444939832927623), H, (5, 3), 21, 155),
    6: (-F(32185011642260912, 30857807938406456805), H, (5, 4), 18, 155),
}
EXPECTED_FINITE_WEAKEST = (
    F(8011730170, 531042645861), H, (5, 0), 6, 2, 151,
    F(48, 2015), F(23192086126, 2655213229305),
)
EXPECTED_GENERIC_G48 = (
    F(2726527595145600839, 13843540499574717938160), H, (5, 4),
    (52, 51, 50, 49, 336, 48),
)
EXPECTED_GENERIC_G47 = (
    -F(48791489743358, 36036573362694791523), H, (5, 4),
    (51, 50, 49, 48, 329, 47),
)

# Frozen after the first complete exact replay.
EXPECTED_PROFILE_DIGEST = "0ef8854e16019b4b20588f4528a91af922af859f378b58d30528124d12c992db"
EXPECTED_FINITE_DIGEST = "eb4bf1b44d21794157ab384547eafe39402fb120c7fb954859d2162700e74450"
EXPECTED_ENGINE_DIGEST = "c0f387eefa34abe7f219f5e3c1ab090e786baa5627f6e8c0828dc25c98650a95"
EXPECTED_TAIL_DIGEST = "02250a5bc50c59335367f831144833422a5fb377e964e39e2de60691da7a0572"
EXPECTED_GENERIC_DIGEST = "c2dcd8c542e6d38fc4b9032398807c0f75999008c32a1a98861c4d1f60125b1f"
EXPECTED_SEMANTIC = "1fb9c1ee2f0c02e2033b940001b58a76232478a2ae94166e8d1f8cb3d54ae6ac"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


for path, expected in (
    (Q23_SOURCE, EXPECTED_Q23_SOURCE_SHA256),
    (Q23_OUTPUT, EXPECTED_Q23_OUTPUT_SHA256),
    (LOW_SOURCE, EXPECTED_LOW_SOURCE_SHA256),
    (LOW_OUTPUT, EXPECTED_LOW_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream artifact changed", path, sha256(path), expected))
require(f"semantic_sha256={EXPECTED_Q23_SEMANTIC}" in Q23_OUTPUT.read_text(),
        "Q23 semantic token missing")
require(f"semantic_sha256={EXPECTED_LOW_SEMANTIC}" in LOW_OUTPUT.read_text(),
        "low-phase semantic token missing")

Q23 = import_module("extreme_pair_q23_base", Q23_SOURCE)
LOW = import_module("extreme_pair_low_phase_base", LOW_SOURCE)
T = Q23.T
B = Q23.B
require(LOW.FIBER_FLOOR == HIGH_FLOOR, (LOW.FIBER_FLOOR, HIGH_FLOOR))


def selected_debt(body: tuple[int, ...], pair: tuple[int, int],
                  minimum: int, maximum: int):
    """Exact worst six-distinct-level debt for fixed physical extremes."""
    require(maximum >= minimum + 5, (minimum, maximum))
    ruler = 14 * lcm(*body)
    i, j = pair
    remaining = tuple(k for k in range(6) if k not in {i, j})
    descending = tuple(sorted(remaining, key=lambda k: body[k], reverse=True))
    levels = [None] * 6
    levels[i] = minimum
    levels[j] = maximum
    for offset, slot in enumerate(descending, 1):
        levels[slot] = minimum + offset
    levels = tuple(levels)
    debt = sum(
        (F(label, 7 * (level * ruler - label))
         for label, level in zip(body, levels)),
        F(0),
    )
    return debt, levels, ruler


def best_primitive_profile(safe_ranges: tuple[tuple[int, int], ...],
                           ruler: int, cross_difference: int, ratio: int):
    """Maximize the exact F_(1,R) floor on one integer body-safe cell."""
    h = abs(cross_difference)
    if h == 0:
        return F(1, 7 * ratio), safe_ranges[0][0], 0, 0
    support = (ratio + 1) * ruler // 14
    plateau = (ratio - 1) * ruler // 14
    require(14 * support == (ratio + 1) * ruler, (ratio, ruler, support))
    require(14 * plateau == (ratio - 1) * ruler, (ratio, ruler, plateau))
    best = None
    for left, right in safe_ranges:
        for lift in range(h * left // ruler - 2, h * right // ruler + 3):
            midpoint_floor = (2 * lift * ruler - h) // (2 * h)
            for cell in (left, right - 1, midpoint_floor - 1,
                         midpoint_floor, midpoint_floor + 1):
                if not left <= cell < right:
                    continue
                endpoint_radius = max(
                    abs(h * cell - lift * ruler),
                    abs(h * (cell + 1) - lift * ruler),
                )
                if endpoint_radius <= plateau:
                    floor = F(1, 7 * ratio)
                elif endpoint_radius < support:
                    floor = F(support - endpoint_radius, ratio * ruler)
                else:
                    floor = F(0)
                row = (floor, cell, lift, endpoint_radius)
                if best is None or row > best:
                    best = row
    require(best is not None, (safe_ranges, ruler, cross_difference, ratio))
    return best


def full_tooth_overlap(ruler: int, a: int, p: int,
                       b: int, q: int, cell: int) -> F:
    """Independent integer two-pointer overlap of two untruncated combs."""
    require(p >= 1 and q >= 1, (p, q))
    r1 = (a * cell) % ruler
    r2 = (b * cell) % ruler
    require(ruler // 14 <= r1 <= 13 * ruler // 14 - a,
            (ruler, a, p, cell, r1))
    require(ruler // 14 <= r2 <= 13 * ruler // 14 - b,
            (ruler, b, q, cell, r2))
    z1 = p * ruler - a
    z2 = q * ruler - b
    d1 = 14 * z1
    d2 = 14 * z2
    left_index = 0
    right_index = 0
    common_numerator = 0
    while left_index < p and right_index < q:
        left1 = 14 * (r1 + left_index * ruler) - ruler
        right1 = left1 + 2 * ruler
        left2 = 14 * (r2 + right_index * ruler) - ruler
        right2 = left2 + 2 * ruler
        common_numerator += max(
            0,
            min(right1 * d2, right2 * d1)
            - max(left1 * d2, left2 * d1),
        )
        if right1 * d2 < right2 * d1:
            left_index += 1
        else:
            right_index += 1
    return F(common_numerator, d1 * d2)


def ray_transport_margin(body: tuple[int, ...], pair: tuple[int, int],
                         ratio: int, scale: int, floor: F):
    debt, levels, ruler = selected_debt(body, pair, scale, ratio * scale)
    a, b = body[pair[0]], body[pair[1]]
    margin = floor - F(4 * (a + b), scale * ruler) - debt
    return margin, debt, levels, ruler


def ray_margin_increment(body: tuple[int, ...], pair: tuple[int, int],
                         ratio: int, scale: int, floor: F) -> F:
    """Positive exact increase of the affine-tail invoice from g to g+1."""
    now = ray_transport_margin(body, pair, ratio, scale, floor)[0]
    later = ray_transport_margin(body, pair, ratio, scale + 1, floor)[0]
    ruler = 14 * lcm(*body)
    levels = ray_transport_margin(body, pair, ratio, scale, floor)[2]
    next_levels = ray_transport_margin(body, pair, ratio, scale + 1, floor)[2]
    transport = F(4 * (body[pair[0]] + body[pair[1]]),
                  ruler * scale * (scale + 1))
    debt_drop = F(0)
    for label, level, next_level in zip(body, levels, next_levels):
        step = next_level - level
        require(step in {1, ratio}, (body, pair, ratio, scale, levels, next_levels))
        debt_drop += F(
            label * step * ruler,
            7 * (level * ruler - label) * (next_level * ruler - label),
        )
    closed = transport + debt_drop
    require(later - now == closed > 0,
            (body, pair, ratio, scale, floor, now, later, closed))
    return closed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = T.residual_bodies()
    require(len(bodies) == BODY_COUNT, len(bodies))
    repeated_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (repeated_exceptions & set(bodies)),
            ("same-level exceptions in residual", repeated_exceptions & set(bodies)))

    low_channels = tuple(
        (p, q)
        for p in range(1, 8)
        for q in range(p + 1, 8)
        if gcd(p, q) == 1 and q > 3 * p and p + q <= 7
    )
    require(low_channels == ((1, 4), (1, 5), (1, 6)), low_channels)

    profiles_by_ray = {}
    profile_digest = hashlib.sha256()
    profile_minima = {}
    for ratio in RAYS:
        profiles = []
        for body in bodies:
            ruler, safe_ranges = T.R.safe_cell_ranges(body)
            for i in range(6):
                for j in range(6):
                    if i == j:
                        continue
                    cross_difference = body[j] - ratio * body[i]
                    floor, cell, lift, endpoint_radius = best_primitive_profile(
                        safe_ranges, ruler, cross_difference, ratio
                    )
                    row = (
                        floor, body, (i, j), cell, lift, endpoint_radius,
                        cross_difference, ruler,
                    )
                    profiles.append(row)
                    profile_digest.update(f"{ratio}|{row}\n".encode())
        require(len(profiles) == BODY_COUNT * ORIENTATION_COUNT,
                (ratio, len(profiles)))
        profiles_by_ray[ratio] = tuple(profiles)
        profile_minima[ratio] = min(profiles)
        require(profile_minima[ratio] == EXPECTED_PROFILE_MINIMA[ratio],
                (ratio, profile_minima[ratio], EXPECTED_PROFILE_MINIMA[ratio]))
    require(sum(map(len, profiles_by_ray.values())) == PROFILE_COUNT, PROFILE_COUNT)

    finite_digest = hashlib.sha256()
    finite_count = 0
    finite_weakest = None
    finite_scale_minima = []
    for ratio in RAYS:
        for scale in range(2, TAIL_START[ratio]):
            scale_weakest = None
            for profile in profiles_by_ray[ratio]:
                floor, body, pair, cell, _, _, _, ruler = profile
                i, j = pair
                overlap = full_tooth_overlap(
                    ruler, body[i], scale, body[j], ratio * scale, cell
                )
                debt, levels, debt_ruler = selected_debt(
                    body, pair, scale, ratio * scale
                )
                require(debt_ruler == ruler, (body, pair, ruler, debt_ruler))
                margin = overlap - debt
                row = (
                    margin, body, pair, ratio, scale, cell,
                    overlap, debt, levels, floor,
                )
                require(margin > 0, ("finite resonance failure", row))
                if finite_weakest is None or row < finite_weakest:
                    finite_weakest = row
                if scale_weakest is None or row < scale_weakest:
                    scale_weakest = row
                finite_digest.update(f"{row}\n".encode())
                finite_count += 1
            require(scale_weakest is not None, (ratio, scale))
            finite_scale_minima.append(scale_weakest)
    require(finite_count == FINITE_ROW_COUNT, finite_count)
    require(finite_weakest[:8] == EXPECTED_FINITE_WEAKEST,
            (finite_weakest, EXPECTED_FINITE_WEAKEST))

    engine_digest = hashlib.sha256()
    engine_count = 0
    for ratio in RAYS:
        for profile in profiles_by_ray[ratio]:
            _, body, pair, cell, _, _, _, ruler = profile
            i, j = pair
            independent = full_tooth_overlap(
                ruler, body[i], 2, body[j], 2 * ratio, cell
            )
            imported = B.intersection_mass(
                T.R.reflected_level_arcs(ruler, body[i], 2, cell),
                T.R.reflected_level_arcs(ruler, body[j], 2 * ratio, cell),
            )
            require(independent == imported,
                    ("engine disagreement", ratio, body, pair, cell,
                     independent, imported))
            row = (ratio, body, pair, cell, independent)
            engine_digest.update(f"{row}\n".encode())
            engine_count += 1
    require(engine_count == ENGINE_CONTROL_COUNT, engine_count)

    tail_digest = hashlib.sha256()
    tail_minima = {}
    tail_hostiles = {}
    for ratio in RAYS:
        start = TAIL_START[ratio]
        rows = []
        prior_rows = []
        for profile in profiles_by_ray[ratio]:
            floor, body, pair, cell, _, _, _, _ = profile
            margin = ray_transport_margin(body, pair, ratio, start, floor)[0]
            prior = ray_transport_margin(body, pair, ratio, start - 1, floor)[0]
            row = (margin, body, pair, start, cell, floor)
            prior_row = (prior, body, pair, start - 1, cell)
            require(margin > 0, ("tail start failure", ratio, row))
            ray_margin_increment(body, pair, ratio, start, floor)
            rows.append(row)
            prior_rows.append(prior_row)
            tail_digest.update(f"{ratio}|{row}|{prior_row}\n".encode())
        tail_minima[ratio] = min(rows)
        tail_hostiles[ratio] = min(prior_rows)
        require(tail_minima[ratio] == EXPECTED_TAIL_MINIMA[ratio],
                (ratio, tail_minima[ratio], EXPECTED_TAIL_MINIMA[ratio]))
        require(tail_hostiles[ratio] == EXPECTED_TAIL_HOSTILES[ratio]
                and tail_hostiles[ratio][0] < 0,
                (ratio, tail_hostiles[ratio], EXPECTED_TAIL_HOSTILES[ratio]))

    generic_digest = hashlib.sha256()
    generic_rows = {}
    generic_nonpositive = {}
    for scale in (GENERIC_TAIL_START - 1, GENERIC_TAIL_START):
        rows = []
        for body in bodies:
            for i in range(6):
                for j in range(6):
                    if i == j:
                        continue
                    pair = (i, j)
                    margin, _, levels, _ = ray_transport_margin(
                        body, pair, 7, scale, HIGH_FLOOR
                    )
                    row = (margin, body, pair, levels)
                    rows.append(row)
                    generic_digest.update(f"{scale}|{row}\n".encode())
                    if scale == GENERIC_TAIL_START:
                        require(margin > 0, ("generic g48 failure", row))
                        ray_margin_increment(body, pair, 7, scale, HIGH_FLOOR)
        generic_rows[scale] = min(rows)
        generic_nonpositive[scale] = sum(row[0] <= 0 for row in rows)
    require(generic_rows[48] == EXPECTED_GENERIC_G48,
            (generic_rows[48], EXPECTED_GENERIC_G48))
    require(generic_rows[47] == EXPECTED_GENERIC_G47,
            (generic_rows[47], EXPECTED_GENERIC_G47))
    require(generic_nonpositive == {47: 1, 48: 0}, generic_nonpositive)

    digests = {
        "profile": profile_digest.hexdigest(),
        "finite": finite_digest.hexdigest(),
        "engine": engine_digest.hexdigest(),
        "tail": tail_digest.hexdigest(),
        "generic": generic_digest.hexdigest(),
    }
    for actual, expected, label in (
        (digests["profile"], EXPECTED_PROFILE_DIGEST, "profile"),
        (digests["finite"], EXPECTED_FINITE_DIGEST, "finite"),
        (digests["engine"], EXPECTED_ENGINE_DIGEST, "engine"),
        (digests["tail"], EXPECTED_TAIL_DIGEST, "tail"),
        (digests["generic"], EXPECTED_GENERIC_DIGEST, "generic"),
    ):
        if expected is not None:
            require(actual == expected, (label, actual, expected))

    semantic_payload = (
        tuple(bodies), tuple(sorted(repeated_exceptions)), low_channels,
        tuple((ratio, profile_minima[ratio]) for ratio in RAYS),
        PROFILE_COUNT, FINITE_ROW_COUNT, finite_weakest,
        tuple(finite_scale_minima), ENGINE_CONTROL_COUNT,
        tuple((ratio, tail_minima[ratio], tail_hostiles[ratio]) for ratio in RAYS),
        GENERIC_TAIL_START, generic_rows, generic_nonpositive, digests,
        "D>=6;m>=2;3<Q/m<23;Q/m not in {4,5,6};gcd(m,Q)<=47",
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, (semantic, EXPECTED_SEMANTIC))

    lines = [
        "LRC14 reflected extreme-pair resonance rays and g48 tail",
        f"universe=residual_bodies:{len(bodies)};physical_orientations:{ORIENTATION_COUNT};profiles:{PROFILE_COUNT};same_level_exceptions_disjoint:{tuple(sorted(repeated_exceptions))}",
        f"ratio_classification=for coprime P<Q with Q>3P,the condition P+Q<=7 gives exactly {low_channels};all other channels have inherited primitive floor 1/105",
        "primitive_optimizer=on each safe integer range,minimize the convex endpoint radius at range endpoints and integers adjacent to kL/h-1/2;h=P*b-Q*a;h=0 uses first safe cell",
    ]
    for ratio in RAYS:
        lines.append(
            f"R{ratio}_profile=minimum:{qtext(profile_minima[ratio][0])};row:{profile_minima[ratio][1:]}"
        )
    lines.extend((
        f"finite_prefixes=R4:g2..16;R5:g2..21;R6:g2..18;rows:{finite_count};weakest_margin:{qtext(finite_weakest[0])};body:{finite_weakest[1]};pair:{finite_weakest[2]};R:{finite_weakest[3]};g:{finite_weakest[4]};cell:{finite_weakest[5]};overlap:{qtext(finite_weakest[6])};debt:{qtext(finite_weakest[7])}",
        f"independent_engine_controls=g2 on every profile;rows:{engine_count};integer_full-tooth_two-pointer_equals_promoted_Fraction_engine",
        "affine_transport=under u=(r+x)/g,the primitive fibre loses at most 4(a+b)/(gL);exact distinct debt and transport both decrease strictly with g",
    ))
    for ratio in RAYS:
        lines.append(
            f"R{ratio}_tail=start_g:{TAIL_START[ratio]};weakest:{qtext(tail_minima[ratio][0])}@body={tail_minima[ratio][1]},pair={tail_minima[ratio][2]},cell={tail_minima[ratio][4]},floor={qtext(tail_minima[ratio][5])};prior_hostile_invoice:{qtext(tail_hostiles[ratio][0])}"
        )
    lines.extend((
        f"generic_high_phase_tail=all nonresonant Q/P>3 channels with gcd scale g>={GENERIC_TAIL_START};g48_weakest:{qtext(generic_rows[48][0])}@body={generic_rows[48][1]},pair={generic_rows[48][2]};g47_hostile:{qtext(generic_rows[47][0])};g47_nonpositive_rows:{generic_nonpositive[47]}",
        "assembly=the physical extreme pair closes the complete rays Q/m=4,5,6;every other reduced channel with Q/m>3 closes when gcd(m,Q)>=48",
        "residual=inside the inherited D>=6 stage:561 bodies;m>=2;3<Q/m<23;Q/m not in {4,5,6};gcd(m,Q)<=47",
        "scope=reflected THM-2941 sufficient family only;not a physical-survivor census and not LRC14",
        f"profile_digest={digests['profile']}",
        f"finite_digest={digests['finite']}",
        f"engine_digest={digests['engine']}",
        f"tail_digest={digests['tail']}",
        f"generic_digest={digests['generic']}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"q23_source_sha256={EXPECTED_Q23_SOURCE_SHA256}",
        f"q23_output_sha256={EXPECTED_Q23_OUTPUT_SHA256}",
        f"q23_semantic_sha256={EXPECTED_Q23_SEMANTIC}",
        f"low_source_sha256={EXPECTED_LOW_SOURCE_SHA256}",
        f"low_output_sha256={EXPECTED_LOW_OUTPUT_SHA256}",
        f"low_semantic_sha256={EXPECTED_LOW_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
