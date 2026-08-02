#!/usr/bin/env python3
r"""Close every THM-2941 reflected midratio channel of gcd scale at least six.

This verifier is scoped to the 561-body residual of THM-2941.  Write the
physical extreme levels as

    m=gP,  Q=gR,  gcd(P,R)=1,  3P<R<6P.

The already closed rays 4 and 5 imply P>=2 and hence R>=7.  Two improvements
replace the inherited generic ``g>=48`` invoice.

First, the exact primitive fibre

    F_(P,R)(z)=[T_((P+R)/14)(z)-T_((R-P)/14)(z)]/(PR)

has minimum at least 1/77 in this ratio cone, with equality only for (2,11)
and (3,11).  The Fourier bound handles PR>=68; eleven channels with PR<68
are checked on the complete fourteen-grid breakpoint set.

Second, the reflected perturbation is one-sided.  For integer N and
epsilon=c/(gL), compare

    ||(N-epsilon)x-theta||<1/14  and  ||Nx-theta||<1/14.

Every boundary moves monotonically right.  Pairing boundaries and summing
their displacements gives the sharp uniform ledger

    mu(A_epsilon triangle A_0) <= c(N+1)/(gLN-c).

For the two clauses the total loss is therefore

    a(P+1)/(gLP-a) + b(R+1)/(gLR-b),

instead of 4(a+b)/(gL).

At g=6, P>=13 reduces monotonically to P=13,R=40 and is positive on all
561*30 body/orientation rows using the Fourier floor.  The complete remaining
bank P=2..12 contains 135 primitive channels and 2,272,050 analytic rows.
Exactly 150 analytic invoices are nonpositive.  An independent integer
two-pointer full-tooth engine searches every body-safe cell on those rows;
all 150 have a positive located physical margin and agree cellwise with the
promoted Fraction interval engine.  Every loss and debt term decreases with
g, so all g>=6 follow.

This closes only a sufficient-certificate branch.  It reduces the current
THM-2941 reflected residual from gcd(m,Q)<=47 to gcd(m,Q)<=5; it is not a
physical-survivor census and does not prove LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
UPSTREAM_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_extreme_pair_resonance_g48_closure_thm2941.py"
)
UPSTREAM_OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_extreme_pair_resonance_g48_closure_thm2941.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_reflected_midratio_primitive_floor_g6_closure_thm3210.out"
)

EXPECTED_UPSTREAM_SOURCE_SHA256 = (
    "5ed1a4235466a801022b67342f6ded07cdea8e1449889ac3cc4b68ed5eb1d87f"
)
EXPECTED_UPSTREAM_OUTPUT_SHA256 = (
    "817047f36413e0ee707bcdb92dd20475b554da740651f387da86587ec281663c"
)
EXPECTED_UPSTREAM_SEMANTIC = (
    "1fb9c1ee2f0c02e2033b940001b58a76232478a2ae94166e8d1f8cb3d54ae6ac"
)

BODY_COUNT = 561
ORIENTATION_COUNT = 30
BASE_SCALE = 6
LARGE_P_START = 13
LARGE_R_BASE = 40
FINITE_CHANNEL_COUNT = 135
FINITE_ROW_COUNT = 2_272_050
ANALYTIC_FAILURE_COUNT = 150
RESTRICTED_FLOOR = F(1, 77)
FOURIER_PRODUCT_START = 68
H = (1, 2, 3, 4, 6, 12)

EXPECTED_CHANNELS_BY_P = (
    (2, 3), (3, 6), (4, 6), (5, 12), (6, 6), (7, 18),
    (8, 12), (9, 18), (10, 12), (11, 30), (12, 12),
)
EXPECTED_SMALL_FLOOR_ROWS = (
    (F(1, 77), 2, 11, (0, 1, 13), 22),
    (F(1, 77), 3, 11, (0,), 33),
    (F(1, 70), 3, 10, (0, 1, 13), 30),
    (F(1, 63), 2, 9, (0, 1, 2, 3, 11, 12, 13), 18),
    (F(2, 119), 3, 17, (6, 7, 8), 51),
    (F(1, 56), 3, 16, (5, 6, 7, 8, 9), 48),
    (F(5, 273), 3, 13, (0, 1, 2, 12, 13), 39),
    (F(2, 105), 4, 15, (5, 6, 7, 8, 9), 60),
    (F(1, 52), 4, 13, (0, 1, 2, 3, 11, 12, 13), 52),
    (F(1, 49), 2, 7, tuple(range(14)), 14),
    (F(1, 49), 3, 14, tuple(range(14)), 42),
)
EXPECTED_FINITE_WEAKEST = (
    -F(446346556060136476, 34286705143599831615),
    H,
    (5, 4),
    2,
    11,
    F(1, 77),
    F(7545, 308449),
    F(4812556518843536, 3116973194872711965),
    (16, 15, 14, 13, 66, 12),
    (0, 1, 13),
)
EXPECTED_LARGE_WEAKEST = (
    F(26392028086719269675209, 104017741322142148178386800),
    H,
    (5, 4),
    F(991, 50960),
    F(138797, 7330429),
    F(18470013214620220016, 71440756402570156715925),
    (82, 81, 80, 79, 240, 78),
)
EXPECTED_DIRECT_WEAKEST = (
    F(12808546596795938, 662437112164402113),
    (1, 2, 3, 4, 6, 12),
    (3, 5),
    2,
    7,
    152,
    F(6084, 295261),
    F(841300116478834, 662437112164402113),
    (16, 15, 14, 12, 13, 42),
)
EXPECTED_FAILURE_BODY_HISTOGRAM = (
    ((1, 2, 3, 4, 6, 12), 136),
    ((1, 2, 3, 4, 8, 12), 2),
    ((1, 2, 3, 6, 8, 12), 3),
    ((1, 2, 4, 6, 8, 12), 3),
    ((1, 3, 4, 6, 8, 12), 3),
    ((2, 3, 4, 6, 8, 12), 3),
)

# Frozen after the first complete exact replay.
EXPECTED_CHANNEL_MINIMA_DIGEST = (
    "355d703d1e6ae7bb854f509d81ddffd53604474318151ef2c373915aeb3ed7bd"
)
EXPECTED_FAILURE_DIGEST = (
    "8290d39d340c24bba6f6f883933ce7a90caf73f3c65bf972731182170140f1cb"
)
EXPECTED_DIRECT_DIGEST = (
    "90aadbacf98fb857657319f6721c273084d088ae3a4fa417d4e228ba9f884a8a"
)
EXPECTED_BOUNDARY_AUDIT_DIGEST = (
    "d2848403cd1414fb2e75a88bb284dbd0f2a7fe0ee58fc82d106fad1e262f286b"
)
EXPECTED_SEMANTIC = (
    "fb4759685fbd559d5526777b230ef94388eb0b72b337debad6e86915b47aeda5"
)


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


require(
    sha256(UPSTREAM_SOURCE) == EXPECTED_UPSTREAM_SOURCE_SHA256,
    ("upstream source changed", sha256(UPSTREAM_SOURCE)),
)
require(
    sha256(UPSTREAM_OUTPUT) == EXPECTED_UPSTREAM_OUTPUT_SHA256,
    ("upstream output changed", sha256(UPSTREAM_OUTPUT)),
)
require(
    f"semantic_sha256={EXPECTED_UPSTREAM_SEMANTIC}" in UPSTREAM_OUTPUT.read_text(),
    "upstream semantic token missing",
)

U = import_module("thm3210_resonance_base", UPSTREAM_SOURCE)
T = U.T
B = U.B


def periodic_tent(radius: F, z: F) -> F:
    bound = radius.numerator // radius.denominator + 3
    return sum(
        (max(F(0), radius - abs(z + shift)) for shift in range(-bound, bound + 1)),
        F(0),
    )


def fiber_density(p: int, r: int, z: F) -> F:
    require(1 <= p < r and gcd(p, r) == 1, (p, r))
    return (
        periodic_tent(F(p + r, 14), z)
        - periodic_tent(F(r - p, 14), z)
    ) / (p * r)


def fiber_minimum(p: int, r: int) -> tuple[F, tuple[int, ...]]:
    # Every tent breakpoint is on the fourteen-grid, so this is complete.
    rows = tuple((fiber_density(p, r, F(k, 14)), k) for k in range(14))
    minimum = min(value for value, _ in rows)
    return minimum, tuple(k for value, k in rows if value == minimum)


def selected_debt(
    body: tuple[int, ...], pair: tuple[int, int], minimum: int, maximum: int
):
    """Worst exact singleton debt for six distinct levels and fixed extremes."""
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
    level_tuple = tuple(levels)
    require(all(level is not None for level in level_tuple), level_tuple)
    debt = sum(
        (
            F(label, 7 * (level * ruler - label))
            for label, level in zip(body, level_tuple)
        ),
        F(0),
    )
    return debt, level_tuple, ruler


def perturbation_loss(
    body: tuple[int, ...], pair: tuple[int, int], g: int, p: int, r: int
) -> F:
    ruler = 14 * lcm(*body)
    a, b = body[pair[0]], body[pair[1]]
    require(g * ruler * p > a and g * ruler * r > b, (body, pair, g, p, r))
    return F(a * (p + 1), g * ruler * p - a) + F(
        b * (r + 1), g * ruler * r - b
    )


def boundary_displacement_sum(n: int, c: int, g: int, ruler: int, theta: F) -> F:
    """Exact paired-root displacement used as a hostile lemma audit."""
    epsilon = F(c, g * ruler)
    require(F(0) < epsilon < n, (n, c, g, ruler))
    roots = set()
    delta = F(1, 14)
    for sign in (-1, 1):
        offset = theta + sign * delta
        for shift in range(-n - 3, n + 4):
            beta = offset + shift
            if F(0) <= beta <= n:
                roots.add((sign, beta))
    displacement = sum(
        (epsilon * beta / (n * (n - epsilon)) for _, beta in roots), F(0)
    )
    bound = F(c * (n + 1), g * ruler * n - c)
    require(displacement <= bound, (n, c, g, ruler, theta, displacement, bound))
    return bound - displacement


def full_tooth_overlap(
    ruler: int, a: int, p: int, b: int, q: int, cell: int
) -> F:
    """Independent integer two-pointer overlap of two untruncated combs."""
    require(p >= 1 and q >= 1, (p, q))
    r1 = (a * cell) % ruler
    r2 = (b * cell) % ruler
    require(
        ruler // 14 <= r1 <= 13 * ruler // 14 - a,
        (ruler, a, p, cell, r1),
    )
    require(
        ruler // 14 <= r2 <= 13 * ruler // 14 - b,
        (ruler, b, q, cell, r2),
    )
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = U.T.residual_bodies()
    require(len(bodies) == BODY_COUNT, len(bodies))
    repeated_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (repeated_exceptions & set(bodies)), repeated_exceptions & set(bodies))

    # Restricted primitive floor: a complete small-product bank plus Fourier.
    small_channels = tuple(
        (p, r)
        for p in range(1, FOURIER_PRODUCT_START)
        for r in range(3 * p + 1, 6 * p)
        if p >= 2 and gcd(p, r) == 1 and p * r < FOURIER_PRODUCT_START
    )
    small_floor_rows = tuple(
        sorted(
            (
                minimum,
                p,
                r,
                minimizers,
                p * r,
            )
            for p, r in small_channels
            for minimum, minimizers in (fiber_minimum(p, r),)
        )
    )
    require(small_floor_rows == EXPECTED_SMALL_FLOOR_ROWS, small_floor_rows)
    equality_channels = tuple(
        (p, r) for value, p, r, _, _ in small_floor_rows if value == RESTRICTED_FLOOR
    )
    require(equality_channels == ((2, 11), (3, 11)), equality_channels)
    fourier_floor_at_68 = F(1, 49) - F(1, 2 * FOURIER_PRODUCT_START)
    require(fourier_floor_at_68 > RESTRICTED_FLOOR, fourier_floor_at_68)

    # Hostile finite audit of the all-real boundary-displacement proof.
    boundary_digest = hashlib.sha256()
    boundary_slacks = []
    for n in range(2, 41):
        for c in range(1, 15):
            for theta_index in range(14):
                slack = boundary_displacement_sum(
                    n, c, BASE_SCALE, 168, F(theta_index, 14)
                )
                row = (slack, n, c, theta_index)
                boundary_slacks.append(row)
                boundary_digest.update(f"{row}\n".encode())

    # P>=13: Fourier, loss, and debt all reduce to g=6,P=13,R=40.
    large_floor = F(1, 49) - F(1, 2 * LARGE_P_START * LARGE_R_BASE)
    require(large_floor == F(991, 50960), large_floor)
    large_rows = []
    for body in bodies:
        for i in range(6):
            for j in range(6):
                if i == j:
                    continue
                pair = (i, j)
                debt, levels, _ = selected_debt(
                    body,
                    pair,
                    BASE_SCALE * LARGE_P_START,
                    BASE_SCALE * LARGE_R_BASE,
                )
                loss = perturbation_loss(
                    body, pair, BASE_SCALE, LARGE_P_START, LARGE_R_BASE
                )
                row = (large_floor - loss - debt, body, pair, large_floor, loss, debt, levels)
                require(row[0] > 0, ("large-P base failure", row))
                large_rows.append(row)
    large_weakest = min(large_rows)
    require(large_weakest == EXPECTED_LARGE_WEAKEST, large_weakest)

    # Exact finite bank P=2..12.
    finite_channels = tuple(
        (p, r)
        for p in range(2, LARGE_P_START)
        for r in range(3 * p + 1, 6 * p)
        if gcd(p, r) == 1
    )
    channel_histogram = tuple(sorted(Counter(p for p, _ in finite_channels).items()))
    require(channel_histogram == EXPECTED_CHANNELS_BY_P, channel_histogram)
    require(len(finite_channels) == FINITE_CHANNEL_COUNT, len(finite_channels))

    channel_minima_digest = hashlib.sha256()
    finite_weakest = None
    failures = []
    finite_count = 0
    channel_profiles = []
    for p, r in finite_channels:
        floor, minimizers = fiber_minimum(p, r)
        imported_floor, imported_minimizers = U.LOW.fiber_minimum(p, r)
        require((floor, minimizers) == (imported_floor, imported_minimizers), (p, r))
        channel_weakest = None
        channel_failure_count = 0
        for body in bodies:
            for i in range(6):
                for j in range(6):
                    if i == j:
                        continue
                    pair = (i, j)
                    debt, levels, ruler = selected_debt(
                        body, pair, BASE_SCALE * p, BASE_SCALE * r
                    )
                    imported_debt, imported_levels, imported_ruler = U.selected_debt(
                        body, pair, BASE_SCALE * p, BASE_SCALE * r
                    )
                    require(
                        (debt, levels, ruler)
                        == (imported_debt, imported_levels, imported_ruler),
                        (body, pair, p, r),
                    )
                    loss = perturbation_loss(body, pair, BASE_SCALE, p, r)
                    row = (
                        floor - loss - debt,
                        body,
                        pair,
                        p,
                        r,
                        floor,
                        loss,
                        debt,
                        levels,
                        minimizers,
                    )
                    if finite_weakest is None or row < finite_weakest:
                        finite_weakest = row
                    if channel_weakest is None or row < channel_weakest:
                        channel_weakest = row
                    if row[0] <= 0:
                        failures.append(row)
                        channel_failure_count += 1
                    finite_count += 1
        profile = (p, r, floor, minimizers, channel_weakest, channel_failure_count)
        channel_profiles.append(profile)
        channel_minima_digest.update(f"{profile}\n".encode())

    require(finite_count == FINITE_ROW_COUNT, finite_count)
    require(finite_count == FINITE_CHANNEL_COUNT * BODY_COUNT * ORIENTATION_COUNT, finite_count)
    require(finite_weakest == EXPECTED_FINITE_WEAKEST, finite_weakest)
    require(len(failures) == ANALYTIC_FAILURE_COUNT, len(failures))
    failure_body_histogram = tuple(sorted(Counter(row[1] for row in failures).items()))
    require(failure_body_histogram == EXPECTED_FAILURE_BODY_HISTOGRAM, failure_body_histogram)

    failure_digest = hashlib.sha256()
    for row in failures:
        failure_digest.update(f"{row}\n".encode())

    # Located repair: independent integer engine and promoted engine on every cell.
    direct_digest = hashlib.sha256()
    direct_rows = []
    engine_comparisons = 0
    for analytic_row in failures:
        _, body, pair, p, r, floor, loss, debt, levels, minimizers = analytic_row
        ruler, safe_ranges = T.R.safe_cell_ranges(body)
        best = None
        for left, right in safe_ranges:
            for cell in range(left, right):
                independent = full_tooth_overlap(
                    ruler,
                    body[pair[0]],
                    BASE_SCALE * p,
                    body[pair[1]],
                    BASE_SCALE * r,
                    cell,
                )
                imported = B.intersection_mass(
                    T.R.reflected_level_arcs(
                        ruler, body[pair[0]], BASE_SCALE * p, cell
                    ),
                    T.R.reflected_level_arcs(
                        ruler, body[pair[1]], BASE_SCALE * r, cell
                    ),
                )
                require(
                    independent == imported,
                    ("engine disagreement", body, pair, p, r, cell, independent, imported),
                )
                candidate = (
                    independent - debt,
                    body,
                    pair,
                    p,
                    r,
                    cell,
                    independent,
                    debt,
                    levels,
                )
                if best is None or candidate > best:
                    best = candidate
                engine_comparisons += 1
        require(best is not None and best[0] > 0, ("located failure", analytic_row, best))
        direct_rows.append(best)
        direct_digest.update(f"{best}|analytic={analytic_row}\n".encode())

    direct_weakest = min(direct_rows)
    require(direct_weakest == EXPECTED_DIRECT_WEAKEST, direct_weakest)

    digests = {
        "channel_minima": channel_minima_digest.hexdigest(),
        "failure": failure_digest.hexdigest(),
        "direct": direct_digest.hexdigest(),
        "boundary_audit": boundary_digest.hexdigest(),
    }
    for label, expected in (
        ("channel_minima", EXPECTED_CHANNEL_MINIMA_DIGEST),
        ("failure", EXPECTED_FAILURE_DIGEST),
        ("direct", EXPECTED_DIRECT_DIGEST),
        ("boundary_audit", EXPECTED_BOUNDARY_AUDIT_DIGEST),
    ):
        if expected is not None:
            require(digests[label] == expected, (label, digests[label], expected))

    semantic_payload = (
        tuple(bodies),
        tuple(sorted(repeated_exceptions)),
        small_floor_rows,
        equality_channels,
        fourier_floor_at_68,
        min(boundary_slacks),
        large_weakest,
        channel_histogram,
        tuple(channel_profiles),
        finite_count,
        finite_weakest,
        failure_body_histogram,
        tuple(failures),
        tuple(direct_rows),
        engine_comparisons,
        digests,
        "561;m>=2;3<Q/m<6;Q/m not in {4,5};gcd(m,Q)<=5",
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, (semantic, EXPECTED_SEMANTIC))

    lines = [
        "LRC14 reflected midratio primitive floor and g6 closure (THM-3210)",
        f"universe=residual_bodies:{len(bodies)};physical_orientations:{ORIENTATION_COUNT};same_level_exceptions_disjoint:{tuple(sorted(repeated_exceptions))}",
        "physical_map=m=gP;Q=gR;g=gcd(m,Q);gcd(P,R)=1;3P<R<6P;closed primitive rays (1,4),(1,5) imply P>=2,R>=7",
        f"restricted_fibre_floor=min_z F_(P,R)>=1/77;equality_channels:{equality_channels};small_product_channels:{len(small_floor_rows)};Fourier_handles_PR>={FOURIER_PRODUCT_START};floor_at_68:{qtext(fourier_floor_at_68)}",
        "one_sided_boundary_lemma=mu(A_epsilon triangle A_0)<=c(N+1)/(gLN-c);two-clause loss=a(P+1)/(gLP-a)+b(R+1)/(gLR-b)",
        f"boundary_hostile_audit=rows:{len(boundary_slacks)};weakest_slack:{min(boundary_slacks)}",
        f"large_P_monotone_base=g:{BASE_SCALE};P:{LARGE_P_START};R:{LARGE_R_BASE};rows:{len(large_rows)};weakest_margin:{qtext(large_weakest[0])};body:{large_weakest[1]};pair:{large_weakest[2]};Fourier_floor:{qtext(large_weakest[3])};loss:{qtext(large_weakest[4])};debt:{qtext(large_weakest[5])}",
        f"finite_bank=P2..12;channels:{len(finite_channels)};per_P:{channel_histogram};rows:{finite_count};formula:{FINITE_CHANNEL_COUNT}*{BODY_COUNT}*{ORIENTATION_COUNT};analytic_failures:{len(failures)};failure_body_histogram:{failure_body_histogram}",
        f"finite_weakest_analytic={qtext(finite_weakest[0])}@body={finite_weakest[1]},pair={finite_weakest[2]},P={finite_weakest[3]},R={finite_weakest[4]},floor={qtext(finite_weakest[5])},loss={qtext(finite_weakest[6])},debt={qtext(finite_weakest[7])}",
        f"located_repairs={len(direct_rows)};all_positive;weakest_margin:{qtext(direct_weakest[0])};body:{direct_weakest[1]};pair:{direct_weakest[2]};P:{direct_weakest[3]};R:{direct_weakest[4]};cell:{direct_weakest[5]};overlap:{qtext(direct_weakest[6])};debt:{qtext(direct_weakest[7])}",
        f"independent_engine=integer full-tooth two-pointer equals promoted Fraction interval engine on every repair cell;comparisons:{engine_comparisons}",
        "monotonicity=Fourier floor increases with PR;each perturbation loss decreases in its primitive coordinate and in g;every exact debt denominator increases with P,R,g",
        "assembly=inside THM-2941 reflected sufficient family every live midratio primitive channel closes for gcd scale g>=6",
        "residual=561 bodies;m>=2;3<Q/m<6;Q/m not in {4,5};gcd(m,Q)<=5",
        "scope=sufficient-certificate residual only;not a physical-survivor census and not LRC14",
        f"channel_minima_digest={digests['channel_minima']}",
        f"failure_digest={digests['failure']}",
        f"direct_digest={digests['direct']}",
        f"boundary_audit_digest={digests['boundary_audit']}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"upstream_source_sha256={EXPECTED_UPSTREAM_SOURCE_SHA256}",
        f"upstream_output_sha256={EXPECTED_UPSTREAM_OUTPUT_SHA256}",
        f"upstream_semantic_sha256={EXPECTED_UPSTREAM_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
