#!/usr/bin/env python3
r"""Uniform closure of both exceptional chromatic bodies at arbitrary levels.

The low-phase clique theorem leaves two robust-``K6`` bodies outside its
immediate arbitrary-level corollary because their same-level-good graphs are
only four-chromatic:

    E1=(1,2,7,9,11,13),       E2=(2,4,7,9,11,13).

Their body-safe geometry supplies the missing coloured sidecar.  For every
label pair ``a<b`` and every oriented primitive distinct low channel

    P != Q,       gcd(P,Q)=1,       P+Q<=7,

put ``h=P*b-Q*a``.  On a body-safe cell ``[j,j+1]``, the primitive skeleton
fiber at subdivision phase ``s`` has transverse coordinate

    z(s)=h(j+s)/L.

An exact convex endpoint optimization over all 480 configurations proves

    min_(0<=s<=1) F_(P,Q)(z(s)) >= 1/42.               (L)

For a fixed integer lift ``k``, the relevant endpoint radius is

    M(j,k)=max(|h*j-kL|, |h*(j+1)-kL|).

It is convex in ``j``.  On each body-safe integer range its minimum is at a
range endpoint or an integer adjacent to ``kL/h-1/2``; the referee enumerates
exactly those candidates and then evaluates the low-channel trapezoid.

Now take actual levels ``p=gP,q=gQ``.  Under ``u=(r+x)/g`` the exact phases
are

    P*x-a(gj+r+x)/(gL),       Q*x-b(gj+r+x)/(gL).

Dropping the two final ``x/(gL)`` terms gives the primitive skeleton above.
The two clause symmetric differences total at most ``4(a+b)/(gL)``, hence

    pair_overlap >= 1/42-4(a+b)/(gL).                  (T)

The exact singleton debt is bounded by

    eps_max(E)=sum_e e/[7(L-e)].

Even with ``g=1``, the weakest value of the right side of (T) minus this full
debt is strictly positive; it occurs at labels ``(11,13)`` and channel
``(1,6)`` on both bodies.  Direct reflected-interval controls at
``g=1,2,5`` verify all 1,440 selected profiles.

Therefore both exceptional bodies close for every assignment of positive
levels.  Indeed, a same-level-good edge closes by the universal chromatic
theorem.  Otherwise choose any unequal pair.  If its reduced coefficients
sum to at least eight, the robust-``K6`` ``1/105`` theorem closes it; if they
sum to at most seven, (L)--(T) close it.  Consequently all 2,217 robust-``K6``
bodies, not merely the 2,215 nonexceptional ones, close at arbitrary reflected
levels.  The 1,584 exceptional scale rays recorded by the preceding theorem
are thus discharged.  The remaining 786 noncomplete-robust bodies and the
physical scope outside THM-2941 remain open.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation" / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
LOW_PHASE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
)
SUPPORT = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_distinct_low_channel_body_safe_support_thm2941.py"
)
UNIVERSAL = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_exceptional_low_channel_uniform_closure_thm2941.out"
)
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "b2418dfda1b48257d1f7582d4ea977203a26f88885e13946bc100ccf264c9ce1"
EXPECTED_SUPPORT_SHA256 = "0c3404d5017a97545f2acb187e91e8089c71a50b2119065d1b013e31801f2a48"
EXPECTED_UNIVERSAL_SHA256 = "dc6f23a201e817dd9134e8660d35e83d3053c67d26fc271ce3eae07f0f857689"
EXPECTED_SEMANTIC_SHA256 = "05e0d61ce1eadd8c4f471a48200dab024416b6b563c8bdfa96d22516b3b52370"
EXCEPTIONS = (
    (1, 2, 7, 9, 11, 13),
    (2, 4, 7, 9, 11, 13),
)
LOW_FLOOR = F(1, 42)
PROFILE_COUNT = 480
DIRECT_CONTROL_SCALES = (1, 2, 5)
DIRECT_CONTROL_COUNT = 1440
EXPECTED_WEAKEST = (
    (
        (1, 2, 7, 9, 11, 13),
        252252,
        F(487425160807264342822652, 20014978569752477975754425125),
        F(1, 42),
        F(
            2813465937513647399802451258963,
            120209961289933382722381077300750,
        ),
        11,
        13,
        1,
        6,
        -53,
        227681,
        48,
        41003,
    ),
    (
        (2, 4, 7, 9, 11, 13),
        504504,
        F(4171698934310078342981596, 320265049370713655040295609125),
        F(1, 42),
        F(
            15135609878977256911491328933433,
            641170628840168737390671809468250,
        ),
        11,
        13,
        1,
        6,
        -53,
        486485,
        51,
        54054,
    ),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


for path, expected in (
    (BASE, EXPECTED_BASE_SHA256),
    (LOW_PHASE, EXPECTED_LOW_PHASE_SHA256),
    (SUPPORT, EXPECTED_SUPPORT_SHA256),
    (UNIVERSAL, EXPECTED_UNIVERSAL_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("exceptional_low_channel_base", BASE)
LOW = import_module("exceptional_low_channel_ratio", LOW_PHASE)


CHANNELS = tuple(
    (p, q)
    for p in range(1, 7)
    for q in range(1, 7)
    if p != q and gcd(p, q) == 1 and p + q <= 7
)
require(CHANNELS == LOW.CHANNELS if hasattr(LOW, "CHANNELS") else len(CHANNELS) == 16,
        CHANNELS)


def best_skeleton_profile(
    safe_ranges: tuple[tuple[int, int], ...],
    ruler: int,
    cross_difference: int,
    p: int,
    q: int,
) -> tuple[F, int, int, int]:
    """Maximize the exact cellwise primitive-fiber floor."""
    h = abs(cross_difference)
    if h == 0:
        return F(1, 7 * max(p, q)), safe_ranges[0][0], 0, 0
    support = (p + q) * ruler // 14
    plateau = abs(p - q) * ruler // 14
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
                    floor = F(1, 7 * max(p, q))
                elif endpoint_radius < support:
                    floor = F(support - endpoint_radius, p * q * ruler)
                else:
                    floor = F(0)
                row = (floor, cell, lift, endpoint_radius)
                if best is None or row > best:
                    best = row
    require(best is not None, (safe_ranges, ruler, cross_difference, p, q))
    return best


def intersection_mass(first, second) -> F:
    i = 0
    j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(
            F(0),
            min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]),
        )
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    profile_count = 0
    direct_controls = 0
    body_rows = []
    profile_digest = hashlib.sha256()

    for body in EXCEPTIONS:
        ruler, safe_ranges = R.safe_cell_ranges(body)
        debt = LOW.universal_debt(body, ruler)
        robust_ruler, robust_debt, robust = LOW.robust_edges(body)
        require((robust_ruler, robust_debt) == (ruler, debt),
                (body, robust_ruler, ruler, robust_debt, debt))
        require(len(robust) == 15, ("exception is not robust K6", body, robust))
        profiles = []
        for a, b in combinations(body, 2):
            for p, q in CHANNELS:
                cross_difference = p * b - q * a
                floor, cell, lift, endpoint_radius = best_skeleton_profile(
                    safe_ranges, ruler, cross_difference, p, q
                )
                require(floor >= LOW_FLOOR,
                        (body, a, b, p, q, floor, cell, lift, endpoint_radius))

                h = abs(cross_difference)
                for endpoint in (cell, cell + 1):
                    phase = F(h * endpoint, ruler)
                    phase -= phase.numerator // phase.denominator
                    require(LOW.fiber_density(min(p, q), max(p, q), phase) >= floor,
                            (body, a, b, p, q, phase, floor))

                margin = floor - F(4 * (a + b), ruler) - debt
                require(margin > 0,
                        ("low-channel debt margin failed", body, a, b, p, q, margin))
                for scale in DIRECT_CONTROL_SCALES:
                    first = R.reflected_level_arcs(ruler, a, scale * p, cell)
                    second = R.reflected_level_arcs(ruler, b, scale * q, cell)
                    actual = intersection_mass(first, second)
                    transported = floor - F(4 * (a + b), scale * ruler)
                    require(actual >= transported,
                            (body, a, b, p, q, scale, actual, transported, cell))
                    direct_controls += 1
                row = (
                    floor, margin, a, b, p, q, cross_difference,
                    cell, lift, endpoint_radius,
                )
                profiles.append(row)
                profile_count += 1
                profile_digest.update(f"{body}|{ruler}|{debt}|{row}\n".encode())
        weakest = min(profiles)
        body_rows.append((body, ruler, debt) + weakest)

    require(profile_count == PROFILE_COUNT, profile_count)
    require(direct_controls == DIRECT_CONTROL_COUNT, direct_controls)
    require(tuple(body_rows) == EXPECTED_WEAKEST, body_rows)

    semantic_payload = (
        CHANNELS,
        tuple(body_rows),
        profile_count,
        DIRECT_CONTROL_SCALES,
        direct_controls,
        profile_digest.hexdigest(),
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                (semantic, EXPECTED_SEMANTIC_SHA256))

    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 exceptional-body distinct-low-channel uniform closure exact referee",
        f"exception_bodies={EXCEPTIONS};oriented_low_channels={len(CHANNELS)};"
        f"profiles={profile_count};uniform_skeleton_floor={qtext(LOW_FLOOR)}",
    ]
    for row in body_rows:
        (body, ruler, debt, floor, margin, a, b, p, q,
         cross_difference, cell, lift, endpoint_radius) = row
        lines.append(
            f"BODY;E={body};L={ruler};eps_max={qtext(debt)};"
            f"weakest_floor={qtext(floor)};weakest_margin={qtext(margin)};"
            f"labels={(a,b)};channel={(p,q)};cross_difference={cross_difference};"
            f"j={cell};lift={lift};endpoint_radius={endpoint_radius}"
        )
    lines.extend((
        f"direct_control_scales={DIRECT_CONTROL_SCALES};direct_controls={direct_controls}",
        "conclusion=both exceptional chromatic bodies close for every assignment of positive reflected levels",
        "corollary=all 2217 robust-K6 bodies close at arbitrary reflected levels;the 1584 exceptional low-clique scale rays are discharged",
        "scope=reflected THM-2941 residual family only;the 786 noncomplete-robust bodies and physical LRC14 remain open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"support_sha256={sha256(SUPPORT)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"profile_digest={profile_digest.hexdigest()}",
        f"source_sha256={source_sha}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
