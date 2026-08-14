#!/usr/bin/env python3
r"""Exact support and endpoint referee for the D=14 complement-clock terminal.

The mathematical contradiction uses the cited lonely-runner theorem for eight
integer speeds.  This companion independently freezes the finite input: the
canonical body's safe-address support modulo 14, the 26 denominator shapes,
and the pointwise identity

    Y_14(alpha) = [1/14,13/14] \ D_(14 alpha).
"""

from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations_with_replacement
from math import lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/lrc14_two_drift_body_projection_support_thm2928.py"
EXPECTED_SUPPORT_SHA256 = "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"
BODY = (1, 2, 3, 4, 6, 12)
EXPECTED_RANGES = (
    (12, 13), (15, 26), (30, 39), (45, 52), (60, 69), (71, 78),
    (90, 97), (99, 108), (116, 123), (129, 138), (142, 153), (155, 156),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_support():
    require(file_hash(SUPPORT_PATH) == EXPECTED_SUPPORT_SHA256, "support dependency changed")
    spec = spec_from_file_location("thm3363_support", SUPPORT_PATH)
    require(spec is not None and spec.loader is not None, "support import")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def circular_distance(value):
    residue = value - value.numerator // value.denominator
    return min(residue, 1 - residue)


def danger(speed, x):
    return circular_distance(speed * x) < Q(1, 14)


def in_projected_carrier(alpha, support, x):
    return any(
        Q(0) <= 14 * x - r <= Q(1) and not danger(alpha, 14 * x - r)
        for r in support
    )


def identity_control(alpha, support):
    points = {Q(0), Q(1)} | {Q(r, 14) for r in range(15)}
    speed = 14 * alpha
    for k in range(speed + 1):
        for sign in (-1, 1):
            x = Q(14 * k + sign, 14 * speed)
            if Q(0) <= x <= Q(1):
                points.add(x)
    ordered = sorted(points)
    probes = set(ordered)
    probes.update((left + right) / 2 for left, right in zip(ordered, ordered[1:]))
    for x in probes:
        left = in_projected_carrier(alpha, support, x)
        right = Q(1, 14) <= x <= Q(13, 14) and not danger(speed, x)
        require(left == right, ("projected identity", alpha, x, left, right))
    return len(probes)


def main():
    support_module = load_support()
    L, ranges = support_module.safe_cell_ranges(BODY)
    require(L == 168, L)
    require(tuple(ranges) == EXPECTED_RANGES, ranges)
    addresses = {j for left, right in ranges for j in range(left, right)}
    support = tuple(sorted({j % 14 for j in addresses}))
    require(len(addresses) == 88, len(addresses))
    require(support == tuple(range(1, 13)), support)

    shapes = tuple(
        values
        for values in combinations_with_replacement((2, 7, 14), 6)
        if lcm(*values) == 14
    )
    require(len(shapes) == 26, len(shapes))
    # Record the more useful composition triples directly.
    compositions = tuple(sorted((row.count(2), row.count(7), row.count(14)) for row in shapes))
    require(len(set(compositions)) == 26, "denominator compositions")

    probe_count = sum(identity_control(alpha, support) for alpha in range(1, 65))
    require(not danger(1, Q(1, 14)), "D_1 endpoint must be safe")
    require(danger(14, Q(1, 14)), "aligned quotient clock must own the endpoint")

    semantic = sha256(
        repr((BODY, L, tuple(ranges), tuple(sorted(addresses)), support, compositions, probe_count)).encode()
    ).hexdigest()
    print("LRC14 D=14 COMPLEMENT-CLOCK TERMINAL REFEREE")
    print("status=VERIFIED-EXACT finite support and pointwise identity; small-LRC contradiction analytic/cited")
    print(f"body={BODY};L={L};safe_addresses={len(addresses)};safe_ranges={tuple(ranges)}")
    print(f"S14={support};denominator_shapes={len(shapes)}")
    print(f"alpha_controls=64;endpoint_and_midpoint_probes={probe_count}")
    print("identity=Y14(alpha)=[1/14,13/14]\\D_(14alpha)")
    print("completion=T=D_1 union D_(14alpha) union six quotient danger combs")
    print(f"support_dependency_sha256={file_hash(SUPPORT_PATH)}")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
