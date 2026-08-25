#!/usr/bin/env python3
"""Exact finite probe for the THM-4030/4032/4041 pack-safe sidecar.

This deliberately separates two universes:

1. the finite exception-profile producer banks declared by the three theorem
   audits, paired with their canonical typed H10/K11 controls; and
2. the actual THM-3818 rank-eleven physical survivor branch.

The declared banks are controls, not physical survivors.  We apply the
THM-3818 primitive-ratio and crossing-height filters explicitly so this
distinction cannot be lost in the output.
"""

from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations
from math import gcd
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")

ONE = F(1)
RHO = F(1, 14)
Q3818 = 91**6
ROOT = Path(__file__).resolve().parents[2]


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def norm(x):
    residue = x % ONE
    return min(residue, ONE - residue)


def pack_safe(H, y):
    """Closed divided-pack safe set: every clearance is >= 1/14."""
    return all(norm(h * y) >= RHO for h in H)


def bad_labels(d, exceptions, y):
    """Strictly dangerous lift labels for each exception."""
    return tuple(
        tuple(
            j for j in range(d)
            if norm(w * (y + j) / d) < RHO
        )
        for w in exceptions
    )


def fully_spoiled(d, exceptions, y):
    return set().union(*(set(mask) for mask in bad_labels(d, exceptions, y))) == set(range(d))


def surviving_labels(d, exceptions, y):
    return tuple(
        j for j in range(d)
        if all(norm(w * (y + j) / d) >= RHO for w in exceptions)
    )


def wall_points(H, d, exceptions):
    """All exact predicate walls in the linear fundamental domain [0,1]."""
    walls = {F(0), F(1)}
    for h in H:
        for integer in range(-1, h + 2):
            for sign in (-1, 1):
                y = (F(integer) + sign * RHO) / h
                if 0 <= y <= 1:
                    walls.add(y)
    for w in exceptions:
        for j in range(d):
            for integer in range(-1, w + 2):
                for sign in (-1, 1):
                    y = d * (F(integer) + sign * RHO) / w - j
                    if 0 <= y <= 1:
                        walls.add(y)
    return tuple(sorted(walls))


def canonical_linear_components(walls, wall_state, cell_state):
    """Intervals (left,left_closed,right,right_closed) plus singleton points.

    This is a linearized [0,1] representation.  The endpoint 1 is the same
    circle point as 0; retaining both makes every endpoint convention visible.
    """
    cell_count = len(walls) - 1
    components = []
    incident_active_wall = set()
    i = 0
    while i < cell_count:
        if not cell_state[i]:
            i += 1
            continue
        start = i
        end = i
        while (
            end + 1 < cell_count
            and cell_state[end + 1]
            and wall_state[end + 1]
        ):
            end += 1
        left_index = start
        right_index = end + 1
        left_closed = wall_state[left_index]
        right_closed = wall_state[right_index]
        if left_closed:
            incident_active_wall.add(left_index)
        if right_closed:
            incident_active_wall.add(right_index)
        for k in range(start + 1, end + 1):
            if wall_state[k]:
                incident_active_wall.add(k)
        components.append((walls[left_index], left_closed, walls[right_index], right_closed))
        i = end + 1

    # Use 0, not the duplicate endpoint 1, for an isolated seam point.
    for k in range(len(walls) - 1):
        if wall_state[k] and k not in incident_active_wall:
            left_active = k > 0 and cell_state[k - 1]
            right_active = k < cell_count and cell_state[k]
            if not left_active and not right_active:
                components.append((walls[k], True, walls[k], True))
    return tuple(sorted(components))


def analyze(H, d, exceptions):
    walls = wall_points(H, d, exceptions)
    predicates = {
        "G": lambda y: pack_safe(H, y),
        "Sigma": lambda y: fully_spoiled(d, exceptions, y),
        "R": lambda y: pack_safe(H, y) and not fully_spoiled(d, exceptions, y),
    }
    result = {"walls": walls}
    for name, predicate in predicates.items():
        wall_state = tuple(predicate(y) for y in walls)
        cell_state = tuple(
            predicate((left + right) / 2)
            for left, right in zip(walls, walls[1:])
        )
        measure = sum(
            (right - left for left, right, active in zip(walls, walls[1:], cell_state) if active),
            F(0),
        )
        components = canonical_linear_components(walls, wall_state, cell_state)
        witness = None
        for k, active in enumerate(wall_state[:-1]):
            if active:
                witness = walls[k]
                break
        if witness is None:
            for left, right, active in zip(walls, walls[1:], cell_state):
                if active:
                    witness = (left + right) / 2
                    break
        result[name] = {
            "measure": measure,
            "components": components,
            "witness": witness,
            "wall_state": wall_state,
            "cell_state": cell_state,
        }

    intersection_measure = sum(
        (right - left for left, right in zip(walls, walls[1:])
         if pack_safe(H, (left + right) / 2)
         and fully_spoiled(d, exceptions, (left + right) / 2)),
        F(0),
    )
    require(
        result["G"]["measure"] == intersection_measure + result["R"]["measure"],
        "G=(G intersect Sigma) disjoint-union R measure bookkeeping",
    )
    return result


def factor(n):
    factors = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            factors.append((p, exponent))
        p += 1
    if n > 1:
        factors.append((n, 1))
    return tuple(factors)


def atlas_pair(p, q):
    if not (0 < p < q and gcd(p, q) == 1 and p + q <= 356):
        return False
    return all(prime % 3 == 2 and exponent <= 2 for prime, exponent in factor(p + q))


def physical_filter(d, exceptions, body_H, pair_H):
    """Apply the explicit primitive/gcd/height gates inherited from THM-3818."""
    s = 1
    t = d
    p, q = pair_H
    body = tuple(d * h for h in body_H) + tuple(exceptions)
    pair = (t * p, t * q)
    row = body + pair
    basic = (
        len(body) == 11
        and len(row) == 13
        and len(set(row)) == 13
        and gcd(*body) == 1
        and gcd(s, t) == 1
        and gcd(p, q) == 1
        and t % d == 0
    )
    crossing_heights = tuple(
        max(pair_speed, body_speed) // gcd(pair_speed, body_speed)
        for pair_speed in pair
        for body_speed in body
    )
    height = all(value > Q3818 for value in crossing_heights)
    return {
        "basic": basic,
        "pair_atlas": atlas_pair(p, q),
        "height": height,
        "max_crossing_height": max(crossing_heights),
        "min_crossing_height": min(crossing_heights),
        "pass": basic and atlas_pair(p, q) and height,
    }


def full_row(d, H, exceptions):
    return tuple(sorted(tuple(d * h for h in H) + tuple(exceptions)))


def verify_reconstructed_witness(H, d, exceptions, analysis):
    y = analysis["R"]["witness"]
    if y is None:
        return None
    labels = surviving_labels(d, exceptions, y)
    require(labels, ("R witness has no safe lift", d, exceptions, y))
    row = full_row(d, H, exceptions)
    candidates = []
    for j in labels:
        x = (y + j) / d
        clearance = min(norm(w * x) for w in row)
        require(clearance >= RHO, ("reconstructed lift failed", d, exceptions, y, j))
        candidates.append((j, x % ONE, clearance))
    return y, tuple(candidates)


def digest(value):
    return sha256(repr(value).encode("utf-8")).hexdigest()


def load_canonical_module(stem):
    path = ROOT / "04-computation" / f"{stem}.py"
    spec = importlib.util.spec_from_file_location(f"scratch_import_{stem}", path)
    require(spec is not None and spec.loader is not None, ("cannot load", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def producer_rows():
    H2 = tuple(range(1, 11)) + (12,)
    H10 = tuple(range(1, 11))
    return {
        2: {
            "H": H2,
            "body_H": tuple(range(2, 11)),
            "pair_H": (1, 12),
            "profiles": tuple(combinations(range(1, 80, 2), 2)),
        },
        3: {
            "H": H10,
            "body_H": tuple(range(2, 10)),
            "pair_H": (1, 10),
            "profiles": tuple(combinations((x for x in range(1, 24) if x % 3), 3)),
        },
        4: {
            "H": H10,
            "body_H": (2, 3, 5, 6, 7, 8, 9, 10),
            "pair_H": (1, 4),
            "profiles": tuple(
                (2 * r, first, second)
                for r in range(1, 12, 2)
                for first, second in combinations(range(1, 20, 2), 2)
            ),
        },
    }


def main():
    banks = producer_rows()
    canonical_modules = {
        3: load_canonical_module("lrc14_d3_affine_defect_lattice_boundary_thm4032"),
        4: load_canonical_module("lrc14_d4_affine_defect_lattice_boundary_thm4030"),
    }
    canonical = {
        2: (1, 11),
        3: (1, 5, 11),
        4: (2, 9, 11),
    }
    canonical_spoiled_phase = {2: F(1, 11), 3: F(2, 11), 4: F(1, 11)}
    canonical_full_witness = {2: F(229, 560), 3: F(1, 14), 4: F(21, 22)}
    # A cheap independent noncontainment certificate for each whole declared
    # bank: if the first closed-safe endpoint is spoiled, the second is not.
    endpoint_certificate = {
        2: (F(1, 14), F(5, 56)),
        3: (F(1, 14), F(13, 140)),
        4: (F(15, 98), F(1, 14)),
    }

    print("LRC14 D2/D3/D4 PACK-SAFE INTERSECTION EXACT PROBE")
    print("status=FINITE-EXACT CONTROL BANK; NOT PHYSICAL-UNIVERSE EXHAUSTION")
    print("endpoint_convention=G_closed(>=1/14);Sigma_open(<1/14);R=G\\Sigma_closed")
    print("circle_representation=linearized_[0,1]_with_0~1_and_explicit_endpoint_flags")
    print(f"THM3818_crossing_height_Q={Q3818}")

    summaries = []
    for d in (2, 3, 4):
        bank = banks[d]
        H = bank["H"]
        profiles = bank["profiles"]
        positives = []
        containments = []
        survivors = []
        survivor_measures = []
        physical_pass = []
        pair_atlas_failures = 0
        height_failures = 0
        maximum_crossing_height = 0
        first_endpoint_spoiled = []
        both_endpoints_spoiled = []

        endpoint_first, endpoint_second = endpoint_certificate[d]
        require(pack_safe(H, endpoint_first), ("first endpoint not in G", d))
        require(pack_safe(H, endpoint_second), ("second endpoint not in G", d))

        for exceptions in profiles:
            analysis = analyze(H, d, exceptions)
            positive = analysis["Sigma"]["measure"] > 0
            if d == 2:
                g = gcd(*exceptions)
                expected = exceptions[0] // g + exceptions[1] // g > 7
                require(positive == expected, ("d2 certificate drift", exceptions))
            elif d == 3:
                module = canonical_modules[d]
                canonical_direct = module.direct_phase(exceptions) is not None
                canonical_certificate = module.defect_certificate(exceptions) is not None
                require(
                    positive == canonical_direct == canonical_certificate,
                    ("d3 canonical certificate drift", exceptions),
                )
            else:
                module = canonical_modules[d]
                canonical_direct = module.direct_covering_phase(*exceptions) is not None
                canonical_certificate = module.defect_certificate(*exceptions) is not None
                require(
                    positive == canonical_direct == canonical_certificate,
                    ("d4 canonical certificate drift", exceptions),
                )
            if not positive:
                continue
            positives.append(exceptions)
            if fully_spoiled(d, exceptions, endpoint_first):
                first_endpoint_spoiled.append(exceptions)
                if fully_spoiled(d, exceptions, endpoint_second):
                    both_endpoints_spoiled.append(exceptions)
            if analysis["R"]["witness"] is None:
                containments.append(exceptions)
            else:
                survivors.append((exceptions, verify_reconstructed_witness(H, d, exceptions, analysis)))
                survivor_measures.append((analysis["R"]["measure"], exceptions))

            filt = physical_filter(d, exceptions, bank["body_H"], bank["pair_H"])
            require(filt["basic"], ("declared control lost basic typing", d, exceptions))
            if not filt["pair_atlas"]:
                pair_atlas_failures += 1
            if not filt["height"]:
                height_failures += 1
            maximum_crossing_height = max(maximum_crossing_height, filt["max_crossing_height"])
            if filt["pass"]:
                physical_pass.append(exceptions)

        control = canonical[d]
        control_analysis = analyze(H, d, control)
        spoiled_y = canonical_spoiled_phase[d]
        require(pack_safe(H, spoiled_y), ("canonical phase not H-safe", d))
        require(fully_spoiled(d, control, spoiled_y), ("canonical phase not spoiled", d))
        x = canonical_full_witness[d]
        mapped_y = (d * x) % ONE
        require(pack_safe(H, mapped_y), ("canonical full witness not pack-safe", d, mapped_y))
        require(not fully_spoiled(d, control, mapped_y),
                ("canonical full witness did not survive spoilage", d, mapped_y))
        row = full_row(d, H, control)
        clearance = min(norm(w * x) for w in row)
        require(clearance >= RHO, ("canonical full-row witness failed", d, x))

        summary = (
            d,
            len(profiles),
            len(positives),
            len(containments),
            len(survivors),
            len(physical_pass),
            pair_atlas_failures,
            height_failures,
            maximum_crossing_height,
            tuple(containments),
        )
        summaries.append(summary)
        minimum_survivor_measure = min(measure for measure, _ in survivor_measures)
        minimum_measure_profiles = tuple(
            exceptions
            for measure, exceptions in survivor_measures
            if measure == minimum_survivor_measure
        )
        require(
            not both_endpoints_spoiled,
            ("two-point G\\Sigma certificate failed", d, both_endpoints_spoiled),
        )
        print(
            f"d={d};H={H};declared_profiles={len(profiles)};"
            f"certificate_positive={len(positives)};G_subset_Sigma={len(containments)};"
            f"G_minus_Sigma_nonempty={len(survivors)}"
        )
        print(
            f"d={d};THM3818_filter_pass={len(physical_pass)};"
            f"pair_atlas_failures={pair_atlas_failures};height_failures={height_failures};"
            f"bank_max_crossing_height={maximum_crossing_height}"
        )
        print(
            f"d={d};minimum_G_minus_Sigma_measure={minimum_survivor_measure};"
            f"attaining_profiles={minimum_measure_profiles}"
        )
        print(
            f"d={d};two_point_G_minus_Sigma_certificate={endpoint_certificate[d]};"
            f"first_endpoint_spoiled_profiles={tuple(first_endpoint_spoiled)};"
            f"both_endpoints_spoiled=()"
        )
        print(
            f"d={d};canonical_exceptions={control};spoiled_H_phase={spoiled_y};"
            f"surviving_H_phase={mapped_y};full_x={x};clearance={clearance};"
            f"G_measure={control_analysis['G']['measure']};"
            f"Sigma_measure={control_analysis['Sigma']['measure']};"
            f"R_measure={control_analysis['R']['measure']}"
        )
        print(f"d={d};canonical_G_components={control_analysis['G']['components']}")
        print(f"d={d};canonical_Sigma_components={control_analysis['Sigma']['components']}")
        print(f"d={d};canonical_R_components={control_analysis['R']['components']}")
        print(f"d={d};positive_profile_digest={digest(tuple(positives))}")
        print(f"d={d};survivor_digest={digest(tuple(survivors))}")

    require(all(summary[3] == 0 for summary in summaries),
            "some declared control-bank profile contains all of G(H)")
    require(all(summary[5] == 0 for summary in summaries),
            "a declared low-height producer unexpectedly entered THM-3818")
    print(f"summary_digest={digest(tuple(summaries))}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
