#!/usr/bin/env python3
"""Componentwise residue-ray law and pure-cone no-go for THM-2941.

Let a positive component of the six-body carrier be

    I_s=[a_s/L,b_s/L],

and put ``delta_s(z)=mu(I_s intersect D_z)-|I_s|/7``.  For the primitive
``P(x)=integral_0^x 1_{D_1}``, periodicity and mean ``1/7`` give
``P(x+m)=P(x)+m/7``.  Substitution in

    z mu(I_s intersect D_z)
      = P(z b_s/L)-P(z a_s/L)

proves, component by component,

    (z+L) delta_s(z+L)=z delta_s(z),
    delta_s(La)=0,
    A_s(L-b)=-A_s(b),       A_s(b)=z delta_s(z), z=b mod L.

If seven normalized tail combs cover ``I_s``, their multiplicity is at
least one pointwise, so the sum of their component excesses is nonnegative.
It is strict when ``|I_s|>=1/105``: every tail label is at least 15, hence
every danger tooth has length at most ``1/105``.  No one relatively open
tooth contains the closed component, and a finite proper open cover of a
connected interval has a positive-measure overlap.

The component law is exact but its *pure positive-cone relaxation* is
universally non-pruning.  If ``e`` is any body speed, the carrier is disjoint
from ``D_e`` up to endpoints, hence

    A_s(e)=-e|I_s|/7,
    A_s(L-e)=e|I_s|/7>0.

Thus every body already has a support-one ray in the strict positive
orthant.  The exact hostile triple below shows the same failure after three
specific directions and their actual reciprocal label weights are retained.
This is only a necessary-cone witness, never a cover; THM-2928 independently
closes the three-drift (aligned k=4) sector.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_component_residue_ray_cone_no_go_thm2941.out"
)
EXPECTED_SOURCE_SHA256 = (
    "5d25a955fe184d6c1a3d8b632b4bbf901dc996ee46ad67c5748836fcc7134404"
)
EXPECTED_SEMANTIC_SHA256 = (
    "1e6b330b2491feab5965aab5b4c157597891880620705eadf49e45797f732525"
)

TAIL_START = 15
LONG_THRESHOLD = F(1, 105)
BODY = (1, 2, 6, 8, 12, 14)
RESIDUES = (378, 1062, 1622)
EXTREMAL_WEIGHTS = (F(536, 2664), F(1001, 2664), F(1127, 2664))
EQUAL_WEIGHTS = (F(1, 3), F(1, 3), F(1, 3))
RECIPROCAL_WEIGHTS = tuple(F(1, z) for z in RESIDUES)
EXPECTED_TYPES = (
    (2, 4, (F(19, 196), F(19, 1372), F(-223, 4116))),
    (9, 4, (F(31, 392), F(-25, 2744), F(-81, 2744))),
    (16, 4, (F(3, 49), F(-11, 343), F(10, 343))),
    (51, 4, (F(45, 392), F(33, 1372), F(89, 1372))),
    (63, 4, (F(-1, 56), F(31, 392), F(31, 392))),
    (72, 6, (F(3, 49), F(24, 343), F(17, 343))),
)
EXPECTED_CERTIFICATE_VALUES = {
    "extremal": (
        F(25, 13986),
        F(0),
        F(235, 18648),
        F(1111, 18648),
        F(1111, 18648),
        F(1111, 18648),
    ),
    "equal": (
        F(233, 12348),
        F(37, 2744),
        F(20, 1029),
        F(559, 8232),
        F(55, 1176),
        F(62, 1029),
    ),
    "reciprocal": (
        F(836945, 3545036712),
        F(1293461, 7090073424),
        F(7373, 49236621),
        F(288593, 787785936),
        F(76957, 1012867632),
        F(38170, 147709863),
    ),
}


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot load", path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(SOURCE) == EXPECTED_SOURCE_SHA256, "source changed")
wall = load_module("component_ray_wall", SOURCE)


def component_amplitude(interval: tuple[int, int], label: int) -> F:
    """Return ``label * delta_s(label)`` exactly."""

    left, right = interval
    primitive = (
        wall.danger_primitive(label * right)
        - wall.danger_primitive(label * left)
    )
    return F(primitive, wall.RULER) - F(
        label * (right - left), 7 * wall.RULER
    )


def component_delta(interval: tuple[int, int], label: int) -> F:
    return component_amplitude(interval, label) / label


def weighted_values(
    columns: tuple[tuple[F, ...], ...], weights: tuple[F, ...]
) -> tuple[F, ...]:
    return tuple(
        sum(
            (
                weight * column[index]
                for weight, column in zip(weights, columns, strict=True)
            ),
            F(),
        )
        for index in range(len(columns[0]))
    )


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    require(
        F(1, 7 * TAIL_START) == LONG_THRESHOLD,
        "tail-tooth strictness threshold changed",
    )
    carrier = wall.carrier_for(BODY)
    L = 14 * math.lcm(*BODY)
    require(wall.RULER % L == 0, "hostile ruler mismatch")
    scale = wall.RULER // L
    require(
        all(left % scale == right % scale == 0 for left, right in carrier),
        "hostile component is off the body ruler",
    )
    grid_components = tuple(
        (left // scale, right // scale) for left, right in carrier
    )

    columns = tuple(
        tuple(component_amplitude(interval, residue) for interval in carrier)
        for residue in RESIDUES
    )
    component_types = tuple(zip(*columns, strict=True))
    type_counts = Counter(component_types)
    lengths_by_type: dict[tuple[F, ...], set[int]] = defaultdict(set)
    for interval, row in zip(grid_components, component_types, strict=True):
        lengths_by_type[row].add(interval[1] - interval[0])
    require(
        all(len(lengths) == 1 for lengths in lengths_by_type.values()),
        "hostile amplitude type is not length-determined",
    )
    type_table = tuple(
        sorted(
            (
                next(iter(lengths_by_type[row])),
                type_counts[row],
                row,
            )
            for row in type_counts
        )
    )
    require(type_table == EXPECTED_TYPES, "hostile amplitude types changed")

    # Full-period exact controls for the analytic component identities.
    recurrence_checks = antipodal_checks = 0
    identity_digest = hashlib.sha256()
    for residue in range(1, L):
        for component, interval in enumerate(carrier):
            delta = component_delta(interval, residue)
            next_delta = component_delta(interval, residue + L)
            amplitude = component_amplitude(interval, residue)
            opposite = component_amplitude(interval, L - residue)
            require(
                (residue + L) * next_delta == residue * delta,
                ("component recurrence failed", interval, residue),
            )
            require(
                opposite == -amplitude,
                ("component antipode failed", interval, residue),
            )
            identity_digest.update(
                repr((residue, component, delta, next_delta, amplitude, opposite)).encode()
            )
            recurrence_checks += 1
            antipodal_checks += 1
    aligned_checks = 0
    for multiple in range(1, 8):
        for component, interval in enumerate(carrier):
            delta = component_delta(interval, multiple * L)
            require(
                delta == 0,
                ("aligned component excess nonzero", interval, multiple),
            )
            identity_digest.update(repr(("aligned", multiple, component, delta)).encode())
            aligned_checks += 1

    long = tuple(
        F(right - left, wall.RULER) >= LONG_THRESHOLD
        for left, right in carrier
    )
    require(sum(long) == 14, "hostile forced-long component count changed")
    certificates: dict[str, tuple[F, ...]] = {}
    weights_by_name = (
        ("extremal", EXTREMAL_WEIGHTS),
        ("equal", EQUAL_WEIGHTS),
        ("reciprocal", RECIPROCAL_WEIGHTS),
    )
    for name, weights in weights_by_name:
        values = weighted_values(columns, weights)
        require(all(value >= 0 for value in values), (name, "negative component"))
        require(
            all(value > 0 for value, is_long in zip(values, long, strict=True) if is_long),
            (name, "non-strict forced-long component"),
        )
        by_type = tuple(
            values[component_types.index(row)]
            for _length, _count, row in type_table
        )
        require(
            by_type == EXPECTED_CERTIFICATE_VALUES[name],
            (name, "certificate values changed", by_type),
        )
        certificates[name] = values
    require(min(certificates["extremal"]) == 0, "extremal face moved")
    require(min(certificates["equal"]) == F(37, 2744), "equal margin moved")
    require(
        min(certificates["reciprocal"]) == F(76957, 1012867632),
        "reciprocal margin moved",
    )
    normalized_reciprocal = tuple(
        weight / sum(RECIPROCAL_WEIGHTS, F()) for weight in RECIPROCAL_WEIGHTS
    )
    require(
        normalized_reciprocal
        == (F(47849, 76031), F(17031, 76031), F(11151, 76031)),
        "normalized reciprocal weights changed",
    )

    # Exhaust every literal body root and every one of its six body speeds.
    # This is the decisive universal pure-cone no-go audit.
    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == 3_003, "body universe changed")
    universal_checks = 0
    minimum_positive = None
    universal_digest = hashlib.sha256()
    for body in bodies:
        body_carrier = wall.carrier_for(body)
        body_L = 14 * math.lcm(*body)
        body_scale = wall.RULER // body_L
        require(
            wall.RULER % body_L == 0
            and all(
                left % body_scale == right % body_scale == 0
                for left, right in body_carrier
            ),
            ("body carrier ruler changed", body),
        )
        for speed in body:
            reflected = body_L - speed
            require(reflected >= TAIL_START, ("reflected ray is not a tail", body, speed))
            for component, interval in enumerate(body_carrier):
                length = F(interval[1] - interval[0], wall.RULER)
                negative = -F(speed, 7) * length
                positive = F(speed, 7) * length
                actual_negative = component_amplitude(interval, speed)
                actual_positive = component_amplitude(interval, reflected)
                require(
                    actual_negative == negative,
                    ("body-safe negative ray failed", body, speed, interval),
                )
                require(
                    actual_positive == positive,
                    ("reflected positive ray failed", body, speed, interval),
                )
                require(positive > 0, ("nonpositive carrier component", body, interval))
                minimum_positive = (
                    positive
                    if minimum_positive is None
                    else min(minimum_positive, positive)
                )
                universal_digest.update(
                    repr(
                        (
                            body,
                            speed,
                            component,
                            (interval[0] // body_scale, interval[1] // body_scale),
                            negative,
                            positive,
                        )
                    ).encode()
                )
                universal_checks += 1
    require(universal_checks == 402_942, "universal audit count changed")
    require(minimum_positive == F(1, 17836), "minimum positive ray changed")

    values_by_type = {
        name: tuple(
            values[component_types.index(row)]
            for _length, _count, row in type_table
        )
        for name, values in certificates.items()
    }
    semantic_payload = (
        EXPECTED_SOURCE_SHA256,
        TAIL_START,
        LONG_THRESHOLD,
        BODY,
        L,
        grid_components,
        RESIDUES,
        type_table,
        long,
        tuple(weights_by_name),
        tuple(sorted(values_by_type.items())),
        normalized_reciprocal,
        recurrence_checks,
        antipodal_checks,
        aligned_checks,
        identity_digest.hexdigest(),
        len(bodies),
        universal_checks,
        minimum_positive,
        universal_digest.hexdigest(),
    )
    semantic_sha256 = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(
            semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            "component-ray semantic digest changed",
        )

    lines = [
        "LRC14 componentwise residue-ray law and pure-cone no-go",
        f"source_sha256={EXPECTED_SOURCE_SHA256}",
        (
            "analytic_ray_law=(z+L)delta_s(z+L)=zdelta_s(z);"
            "aligned_delta_s(La)=0;A_s(L-b)=-A_s(b)"
        ),
        (
            "analytic_cover_law=seven-tail cover forces sum(delta_s)>=0;"
            "strict for |I_s|>=1/105 because tail teeth are <=1/105"
        ),
        f"body={BODY};L={L};components={len(carrier)};forced_long={sum(long)}",
        f"component_grid={grid_components}",
        f"amplitude_types={type_table}",
        (
            f"ray_checks=recurrence:{recurrence_checks};"
            f"antipodal:{antipodal_checks};aligned:{aligned_checks};"
            f"identity_sha256={identity_digest.hexdigest()}"
        ),
        f"extremal_weights={EXTREMAL_WEIGHTS};values_by_type={values_by_type['extremal']}",
        f"equal_weights={EQUAL_WEIGHTS};values_by_type={values_by_type['equal']}",
        (
            f"reciprocal_weights={RECIPROCAL_WEIGHTS};"
            f"normalized={normalized_reciprocal};"
            f"values_by_type={values_by_type['reciprocal']};"
            f"minimum={ftext(min(certificates['reciprocal']))}"
        ),
        (
            "hostile_scope=necessary component cone only;not a cover;"
            "THM-2928 closes aligned k=4 independently"
        ),
        (
            f"all_body_roots={len(bodies)};"
            f"universal_reflected_ray_checks={universal_checks};"
            f"minimum_positive_amplitude={ftext(minimum_positive)};"
            f"universal_sha256={universal_digest.hexdigest()}"
        ),
        (
            "universal_formula=for every body speed e and carrier component I_s,"
            "A_s(e)=-e|I_s|/7 and A_s(L-e)=e|I_s|/7>0"
        ),
        (
            "conclusion=pure component-amplitude cone prunes zero aligned sectors;"
            "fixed reciprocal magnitudes/order/unit directions/overlap sidecars are required"
        ),
        f"semantic_sha256={semantic_sha256}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
