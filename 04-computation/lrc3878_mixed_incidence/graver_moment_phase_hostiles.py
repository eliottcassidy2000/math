#!/usr/bin/env python3
"""Scratch hostiles for THM-3743 and phase-free packet moments.

The first control enumerates the complete minimum-l1 relation fibre of the
two t=11 AP endpoint-liar rows.  The second is an abstract information-loss
control: a half-translate preserves inversion, all component lengths, and the
entire Fourier power spectrum, but changes containment in a fixed packet.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
import json
from pathlib import Path


CHECKS = 0
EXPECTED_SEMANTIC_SHA256 = "9124631b862e70b6f2aeb2c150ee363d3bbdfc527bd058653d397d39bc16f966"


def require(ok: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


HERE = Path(__file__).resolve().parent
SPEC = spec_from_file_location("packet", HERE / "packet_covariance_probe.py")
require(SPEC is not None and SPEC.loader is not None, "load packet module")
packet = module_from_spec(SPEC)
SPEC.loader.exec_module(packet)


def positive_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in positive_compositions(total - first, parts - 1):
            yield (first,) + rest


def canonical_sign(row: tuple[int, ...]) -> tuple[int, ...]:
    first = next(x for x in row if x)
    return row if first > 0 else tuple(-x for x in row)


def relation_fibre(speeds: tuple[int, ...], max_norm: int = 3):
    """Return the first nonempty l1 shell, modulo overall sign."""
    n = len(speeds)
    for norm in range(1, max_norm + 1):
        relations = set()
        for support_size in range(1, min(n, norm) + 1):
            for support in combinations(range(n), support_size):
                for magnitudes in positive_compositions(norm, support_size):
                    for signs in product((-1, 1), repeat=support_size):
                        row = [0] * n
                        for i, magnitude, sign in zip(support, magnitudes, signs):
                            row[i] = sign * magnitude
                        if sum(a * v for a, v in zip(row, speeds)) == 0:
                            relations.add(canonical_sign(tuple(row)))
        if relations:
            return norm, tuple(sorted(relations))
    raise RuntimeError("empty bounded relation search")


def circle_shift(pieces: list[tuple[Q, Q]], delta: Q):
    out = []
    for a, b in pieces:
        a += delta
        b += delta
        while a >= 1:
            a -= 1
            b -= 1
        if b <= 1:
            out.append((a, b))
        else:
            out.extend(((a, Q(1)), (Q(0), b - 1)))
    return packet.merge(out)


def circle_reflect(pieces: list[tuple[Q, Q]]):
    return circle_shift([((1 - b) % 1, (1 - a) % 1)
                         for a, b in pieces if a != 0 and b != 1], Q(0)) \
        if all(a != 0 and b != 1 for a, b in pieces) else packet.merge(
            [(Q(0), 1 - a) if b == 1 else
             (1 - b, Q(1)) if a == 0 else
             (1 - b, 1 - a) for a, b in pieces]
        )


def main() -> None:
    body = tuple(range(1, 12))
    rows = (
        body + (33, 154),
        body + (77, 110),
    )
    fibres = [relation_fibre(row) for row in rows]
    require(fibres[0] == fibres[1], "identical minimum relation fibres")
    norm, fibre = fibres[0]
    require(norm == 3 and len(fibre) == 30, "minimum norm/count")
    support_histogram = {
        size: sum(sum(x != 0 for x in relation) == size for relation in fibre)
        for size in (1, 2, 3)
    }
    require(support_histogram == {1: 0, 2: 5, 3: 25},
            "doubling/Schur split")

    tail_relations = ((14, -3), (10, -7))
    for row, coefficients in zip(rows, tail_relations):
        require(coefficients[0] * row[-2] + coefficients[1] * row[-1] == 0,
                "tail relation")
        require(sum(map(abs, coefficients)) == 17, "tail relation norm")

    walls = tuple(Q(r, 14) for r in (3, 5, 9, 11))
    incidences = []
    for p, q in ((3, 14), (7, 10)):
        statuses = tuple(
            packet.circle_distance(11 * p * y) >= packet.DELTA
            and packet.circle_distance(11 * q * y) >= packet.DELTA
            for y in walls
        )
        incidences.append(statuses)
    require(incidences == [(False,) * 4, (True,) * 4],
            "opposite isolated-wall incidence")

    # Information-theoretic hostile to every phase-free autocorrelation
    # statistic.  Both sets are inversion symmetric about zero; translation
    # multiplies Fourier coefficients by phases and preserves all magnitudes.
    obstruction = packet.union_danger((1, 3))
    shifted = circle_shift(obstruction, Q(1, 2))
    require(circle_reflect(obstruction) == obstruction, "packet inversion")
    require(circle_reflect(shifted) == shifted, "half-translate inversion")
    require(packet.measure(shifted) == packet.measure(obstruction) == Q(5, 21),
            "equal packet masses")
    require(packet.circle_component_lengths(shifted)
            == packet.circle_component_lengths(obstruction),
            "equal component lengths")
    for t in range(1, 65):
        require(packet.grid_energy(shifted, t)
                == packet.grid_energy(obstruction, t),
                f"equal full grid-energy spectrum t={t}")
    contained_overlap = packet.intersection_measure(obstruction, obstruction)
    shifted_overlap = packet.intersection_measure(shifted, obstruction)
    require(contained_overlap == Q(5, 21), "identity containment")
    require(shifted_overlap < packet.measure(shifted),
            "half-translate breaks containment")

    semantic = {
        "rows": rows,
        "minimum_norm": norm,
        "minimum_fibre": fibre,
        "support_histogram": support_histogram,
        "tail_relations": tail_relations,
        "wall_incidence": incidences,
        "phase_hostile": {
            "packet": (1, 3),
            "mass": str(packet.measure(obstruction)),
            "component_lengths": tuple(map(str, packet.circle_component_lengths(obstruction))),
            "identity_overlap": str(contained_overlap),
            "half_shift_overlap": str(shifted_overlap),
            "checked_grid_energies": 64,
        },
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()
    if digest != EXPECTED_SEMANTIC_SHA256:
        raise RuntimeError("frozen semantic transcript")

    print("THM3878_GRAVER_AND_MOMENT_PHASE_HOSTILES_20260823")
    print("scope=information_loss_controls;not_LRC_counterexamples;LRC14=OPEN")
    print(f"t11_rows={rows}")
    print(f"minimum_l1={norm};common_fibre_size={len(fibre)};support_histogram={support_histogram}")
    print("minimum_fibre_uses_only_AP11_body=yes;tails_invisible_to_THM3743_minimizer=yes")
    print("tail_pair_relations=((14,-3),(10,-7));both_l1=17;THM778_pair_words_distinguish_ratios")
    print(f"isolated_wall_pair_safe={incidences}")
    print("graver_result=same_complete_minimum_relation_fibre_but_opposite_mixed_endpoint_incidence")
    print(f"autocorrelation_hostile_mass={packet.measure(obstruction)};components={packet.circle_component_lengths(obstruction)}")
    print(f"identity_overlap={contained_overlap};half_shift_overlap={shifted_overlap}")
    print("autocorrelation_result=same_inversion_components_and_all_Fourier_magnitudes_but_different_containment")
    print(f"semantic_sha256={digest}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
