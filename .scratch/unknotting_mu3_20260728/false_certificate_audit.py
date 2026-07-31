#!/usr/bin/env python3
"""Independent finite audit of the repaired engine's FALSE lane.

This does not replace a proof of the Gordon--Litherland implementation.
It checks the algebra behind the finite-abelian-group decisions against
an independent Smith-normal-form computation, verifies the classic 7_4
Lickorish obstruction exactly, and probes PD relabel/order invariance.
"""

import importlib.util
import random
import sys
from math import gcd
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
ENGINE_PATH = (
    Path(sys.argv[1]).resolve()
    if len(sys.argv) > 1
    else REPO_ROOT / "04-computation" / "unknot1_decider.py"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_engine():
    spec = importlib.util.spec_from_file_location("unknot1_repaired", ENGINE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def independent_cyclic_snf(matrix):
    snf = smith_normal_form(sp.Matrix(matrix), domain=sp.ZZ)
    diagonal = [
        abs(int(snf[i, i]))
        for i in range(min(snf.rows, snf.cols))
        if snf[i, i] != 0
    ]
    require(len(diagonal) == len(matrix), "matrix was unexpectedly singular")
    return sum(value > 1 for value in diagonal) <= 1, diagonal


def random_symmetric_matrix(rng, dimension):
    matrix = [[0] * dimension for _ in range(dimension)]
    for i in range(dimension):
        for j in range(i, dimension):
            value = rng.randint(-5, 5)
            matrix[i][j] = matrix[j][i] = value
    return matrix


def relabel_and_reorder(rng, pd):
    labels = list(range(1, 2 * len(pd) + 1))
    shuffled = list(labels)
    rng.shuffle(shuffled)
    rename = dict(zip(labels, shuffled))
    transformed = [[rename[value] for value in row] for row in pd]
    rng.shuffle(transformed)
    return transformed


def main():
    engine = load_engine()

    # Named exact controls used by the FALSE lane.
    trefoil = engine.goeritz_invariants(engine.parse_pd(engine.TREFOIL))
    figure8 = engine.goeritz_invariants(engine.parse_pd(engine.FIG8))
    knot74 = engine.goeritz_invariants(engine.parse_pd(engine.K7_4))
    require((trefoil["det"], trefoil["sigma"]) == (3, -2),
            "trefoil calibration changed")
    require((figure8["det"], figure8["sigma"]) == (5, 0),
            "figure-eight calibration changed")
    require(knot74["det"] == 15 and abs(knot74["sigma"]) == 2,
            "7_4 calibration changed")

    lk74 = engine.linking_form_analysis(knot74["G"])
    units15 = [unit for unit in range(1, 15) if gcd(unit, 15) == 1]
    orbit74 = {(lk74["q"] * unit * unit) % 15 for unit in units15}
    require(lk74["q"] == 11, "unexpected chosen 7_4 generator form")
    require(orbit74 == {11, 14}, "wrong 7_4 generator square orbit")
    require(orbit74.isdisjoint({2, 13}) and lk74["obstructed"],
            "7_4 should miss both +2 and -2 modulo 15")
    require(engine.decide(engine.K7_4)["verdict"] == engine.VERDICT_FALSE,
            "7_4 no longer receives its exact Lickorish false certificate")

    # The engine intentionally searches only basis vectors and two-coordinate
    # +/- combinations for a cyclic generator.  This symmetric presentation
    # has coker(G)=Z/105; the coordinate images are 15,21,35, no supported
    # pair is a generator, while (1,1,1) is.  The implementation must fail
    # safe (no obstruction), not manufacture a FALSE verdict.
    incomplete_matrix = [
        [998137, -666540, -27741],
        [-666540, 445105, 18525],
        [-27741, 18525, 771],
    ]
    incomplete = engine.linking_form_analysis(incomplete_matrix)
    incomplete_cyclic, incomplete_snf = independent_cyclic_snf(
        incomplete_matrix)
    require(incomplete_cyclic and incomplete_snf == [1, 1, 105],
            "hostile cyclic presentation changed")
    require(incomplete["q"] is None and not incomplete["obstructed"],
            "incomplete generator search did not fail safe")

    # Exact PD label/crossing-order invariance checks.
    rng = random.Random(20260728)
    named = [
        (engine.TREFOIL, (3, -2)),
        (engine.FIG8, (5, 0)),
        (engine.K7_4, (15, knot74["sigma"])),
        (engine.EXAMPLE11, (43, -2)),
    ]
    relabel_checks = 0
    for pd, expected in named:
        for _ in range(100):
            variant = relabel_and_reorder(rng, pd)
            invariant = engine.goeritz_invariants(engine.parse_pd(variant))
            require((invariant["det"], invariant["sigma"]) == expected,
                    "PD relabel/order invariance failed")
            relabel_checks += 1

    # Compare cyclicity and the final square-orbit decision against an
    # independent SNF computation on arbitrary symmetric presentations.
    matrix_checks = 0
    cyclic_checks = 0
    noncyclic_checks = 0
    safe_incomplete_generator_searches = 0
    while matrix_checks < 1000:
        dimension = rng.randint(1, 5)
        matrix = random_symmetric_matrix(rng, dimension)
        determinant = abs(int(sp.Matrix(matrix).det()))
        if determinant <= 1 or determinant % 2 == 0 \
                or determinant > engine.LICKORISH_D_CAP:
            continue
        expected_cyclic, diagonal = independent_cyclic_snf(matrix)
        analysis = engine.linking_form_analysis(matrix)
        require(analysis["applicable"], "odd finite presentation was skipped")
        require(analysis["cyclic"] == expected_cyclic,
                "adjugate cyclicity disagrees with SNF %s" % diagonal)
        if expected_cyclic:
            cyclic_checks += 1
            if analysis["q"] is None:
                safe_incomplete_generator_searches += 1
                require(not analysis["obstructed"],
                        "missing generator candidate must never obstruct")
            else:
                d = analysis["d"]
                targets = {2 % d, (-2) % d}
                orbit = {
                    (analysis["q"] * unit * unit) % d
                    for unit in range(1, d) if gcd(unit, d) == 1
                }
                require(analysis["obstructed"] == orbit.isdisjoint(targets),
                        "reported Lickorish result disagrees with full orbit")
        else:
            noncyclic_checks += 1
            require(analysis["obstructed"],
                    "a noncyclic surgery homology group must obstruct u=1")
        matrix_checks += 1

    print("named controls = trefoil, figure-eight, 7_4")
    print("7_4 form orbit =", sorted(orbit74), "; targets =", [2, 13])
    print("explicit incomplete-generator SNF =", incomplete_snf,
          "; obstructed =", incomplete["obstructed"])
    print("PD relabel/order checks =", relabel_checks)
    print("independent SNF checks =", matrix_checks)
    print("cyclic presentations =", cyclic_checks)
    print("noncyclic presentations =", noncyclic_checks)
    print("safe incomplete generator searches =",
          safe_incomplete_generator_searches)
    print("FALSE-lane finite audit = PASS")
    print("scope = algebra/implementation probes, not a replacement for the "
          "Gordon--Litherland and Montesinos--Lickorish theorems")


if __name__ == "__main__":
    main()
