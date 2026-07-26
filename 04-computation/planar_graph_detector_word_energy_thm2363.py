#!/usr/bin/env python3
"""Exact finite-field companion for THM-2363.

All coefficient arithmetic is over the Gaussian integers.  Field elements
of F_169 are pairs a+b*t with t^2=2 over F_13.
"""

from __future__ import annotations

from typing import Dict, Iterable, Tuple


P = 13
N = P * P
INV2 = 7

Field = Tuple[int, int]
Gaussian = Tuple[int, int]
Key = Tuple[Field, Field]
Array = Dict[Key, Gaussian]

ZERO_F: Field = (0, 0)
ZERO_G: Gaussian = (0, 0)
ELEMENTS = tuple((a, b) for a in range(P) for b in range(P))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fadd(x: Field, y: Field) -> Field:
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def fneg(x: Field) -> Field:
    return ((-x[0]) % P, (-x[1]) % P)


def fsub(x: Field, y: Field) -> Field:
    return fadd(x, fneg(y))


def fmul(x: Field, y: Field) -> Field:
    # (a+b*t)(c+d*t)=(ac+2bd)+(ad+bc)t, t^2=2.
    return (
        (x[0] * y[0] + 2 * x[1] * y[1]) % P,
        (x[0] * y[1] + x[1] * y[0]) % P,
    )


def fhalf_square(x: Field) -> Field:
    square = fmul(x, x)
    return ((INV2 * square[0]) % P, (INV2 * square[1]) % P)


def gadd(x: Gaussian, y: Gaussian) -> Gaussian:
    return (x[0] + y[0], x[1] + y[1])


def gneg(x: Gaussian) -> Gaussian:
    return (-x[0], -x[1])


def gabs2(x: Gaussian) -> int:
    return x[0] * x[0] + x[1] * x[1]


def asum(values: Iterable[Gaussian]) -> Gaussian:
    total = ZERO_G
    for value in values:
        total = gadd(total, value)
    return total


def coefficient(array: Array, q: Field, z: Field) -> Gaussian:
    return array.get((q, z), ZERO_G)


def aggregate(array: Array, q: Field) -> Gaussian:
    return asum(coefficient(array, q, z) for z in ELEMENTS)


def graph_value(array: Array, c: Field, q: Field) -> Gaussian:
    return coefficient(array, q, fadd(fhalf_square(q), c))


def detector_terms(array: Array) -> Tuple[int, int, int]:
    row_energy = 0
    graph_sum_energy = 0
    for c in ELEMENTS:
        nonzero_values = [
            graph_value(array, c, q) for q in ELEMENTS if q != ZERO_F
        ]
        row_energy += sum(gabs2(value) for value in nonzero_values)
        graph_sum_energy += gabs2(asum(nonzero_values))
    return row_energy + graph_sum_energy, row_energy, graph_sum_energy


def joint_nonzero_energy(array: Array) -> int:
    return sum(
        gabs2(coefficient(array, q, z))
        for q in ELEMENTS
        if q != ZERO_F
        for z in ELEMENTS
    )


def aggregate_energy(array: Array) -> int:
    return sum(
        gabs2(aggregate(array, q)) for q in ELEMENTS if q != ZERO_F
    )


def in_mask(name: str, q: Field) -> bool:
    if name == "a":
        return q[0] != 0 and q[1] == 0
    if name == "b":
        return q[0] == 0 and q[1] != 0
    if name == "ab":
        return q[0] != 0 and q[1] != 0
    raise RuntimeError(f"unknown mask {name}")


def masked_energy(array: Array, name: str) -> int:
    return sum(
        gabs2(aggregate(array, q)) for q in ELEMENTS if in_mask(name, q)
    )


def check_array(array: Array, label: str) -> Tuple[int, int]:
    detector, graph_row, graph_sum = detector_terms(array)
    joint = joint_nonzero_energy(array)
    total = aggregate_energy(array)
    require(
        graph_row == joint,
        f"{label}: graph foliation stopped partitioning the joint array",
    )
    require(
        detector == joint + graph_sum,
        f"{label}: summed detector identity changed",
    )
    require(
        N * joint >= total,
        f"{label}: rowwise Cauchy domination failed",
    )
    for name in ("a", "b", "ab"):
        energy = masked_energy(array, name)
        require(
            total >= energy and N * detector >= energy,
            f"{label}: word-support domination failed for {name}",
        )
    return detector, total


def deterministic_control(seed: int) -> Array:
    """Make an exact sparse Gaussian-integer array without randomness."""
    array: Array = {}
    state = seed + 1
    for step in range(31):
        state = (1103515245 * state + 12345) % (2**31)
        qi = state % N
        state = (1103515245 * state + 12345) % (2**31)
        zi = state % N
        state = (1103515245 * state + 12345) % (2**31)
        real = state % 9 - 4
        state = (1103515245 * state + 12345) % (2**31)
        imag = state % 9 - 4
        q = ELEMENTS[qi]
        z = ELEMENTS[zi]
        old = array.get((q, z), ZERO_G)
        array[(q, z)] = gadd(old, (real + (step % 3 == 0), imag))
    return {key: value for key, value in array.items() if value != ZERO_G}


def sharp_array(q1: Field, q2: Field, alpha: Gaussian) -> Array:
    array: Array = {}
    for z in ELEMENTS:
        array[(q1, z)] = alpha
        array[(q2, z)] = gneg(alpha)
    return array


def main() -> None:
    # The graph coordinate c=z-q^2/2 is a bijection for every fixed q.
    for q in ELEMENTS:
        seen = {fsub(z, fhalf_square(q)) for z in ELEMENTS}
        require(len(seen) == N, "planar graph coordinate stopped being bijective")

    for seed in range(80):
        check_array(deterministic_control(seed), f"control-{seed}")

    alpha = (2, -3)
    sharp_pairs = {
        "a": ((1, 0), (2, 0)),
        "b": ((0, 1), (0, 2)),
        "ab": ((1, 1), (2, 1)),
    }
    sharp_detector = None
    sharp_energy = None
    for name, (q1, q2) in sharp_pairs.items():
        array = sharp_array(q1, q2, alpha)
        detector, total = check_array(array, f"sharp-{name}")
        energy = masked_energy(array, name)
        _, joint, graph_sum = detector_terms(array)
        require(
            graph_sum == 0
            and total == energy
            and N * detector == energy
            and joint == detector,
            f"sharp equality control changed for {name}",
        )
        sharp_detector = detector
        sharp_energy = energy

    # A nonzero target row may cancel after summing over the jet coordinate.
    # The graph detector still sees it through its positive row term.
    q = (1, 1)
    z1 = (0, 0)
    z2 = (1, 0)
    jet_cancel: Array = {(q, z1): alpha, (q, z2): gneg(alpha)}
    jet_detector, jet_total = check_array(jet_cancel, "jet-cancellation")
    require(
        jet_total == 0 and jet_detector > 0,
        "jet-resolved strictness control stopped separating the two energies",
    )

    require(
        sharp_detector is not None and sharp_energy is not None,
        "sharp controls did not run",
    )
    print("THM-2363 planar-graph detector / word-energy exact companion")
    print("field: F_169=F_13[t]/(t^2-2)")
    print("graph partition cells: 28561")
    print("deterministic Gaussian-integer controls: 80")
    print("word masks checked: a, b, ab")
    print(f"sharp detector: {sharp_detector}")
    print(f"sharp masked energy: {sharp_energy}")
    print("sharp ratio: D_graph/E_sigma=1/169")
    print(f"jet-cancellation aggregate energy: {jet_total}")
    print(f"jet-cancellation graph detector: {jet_detector}")
    print("VERDICT: D_graph >= E_sigma/169, and 1/169 is sharp")


if __name__ == "__main__":
    main()
