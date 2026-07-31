#!/usr/bin/env python3
"""Independent finite-field audit of the THM-2790 response candidate.

This intentionally reuses only the proved THM-2779 reconstruction of the
THM-2625 endpoint factors.  It does not use the candidate's orbit partition
or analysis routine.  Instead it chooses a canonical transverse vector for
each nonzero step and parametrizes every determinant fibre directly by

    R(w, delta) = w*s + delta*t,  det(s,t)=1.

It then checks the spatial forward difference and all twelve nontrivial
Fourier evaluations in both certified finite fields.

Script:
  04-computation/lrc14_universal_central_response_independent_audit_thm2790.py
Output:
  05-knowledge/results/lrc14_universal_central_response_independent_audit_thm2790.out
"""

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_bockstein_symplectic_endpoint_gate_thm2779.py"
)
SPEC = importlib.util.spec_from_file_location("thm2779_endpoint_gate_audit", SOURCE)
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)

P = 13
POINTS = GATE.POINTS
NONZERO = GATE.NONZERO


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def scale(a, x):
    return (a * x[0] % P, a * x[1] % P)


def transverse(step):
    a, b = step
    if a:
        return (0, pow(a, -1, P))
    return (-pow(b, -1, P) % P, 0)


def main():
    fields = GATE.build_endpoint_factors()
    require(
        tuple(field[0] for field in fields)
        == (352341050142921841, 956354278959359281),
        "certified field universe drift",
    )
    total_spatial = 0
    total_modes = 0
    total_carry_spatial = 0
    total_carry_modes = 0

    print("THM-2790 INDEPENDENT DETERMINANT-COORDINATE RESPONSE AUDIT")
    print("status=VERIFIED-EXACT independent dual-field companion")
    for field_index, (prime, left, right) in enumerate(fields):
        primitive_root = GATE.BASE.MODS[field_index][1]
        zeta = pow(primitive_root, GATE.BASE.NN // P, prime)
        field_spatial = 0
        field_modes = 0
        field_carry_spatial = 0
        field_carry_modes = 0

        for step in NONZERO:
            across = transverse(step)
            require(GATE.det(step, across) == 1, "bad transverse vector")
            seen = set()
            step_spatial = 0
            step_modes = 0
            step_carry_spatial = 0
            step_carry_modes = 0

            for delta in range(P):
                cycle = tuple(
                    add(scale(w, step), scale(delta, across))
                    for w in range(P)
                )
                require(len(set(cycle)) == P, "cycle collision")
                require(
                    all(GATE.det(step, point) == delta for point in cycle),
                    "determinant fibre mismatch",
                )
                seen.update(cycle)

                values = tuple(
                    left[add(point, step)] * right[point] % prime
                    for point in cycle
                )
                differences = tuple(
                    (values[(w + 1) % P] - values[w]) % prime
                    for w in range(P)
                )
                spatial = sum(value != 0 for value in differences)
                step_spatial += spatial
                if delta == P - 1:
                    step_carry_spatial += spatial

                for character in range(1, P):
                    j_mode = sum(
                        values[w] * pow(zeta, -character * w, prime)
                        for w in range(P)
                    ) % prime
                    response_mode = sum(
                        differences[w] * pow(zeta, -character * w, prime)
                        for w in range(P)
                    ) % prime
                    require(
                        response_mode
                        == (pow(zeta, character, prime) - 1) * j_mode % prime,
                        "Fourier multiplier sign mismatch",
                    )
                    require(j_mode != 0, "central character vanished")
                    step_modes += 1
                    if delta == P - 1:
                        step_carry_modes += 1

            require(seen == set(POINTS), "cycles do not partition endpoint plane")
            require(step_spatial == P * P, "spatial response hole")
            require(step_modes == P * (P - 1), "central mode hole")
            require(step_carry_spatial == P, "carry-wall spatial hole")
            require(step_carry_modes == P - 1, "carry-wall mode hole")
            field_spatial += step_spatial
            field_modes += step_modes
            field_carry_spatial += step_carry_spatial
            field_carry_modes += step_carry_modes

        require(field_spatial == len(NONZERO) * P * P, "field spatial census")
        require(field_modes == len(NONZERO) * P * (P - 1), "field mode census")
        require(
            field_carry_spatial == len(NONZERO) * P,
            "field carry spatial census",
        )
        require(
            field_carry_modes == len(NONZERO) * (P - 1),
            "field carry mode census",
        )
        total_spatial += field_spatial
        total_modes += field_modes
        total_carry_spatial += field_carry_spatial
        total_carry_modes += field_carry_modes
        print(
            f"p={prime} spatial={field_spatial} modes={field_modes} "
            f"carry_spatial={field_carry_spatial} "
            f"carry_modes={field_carry_modes}"
        )

    print(
        f"both_fields spatial={total_spatial} modes={total_modes} "
        f"carry_spatial={total_carry_spatial} carry_modes={total_carry_modes}"
    )
    print("INDEPENDENT CENTRAL-RESPONSE AUDIT PASSED")


if __name__ == "__main__":
    main()
