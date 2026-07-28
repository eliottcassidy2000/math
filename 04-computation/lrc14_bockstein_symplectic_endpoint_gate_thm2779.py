#!/usr/bin/env python3
"""Exact universal endpoint-square certificate for THM-2779.

This companion reconstructs the canonical THM-2625 endpoint factors in both
of that theorem's certified finite-field specializations.  It proves that
every nonzero endpoint increment changes the corresponding factor by neither
sign, and therefore every one of the

    2,184 * 169 * 169 = 62,377,224

normalized symplectic endpoint squares passes the full THM-2772
Segre--Pluecker gate.  The four corners are coefficients of one factorized
marked triangle.  They are not four Boolean allocation states on one physical
support.

Script:
  04-computation/lrc14_bockstein_symplectic_endpoint_gate_thm2779.py
Output:
  05-knowledge/results/lrc14_bockstein_symplectic_endpoint_gate_thm2779.out
"""

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "04-computation/lrc14_canonical_endpoint_current_thm2625.py"
SPEC = importlib.util.spec_from_file_location("thm2625_endpoint_control", SOURCE)
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

P = 13
POINTS = tuple((x, y) for x in range(P) for y in range(P))
NONZERO = tuple(point for point in POINTS if point != (0, 0))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def det(x, y):
    return (x[0] * y[1] - x[1] * y[0]) % P


def build_endpoint_factors():
    """Rebuild the separate THM-2625 left/right endpoint transforms."""
    q_intervals = BASE.build_set(BASE.PAT_QA, BASE.ZERO_ELL)
    q_starts = [start for start, _ in q_intervals]
    tabs = BASE.make_tabs(q_intervals, BASE.X, BASE.MODS)
    p_banks = [dict() for _ in BASE.MODS]
    q_banks = [dict() for _ in BASE.MODS]
    z13 = [pow(root, BASE.NN // P, prime) for prime, root in BASE.MODS]

    for address, ell in BASE.REPS.items():
        e_intervals = BASE.build_set(BASE.PAT_E, ell)
        a_values, _overlap = BASE.x_sweep(
            e_intervals,
            q_intervals,
            q_starts,
            BASE.X,
            BASE.MODS,
            tabs,
        )
        b_values = BASE.endpoint_sum(e_intervals, -BASE.Y, BASE.MODS)
        for field_index, (prime, _root) in enumerate(BASE.MODS):
            phase = pow(
                z13[field_index],
                BASE.M_DEEP * ell[BASE.TB] % P,
                prime,
            )
            p_banks[field_index][address] = (
                phase * a_values[field_index] % prime
            )
            q_banks[field_index][address] = b_values[field_index] % prime

    fields = []
    for field_index, (prime, _root) in enumerate(BASE.MODS):
        powers = tuple(
            pow(z13[field_index], exponent, prime) for exponent in range(P)
        )
        left = {}
        right = {}
        for point in POINTS:
            left_value = 0
            right_value = 0
            for address in POINTS:
                tau0, tau1 = BASE.TAU[address]
                pairing = (tau0 * point[0] + tau1 * point[1]) % P
                left_value = (
                    left_value
                    + p_banks[field_index][address]
                    * powers[-pairing % P]
                ) % prime
                right_value = (
                    right_value
                    + q_banks[field_index][address]
                    * powers[pairing]
                ) % prime
            left[point] = left_value
            right[point] = right_value
        require(all(left.values()), f"left support hole modulo {prime}")
        require(all(right.values()), f"right support hole modulo {prime}")
        fields.append((prime, left, right))
    return tuple(fields)


def hadamard_gate(field, left_origin, right_origin, source_step, target_step):
    prime, left, right = field
    p0 = left[left_origin]
    p1 = left[add(left_origin, source_step)]
    q0 = right[right_origin]
    q1 = right[add(right_origin, target_step)]
    values = (
        p0 * q0 % prime,
        p1 * q0 % prime,
        p0 * q1 % prime,
        p1 * q1 % prime,
    )
    d0 = sum(values) % prime
    d1 = (values[0] + values[1] - values[2] - values[3]) % prime
    d2 = (values[0] - values[1] + values[2] - values[3]) % prime
    d3 = (values[0] - values[1] - values[2] + values[3]) % prime
    require(all((p0, p1, q0, q1)), "zero endpoint factor in gate")
    require(all((d0, d1, d2, d3)), "zero Hadamard gate coordinate")
    require(d0 * d3 % prime == d1 * d2 % prime, "Pluecker identity")
    return (p0, p1, q0, q1), (d0, d1, d2, d3)


EXPECTED_WITNESSES = (
    (
        352341050142921841,
        (
            153528379476679734,
            205065739500452496,
            153829297094820499,
            338655797145668320,
        ),
        (
            340901553381135615,
            319769238334939305,
            60314403185657026,
            233858556407273320,
        ),
    ),
    (
        956354278959359281,
        (
            955153878005214586,
            313954411340872907,
            902111039187286496,
            3971292277140798,
        ),
        (
            511989854700512392,
            62552333514892978,
            86125904489159431,
            419745850041132014,
        ),
    ),
)


def main():
    fields = build_endpoint_factors()
    require(tuple(field[0] for field in fields)
            == tuple(row[0] for row in EXPECTED_WITNESSES),
            "embedding primes changed")

    normalized_frames = tuple(
        (source_step, target_step)
        for source_step in NONZERO
        for target_step in NONZERO
        if det(source_step, target_step) == 1
    )
    require(len(normalized_frames) == 2184, "normalized-frame census")

    source_edge_counts = []
    target_edge_counts = []
    for prime, left, right in fields:
        source_edges = 0
        target_edges = 0
        for step in NONZERO:
            for origin in POINTS:
                left0 = left[origin]
                left1 = left[add(origin, step)]
                right0 = right[origin]
                right1 = right[add(origin, step)]
                require((left0 - left1) % prime != 0,
                        f"left difference zero at {(prime, origin, step)}")
                require((left0 + left1) % prime != 0,
                        f"left sum zero at {(prime, origin, step)}")
                require((right0 - right1) % prime != 0,
                        f"right difference zero at {(prime, origin, step)}")
                require((right0 + right1) % prime != 0,
                        f"right sum zero at {(prime, origin, step)}")
                source_edges += 1
                target_edges += 1
        require(source_edges == 168 * 169, "source edge census")
        require(target_edges == 168 * 169, "target edge census")
        source_edge_counts.append(source_edges)
        target_edge_counts.append(target_edges)

    left_origin = (0, 0)
    right_origin = (0, 0)
    source_step = (0, 1)
    target_step = (12, 0)
    require(det(source_step, target_step) == 1, "witness orientation")
    witnesses = tuple(
        hadamard_gate(
            field,
            left_origin,
            right_origin,
            source_step,
            target_step,
        )
        for field in fields
    )
    require(
        tuple((field[0], factors, hadamard)
              for field, (factors, hadamard) in zip(fields, witnesses))
        == EXPECTED_WITNESSES,
        "dual-field witness changed",
    )

    square_count = len(normalized_frames) * len(POINTS) * len(POINTS)
    require(square_count == 62377224, "universal square count")

    print("THM-2779 UNIVERSAL BOCKSTEIN/SYMPLECTIC ENDPOINT GATE")
    print("status=VERIFIED-EXACT coefficient-side secondary companion")
    print("validity_gate=THM2625 dual Lucas/order embeddings PASS")
    print(f"normalized_frames={len(normalized_frames)}")
    print(f"source_edges_per_field={tuple(source_edge_counts)}")
    print(f"target_edges_per_field={tuple(target_edge_counts)}")
    print(f"normalized_endpoint_squares={square_count}")
    print(
        "witness="
        f"L={left_origin},R={right_origin},s={source_step},t={target_step},"
        f"det={det(source_step, target_step)}"
    )
    for field, (factors, hadamard) in zip(fields, witnesses):
        print(
            f"p={field[0]} factors={factors} "
            f"Hadamard={hadamard} Pluecker=0"
        )
    print("scope=all four corners are endpoint coefficients of one marked triangle, not Boolean allocation states on one physical support")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
