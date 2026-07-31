#!/usr/bin/env python3
"""Exact companion for the THM-2790 depth-two central response theorem.

It imports the immutable secondary THM-2779 endpoint-gate companion,
reconstructs the two certified finite-field THM-2625 endpoint factors, and
tests the unnormalized current

    J_s(R) = P(R+s) Q(R),
    ((Z-1)J)_s(R) = J_s(R+s) - J_s(R)

for every nonzero s in F_13^2 and every R in F_13^2.

Digest basis: J and response use lexicographic s and then lexicographic R.
Modes use lexicographic s, the central cycles ordered by their lexicographically
least representative, and characters k=1,...,12.  Each field element is
serialized as one unsigned eight-byte big-endian integer.

Script:
  04-computation/lrc14_universal_central_response_thm2790.py
Output:
  05-knowledge/results/lrc14_universal_central_response_thm2790.out
"""

import hashlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_bockstein_symplectic_endpoint_gate_thm2779.py"
)
SPEC = importlib.util.spec_from_file_location("thm2779_endpoint_gate", SOURCE)
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)

P = 13
N = P * P
POINTS = GATE.POINTS
NONZERO = GATE.NONZERO

EXPECTED_FIELDS = (
    (
        352341050142921841,
        "fb738f7aca8dc6c08e0a6178e0fd9603f60deffec4bb26e16c7b9f4abae33d3f",
        "4bf698c20fe9099ade588112c44f8d0506a4a23b7261a2db3d2c744c240476da",
        "29916cc2fb6a66c939cda174be1623e57643e7d4a0707702e2b4f3a66189b4c2",
    ),
    (
        956354278959359281,
        "af3c9c215bf56b9ecb80f0f6d83ac3d746aed7cf0cf65100d94ce3af0fe5fe71",
        "7635f9584e36aa3227432125749393f35169ad7bfba1bc3220087625db44cb72",
        "e07cd3af43dd3969c2d4afb5706debc178375e04e76440f0a624d1aebe840ff4",
    ),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(values):
    payload = b"".join(int(value).to_bytes(8, "big") for value in values)
    return hashlib.sha256(payload).hexdigest()


def point_add(x, y):
    return ((x[0] + y[0]) % P, (x[1] + y[1]) % P)


def digit_pair(n):
    return n % P, (n // P) % P


def x_action(n):
    return 14 * n % N


def y_action(n):
    v, w = digit_pair(n)
    return (v + 1) % P + P * w


def z_action(n):
    return (n + P) % N


def compose(f, g, n):
    return f(g(n))


def verify_heisenberg_carrier():
    """Verify XY=ZYX and isolate the physical +1 carry."""
    carry_count = 0
    for n in range(N):
        v, w = digit_pair(n)
        require(digit_pair(x_action(n)) == (v, (w + v) % P), "X digits")
        require(compose(x_action, y_action, n)
                == compose(z_action, lambda m: compose(y_action, x_action, m), n),
                "XY=ZYX")
        physical_successor = (n + 1) % N
        if v == P - 1:
            require(physical_successor == compose(z_action, y_action, n),
                    "carry successor")
            carry_count += 1
        else:
            require(physical_successor == y_action(n), "noncarry successor")
    require(carry_count == P, "carry-wall census")
    return carry_count


def verify_edge_gate_hostile():
    """Show that separate THM-2779 edge gates do not imply this theorem.

    On a C_13 cycle, P_j=zeta^j and Q_j=zeta^-j have nonzero adjacent
    sums and differences, but P_(j+1)Q_j=zeta is constant.  Thus every
    product-current response and every nontrivial product-current mode
    vanishes.
    """
    checked = 0
    for field_index, (prime, root) in enumerate(GATE.BASE.MODS):
        zeta = pow(root, GATE.BASE.NN // P, prime)
        values = tuple(pow(zeta, j, prime) for j in range(P))
        inverse_values = tuple(pow(zeta, -j, prime) for j in range(P))
        require(
            all(
                (values[j] - values[(j + 1) % P]) % prime
                and (values[j] + values[(j + 1) % P]) % prime
                and (inverse_values[j] - inverse_values[(j + 1) % P]) % prime
                and (inverse_values[j] + inverse_values[(j + 1) % P]) % prime
                for j in range(P)
            ),
            f"edge hostile failed separate gate in field {field_index}",
        )
        current = tuple(
            values[(j + 1) % P] * inverse_values[j] % prime
            for j in range(P)
        )
        require(len(set(current)) == 1, "edge hostile current is not constant")
        require(
            all(
                sum(
                    current[j] * pow(zeta, -character * j, prime)
                    for j in range(P)
                )
                % prime
                == 0
                for character in range(1, P)
            ),
            "edge hostile acquired a nontrivial central mode",
        )
        checked += P
    return checked


def orbit_partition(step):
    """Return the 13 translation-by-step cycles in deterministic order."""
    unseen = set(POINTS)
    cycles = []
    while unseen:
        start = min(unseen)
        orbit = []
        point = start
        for _ in range(P):
            orbit.append(point)
            unseen.remove(point)
            point = point_add(point, step)
        require(point == start, "central orbit failed to close")
        cycles.append(tuple(orbit))
    require(len(cycles) == P, "central cycle census")
    return tuple(cycles)


def analyze_field(field, field_index):
    prime, left, right = field
    primitive_root = GATE.BASE.MODS[field_index][1]
    z13 = pow(primitive_root, GATE.BASE.NN // P, prime)
    z13_powers = tuple(pow(z13, exponent, prime) for exponent in range(P))
    j_bank = []
    response_bank = []
    central_mode_bank = []
    support_by_step = []
    active_cycles_by_step = []
    carry_wall_support_by_step = []
    central_mode_support_by_step = []
    carry_wall_mode_support_by_step = []

    for step in NONZERO:
        j = {}
        response = {}
        for origin in POINTS:
            value = left[point_add(origin, step)] * right[origin] % prime
            j[origin] = value
            j_bank.append(value)
        for origin in POINTS:
            shifted = point_add(origin, step)
            value = (j[shifted] - j[origin]) % prime
            response[origin] = value
            response_bank.append(value)

        support_by_step.append(sum(value != 0 for value in response.values()))
        carry_wall_support_by_step.append(
            sum(
                response[origin] != 0
                for origin in POINTS
                if GATE.det(step, origin) == P - 1
            )
        )
        active_cycles = 0
        central_mode_support = 0
        carry_wall_mode_support = 0
        for cycle in orbit_partition(step):
            require(
                sum(response[origin] for origin in cycle) % prime == 0,
                "central telescoping identity",
            )
            active_cycles += any(response[origin] != 0 for origin in cycle)
            is_carry_wall = GATE.det(step, cycle[0]) == P - 1
            require(
                all(
                    GATE.det(step, origin) == GATE.det(step, cycle[0])
                    for origin in cycle
                ),
                "determinant changed on central cycle",
            )
            for character in range(1, P):
                j_mode = sum(
                    j[origin] * z13_powers[-character * index % P]
                    for index, origin in enumerate(cycle)
                ) % prime
                response_mode = sum(
                    response[origin] * z13_powers[-character * index % P]
                    for index, origin in enumerate(cycle)
                ) % prime
                require(
                    response_mode
                    == (z13_powers[character] - 1) * j_mode % prime,
                    "central response multiplier identity",
                )
                central_mode_bank.append(j_mode)
                central_mode_support += j_mode != 0
                if is_carry_wall:
                    carry_wall_mode_support += j_mode != 0
        active_cycles_by_step.append(active_cycles)
        central_mode_support_by_step.append(central_mode_support)
        carry_wall_mode_support_by_step.append(carry_wall_mode_support)

    expected_spatial = len(NONZERO) * len(POINTS)
    expected_modes = len(NONZERO) * P * (P - 1)
    expected_carry_spatial = len(NONZERO) * P
    expected_carry_modes = len(NONZERO) * (P - 1)
    require(len(j_bank) == expected_spatial, "endpoint current universe drift")
    require(all(j_bank), "endpoint current support hole")
    require(len(response_bank) == expected_spatial, "response universe drift")
    require(all(response_bank), "central response support hole")
    require(
        all(value == len(POINTS) for value in support_by_step),
        "per-step central response support hole",
    )
    require(
        all(value == P for value in active_cycles_by_step),
        "inactive determinant cycle",
    )
    require(
        sum(carry_wall_support_by_step) == expected_carry_spatial
        and all(value == P for value in carry_wall_support_by_step),
        "carry-wall response support hole",
    )
    require(
        sum(central_mode_support_by_step) == expected_modes
        and all(value == P * (P - 1) for value in central_mode_support_by_step),
        "nontrivial central mode support hole",
    )
    require(
        sum(carry_wall_mode_support_by_step) == expected_carry_modes
        and all(value == P - 1 for value in carry_wall_mode_support_by_step),
        "carry-wall central mode support hole",
    )
    result = {
        "prime": prime,
        "j_support": sum(value != 0 for value in j_bank),
        "response_support": sum(value != 0 for value in response_bank),
        "response_support_min": min(support_by_step),
        "response_support_max": max(support_by_step),
        "full_response_steps": sum(
            support == len(POINTS) for support in support_by_step
        ),
        "active_cycles_min": min(active_cycles_by_step),
        "active_cycles_max": max(active_cycles_by_step),
        "carry_wall_support": sum(carry_wall_support_by_step),
        "carry_wall_support_min": min(carry_wall_support_by_step),
        "carry_wall_support_max": max(carry_wall_support_by_step),
        "central_mode_support": sum(central_mode_support_by_step),
        "central_mode_support_min": min(central_mode_support_by_step),
        "central_mode_support_max": max(central_mode_support_by_step),
        "carry_wall_mode_support": sum(carry_wall_mode_support_by_step),
        "carry_wall_mode_support_min": min(carry_wall_mode_support_by_step),
        "carry_wall_mode_support_max": max(carry_wall_mode_support_by_step),
        "j_digest": digest(j_bank),
        "response_digest": digest(response_bank),
        "central_mode_digest": digest(central_mode_bank),
    }
    expected = EXPECTED_FIELDS[field_index]
    require(
        (
            result["prime"],
            result["j_digest"],
            result["response_digest"],
            result["central_mode_digest"],
        )
        == expected,
        "field or deterministic digest drift",
    )
    return result


def main():
    carry_count = verify_heisenberg_carrier()
    hostile_edges = verify_edge_gate_hostile()
    fields = GATE.build_endpoint_factors()
    results = tuple(
        analyze_field(field, field_index)
        for field_index, field in enumerate(fields)
    )

    print("THM-2790 UNIVERSAL DEPTH-TWO CENTRAL RESPONSE")
    print("status=VERIFIED-EXACT dual-field coefficient-side companion")
    print(f"carrier_points={N} carry_wall_points={carry_count}")
    print("carrier_law=XY=ZYX; [X,Y]=Z under XYX^-1Y^-1")
    print(f"edge_gate_constant_product_hostile={hostile_edges}/{2 * P}")
    print("current_normalization=unnormalized_P_times_Q; full_current_differs_by_one_common_nonzero_scalar")
    print("digest_basis=J_response:lex_s_then_R;modes:lex_s_then_lex_min_cycle_then_k;unsigned_8_byte_big_endian")
    for result in results:
        print(
            f"p={result['prime']} "
            f"J_support={result['j_support']}/{len(NONZERO) * len(POINTS)} "
            f"response_support={result['response_support']}/"
            f"{len(NONZERO) * len(POINTS)} "
            f"per_s={result['response_support_min']}.."
            f"{result['response_support_max']} "
            f"full_s={result['full_response_steps']}/{len(NONZERO)} "
            f"active_cycles_per_s={result['active_cycles_min']}.."
            f"{result['active_cycles_max']} "
            f"carry_wall_support={result['carry_wall_support']}/"
            f"{len(NONZERO) * P} "
            f"per_s_carry={result['carry_wall_support_min']}.."
            f"{result['carry_wall_support_max']}"
        )
        print(
            f"p={result['prime']} "
            f"central_nonzero_modes={result['central_mode_support']}/"
            f"{len(NONZERO) * P * (P - 1)} "
            f"per_s_modes={result['central_mode_support_min']}.."
            f"{result['central_mode_support_max']} "
            f"carry_wall_modes={result['carry_wall_mode_support']}/"
            f"{len(NONZERO) * (P - 1)} "
            f"per_s_carry_modes={result['carry_wall_mode_support_min']}.."
            f"{result['carry_wall_mode_support_max']}"
        )
        print(
            f"p={result['prime']} J_sha256={result['j_digest']} "
            f"response_sha256={result['response_digest']} "
            f"central_modes_sha256={result['central_mode_digest']}"
        )
    print("scope=coefficient endpoint response to the depth-two central Z quotient; no physical high-sheet descent, Boolean carrier, or wing-current identification")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
