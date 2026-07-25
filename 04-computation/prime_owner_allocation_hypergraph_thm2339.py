#!/usr/bin/env python3
"""Optimization-safe exact controls for THM-2339.

This is not a knot-distance computation.  It verifies the finite N^4
word-metric hostile, the two-token allocation DP and Boolean Mobius energy,
and exhaustive integer instances of the deletion/collision identities.
Every load-bearing check uses ``require`` so ``python3`` and ``python3 -O``
execute the same verification.
"""

from __future__ import annotations

from itertools import product


Vector = tuple[int, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def add(*vectors: Vector) -> Vector:
    return tuple(sum(coordinates) for coordinates in zip(*vectors))


def sub(a: Vector, *vectors: Vector) -> Vector:
    answer = a
    for vector in vectors:
        answer = tuple(x - y for x, y in zip(answer, vector))
    return answer


def scale(k: int, vector: Vector) -> Vector:
    return tuple(k * x for x in vector)


def l1(vector: Vector) -> int:
    return sum(abs(x) for x in vector)


p = (1, 0, 0, 0)
q = (0, 1, 0, 0)
a = (0, 0, 1, 0)
b = (0, 0, 0, 1)
g = (1, 1, -1, -1)
v = (1, 1, -1, 0)
zero = (0, 0, 0, 0)


def norm0(vector: Vector) -> int:
    """Word norm with coordinate generators and g, all of cost one."""
    upper = l1(vector)
    return min(
        abs(k) + l1(sub(vector, scale(k, g)))
        for k in range(-upper, upper + 1)
    )


def norm1(vector: Vector) -> int:
    """Word norm with coordinate generators, g, and v, all of cost one."""
    upper = l1(vector)
    return min(
        abs(k) + abs(j) + l1(sub(vector, scale(k, g), scale(j, v)))
        for k in range(-upper, upper + 1)
        for j in range(-upper, upper + 1)
        if abs(k) + abs(j) <= upper
    )


def subset_vector(mask: int) -> Vector:
    answer = zero
    if mask & 1:
        answer = add(answer, p)
    if mask & 2:
        answer = add(answer, q)
    return answer


def subset_mobius(values: tuple[int, ...]) -> tuple[int, ...]:
    answer = []
    for mask in range(4):
        total = 0
        submask = mask
        while True:
            total += (-1) ** (
                mask.bit_count() - submask.bit_count()
            ) * values[submask]
            if submask == 0:
                break
            submask = (submask - 1) & mask
        answer.append(total)
    return tuple(answer)


def invert_mobius(coefficients: tuple[int, ...], mask: int) -> int:
    return sum(
        coefficients[submask]
        for submask in range(4)
        if submask & ~mask == 0
    )


def allocation_costs(norm) -> tuple[int, ...]:
    costs = []
    # Bits record the owner block of p and q: 0 for a, 1 for b.
    for owner_p, owner_q in product(range(2), repeat=2):
        mask_a = (1 if owner_p == 0 else 0) | (2 if owner_q == 0 else 0)
        mask_b = 3 ^ mask_a
        costs.append(
            norm(sub(a, subset_vector(mask_a)))
            + norm(sub(b, subset_vector(mask_b)))
        )
    return tuple(costs)


def subset_dp(w_a: tuple[int, ...], w_b: tuple[int, ...]) -> int:
    infinity = 10**9
    previous = (0, infinity, infinity, infinity)
    for weights in (w_a, w_b):
        current = []
        for mask in range(4):
            candidates = []
            assigned = mask
            while True:
                candidates.append(previous[mask ^ assigned] + weights[assigned])
                if assigned == 0:
                    break
                assigned = (assigned - 1) & mask
            current.append(min(candidates))
        previous = tuple(current)
    return previous[3]


def hypergraph_energy(
    coefficients: tuple[tuple[int, ...], tuple[int, ...]],
    owner_p: int,
    owner_q: int,
) -> int:
    owners = (owner_p, owner_q)
    energy = 0
    for mask in (1, 2, 3):
        token_indices = [i for i in range(2) if mask & (1 << i)]
        colours = {owners[i] for i in token_indices}
        if len(colours) == 1:
            colour = next(iter(colours))
            energy += coefficients[colour][mask]
    return energy


def check_hostile() -> None:
    vectors = {
        "a": a,
        "b": b,
        "a+b": add(a, b),
        "p": p,
        "q": q,
        "p+q": add(p, q),
        "p-a": sub(p, a),
        "q-a": sub(q, a),
        "p-b": sub(p, b),
        "q-b": sub(q, b),
        "p-a-b": sub(p, a, b),
        "q-a-b": sub(q, a, b),
        "p+q-a": sub(add(p, q), a),
        "p+q-b": sub(add(p, q), b),
        "p+q-a-b": sub(add(p, q), a, b),
    }
    expected = {
        "a": (1, 1),
        "b": (1, 1),
        "a+b": (2, 2),
        "p": (1, 1),
        "q": (1, 1),
        "p+q": (2, 2),
        "p-a": (2, 2),
        "q-a": (2, 2),
        "p-b": (2, 2),
        "q-b": (2, 2),
        "p-a-b": (2, 2),
        "q-a-b": (2, 2),
        "p+q-a": (2, 1),
        "p+q-b": (2, 2),
        "p+q-a-b": (1, 1),
    }
    observed = {
        name: (norm0(vector), norm1(vector))
        for name, vector in vectors.items()
    }
    require(observed == expected, f"word-length table mismatch: {observed}")

    print("hostile word-length table")
    for name in vectors:
        print(f"{name}: d0={observed[name][0]} d1={observed[name][1]}")

    for metric_name, norm, expected_costs, expected_lambda in (
        ("d0", norm0, (3, 4, 4, 3), 3),
        ("d1", norm1, (2, 4, 4, 3), 2),
    ):
        costs = allocation_costs(norm)
        require(costs == expected_costs, f"{metric_name} allocation mismatch")
        lift = min(costs)
        direct = norm(sub(add(a, b), add(p, q)))
        require(lift == expected_lambda, f"{metric_name} lift mismatch")

        w_a = tuple(
            norm(sub(a, subset_vector(mask))) - norm(a)
            for mask in range(4)
        )
        w_b = tuple(
            norm(sub(b, subset_vector(mask))) - norm(b)
            for mask in range(4)
        )
        mu_a = subset_mobius(w_a)
        mu_b = subset_mobius(w_b)

        require(
            all(invert_mobius(mu_a, mask) == w_a[mask] for mask in range(4)),
            f"{metric_name} Mobius inversion failed at a",
        )
        require(
            all(invert_mobius(mu_b, mask) == w_b[mask] for mask in range(4)),
            f"{metric_name} Mobius inversion failed at b",
        )

        baseline = norm(a) + norm(b)
        dp_energy = subset_dp(w_a, w_b)
        require(
            baseline + dp_energy == lift,
            f"{metric_name} subset DP failed",
        )

        for owner_p, owner_q in product(range(2), repeat=2):
            mask_a = (1 if owner_p == 0 else 0) | (
                2 if owner_q == 0 else 0
            )
            mask_b = 3 ^ mask_a
            direct_energy = w_a[mask_a] + w_b[mask_b]
            mobius_energy = hypergraph_energy(
                (mu_a, mu_b), owner_p, owner_q
            )
            require(
                direct_energy == mobius_energy,
                f"{metric_name} hypergraph energy failed",
            )

        print(
            f"{metric_name}: allocations={costs} "
            f"w_a={w_a} mu_a={mu_a} w_b={w_b} mu_b={mu_b} "
            f"Lambda={lift} Omega={lift-direct} dp_energy={dp_energy}"
        )


def check_deletion_collision_banks() -> None:
    deletion_cases = 0
    for delta_fine in range(5):
        for delta_coarse in range(delta_fine + 1):
            for translation in range(5):
                for bypass in range(5):
                    omega_j_fine = delta_fine + translation
                    omega_j_coarse = delta_coarse + translation
                    omega_u_fine = delta_fine + translation + bypass
                    omega_u_coarse = delta_coarse + translation + bypass
                    target_drop = omega_j_fine - omega_j_coarse
                    root_drop = omega_u_fine - omega_u_coarse
                    require(
                        target_drop == root_drop == delta_fine - delta_coarse,
                        "prime-owner deletion cocycle mismatch",
                    )
                    require(
                        omega_u_fine - omega_j_fine == bypass,
                        "prime-owner bypass gap mismatch",
                    )
                    deletion_cases += 1

    collision_cases = 0
    for translation in range(5):
        for bypass in range(5):
            sigma = translation + bypass
            for competitor_count in range(4):
                for excesses in product(range(5), repeat=competitor_count):
                    q_value = min((bypass, *excesses))
                    target_drop = sigma - q_value
                    require(
                        target_drop
                        == translation + (bypass - q_value),
                        "collision split mismatch",
                    )
                    require(
                        translation <= target_drop <= sigma,
                        "collision bounds mismatch",
                    )
                    require(
                        sigma - target_drop == q_value,
                        "root-target collision gap mismatch",
                    )
                    merged_owner = bypass == q_value
                    require(
                        merged_owner
                        == all(bypass <= excess for excess in excesses),
                        "merged-owner criterion mismatch",
                    )
                    for excess in excesses:
                        untouched_owner = excess == q_value
                        require(
                            untouched_owner
                            == (
                                excess <= bypass
                                and all(excess <= other for other in excesses)
                            ),
                            "untouched-owner criterion mismatch",
                        )
                    collision_cases += 1

    print("symbolic integer identity banks")
    print(f"deletion_cases={deletion_cases}")
    print(f"collision_cases={collision_cases}")


def main() -> None:
    print("THM-2339 exact allocation and prime-owner controls")
    check_hostile()
    check_deletion_collision_banks()
    print("all exact controls passed")


if __name__ == "__main__":
    main()
