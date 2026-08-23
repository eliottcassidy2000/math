#!/usr/bin/env python3
"""Exact hostile companion for the THM-3825 inert cube decoder.

The checker verifies the finite carrier directly, locates the first losses
when one scope gate is removed, and exhibits and repairs the orientation
collision in the 78-channel label code.  ``require`` remains active under
``python -O``.
"""

from __future__ import annotations

from collections import defaultdict
from math import gcd, isqrt
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CAP = 356


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factor(n: int) -> dict[int, int]:
    require(n >= 1, "factor domain")
    factors: dict[int, int] = {}
    divisor = 2
    while divisor * divisor <= n:
        while n % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            n //= divisor
        divisor = 3 if divisor == 2 else divisor + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def inert_only(n: int) -> bool:
    return n > 1 and all(prime % 3 == 2 for prime in factor(n))


def admissible_shell(shell: int) -> bool:
    shell_factors = factor(shell)
    return (
        shell >= 3
        and all(prime % 3 == 2 for prime in shell_factors)
        and all(exponent <= 2 for exponent in shell_factors.values())
    )


def decode(value: int) -> tuple[int, int, int, int] | None:
    """Return (scale, shell, a, b), or reject outside the restricted image."""
    value_factors = factor(value)
    if 3 in value_factors:
        return None
    scale = 1
    shell = 1
    cofactor = 1
    for prime, exponent in value_factors.items():
        if prime % 3 == 2:
            carry, digit = divmod(exponent, 3)
            scale *= prime**carry
            shell *= prime**digit
        elif prime % 3 == 1:
            cofactor *= prime**exponent
        else:
            return None
    numerator = 4 * cofactor - shell * shell
    if numerator <= 0 or numerator % 3 != 0:
        return None
    difference_squared = numerator // 3
    difference = isqrt(difference_squared)
    if difference * difference != difference_squared:
        return None
    if not 0 < difference < shell or (shell - difference) % 2 != 0:
        return None
    first = (shell - difference) // 2
    second = (shell + difference) // 2
    if gcd(first, second) != 1:
        return None
    return scale, shell, first, second


def label_index(i: int, j: int) -> int:
    require(1 <= i < j <= 13, "label pair")
    return (j - 1) * (j - 2) // 2 + i - 1


def main() -> None:
    checks = 0
    selected: list[tuple[int, int]] = []
    admissible_shells = []
    primitive_values: dict[int, tuple[int, int]] = {}
    for shell in range(3, CAP + 1):
        if admissible_shell(shell):
            admissible_shells.append(shell)
        for first in range(1, (shell + 1) // 2):
            second = shell - first
            if gcd(first, second) != 1 or not admissible_shell(shell):
                continue
            selected.append((first, second))
            value = first**3 + second**3
            require(value not in primitive_values, "restricted primitive collision")
            primitive_values[value] = (first, second)
            require(decode(value) == (1, shell, first, second), "primitive decode")
            checks += 3

    require(len(admissible_shells) == 94, "shell census")
    require(len(selected) == 5_855, "pair census")
    require(len(primitive_values) == 5_855, "value census")

    inert_scales = [
        scale
        for scale in range(1, 101)
        if scale == 1 or all(prime % 3 == 2 for prime in factor(scale))
    ]
    for first, second in selected:
        shell = first + second
        base = first**3 + second**3
        for scale in inert_scales:
            value = scale**3 * base
            require(decode(value) == (scale, shell, first, second), "scaled decode")
            checks += 1

    label_indices = {
        label_index(i, j)
        for j in range(2, 14)
        for i in range(1, j)
    }
    require(label_indices == set(range(78)), "label-index interval")

    smallest_base, smallest_pair = min(
        (first**3 + second**3, (first, second)) for first, second in selected
    )
    require((smallest_base, smallest_pair) == (28, (1, 3)), "smallest carrier state")

    # Ramified prime 3: already the shell (1,2) is rejected.  If 3 is admitted
    # into the scale, its cubic valuation aliases exactly three label channels.
    ramified_shell_value = 1**3 + 2**3
    require(ramified_shell_value == 9 and decode(ramified_shell_value) is None,
            "ramified shell hostile")
    ramified_scale_value = 3**3 * smallest_base
    require(ramified_scale_value == 756 and decode(ramified_scale_value) is None,
            "ramified scale hostile")
    require(3**0 * ramified_scale_value == 3**3 * smallest_base,
            "ramified scale/label alias")

    # Seven is the least split prime.  Keeping every other carrier gate fixed,
    # the least split-scale state is the least primitive carrier value times 7^3.
    split_scale_value = 7**3 * smallest_base
    require(split_scale_value == 9_604 and decode(split_scale_value) is None,
            "least split-scale hostile")

    # Decoder completeness first fails at inert shell 8=2^3.  The stronger
    # singleton conclusion first fails (within the bounded exhaustive search)
    # at shell 64; exact exponent three first fails at shell 125.
    high_exponent_states = []
    for shell in range(3, CAP + 1):
        shell_factors = factor(shell)
        if not inert_only(shell) or max(shell_factors.values()) <= 2:
            continue
        for first in range(1, (shell + 1) // 2):
            second = shell - first
            if gcd(first, second) == 1:
                high_exponent_states.append((first**3 + second**3, first, second, shell))
    first_decoder_loss = min(high_exponent_states)
    require(first_decoder_loss == (152, 3, 5, 8), "least exponent-cap decoder loss")
    require(decode(first_decoder_loss[0]) is None, "exponent-cap hostile rejection")

    fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for first in range(1, CAP):
        for second in range(first + 1, CAP):
            fibres[first**3 + second**3].append((first, second))

    bad_collisions = []
    exact_three_collisions = []
    for value, first, second, shell in high_exponent_states:
        if len(fibres[value]) < 2:
            continue
        state = (value, first, second, shell, tuple(fibres[value]))
        bad_collisions.append(state)
        exponents = factor(shell).values()
        if max(exponents) == 3:
            exact_three_collisions.append(state)
    least_bad_collision = min(bad_collisions)
    least_exact_three_collision = min(exact_three_collisions)
    require(
        least_bad_collision
        == (65_728, 31, 33, 64, ((12, 40), (31, 33))),
        "least high-exponent singleton failure",
    )
    require(
        least_exact_three_collision
        == (515_375, 54, 71, 125, ((15, 80), (54, 71))),
        "least exact-exponent-three singleton failure",
    )

    ramified_shell_collisions = []
    split_scale_collisions = []
    for value, pairs in fibres.items():
        if len(pairs) < 2:
            continue
        for first, second in pairs:
            scale = gcd(first, second)
            shell = (first + second) // scale
            shell_factors = factor(shell)
            if (
                3 in shell_factors
                and all(
                    prime == 3 or (prime % 3 == 2 and exponent <= 2)
                    for prime, exponent in shell_factors.items()
                )
                and (scale == 1 or all(prime % 3 == 2 for prime in factor(scale)))
            ):
                ramified_shell_collisions.append((value, (first, second), tuple(pairs)))
            if admissible_shell(shell) and any(prime % 3 == 1 for prime in factor(scale)):
                split_scale_collisions.append(
                    (value, (first, second), scale, shell, tuple(pairs))
                )
    least_ramified_shell_collision = min(ramified_shell_collisions)
    least_split_scale_collision = min(split_scale_collisions)
    require(
        least_ramified_shell_collision
        == (4_104, (2, 16), ((2, 16), (9, 15))),
        "least ramified-shell singleton failure",
    )
    require(
        least_split_scale_collision
        == (41_343_640, (86, 344), 86, 5, ((86, 344), (197, 323))),
        "least split-containing-scale singleton failure",
    )

    # The proposed 78-channel address remembers the unordered coordinate pair
    # but not which coordinate receives the smaller primitive speed.  These two
    # labelled states have distinct projective covectors and the same address.
    current_address = 3 ** label_index(1, 2) * smallest_base
    covector_smaller_at_1 = (3, -1)
    covector_smaller_at_2 = (1, -3)
    require(current_address == 28, "orientation hostile address")
    require(
        covector_smaller_at_1 != covector_smaller_at_2
        and covector_smaller_at_1 != tuple(-entry for entry in covector_smaller_at_2),
        "distinct covector sign classes",
    )

    # A minimal scalar repair is lambda=2*kappa+epsilon.  This is still a fixed
    # coordinate-order gauge, but it retains the assignment orientation.
    repaired_checks = 0
    channel_powers = [3**channel for channel in range(156)]
    for base in primitive_values:
        require(base % 3 != 0, "label base must be 3-free")
        for unordered_label in range(78):
            for orientation in range(2):
                channel = 2 * unordered_label + orientation
                channel_power = channel_powers[channel]
                address = channel_power * base
                require(
                    (channel // 2, channel % 2, address // channel_power)
                    == (unordered_label, orientation, base),
                    "oriented label decoder",
                )
                require(address % (3 * channel_power) != 0, "exact oriented channel valuation")
                repaired_checks += 1

    gauge_before = 3 ** label_index(1, 2) * smallest_base
    gauge_after = 3 ** label_index(2, 3) * smallest_base
    require((gauge_before, gauge_after) == (28, 252), "permutation gauge hostile")

    print("INERT CUBE DECODER HOSTILE AUDIT")
    print(f"universe=coprime 1<=a<b, a+b<={CAP}")
    print(f"admissible_shells={len(admissible_shells)}; selected_pairs={len(selected)}")
    print(f"primitive_values={len(primitive_values)}; inert_scales_through_100={len(inert_scales)}")
    print(f"least_carrier_value={smallest_base} from pair={smallest_pair}")
    print("ramified_shell_loss=9=(1^3+2^3); decoder rejects v3>0")
    print("ramified_scale_label_alias=756: (label=0,g=3) equals (label=3,g=1) on pair (1,3)")
    print("least_ramified_shell_singleton_failure=4104=(2,16)=(9,15); shell=9=3^2")
    print("least_split_scale_loss=9604=7^3*(1^3+3^3); decoder rejects")
    print("least_split_scale_singleton_failure=41343640=(86,344)=(197,323); g=86=2*43, shell=5")
    print("least_exponent_cap_decoder_loss=152=3^3+5^3; primitive shell=8=2^3")
    print("least_high_exponent_singleton_failure=65728=(12,40)=(31,33); shell=64=2^6")
    print("least_exact_exponent_three_failure=515375=(15,80)=(54,71); shell=125=5^3")
    print("orientation_collision=address 28 for labelled speeds (1,3) and (3,1) on labels {1,2}")
    print("covector_sign_classes=(3,-1) versus (1,-3)")
    print(f"unoriented_labelled_count={78 * len(selected)}")
    print(f"oriented_labelled_count={156 * len(selected)}; repaired_channel=2*kappa+epsilon")
    print("coordinate_gauge_example=28 on {1,2} maps to 252 on {2,3} under a relabelling")
    print(f"active_checks={checks + 2 * repaired_checks + 19}")
    print("RESULT PASS")


if __name__ == "__main__":
    main()
