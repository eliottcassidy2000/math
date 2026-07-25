#!/usr/bin/env python3
"""Exact finite referee for THM-2231.

The paper proof identifies infinite integral relation-carry paths with the
kernel of the relation map on the q-adic completion.  This companion checks:

* finite-cylinder carry/path equivalence for prime and composite radices;
* the sharp mixed-sign carry bounds on every accepted cylinder;
* the primitive thirteen-speed base-13 state-plateau family; and
* the exact one-step decrease of the remaining radix-length clock.

Only exact integer arithmetic is used.  Every check raises explicitly, so
``python`` and ``python -O`` execute the same verification.
"""

from functools import reduce
from hashlib import sha256
from itertools import product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def digit_layers(values, q, depth):
    """Return the first ``depth`` coordinatewise base-q digit vectors."""
    return tuple(
        tuple((value // (q**j)) % q for value in values)
        for j in range(depth)
    )


def carry_path(digits, relation, q):
    """Return the integral carries, or ``None`` if one transition fails."""
    carry = 0
    carries = [carry]
    for layer in digits:
        numerator = carry + dot(relation, layer)
        if numerator % q:
            return None
        carry = numerator // q
        carries.append(carry)
    return tuple(carries)


def radix_length(value, q):
    """Number of base-q digits of a nonnegative integer, with length(0)=0."""
    length = 0
    while value:
        value //= q
        length += 1
    return length


def quotient_state(values, relation, q, level):
    modulus = q**level
    quotients = tuple(value // modulus for value in values)
    remainders = tuple(value % modulus for value in values)
    numerator = dot(relation, remainders)
    require(
        numerator % modulus == 0,
        f"nonintegral carry at level {level}",
    )
    carry = numerator // modulus
    digits = tuple(quotient % q for quotient in quotients)
    owners = tuple(i + 1 for i, quotient in enumerate(quotients) if quotient)
    height = radix_length(max(quotients, default=0), q)
    return carry, owners, digits, quotients, height


def audit_finite_cylinders():
    # The q=4 case is deliberately composite.
    cases = (
        (2, 5, (1, 1, -1)),
        (3, 4, (2, -3, 1)),
        (4, 3, (1, -2, 1)),
    )
    vectors_checked = 0
    accepted = []
    digest = sha256()

    for q, depth, relation in cases:
        modulus = q**depth
        positive_mass = sum(c for c in relation if c > 0)
        negative_mass = sum(-c for c in relation if c < 0)
        require(positive_mass > 0 and negative_mass > 0, "mixed-sign control")
        accepted_here = 0

        for values in product(range(modulus), repeat=len(relation)):
            vectors_checked += 1
            digits = digit_layers(values, q, depth)
            carries = carry_path(digits, relation, q)
            direct = dot(relation, values) % modulus == 0
            require(
                (carries is not None) == direct,
                f"cylinder/path mismatch q={q}, depth={depth}, values={values}",
            )
            if carries is None:
                continue

            accepted_here += 1
            prefix = [0] * len(relation)
            for j, layer in enumerate(digits, start=1):
                place = q ** (j - 1)
                prefix = [
                    prefix[i] + place * layer[i]
                    for i in range(len(relation))
                ]
                require(
                    dot(relation, prefix) == (q**j) * carries[j],
                    "telescoping identity failed",
                )
                require(
                    -negative_mass < carries[j] < positive_mass,
                    "sharp mixed-sign carry bound failed",
                )

            digest.update(
                (
                    f"{q}:{depth}:{relation}:{values}:{carries[-1]}\n"
                ).encode("ascii")
            )

        accepted.append((q, depth, relation, accepted_here))

    return vectors_checked, tuple(accepted), digest.hexdigest()


def audit_primitive_plateau():
    relation = (1,) + (0,) * 10 + (1, -1)
    require(len(relation) == 13, "relation arity drift")
    require(sum(abs(c) for c in relation) == 3, "relation norm drift")
    require(max(abs(c) for c in relation) == 1, "relation height drift")

    carry_owner_state_count = 0
    full_state_count = 0
    clock_transitions = 0
    digest = sha256()

    for exponent in range(2, 129):
        power = 13**exponent
        values = tuple(range(1, 12)) + (power, power + 1)
        require(len(values) == 13, "row arity drift")
        require(len(set(values)) == 13, "row distinctness failed")
        require(reduce(gcd, values) == 1, "row primitivity failed")
        require(dot(relation, values) == 0, "relation identity failed")

        zero_digits = (0,) * 13
        expected_owners = (12, 13)
        for level in range(1, exponent):
            carry, owners, digits, _, _ = quotient_state(
                values, relation, 13, level
            )
            require(
                (carry, owners, digits)
                == (0, expected_owners, zero_digits),
                f"full plateau drift at exponent={exponent}, level={level}",
            )
            full_state_count += 1

        carry, owners, digits, _, _ = quotient_state(
            values, relation, 13, exponent
        )
        expected_top_digits = (0,) * 11 + (1, 1)
        require(
            (carry, owners, digits)
            == (0, expected_owners, expected_top_digits),
            f"top digit drift at exponent={exponent}",
        )
        carry_owner_state_count += exponent

        terminal = quotient_state(values, relation, 13, exponent + 1)
        require(
            terminal[0] == 0 and terminal[1] == () and terminal[4] == 0,
            f"terminal state drift at exponent={exponent}",
        )

        for level in range(0, exponent + 2):
            current = quotient_state(values, relation, 13, level)
            following = quotient_state(values, relation, 13, level + 1)
            require(
                following[4] == max(current[4] - 1, 0),
                f"radix clock drift at exponent={exponent}, level={level}",
            )
            require(
                (current[1] == ()) == (current[4] == 0),
                "owner/clock termination mismatch",
            )
            clock_transitions += 1

        digest.update(
            (
                f"{exponent}:{values[-2]}:{values[-1]}:"
                f"{terminal[0]}:{terminal[4]}\n"
            ).encode("ascii")
        )

    # The carry/owner count includes levels 1,...,N.  The full-state count
    # includes levels 1,...,N-1.
    require(carry_owner_state_count == sum(range(2, 129)), "owner count drift")
    require(full_state_count == sum(range(1, 128)), "full-state count drift")
    return (
        carry_owner_state_count,
        full_state_count,
        clock_transitions,
        digest.hexdigest(),
    )


def audit_generic_radix_clocks():
    transitions = 0
    for q in range(2, 18):
        for seed in range(257):
            values = (seed, seed * q + 1, seed * seed + 3)
            level = 0
            while True:
                quotients = tuple(value // (q**level) for value in values)
                next_quotients = tuple(
                    value // (q ** (level + 1)) for value in values
                )
                height = radix_length(max(quotients), q)
                next_height = radix_length(max(next_quotients), q)
                require(
                    next_height == max(height - 1, 0),
                    f"generic clock drift q={q}, seed={seed}, level={level}",
                )
                transitions += 1
                if height == 0:
                    break
                level += 1
    return transitions


def main():
    vectors_checked, accepted, cylinder_digest = audit_finite_cylinders()
    (
        carry_owner_states,
        full_states,
        plateau_clock_transitions,
        plateau_digest,
    ) = audit_primitive_plateau()
    generic_clock_transitions = audit_generic_radix_clocks()

    require(vectors_checked == 826_353, "finite-cylinder universe drift")
    require(
        accepted
        == (
            (2, 5, (1, 1, -1), 1_024),
            (3, 4, (2, -3, 1), 6_561),
            (4, 3, (1, -2, 1), 4_096),
        ),
        "accepted-cylinder census drift",
    )
    require(
        cylinder_digest
        == "08e26f7d17ac9eaad9ee451fb2d42041c21c1b83cdfc0907b5f0dff7c25287bf",
        "accepted-cylinder digest drift",
    )
    require(carry_owner_states == 8_255, "carry/owner plateau count drift")
    require(full_states == 8_128, "full-state plateau count drift")
    require(
        plateau_clock_transitions == 8_509,
        "plateau clock-transition count drift",
    )
    require(
        plateau_digest
        == "3566f07b001b863dcc5c65480bdeba0eda00d295aab4995efbb2084d4e734973",
        "plateau-family digest drift",
    )
    require(
        generic_clock_transitions == 26_751,
        "generic clock-transition count drift",
    )

    print(f"finite_cylinder_vectors_checked={vectors_checked}")
    print(f"accepted_cylinder_counts={accepted}")
    print(f"accepted_cylinder_digest={cylinder_digest}")
    print("composite_radix_control=q4_depth3_PASS")
    print("qadic_kernel_cylinder_equivalence=PASS")
    print("sharp_mixed_sign_carry_bounds=PASS")
    print("primitive_lrc_family=V_N=(1,...,11,13^N,13^N+1)")
    print("primitive_lrc_exponents_checked=2..128")
    print(f"carry_owner_plateau_states_checked={carry_owner_states}")
    print(f"carry_owner_digit_plateau_states_checked={full_states}")
    print(f"plateau_clock_transitions_checked={plateau_clock_transitions}")
    print(f"plateau_family_digest={plateau_digest}")
    print(f"generic_radix_clock_transitions_checked={generic_clock_transitions}")
    print("exact_radix_clock=PASS")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
