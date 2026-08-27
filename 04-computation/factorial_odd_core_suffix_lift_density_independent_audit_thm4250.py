#!/usr/bin/env python3
"""Independent hostile audit of THM-4250 odd-core binary suffix closure.

Derived without inspecting odd_core_suffix_closure_probe.py.  It uses two
representations: explicit occupied-bit sets and a low-to-high schoolbook
multiplication automaton.  Every gate remains live under python -O.
"""

from hashlib import sha256
from math import ceil, log2


MAX_ODD_CORE = 63
DIRECT_MAX_LENGTH = 13
LONG_SEQUENCE_LENGTH = 30


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def positions(integer, length):
    return frozenset(index for index in range(length) if (integer >> index) & 1)


def pair_indices(odd_core):
    require(odd_core >= 3 and odd_core % 2 == 1, odd_core)
    return tuple(range(1, (odd_core - 1) // 2 + 1))


def explicit_closed(odd_core, length, residue):
    require(1 <= residue < (1 << length) and residue % 2 == 1, (length, residue))
    mask = (1 << length) - 1
    for left in pair_indices(odd_core):
        right = odd_core - left
        left_word = positions((left * residue) & mask, length)
        right_word = positions((right * residue) & mask, length)
        if left_word.isdisjoint(right_word):
            return False
    return True


def initial_automaton_state(odd_core):
    # carries[m-1] is the multiplication carry for multiplier m.
    return (tuple(0 for _ in range(1, odd_core)), 0)


def automaton_step(odd_core, state, input_bit):
    carries, hit_mask = state
    outputs = []
    next_carries = []
    for multiplier, carry in enumerate(carries, start=1):
        require(0 <= carry <= multiplier - 1, (multiplier, carry))
        total = carry + multiplier * input_bit
        outputs.append(total & 1)
        next_carry = total >> 1
        require(0 <= next_carry <= multiplier - 1, (multiplier, next_carry))
        next_carries.append(next_carry)
    for left in pair_indices(odd_core):
        right = odd_core - left
        if outputs[left - 1] and outputs[right - 1]:
            hit_mask |= 1 << (left - 1)
    return (tuple(next_carries), hit_mask)


def automaton_closed(odd_core, state):
    _, hit_mask = state
    target = sum(1 << (left - 1) for left in pair_indices(odd_core))
    return (hit_mask & target) == target


def direct_counts(odd_core, max_length):
    return tuple(
        sum(
            explicit_closed(odd_core, length, residue)
            for residue in range(1, 1 << length, 2)
        )
        for length in range(1, max_length + 1)
    )


def automaton_counts(odd_core, max_length):
    # Odd residues have low input bit one. Aggregate equal finite states.
    state_counts = {
        automaton_step(odd_core, initial_automaton_state(odd_core), 1): 1
    }
    answer = []
    for _length in range(1, max_length + 1):
        answer.append(sum(
            multiplicity
            for state, multiplicity in state_counts.items()
            if automaton_closed(odd_core, state)
        ))
        next_counts = {}
        for state, multiplicity in state_counts.items():
            for bit in (0, 1):
                child = automaton_step(odd_core, state, bit)
                next_counts[child] = next_counts.get(child, 0) + multiplicity
        state_counts = next_counts
    return tuple(answer)


def unresolved_pairs_and_carries(odd_core, length, residue):
    mask = (1 << length) - 1
    carry_parities = {
        multiplier: ((multiplier * residue) >> length) & 1
        for multiplier in range(1, odd_core)
    }
    unresolved = tuple(
        (left, odd_core - left)
        for left in pair_indices(odd_core)
        if ((left * residue) & mask) & (((odd_core - left) * residue) & mask) == 0
    )
    return unresolved, carry_parities


def predicted_closed_lifts(odd_core, length, residue):
    unresolved, carry_parities = unresolved_pairs_and_carries(
        odd_core, length, residue
    )
    if not unresolved:
        return (0, 1)
    required_lifts = []
    for left, right in unresolved:
        even = left if left % 2 == 0 else right
        odd = right if left % 2 == 0 else left
        if carry_parities[even] != 1:
            return ()
        required_lifts.append(1 ^ carry_parities[odd])
    if len(set(required_lifts)) != 1:
        return ()
    return (required_lifts[0],)


def two_lift_audit():
    digest = sha256()
    cells = 0
    child_histogram = {0: 0, 1: 0, 2: 0}
    for odd_core in range(3, MAX_ODD_CORE + 1, 2):
        for length in range(1, DIRECT_MAX_LENGTH + 1):
            modulus = 1 << length
            for residue in range(1, modulus, 2):
                predicted = predicted_closed_lifts(odd_core, length, residue)
                actual = tuple(
                    lift
                    for lift in (0, 1)
                    if explicit_closed(
                        odd_core, length + 1, residue + lift * modulus
                    )
                )
                require(predicted == actual, (
                    odd_core, length, residue, predicted, actual
                ))
                parent = explicit_closed(odd_core, length, residue)
                require((len(actual) == 2) == parent, (
                    odd_core, length, residue, parent, actual
                ))
                child_histogram[len(actual)] += 1
                digest.update(
                    f"{odd_core}:{length}:{residue}:{actual}\n".encode()
                )
                cells += 1
    return cells, child_histogram, digest.hexdigest()


def max_one_run(residue, length):
    current = 0
    maximum = 0
    for index in range(length):
        if (residue >> index) & 1:
            current += 1
            maximum = max(maximum, current)
        else:
            current = 0
    return maximum


def run_and_density_audit():
    digest = sha256()
    cells = 0
    density_checks = 0
    shorter_hostiles = []
    for odd_core in range(3, MAX_ODD_CORE + 1, 2):
        run_length = ceil(log2(odd_core - 1)) + 1
        found_shorter_hostile = None
        counts = direct_counts(odd_core, DIRECT_MAX_LENGTH)
        for length in range(1, DIRECT_MAX_LENGTH + 1):
            total = 1 << (length - 1)
            nonclosed = total - counts[length - 1]
            blocks = (length - 1) // run_length
            # nonclosed/total <= (1-2^-L)^blocks, checked integrally.
            require(
                nonclosed * (1 << (run_length * blocks))
                <= total * ((1 << run_length) - 1) ** blocks,
                (odd_core, length, nonclosed, total, run_length, blocks),
            )
            density_checks += 1
            for residue in range(1, 1 << length, 2):
                longest = max_one_run(residue, length)
                is_closed = explicit_closed(odd_core, length, residue)
                if longest >= run_length:
                    require(is_closed, (
                        odd_core, length, residue, longest, run_length
                    ))
                if (
                    found_shorter_hostile is None
                    and longest >= run_length - 1
                    and not is_closed
                ):
                    found_shorter_hostile = (length, residue)
                cells += 1
        if found_shorter_hostile is not None:
            shorter_hostiles.append((odd_core,) + found_shorter_hostile)
        digest.update(
            f"{odd_core}:{run_length}:{counts}:{found_shorter_hostile}\n".encode()
        )
    # The sufficient length cannot be uniformly reduced.
    for expected in ((3, 1, 1), (5, 2, 3), (9, 3, 7), (17, 4, 15), (33, 5, 31)):
        require(expected in shorter_hostiles, (expected, shorter_hostiles))
    return cells, density_checks, tuple(shorter_hostiles), digest.hexdigest()


def negative_interval_audit():
    digest = sha256()
    certified = 0
    outside_closed = None
    outside_failed = None
    for odd_core in range(3, MAX_ODD_CORE + 1, 2):
        for length in range(1, DIRECT_MAX_LENGTH + 1):
            modulus = 1 << length
            upper = (1 << (length - 1)) // (odd_core - 1)
            for odd_u in range(1, upper + 1, 2):
                residue = (-odd_u) % modulus
                require(explicit_closed(odd_core, length, residue), (
                    odd_core, length, odd_u, residue, upper
                ))
                top_bit = 1 << (length - 1)
                for left in pair_indices(odd_core):
                    right = odd_core - left
                    left_residue = (left * residue) % modulus
                    right_residue = (right * residue) % modulus
                    require(left_residue & top_bit and right_residue & top_bit, (
                        odd_core, length, odd_u, left, right,
                        left_residue, right_residue,
                    ))
                digest.update(f"{odd_core}:{length}:{odd_u}:{residue}\n".encode())
                certified += 1

            first_outside_odd = upper + 1
            if first_outside_odd % 2 == 0:
                first_outside_odd += 1
            for odd_u in range(first_outside_odd, modulus, 2):
                row = (odd_core, length, odd_u, (-odd_u) % modulus)
                if explicit_closed(odd_core, length, row[3]):
                    if outside_closed is None:
                        outside_closed = row
                elif outside_failed is None:
                    outside_failed = row
                if outside_closed is not None and outside_failed is not None:
                    break
    require(outside_closed is not None, "no outside closed hostile")
    require(outside_failed is not None, "no outside failed hostile")
    return certified, outside_closed, outside_failed, digest.hexdigest()


def fibonacci_numbers(count):
    values = []
    previous, current = 0, 1
    for _ in range(count):
        previous, current = current, previous + current
        values.append(previous)
    return tuple(values)


def sequence_and_automaton_audit():
    digest = sha256()
    selected = {}
    for odd_core in range(3, MAX_ODD_CORE + 1, 2):
        direct = direct_counts(odd_core, DIRECT_MAX_LENGTH)
        automatic = automaton_counts(odd_core, DIRECT_MAX_LENGTH)
        require(direct == automatic, (odd_core, direct, automatic))
        if odd_core in (3, 5, 7):
            selected[odd_core] = automaton_counts(odd_core, LONG_SEQUENCE_LENGTH)
        digest.update(f"{odd_core}:{direct}\n".encode())

    fib = fibonacci_numbers(LONG_SEQUENCE_LENGTH)
    b3 = selected[3]
    for length, (closed_count, fib_count) in enumerate(zip(b3, fib), start=1):
        require(
            closed_count == (1 << (length - 1)) - fib_count,
            (length, closed_count, fib_count),
        )
    return selected, digest.hexdigest()


def explicit_hostile_controls():
    controls = {}
    controls["carry_parity_load_b3_l1"] = {
        "parent": explicit_closed(3, 1, 1),
        "children": tuple(
            explicit_closed(3, 2, 1 + lift * 2) for lift in (0, 1)
        ),
        "predicted_lifts": predicted_closed_lifts(3, 1, 1),
    }
    require(
        controls["carry_parity_load_b3_l1"]
        == {"parent": False, "children": (False, True), "predicted_lifts": (1,)},
        controls,
    )
    controls["zero_child_b5_l1"] = predicted_closed_lifts(5, 1, 1)
    require(controls["zero_child_b5_l1"] == (), controls)
    controls["short_run_b5"] = explicit_closed(5, 2, 3)
    require(controls["short_run_b5"] is False, controls)
    controls["negative_outside_close"] = explicit_closed(5, 4, (-3) % 16)
    controls["negative_outside_fail"] = explicit_closed(5, 4, (-5) % 16)
    require(controls["negative_outside_close"] is True, controls)
    require(controls["negative_outside_fail"] is False, controls)
    for length in range(1, DIRECT_MAX_LENGTH + 1):
        for residue in range(1, 1 << length, 2):
            require(
                explicit_closed(3, length, residue)
                == (max_one_run(residue, length) >= 2),
                (length, residue),
            )
    return controls


def main():
    two_lift = two_lift_audit()
    run_density = run_and_density_audit()
    negative = negative_interval_audit()
    sequence = sequence_and_automaton_audit()
    controls = explicit_hostile_controls()

    print("THM-4250 INDEPENDENT ODD-CORE SUFFIX CLOSURE AUDIT")
    print(f"odd cores: 3..{MAX_ODD_CORE}, odd")
    print(f"direct census lengths: 1..{DIRECT_MAX_LENGTH}")
    print("two-lift cells:", two_lift[0])
    print("two-lift child histogram:", two_lift[1])
    print("two-lift digest:", two_lift[2])
    print("run cells:", run_density[0])
    print("density inequalities:", run_density[1])
    print("L-1 hostile count:", len(run_density[2]))
    print("first L-1 hostiles:", run_density[2][:10])
    print("run/density digest:", run_density[3])
    print("negative-interval certified residues:", negative[0])
    print("first outside closed row:", negative[1])
    print("first outside failed row:", negative[2])
    print("negative-interval digest:", negative[3])
    for odd_core in (3, 5, 7):
        print(f"b={odd_core} closed counts l=1..30:", sequence[0][odd_core])
    print("sequence/automaton digest:", sequence[1])
    print("hostile controls:", controls)
    print("b=3 formula: closed(l)=2^(l-1)-F_l, F_1=F_2=1")
    print("PASS")


if __name__ == "__main__":
    main()
