#!/usr/bin/env python3
"""Exact audit for THM-2228, the Mahler 3/2 carry/stabilization theorem."""

from fractions import Fraction
from itertools import product
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def carry_integer(word):
    """C_m(word)=sum c_j 2^j 3^(m-1-j)."""
    m = len(word)
    return sum(c * 2**j * 3 ** (m - 1 - j) for j, c in enumerate(word))


def canonical_residue(word):
    """The unique A mod 2^m producing the finite parity itinerary word."""
    m = len(word)
    if m == 0:
        return 0
    modulus = 2**m
    return (-pow(3, -m, modulus) * carry_integer(word)) % modulus


def step_integer(a):
    c = a & 1
    return (3 * a + c) // 2, c


def recover_word(a, m):
    bits = []
    for _ in range(m):
        a, c = step_integer(a)
        bits.append(c)
    return tuple(bits)


def periodic_tail(block, phase=0):
    block = tuple(block)
    length = len(block)
    rotated = block[phase:] + block[:phase]
    numerator = sum(
        c * Fraction(2, 3) ** (j + 1) for j, c in enumerate(rotated)
    )
    return numerator / (1 - Fraction(2, 3) ** length)


def canonical_cycle(cycle):
    cycle = tuple(cycle)
    rotations = [cycle[j:] + cycle[:j] for j in range(len(cycle))]
    return min(rotations)


def lower_half_cycles(q):
    multiplier = 3 * pow(2, -1, q) % q
    seen = set()
    cycles = set()
    for start in range(1, q):
        if start in seen:
            continue
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            seen.add(current)
            current = multiplier * current % q
        cycle = tuple(orbit[orbit.index(current) :])
        if all(2 * residue < q for residue in cycle):
            cycles.add(canonical_cycle(cycle))
    return tuple(sorted(cycles))


def audit_cylinder_bijection():
    total_words = 0
    final_count = 0
    for m in range(1, 17):
        residues = {}
        modulus = 2**m
        for word in product((0, 1), repeat=m):
            residue = canonical_residue(word)
            require(0 <= residue < modulus, f"residue range drift at m={m}")
            require(residue not in residues, f"cylinder collision at m={m}")
            require(recover_word(residue, m) == word, f"itinerary drift at m={m}")
            if m >= 2:
                shifted = canonical_residue(word[1:])
                stepped, _ = step_integer(residue)
                require(
                    stepped % (2 ** (m - 1)) == shifted,
                    f"finite conjugacy drift at m={m}",
                )
            residues[residue] = word
            total_words += 1
        require(len(residues) == modulus, f"cylinder surjectivity drift at m={m}")
        final_count = len(residues)
    return 16, total_words, final_count


def audit_periodic_hostile_control():
    block = (1, 0, 0)
    phases = tuple(periodic_tail(block, phase) for phase in range(3))
    require(phases == (Fraction(18, 19), Fraction(8, 19), Fraction(12, 19)),
            "periodic safe-tail phases drift")
    require(all(phase < 1 for phase in phases), "formal safe cycle drift")
    require(carry_integer(block) == 9, "period carry drift")
    formal_state = -Fraction(9, 19)
    for repetitions in range(1, 9):
        word = block * repetitions
        modulus = 2 ** len(word)
        rational_residue = (
            formal_state.numerator
            * pow(formal_state.denominator, -1, modulus)
        ) % modulus
        require(
            canonical_residue(word) == rational_residue,
            f"periodic 2-adic residue drift at repetitions={repetitions}",
        )
    return phases, formal_state


def audit_stabilizing_bad_tail():
    a = 1
    word = []
    for _ in range(4):
        a, c = step_integer(a)
        word.append(c)
    require(tuple(word) == (1, 0, 1, 1), "A=1 carry prefix drift")
    partial = sum(
        c * Fraction(2, 3) ** (j + 1) for j, c in enumerate(word)
    )
    require(partial == Fraction(94, 81) > 1, "A=1 tail obstruction drift")
    for m in range(1, 13):
        require(
            canonical_residue(recover_word(1, m)) == 1,
            f"A=1 stabilization drift at m={m}",
        )
    return tuple(word), partial


def audit_greedy_boundary():
    y = Fraction(1)
    bits = []
    partial = Fraction(0)
    states = []
    for n in range(256):
        states.append(y)
        digit = (3 * y.numerator) // (2 * y.denominator)
        require(digit in (0, 1), f"greedy digit drift at n={n}")
        bits.append(digit)
        partial += digit * Fraction(2, 3) ** (n + 1)
        y = Fraction(3, 2) * y - digit
        require(y > 0, f"greedy boundary terminated at n={n}")
        require(y < 1, f"greedy boundary escaped at n={n}")
        require(
            Fraction(1) == partial + Fraction(2, 3) ** (n + 1) * y,
            f"greedy expansion identity drift at n={n}",
        )
        require(partial < 1, f"finite boundary prefix ceased to be strict at n={n}")
        prefix_carry = carry_integer(bits)
        require(2 * prefix_carry < 3 ** (n + 1),
                f"raw boundary cylinder drift at n={n}")
        require(2 * prefix_carry + 2 ** (n + 2) >= 3 ** (n + 1),
                f"boundary falsely gained an open-tail certificate at n={n}")
    require(len(set(states)) == len(states), "boundary state repetition drift")
    return "".join(map(str, bits[:64])), y


def audit_periodic_nontermination():
    checked = 0
    for length in range(1, 13):
        denominator = 3**length - 2**length
        require(denominator > 0, f"period denominator sign drift at L={length}")
        for block in product((0, 1), repeat=length):
            carry = carry_integer(block)
            state = -Fraction(carry, denominator)
            if any(block):
                require(state < 0, f"nonzero periodic state sign drift at L={length}")
            else:
                require(state == 0, f"zero periodic state drift at L={length}")
            checked += 1
    return checked


def infinite_prefix_carry(block, phase, length):
    rotated = block[phase:] + block[:phase]
    word = tuple(rotated[j % len(block)] for j in range(length))
    return carry_integer(word)


def audit_open_tail_certificates():
    checked = 0
    safe = 0
    certified = 0
    for length in range(1, 9):
        for block in product((0, 1), repeat=length):
            for phase in range(length):
                tail = periodic_tail(block, phase)
                witness = None
                for m in range(1, 201):
                    carry = infinite_prefix_carry(block, phase, m)
                    if 2 * carry + 2 ** (m + 1) < 3**m:
                        witness = m
                        break
                if tail < 1:
                    safe += 1
                    require(witness is not None,
                            f"open-tail certificate missing at L={length}, phase={phase}")
                    certified += 1
                else:
                    require(witness is None,
                            f"false open-tail certificate at L={length}, phase={phase}")
                checked += 1
    return checked, safe, certified


def audit_lower_half_cycle():
    scanned = []
    for q in range(5, 20, 2):
        if gcd(q, 6) != 1:
            continue
        cycles = lower_half_cycles(q)
        scanned.append((q, cycles))
    first = next((q, cycles) for q, cycles in scanned if cycles)
    require(first == (19, ((4, 6, 9),)), "first lower-half cycle drift")
    return tuple(scanned), first


def audit_pseudo_z_family():
    samples = []
    points_checked = 0
    sample_horizons = {1, 2, 3, 4, 8, 16, 32, 64, 128}
    for horizon in range(1, 129):
        k = 9 * pow(2, -horizon, 19) % 19
        require(1 <= k <= 18, f"pseudo-Z multiplier drift at M={horizon}")
        numerator = k * 2**horizon
        require(numerator % 19 == 9, f"pseudo-Z initial phase drift at M={horizon}")
        integer_part = (numerator - 9) // 19
        require(
            0 <= integer_part < 2**horizon,
            f"pseudo-Z canonical range drift at M={horizon}",
        )
        formal_residue = (-9 * pow(19, -1, 2**horizon)) % 2**horizon
        require(
            integer_part == formal_residue,
            f"pseudo-Z 2-adic convergence drift at M={horizon}",
        )
        phases = []
        for n in range(horizon + 1):
            residue = (k * 3**n * 2 ** (horizon - n)) % 19
            require(residue in (4, 6, 9), f"pseudo-Z phase cycle drift at M={horizon}")
            require(2 * residue < 19, f"pseudo-Z safe inequality drift at M={horizon}")
            phases.append(residue)
            points_checked += 1
        if horizon in sample_horizons:
            samples.append((horizon, k, integer_part, tuple(phases[:6])))
    require(points_checked == 8384, "pseudo-Z point count drift")
    return tuple(samples), points_checked


def main():
    max_bijection_depth, prefix_words, final_residues = audit_cylinder_bijection()
    phases, formal_state = audit_periodic_hostile_control()
    bad_word, bad_partial = audit_stabilizing_bad_tail()
    boundary_prefix, boundary_remainder = audit_greedy_boundary()
    periodic_blocks = audit_periodic_nontermination()
    tail_rows, safe_tail_rows, certified_tail_rows = audit_open_tail_certificates()
    scanned_cycles, first_cycle = audit_lower_half_cycle()
    pseudo_samples, pseudo_points = audit_pseudo_z_family()

    print(f"cylinder_bijection_depth={max_bijection_depth}")
    print(f"cylinder_words_checked={prefix_words}")
    print(f"depth16_distinct_residues={final_residues}")
    print(f"periodic_safe_word=100")
    print(f"periodic_Y_phases={','.join(map(str, phases))}")
    print(f"periodic_fractional_phases={','.join(str(x / 2) for x in phases)}")
    print(f"periodic_2adic_state={formal_state}")
    print(f"stabilizing_A1_prefix={''.join(map(str, bad_word))}")
    print(f"stabilizing_A1_partial_Y={bad_partial}")
    print(f"greedy_boundary_steps=256")
    print(f"greedy_boundary_prefix64={boundary_prefix}")
    print(f"greedy_boundary_remainder256={boundary_remainder}")
    print(f"periodic_blocks_checked={periodic_blocks}")
    print(f"periodic_tail_rows_checked={tail_rows}")
    print(f"safe_tail_rows={safe_tail_rows}")
    print(f"open_tail_certificates={certified_tail_rows}")
    print(
        "lower_half_cycle_scan="
        + ";".join(f"{q}:{cycles}" for q, cycles in scanned_cycles)
    )
    print(f"first_lower_half_cycle={first_cycle}")
    print(f"pseudo_Z_horizons_checked=1..128")
    print(f"pseudo_Z_points_checked={pseudo_points}")
    print(f"pseudo_Z_samples={pseudo_samples}")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
