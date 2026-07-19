#!/usr/bin/env python3
"""Dependency-free exact referee for THM-1252.

The strict-cover topology (common deletion-minimal subcover, adjacency of
overlapping selected teeth, and wall-provider continuation) is the paper
layer.  This script audits the finite functional-cycle fact, exact wall
address algebra, detuned toothpick law, disjoint quantum packing, compact
rung ceiling, and the sharp centered-spoke guardrail.
"""

from fractions import Fraction
from itertools import product
from math import gcd, lcm


def selected_cycle(mapping: tuple[int, ...], start: int) -> tuple[int, ...]:
    seen: dict[int, int] = {}
    orbit: list[int] = []
    vertex = start
    while vertex not in seen:
        seen[vertex] = len(orbit)
        orbit.append(vertex)
        vertex = mapping[vertex]
    return tuple(orbit[seen[vertex] :])


def blocker_cycle_audit() -> tuple[int, int]:
    maps = 0
    cycles: set[tuple[int, ...]] = set()
    choices = [tuple(j for j in range(6) if j != i) for i in range(6)]
    for mapping in product(*choices):
        maps += 1
        for start in range(6):
            cycle = selected_cycle(mapping, start)
            # Canonical rotation, only to avoid recounting the same cycle.
            pivot = min(range(len(cycle)), key=cycle.__getitem__)
            cycle = cycle[pivot:] + cycle[:pivot]
            cycles.add(cycle)
            assert any(
                cycle[(i + 1) % len(cycle)] > cycle[i]
                for i in range(len(cycle))
            )
    assert maps == 5**6
    return maps, len(cycles)


def all_edge_carrier_margin_audit() -> int:
    rows = 0
    for c in range(1, 81):
        for source in range(c + 1, 3 * c + 1):
            for target in range(c + 1, 3 * c + 1):
                tooth_width = Fraction(1, 7 * target)
                rounding_error = Fraction(1, 2 * (c + source))
                carrier_half_gap = Fraction(3, 7 * c)
                assert tooth_width + rounding_error < carrier_half_gap
                assert tooth_width < Fraction(4, 28 * c)
                rows += 1
    return rows


def containing_tooth(speed: int, wall_numerator: int, wall_denominator: int):
    """Return the unique strict-danger tooth address containing a wall."""
    scaled = speed * wall_numerator
    base = scaled // wall_denominator
    hits = []
    for address in range(base - 1, base + 3):
        if 14 * abs(scaled - address * wall_denominator) < wall_denominator:
            hits.append(address)
    assert len(hits) <= 1
    return hits[0] if hits else None


def wall_address_audit() -> tuple[int, int, int]:
    rows = 0
    spanning_rows = 0
    proper_return_rows = 0
    for j in range(2, 81):
        for h in range(1, 12 * j + 1):
            if h == j:
                continue
            for target_address in range(1, j):
                denominator = 14 * j
                left_num = 14 * target_address - 1
                right_num = 14 * target_address + 1
                m_left = containing_tooth(h, left_num, denominator)
                m_right = containing_tooth(h, right_num, denominator)
                if m_left is None or m_right is None:
                    continue
                rows += 1
                if m_left == m_right:
                    # One provider tooth spans the target tooth.  This is
                    # exactly the branch excluded by irredundancy (or by the
                    # minimum-speed gauge), and width forces h<j.
                    spanning_rows += 1
                    assert h < j
                    continue

                assert m_left < m_right
                rung = m_right - m_left
                left = Fraction(left_num, denominator)
                right = Fraction(right_num, denominator)
                seam_left = Fraction(14 * m_left + 1, 14 * h) - left
                seam_right = right - Fraction(14 * m_right - 1, 14 * h)
                assert seam_left > 0 and seam_right > 0
                detuning = h - (7 * rung - 1) * j
                assert detuning > 0
                assert seam_left + seam_right == Fraction(detuning, 7 * j * h)
                common = gcd(j, h)
                assert detuning % common == 0
                assert detuning >= common
                assert seam_left >= Fraction(common, 14 * j * h)
                assert seam_right >= Fraction(common, 14 * j * h)
                proper_return_rows += 1
    return rows, spanning_rows, proper_return_rows


def compact_rung_audit() -> int:
    rows = 0
    for j in range(1, 151):
        for rung in range(1, 336):
            # Sample every small positive gcd-compatible detuning.  The
            # identity is symbolic; this lane checks indexing and constants.
            for detuning in range(1, 25):
                h = (7 * rung - 1) * j + detuning
                if not h < 2345 * j:
                    continue
                common = gcd(j, h)
                if detuning % common:
                    continue
                seam_sum = Fraction(1, 7 * j) - Fraction(7 * rung - 1, 7 * h)
                assert seam_sum == Fraction(detuning, 7 * j * h)
                assert h > (7 * rung - 1) * j
                assert h > 6 * j
                assert rung <= 335
                assert detuning >= common
                rows += 1

    # The next rung is incompatible with the compact ratio h/j<2345.
    for j in range(1, 151):
        assert (7 * 336 - 1) * j >= 2351 * j > 2345 * j
    return rows


def two_quantum_packing_audit() -> int:
    rows = 0
    for j in range(1, 61):
        for h_left in range(1, 61):
            for h_right in range(1, 61):
                q_left = Fraction(1, 14 * lcm(j, h_left))
                q_right = Fraction(1, 14 * lcm(j, h_right))
                assert q_left <= Fraction(1, 14 * j)
                assert q_right <= Fraction(1, 14 * j)
                assert q_left + q_right <= Fraction(1, 7 * j)
                rows += 1
    return rows


def sharp_guardrail() -> None:
    c, k, source, target, provider = 5, 2, 7, 8, 49
    packet = (7, 8, 49, 50, 51, 52)
    spoke = Fraction(1, 2)
    gap_left = Fraction(14 * k + 1, 14 * c)
    gap_right = Fraction(14 * k + 13, 14 * c)
    assert gap_left < spoke < gap_right

    def distance(speed: int, time: Fraction) -> Fraction:
        value = speed * time
        floor = value.numerator // value.denominator
        frac = value - floor
        return min(frac, 1 - frac)

    assert distance(c, spoke) > Fraction(1, 4)
    assert distance(source, spoke) > Fraction(1, 4)
    dangerous = [v for v in packet if v != source and distance(v, spoke) < Fraction(1, 14)]
    assert min(dangerous) == target

    left, right = Fraction(55, 112), Fraction(57, 112)
    assert gap_left < left < spoke < right < gap_right
    assert containing_tooth(provider, 55, 112) == 24
    assert containing_tooth(provider, 57, 112) == 25
    seam_left = Fraction(14 * 24 + 1, 14 * provider) - left
    seam_right = right - Fraction(14 * 25 - 1, 14 * provider)
    assert seam_left == seam_right == Fraction(1, 5488)
    assert seam_left + seam_right == Fraction(1, 2744)
    assert provider - 6 * target == 1


def binary_rung_family_audit() -> int:
    """Binary centered digits coexist with every compact toothpick rung."""
    c = 10001
    source = 10003
    target = 10002
    k = (c - 1) // 2
    target_address = target // 2
    assert c % 2 == source % 2 == 1 and target % 2 == 0
    assert c < target < source
    middle = k + target_address
    target_clock = c + target
    assert 2 * middle == target_clock - 1
    assert ((target_clock - 1) // 2) - middle == 0
    assert ((target_clock + 1) // 2) - middle == 1

    rows = 0
    for rung in range(1, 336):
        detuning = 1 if rung % 2 else 2
        provider = (7 * rung - 1) * target + detuning
        assert provider < 2345 * c
        assert gcd(target, provider) == detuning
        left_address = (provider - rung) // 2
        right_address = (provider + rung) // 2
        assert right_address - left_address == rung

        left = Fraction(14 * target_address - 1, 14 * target)
        right = Fraction(14 * target_address + 1, 14 * target)
        seam_left = Fraction(14 * left_address + 1, 14 * provider) - left
        seam_right = right - Fraction(14 * right_address - 1, 14 * provider)
        quantum = Fraction(detuning, 14 * target * provider)
        assert seam_left == seam_right == quantum
        assert seam_left + seam_right == Fraction(
            detuning, 7 * target * provider
        )
        rows += 1
    return rows


def main() -> None:
    maps, cycles = blocker_cycle_audit()
    margin_rows = all_edge_carrier_margin_audit()
    wall_rows, spanning_rows, return_rows = wall_address_audit()
    rung_rows = compact_rung_audit()
    packing_rows = two_quantum_packing_audit()
    sharp_guardrail()
    binary_rungs = binary_rung_family_audit()

    print("THM-1252 COHERENT BLOCKER TWO-WALL FORK EXACT AUDIT")
    print(f"loopless blocker maps = {maps}")
    print(f"distinct directed cycles checked = {cycles}")
    print(f"all-edge centered target-tooth margin rows = {margin_rows}")
    print(f"two-wall same-label address rows = {wall_rows}")
    print(f"single-tooth spanning rows excluded by irredundancy = {spanning_rows}")
    print(f"proper detuned address-return rows = {return_rows}")
    print(f"compact gcd-compatible rung rows = {rung_rows}")
    print(f"two-quantum disjoint-packing rows = {packing_rows}")
    print(f"binary-digit-compatible sharp rungs = {binary_rungs}")
    print("finite rung alphabet = 1 <= r <= 335")
    print("sharp local row = (c,k,i,j,h,r,E)=(5,2,7,8,49,1,1)")
    print("sharp wall seams = 1/5488 + 1/5488 = 1/2744")
    print("full-word debt = H >= 1/c + (7/12) sum_a 1/lcm(s_a,s_(a+1))")
    print("RESULT: PASS")


if __name__ == "__main__":
    main()
