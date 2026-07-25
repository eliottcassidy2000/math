#!/usr/bin/env python3
"""Exact checks for THM-2276's shallow pair and phase-cone boundary."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd


P = 13
CUTOFF = 757


def require(condition: bool, message: str) -> None:
    """Optimization-safe exact check."""

    if not condition:
        raise AssertionError(message)


def canonical_epsilon(m: int) -> Fraction:
    """A positive rational half-width satisfying both strict bounds."""

    center = Fraction(1, 4 * m)
    owner_margin = Fraction(1, 14) - center
    return min(center, owner_margin) / 2


def circle_norm(x: Fraction) -> Fraction:
    """Distance of an exact rational circle point to the nearest integer."""

    residue = x % 1
    return min(residue, 1 - residue)


def local_private_epsilon(
    *,
    center: Fraction,
    root: int,
    guard: int,
    units: tuple[int, ...],
    blockers: tuple[int, int, int],
) -> Fraction:
    """Half the exact Lipschitz margin around the two lifted centers."""

    margins_in_y = [center, Fraction(1, 14) - center]
    c1, c2, c3 = blockers

    for sign in (-1, 1):
        y = sign * center
        x = (root + y) / 13

        owner_margin = Fraction(1, 14) - circle_norm(c1 * x)
        require(owner_margin > 0, f"nonpositive owner margin at sign={sign}")
        margins_in_y.append(owner_margin / Fraction(c1, 13))

        guard_margin = circle_norm(guard * x) - Fraction(1, 7)
        require(guard_margin > 0, f"nonpositive guard margin at sign={sign}")
        margins_in_y.append(guard_margin / Fraction(guard, 13))

        for q in units:
            unit_margin = circle_norm(q * x) - Fraction(1, 14)
            require(
                unit_margin > 0,
                f"nonpositive unit margin for q={q} at sign={sign}",
            )
            margins_in_y.append(unit_margin / Fraction(q, 13))

        for blocker in (c2, c3):
            deep_margin = circle_norm(blocker * x) - Fraction(1, 14)
            require(
                deep_margin > 0,
                f"nonpositive deep margin for c={blocker} at sign={sign}",
            )
            margins_in_y.append(deep_margin / Fraction(blocker, 13))

    epsilon = min(margins_in_y) / 2
    require(epsilon > 0, "constructed local epsilon is nonpositive")
    return epsilon


def main() -> None:
    pairs = [
        (a, d)
        for a in range(1, CUTOFF + 1)
        for d in range(1, CUTOFF // a + 1)
        if gcd(a, d) == 1 and a % P != 0 and d % P != 0
    ]
    require(len(pairs) == 3279, "thirteen-unit ordered pair census changed")
    require(
        len({tuple(sorted(pair)) for pair in pairs}) == 1640,
        "thirteen-unit unordered pair census changed",
    )

    guard_owner_height = 0
    owner_unit_height = 0
    guard_multipliers: set[int] = set()
    unit_multipliers: set[int] = set()

    for a, d in pairs:
        # H=g*a, u_1=g*d, c_1=13*g*d.  The common factor g cancels.
        require(
            13 * d * a - a * (13 * d) == 0,
            f"guard-owner lift fails for pair {(a, d)}",
        )
        require(
            gcd(13 * d, a) == 1,
            f"guard-owner lift is imprimitive for pair {(a, d)}",
        )
        require(a % 13 != 0, f"guard multiplier a={a} is not a 13-unit")
        guard_owner_height = max(guard_owner_height, 13 * d, a)
        guard_multipliers.add(a)

        # u_1=g*a, q_i=g*d, c_1=13*g*a.
        require(
            d * (13 * a) - 13 * a * d == 0,
            f"owner-unit lift fails for pair {(a, d)}",
        )
        require(
            gcd(d, 13 * a) == 1,
            f"owner-unit lift is imprimitive for pair {(a, d)}",
        )
        require(d % 13 != 0, f"unit multiplier d={d} is not a 13-unit")
        owner_unit_height = max(owner_unit_height, d, 13 * a)
        unit_multipliers.add(d)

    require(guard_owner_height == 9841, "guard-owner maximum height changed")
    require(owner_unit_height == 9841, "owner-unit maximum height changed")

    possible_multipliers = {
        m for m in range(1, CUTOFF + 1) if m % P != 0
    }
    require(
        guard_multipliers == possible_multipliers,
        "guard-owner multiplier projection is incomplete",
    )
    require(
        unit_multipliers == possible_multipliers,
        "owner-unit multiplier projection is incomplete",
    )
    require(
        len(possible_multipliers) == 699,
        "possible multiplier census changed",
    )

    low_bank = {1, 2, 3}
    hard_bank = {
        m for m in range(4, CUTOFF + 1) if m % P != 0
    }
    require(
        possible_multipliers == low_bank | hard_bank,
        "low/hard banks do not partition the multipliers",
    )
    require(low_bank.isdisjoint(hard_bank), "low and hard banks overlap")
    require(len(hard_bank) == 696, "hard-bank census changed")
    require(
        len(range(4, CUTOFF + 1)) == 754,
        "ambient hard interval census changed",
    )
    require(CUTOFF // P == 58, "multiple-of-thirteen census changed")

    # The low phase cone has half-angle pi*m/7 < pi/2 exactly for m<=3.
    require(
        all(2 * m < 7 for m in low_bank),
        "a low multiplier leaves the positive phase cone",
    )
    require(2 * 4 > 7, "multiplier four did not cross the phase boundary")

    # For every hard multiplier, freeze a rational two-interval carrier.
    # Its centers have phases -i and +i at frequency m, so equal interval
    # amplitudes cancel exactly.
    cancellation_rows: list[str] = []
    for m in sorted(hard_bank):
        center = Fraction(1, 4 * m)
        epsilon = canonical_epsilon(m)
        owner_margin = Fraction(1, 14) - center

        require(epsilon > 0, f"nonpositive abstract epsilon for m={m}")
        require(epsilon < center, f"abstract intervals overlap for m={m}")
        require(
            epsilon < owner_margin,
            f"abstract carrier exits D_1 for m={m}",
        )
        require(
            center + epsilon < Fraction(1, 14),
            f"positive abstract interval exits D_1 for m={m}",
        )
        require(
            -center - epsilon > -Fraction(1, 14),
            f"negative abstract interval exits D_1 for m={m}",
        )

        # Exact quarter turns: m*(+center)=+1/4 and m*(-center)=-1/4.
        require(
            m * center == Fraction(1, 4),
            f"positive center is not a quarter turn for m={m}",
        )
        require(
            m * (-center) == Fraction(-1, 4),
            f"negative center is not a quarter turn for m={m}",
        )

        # Encode exp(-2*pi*i*(+1/4))=-i and
        # exp(-2*pi*i*(-1/4))=+i as exact Gaussian integer pairs.
        phase_plus = (0, -1)
        phase_minus = (0, 1)
        require(
            (
                phase_plus[0] + phase_minus[0],
                phase_plus[1] + phase_minus[1],
            )
            == (0, 0),
            f"quarter-turn phases do not cancel for m={m}",
        )

        cancellation_rows.append(f"{m}:{center}:{epsilon}")

    # D_1 has length 1/7<2/13, hence any subset has root occupancy <=2.
    require(
        Fraction(1, 7) < Fraction(2, 13),
        "D_1 no longer has the two-root fibre cap",
    )

    # A two-sheet nonzero root character cannot cancel: if zeta^j=-1,
    # then zeta^(2j)=1, but 2j is nonzero modulo 13.
    character_checks = 0
    for k in range(1, 13):
        for sheet_gap in range(1, 13):
            j = (k * sheet_gap) % 13
            require(
                j != 0,
                f"nonzero character killed nonzero gap {(k, sheet_gap)}",
            )
            require(
                (2 * j) % 13 != 0,
                f"two-sheet cancellation appeared at {(k, sheet_gap)}",
            )
            character_checks += 1
    require(character_checks == 144, "root-character check census changed")

    # Strengthen the abstract carrier to an exact local c_1-private flow.
    # The profile universe is the 150 strict rows 2<=b<c, 5<=c<=19.
    strict_profiles = [
        (b, c)
        for c in range(5, 20)
        for b in range(2, c)
    ]
    require(len(strict_profiles) == 150, "strict-profile census changed")

    local_rows: list[str] = []
    local_witness_checks = 0
    smallest_local_epsilon: Fraction | None = None
    smallest_local_state: tuple[int, int, int] | None = None

    for m in sorted(hard_bank):
        units = tuple(m + 13 * k for k in range(5))
        guard = m if m % 2 else m + 13
        root = (6 * pow(m, -1, 13)) % 13
        center = Fraction(1, 4 * m)

        require(guard % 2 == 1, f"guard is even for m={m}")
        require(guard % 13 != 0, f"guard is not a 13-unit for m={m}")
        require(len(set(units)) == 5, f"unit labels repeat for m={m}")
        require(
            all(q > 0 and q % 13 != 0 for q in units),
            f"a unit speed is invalid for m={m}",
        )
        require(
            (m * root - 6) % 13 == 0,
            f"fixed-root address is wrong for m={m}",
        )

        for b, c in strict_profiles:
            blockers = (13, 2 * m * 13**b, 2 * m * 13**c)
            c1, c2, c3 = blockers
            state = (m, b, c)
            require(c1 < c2 < c3, f"blockers are not distinct at {state}")
            require(c1 % 13 == 0, f"shallow blocker is invalid at {state}")
            require(
                c2 % 13**b == 0 and c2 % 13 ** (b + 1) != 0,
                f"middle blocker valuation is wrong at {state}",
            )
            require(
                c3 % 13**c == 0 and c3 % 13 ** (c + 1) != 0,
                f"deep blocker valuation is wrong at {state}",
            )

            # The pair relation 13*q_1-m*c_1=0 is primitive and has
            # carry +m*c_1=13m.
            require(
                13 * units[0] - m * c1 == 0,
                f"local pair relation fails at {state}",
            )
            require(gcd(13, m) == 1, f"local pair is imprimitive at {state}")
            require(
                max(13, m) <= 757,
                f"local pair height exceeds 757 at {state}",
            )

            for sign in (-1, 1):
                y = sign * center
                x = (root + y) / 13
                require(
                    circle_norm(c1 * x) == center < Fraction(1, 14),
                    f"shallow owner misses center {(state, sign)}",
                )
                require(
                    circle_norm(guard * x) > Fraction(1, 7),
                    f"guard fails at center {(state, sign)}",
                )
                require(
                    all(
                        circle_norm(q * x) > Fraction(1, 14) for q in units
                    ),
                    f"a unit is dangerous at center {(state, sign)}",
                )
                require(
                    circle_norm(c2 * x) == Fraction(1, 2),
                    f"middle blocker is not antipodal at {(state, sign)}",
                )
                require(
                    circle_norm(c3 * x) == Fraction(1, 2),
                    f"deep blocker is not antipodal at {(state, sign)}",
                )

            epsilon = local_private_epsilon(
                center=center,
                root=root,
                guard=guard,
                units=units,
                blockers=blockers,
            )
            require(epsilon < center, f"local intervals overlap at {state}")
            require(
                center + epsilon < Fraction(1, 14),
                f"local carrier exits D_1 at {state}",
            )

            # Hostile endpoint checks. The norm is 1-Lipschitz, and the
            # margin routine already proves the full open intervals.
            for sign in (-1, 1):
                for delta_sign in (-1, 1):
                    y = sign * center + delta_sign * epsilon
                    x = (root + y) / 13
                    endpoint = (state, sign, delta_sign)
                    require(
                        circle_norm(c1 * x) < Fraction(1, 14),
                        f"shallow owner fails at endpoint {endpoint}",
                    )
                    require(
                        circle_norm(guard * x) > Fraction(1, 7),
                        f"guard fails at endpoint {endpoint}",
                    )
                    require(
                        all(
                            circle_norm(q * x) > Fraction(1, 14)
                            for q in units
                        ),
                        f"a unit is dangerous at endpoint {endpoint}",
                    )
                    require(
                        circle_norm(c2 * x) > Fraction(1, 14),
                        f"middle blocker is dangerous at endpoint {endpoint}",
                    )
                    require(
                        circle_norm(c3 * x) > Fraction(1, 14),
                        f"deep blocker is dangerous at endpoint {endpoint}",
                    )

            # On the fixed root branch T is one-to-one, hence the ancestry
            # multiplicity equals the image indicator. Quarter-turn
            # cancellation at m then scales back to cancellation at 13m.
            require(
                m * center == Fraction(1, 4),
                f"image phase is not a quarter turn at {state}",
            )
            require(
                13 * m * (center / 13) == Fraction(1, 4),
                f"ancestry phase does not scale at {state}",
            )

            local_rows.append(f"{m}:{b}:{c}:{root}:{epsilon}")
            local_witness_checks += 1
            if smallest_local_epsilon is None or epsilon < smallest_local_epsilon:
                smallest_local_epsilon = epsilon
                smallest_local_state = (m, b, c)

    require(
        local_witness_checks == 696 * 150 == 104400,
        "local witness census changed",
    )
    require(
        smallest_local_epsilon is not None,
        "smallest local epsilon was not recorded",
    )
    require(
        smallest_local_state is not None,
        "smallest local epsilon state was not recorded",
    )

    hard_csv = ",".join(map(str, sorted(hard_bank)))
    carrier_digest = sha256("\n".join(cancellation_rows).encode()).hexdigest()
    hard_digest = sha256(hard_csv.encode()).hexdigest()
    local_digest = sha256("\n".join(local_rows).encode()).hexdigest()

    print("THM-2276 shallow-owner phase-cone exact verification")
    print(f"ordered_thirteen_unit_pair_shapes={len(pairs)}")
    print("unordered_thirteen_unit_pair_shapes=1640")
    print(f"guard_owner_max_height={guard_owner_height}")
    print(f"owner_unit_max_height={owner_unit_height}")
    print(f"possible_multiplier_count={len(possible_multipliers)}")
    print("low_phase_cone_multipliers=1,2,3")
    print(f"hard_bank_count={len(hard_bank)}")
    print(f"hard_bank_sha256={hard_digest}")
    print(f"cancellation_carrier_sha256={carrier_digest}")
    print(f"two_sheet_character_checks={character_checks}")
    print(f"strict_profile_count={len(strict_profiles)}")
    print(f"local_private_witness_checks={local_witness_checks}")
    print(f"local_private_witness_sha256={local_digest}")
    print(
        "smallest_local_epsilon="
        f"{smallest_local_epsilon} at_m_b_c={smallest_local_state}"
    )
    for m in (4, 5, 12, 14, 757):
        if m in hard_bank:
            center = Fraction(1, 4 * m)
            epsilon = canonical_epsilon(m)
            print(
                f"carrier_m={m} center={center} epsilon={epsilon} "
                f"left={-center-epsilon} right={center+epsilon}"
            )
    print("PASS")


if __name__ == "__main__":
    main()
