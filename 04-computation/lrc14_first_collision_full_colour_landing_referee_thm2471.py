#!/usr/bin/env python3
"""Independent exact landing referee for THM-2471.

The control is deliberately hostile to the tempting strengthening "root colour
k lands already at ordinary frequency q=k".  Two complementary Boolean step
functions on a 65-cell circle collide after one P_13 step.  Their ordinary
Fourier products vanish unless 5 divides q, so ten of the twelve base
frequencies vanish.  Nevertheless every progression q=k+13h contains exactly
one landing in every five consecutive h, while the common-root correlation
fires all twelve colours exactly over Q(zeta_13).
"""

from fractions import Fraction as F


P = 13
AUXILIARY_PERIOD = 5
GRID = P * AUXILIARY_PERIOD
SUPPORT = frozenset(range(6))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def phi13_reduce(values: list[F], colour: int) -> tuple[F, ...]:
    """Reduce sum_s values[s] zeta^(-colour*s) in the power basis."""

    raw = [F(0) for _ in range(P)]
    for shift, value in enumerate(values):
        raw[(-colour * shift) % P] += value
    return tuple(raw[index] - raw[P - 1] for index in range(P - 1))


def circular_jump_count(values: list[int]) -> int:
    return sum(
        values[index] != values[(index + 1) % len(values)]
        for index in range(len(values))
    )


def ordinary_coefficient_is_nonzero(frequency: int) -> bool:
    """Exact nonvanishing law for the 65-cell support used below.

    The cell word has period 13 cells, so its discrete endpoint sum contains
    sum_(t=0)^4 exp(-2*pi*i*q*t/5), which vanishes unless 5|q.  When 5|q but
    65 does not divide q, the remaining six-term rational mask is nonconstant
    on C_13 and hence is nonzero by Phi_13 irreducibility.  At nonzero
    multiples of 65 the common cell-integral factor vanishes.
    """

    if frequency == 0:
        return True
    return frequency % AUXILIARY_PERIOD == 0 and frequency % GRID != 0


def main() -> None:
    # U is periodic with period 1/5; V is its complement.  Both are unions of
    # rational half-open grid cells, so endpoint seams are null.
    u_cells = [int(index % P in SUPPORT) for index in range(GRID)]
    v_cells = [1 - value for value in u_cells]
    require(
        all(u_cells[index] == u_cells[(index + P) % GRID] for index in range(GRID)),
        "period-13 cell word",
    )
    require(all(left * right == 0 for left, right in zip(u_cells, v_cells)), "U,V disjoint")
    require(all(left + right == 1 for left, right in zip(u_cells, v_cells)), "U,V complementary")

    jumps_u = circular_jump_count(u_cells)
    jumps_v = circular_jump_count(v_cells)
    require((jumps_u, jumps_v) == (10, 10), "jump census")

    # On base cell b, the roots x=(y+r)/13 occupy fine cells 5r+b.
    root_u = [
        [u_cells[AUXILIARY_PERIOD * root + base] for root in range(P)]
        for base in range(AUXILIARY_PERIOD)
    ]
    root_v = [
        [v_cells[AUXILIARY_PERIOD * root + base] for root in range(P)]
        for base in range(AUXILIARY_PERIOD)
    ]
    require(all(sum(row) == 6 for row in root_u), "P_13 U root count")
    require(all(sum(row) == 7 for row in root_v), "P_13 V root count")
    require(
        all(root_u[base][root] * root_v[base][root] == 0 for base in range(5) for root in range(P)),
        "same-root disjointness",
    )

    correlation = [
        F(
            sum(
                root_u[base][(root + shift) % P] * root_v[base][root]
                for base in range(AUXILIARY_PERIOD)
                for root in range(P)
            ),
            AUXILIARY_PERIOD,
        )
        for shift in range(P)
    ]
    expected_correlation = [F(value) for value in (0, 5, 3, 2, 6, 1, 4, 4, 1, 6, 2, 3, 5)]
    require(correlation == expected_correlation, "root correlation profile")

    service = sum(correlation)
    collision = service / (P * P)
    require(service == 42 and collision == F(42, 169), "service/collision normalization")
    for colour in range(1, P):
        require(any(phi13_reduce(correlation, colour)), f"root colour {colour} vanished")

    # Root orthogonality and Parseval, checked over Q without evaluating any
    # transcendental root of unity.
    root_currents = [
        tuple(entry / (P * P) for entry in phi13_reduce(correlation, colour))
        for colour in range(1, P)
    ]
    signed_sum = tuple(sum(current[index] for current in root_currents) for index in range(P - 1))
    require(signed_sum == (-collision,) + (F(0),) * (P - 2), "signed root-current ledger")

    current_energy = (
        sum(value * value for value in correlation) / P
        - (service / P) ** 2
    ) / (P * P)
    require(current_energy == F(602, 28561), "root-current energy")
    require(current_energy >= collision * collision / 12, "THM-2471 energy floor")

    # Exact ordinary-frequency hostile.  Since V=1-U, every nonzero product
    # is -|U_hat(q)|^2.  It is nonzero precisely on q=5 mod the support law.
    base_landings = []
    first_h = []
    first_q = []
    for colour in range(1, P):
        base_landings.append(ordinary_coefficient_is_nonzero(colour))
        landing_h = [
            h
            for h in range(AUXILIARY_PERIOD)
            if ordinary_coefficient_is_nonzero(colour + P * h)
        ]
        require(len(landing_h) == 1, "one landing per five gauges")
        first_h.append(landing_h[0])
        first_q.append(colour + P * landing_h[0])

        # The pattern is periodic in h with period five.  Every hostile block
        # tested here has exactly one landing; this directly audits the
        # syndetic, rather than merely existential, conclusion.
        for start in range(-50, 51):
            block = [
                ordinary_coefficient_is_nonzero(colour + P * h)
                for h in range(start, start + AUXILIARY_PERIOD)
            ]
            require(sum(block) == 1, "syndetic five-gauge block")

    require(sum(base_landings) == 2, "ten base frequencies must vanish")
    require(
        first_h == [3, 1, 4, 2, 0, 3, 1, 4, 2, 0, 3, 1],
        "first-gauge profile",
    )
    require(
        first_q == [40, 15, 55, 30, 5, 45, 20, 60, 35, 10, 50, 25],
        "first-frequency profile",
    )
    require(max(first_q) == 60, "landing maximum")
    require(
        max(first_q) <= P * jumps_u * jumps_v - 1,
        "endpoint-Prony jump-product bound",
    )

    # Independent sharp equality control for the universal I/12 invoice.
    sharp_correlation = [F(0)] + [F(1, 12)] * 12
    sharp_collision = sum(sharp_correlation) / (P * P)
    sharp_current = -F(1, 12 * P * P)
    require(sharp_collision == F(1, 169), "sharp collision")
    require(abs(sharp_current) == sharp_collision / 12, "sharp maximum invoice")
    for colour in range(1, P):
        reduced = phi13_reduce(sharp_correlation, colour)
        require(reduced == (-F(1, 12),) + (F(0),) * 11, "sharp colour profile")

    print("THM-2471 INDEPENDENT FULL-COLOUR LANDING REFEREE")
    print(f"grid={GRID};support_size={len(SUPPORT)};jumps={jumps_u},{jumps_v}")
    print("correlation=" + ",".join(map(str, correlation)) + f";service={service};I={collision}")
    print(f"nonzero_root_colours=12;root_energy={current_energy};floor={collision * collision / 12}")
    print(f"base_ordinary_landings={sum(base_landings)};base_vanishing_colours={12 - sum(base_landings)}")
    print("first_landing_h=" + ",".join(map(str, first_h)))
    print("first_landing_q=" + ",".join(map(str, first_q)) + f";max={max(first_q)}")
    print("syndetic_gap=5;prony_bound=1299;sharp_max=I/12")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
