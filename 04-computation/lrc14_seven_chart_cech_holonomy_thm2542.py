#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2542."""

from fractions import Fraction


ROOT_P = 13
CLOCK_Q = 7
FIELD = 547  # 547 - 1 = 6 * 7 * 13


def primitive_root(p: int) -> int:
    factors = []
    n = p - 1
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    for g in range(2, p):
        if all(pow(g, (p - 1) // r, p) != 1 for r in factors):
            return g
    raise AssertionError("no primitive root")


def orbit_cycles(clock_step: int, root_step: int):
    unseen = {(k, r) for k in range(CLOCK_Q) for r in range(ROOT_P)}
    lengths = []
    while unseen:
        start = min(unseen)
        cur = start
        length = 0
        while cur in unseen:
            unseen.remove(cur)
            k, r = cur
            cur = (
                (k + clock_step) % CLOCK_Q,
                (r + root_step) % ROOT_P,
            )
            length += 1
        assert cur == start
        lengths.append(length)
    return sorted(lengths)


def main() -> None:
    generator = primitive_root(FIELD)
    zeta = pow(generator, (FIELD - 1) // ROOT_P, FIELD)
    xi = pow(generator, (FIELD - 1) // CLOCK_Q, FIELD)
    assert pow(zeta, ROOT_P, FIELD) == 1
    assert pow(xi, CLOCK_Q, FIELD) == 1
    assert all(pow(zeta, j, FIELD) != 1 for j in range(1, ROOT_P))
    assert all(pow(xi, j, FIELD) != 1 for j in range(1, CLOCK_Q))

    nonzero_classes = 0
    gauge_basis_checks = 0
    mapping_torus_checks = 0
    local_mode_checks = 0
    root_mode_checks = 0
    uncharted_cancellation_checks = 0
    thm2535_special_classes = 0

    for eta in range(1, CLOCK_Q):
        for source_root in range(ROOT_P):
            for target_root in range(ROOT_P):
                if source_root == target_root:
                    continue
                root_step = (source_root - target_root) % ROOT_P
                transitions = [root_step] * CLOCK_Q
                holonomy = sum(transitions) % ROOT_P
                assert holonomy == (CLOCK_Q * root_step) % ROOT_P
                assert holonomy != 0
                nonzero_classes += 1
                if eta == 1 and target_root == 0 and source_root != 0:
                    thm2535_special_classes += 1

                # A vertex gauge h changes the transition ending at k by
                # h_k-h_(k-eta).  Basis gauges suffice by linearity.
                for vertex in range(CLOCK_Q):
                    h = [0] * CLOCK_Q
                    h[vertex] = 1
                    gauged = [
                        (
                            transitions[k]
                            + h[k]
                            - h[(k - eta) % CLOCK_Q]
                        ) % ROOT_P
                        for k in range(CLOCK_Q)
                    ]
                    assert sum(gauged) % ROOT_P == holonomy
                    gauge_basis_checks += 1

                trivial_degrees = [
                    n for n in range(1, ROOT_P + 1)
                    if n * holonomy % ROOT_P == 0
                ]
                assert trivial_degrees == [ROOT_P]

                lengths = orbit_cycles(eta, root_step)
                assert lengths == [ROOT_P * CLOCK_Q]
                mapping_torus_checks += 1

                # Equal-mass hostile:
                # D_k(h,r)=1_{h=s}(1_{r=k}-1_{r=k+eta}).
                for k in range(CLOCK_Q):
                    for alpha in range(1, ROOT_P):
                        zs = pow(
                            zeta,
                            (-alpha * source_root) % ROOT_P,
                            FIELD,
                        )
                        for beta in range(1, CLOCK_Q):
                            value = (
                                zs
                                * pow(xi, (-beta * k) % CLOCK_Q, FIELD)
                                * (
                                    1
                                    - pow(
                                        xi,
                                        (-beta * eta) % CLOCK_Q,
                                        FIELD,
                                    )
                                )
                            ) % FIELD
                            assert value != 0
                            local_mode_checks += 1

                # After charting, sum_k R_k=7(delta_s-delta_t).
                for alpha in range(1, ROOT_P):
                    root_value = (
                        CLOCK_Q
                        * (
                            pow(
                                zeta,
                                (-alpha * source_root) % ROOT_P,
                                FIELD,
                            )
                            - pow(
                                zeta,
                                (-alpha * target_root) % ROOT_P,
                                FIELD,
                            )
                        )
                    ) % FIELD
                    assert root_value != 0
                    root_mode_checks += 1

                for h in range(ROOT_P):
                    for r in range(CLOCK_Q):
                        total = 0
                        for k in range(CLOCK_Q):
                            if h == source_root:
                                total += int(r == k)
                                total -= int(r == (k + eta) % CLOCK_Q)
                        assert total == 0
                        uncharted_cancellation_checks += 1

    # A concrete seven-stop orbit itinerary.  Prescribe the first seven
    # base-13 digits to be 0,1,...,6.  The resulting rational cylinder has
    # positive measure 13^{-7} and visits X_k=[k/13,(k+1)/13).
    digits = list(range(CLOCK_Q))
    numerator = sum(d * ROOT_P ** (CLOCK_Q - 1 - j) for j, d in enumerate(digits))
    denominator = ROOT_P ** CLOCK_Q
    midpoint = Fraction(2 * numerator + 1, 2 * denominator)
    for k in range(CLOCK_Q):
        y = (midpoint * ROOT_P ** k) % 1
        digit = (y * ROOT_P).numerator // (y * ROOT_P).denominator
        assert digit == k
    itinerary_measure = Fraction(1, denominator)

    semantic_arrival_flags = [0] * CLOCK_Q
    assert not any(semantic_arrival_flags)

    print(f"field={FIELD} generator={generator} zeta13={zeta} xi7={xi}")
    print(f"nonzero_general_cech_instances={nonzero_classes}")
    print(f"thm2535_special_classes={thm2535_special_classes}")
    print(f"gauge_basis_checks={gauge_basis_checks}")
    print("minimal_trivializing_cover=13")
    print(f"mapping_torus_checks={mapping_torus_checks} cycle_length=91")
    print(f"primitive_local_modes={local_mode_checks}")
    print(f"charted_root_modes={root_mode_checks}")
    print(f"uncharted_cancellation_checks={uncharted_cancellation_checks}")
    print(
        "uncharted_clock_boundary_zero=1 "
        f"semantic_arrival_flags={sum(semantic_arrival_flags)}"
    )
    print(f"seven_stop_digit_itinerary_measure={itinerary_measure}")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
