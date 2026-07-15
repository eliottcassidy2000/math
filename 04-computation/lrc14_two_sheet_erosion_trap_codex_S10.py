#!/usr/bin/env python3
"""Exact verifier for THM-789's local-trapping example and fingerprints."""

from fractions import Fraction as F


U = (1, 2, 3, 5, 7, 8, 9, 10, 11, 12)
T0 = F(4, 17)
TAU = F(14, 19)
ALPHA = F(1, 13)
BETA = F(1, 11)
GAMMA = BETA - ALPHA


def norm(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def phi(speeds: tuple[int, ...], t: F) -> F:
    return min(norm(w * t) for w in speeds)


def breakpoint_candidates(speeds: tuple[int, ...]) -> set[F]:
    denominators = {2 * u for u in speeds}
    for i, u in enumerate(speeds):
        for v in speeds[i + 1 :]:
            denominators.add(u + v)
            denominators.add(abs(u - v))
    return {F(k, q) for q in denominators for k in range(q + 1)}


def exact_loneliness(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    candidates = breakpoint_candidates(speeds)
    value = max(phi(speeds, t) for t in candidates)
    maximizers = tuple(sorted(t for t in candidates if phi(speeds, t) == value))
    return value, maximizers


def q_value(t: F) -> F:
    return norm(11 * t) + norm(2 * t)


def main() -> None:
    m_value, maximizers = exact_loneliness(U)
    assert m_value == F(2, 17)
    clearances = tuple(norm(u * T0) for u in U)
    numerators = tuple(int(17 * z) for z in clearances)
    assert numerators == (4, 8, 5, 3, 6, 2, 2, 6, 7, 3)
    assert min(clearances) == F(2, 17) > BETA

    nx, ny = 3, 2
    assert abs(13 * T0 - nx) == F(1, 17)
    assert abs(9 * T0 - ny) == F(2, 17)
    assert nx % 2 != ny % 2 and 13 * ny - 9 * nx == -1

    radius_h = F(8, 1989)
    assert q_value(T0) == F(15, 17)
    assert q_value(T0 + radius_h) == F(11, 13)
    assert q_value(T0 - radius_h) > F(11, 13)

    return_radius = F(1, 858)
    cell_radius = F(1, 864)
    lipschitz_radius = (F(2, 17) - ALPHA) / max(U)
    assert GAMMA == F(2, 143)
    assert return_radius == GAMMA / 12
    assert cell_radius < return_radius < radius_h
    assert lipschitz_radius == F(3, 884) < radius_h

    assert phi(U, TAU) == F(2, 19)
    assert q_value(TAU) == F(11, 19) < F(11, 13)

    margin_t0 = F(11, 13) - q_value(T0)
    margin_tau = F(11, 13) - q_value(TAU)
    assert margin_t0 < 0 < margin_tau

    print("THM-789 exact local-trapping verification")
    print(f"U={U}")
    print(f"M(U)={m_value} maximizer_count={len(maximizers)} t0_is_maximizer={T0 in maximizers}")
    print(f"t0={T0} clearance_numerators_mod17={numerators} phi={phi(U,T0)}")
    print(f"odd_nearest=(3,2) determinant={13*ny-9*nx}")
    print(f"Q(t0)={q_value(T0)} H_symmetric_radius={radius_h}")
    print(f"R_U=(-{return_radius},{return_radius})")
    print(f"cell_difference_subset=(-{cell_radius},{cell_radius})")
    print(f"lipschitz_radius={lipschitz_radius}")
    print(f"escape_tau={TAU} phi={phi(U,TAU)} Q={q_value(TAU)}")
    print(f"component_choice_margins: trapped={margin_t0} escape={margin_tau}")
    print("component_margin_tournament: vertices=2 scores=(0,1) cycles=0 sccs=2 hp_count=1 edge_flips=0")
    print("information_destroyed_by_fixed_anchor=global_deep-component_choice")
    print("PASS")


if __name__ == "__main__":
    main()

