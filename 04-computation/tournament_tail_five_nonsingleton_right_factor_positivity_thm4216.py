#!/usr/bin/env python3
"""Exact symbolic and transferred-data audit for THM-4216.

The proof starts from THM-4208's unary response at
Q_5=C3 triangleright P_5.  It changes the right rooted coordinates from
(V^0,V^1) to endpoint mass v=V^0+V^1 and Hamilton-start count
t=V^1-V^0, then uses the exact fan identity r=v-t.  The apparently
dangerous current term becomes a nonnegative edge-avoidance term plus three
controlled moments.
"""

from __future__ import annotations

import hashlib

import sympy as sp

import tournament_cycle_prefix_arbitrary_context_thm4208 as response
import tournament_ordinal_cocycle_parity_thm4184 as base


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ordinal(*factors: base.TournamentData) -> base.TournamentData:
    need(bool(factors), "ordinal product needs a factor")
    value = factors[0]
    for factor in factors[1:]:
        value = base.ordinal_data(value, factor)
    return value


def exact_decomposition(data: base.TournamentData) -> tuple[int, int, int]:
    h = data.hamilton
    w = data.mass
    degree, _, incoming = response.directed_masses(data)
    v = [row[0] + row[1] for row in data.ends]
    t = [row[1] - row[0] for row in data.ends]
    need(all(t_i >= 0 for t_i in t), "Hamilton-start coordinates")
    need(
        all(v[i] == incoming[i] + t[i] for i in range(len(v))),
        "endpoint fan identity v=r+t",
    )
    moment = sum(value * value for value in v)
    mixed = sum(v_i * t_i for v_i, t_i in zip(v, t))
    endpoint_square = sum(t_i * t_i for t_i in t)
    avoidance = sum(
        (1143 * v[i] + 9 * t[i]) * (w - degree[i])
        for i in range(len(v))
    )
    need(avoidance >= 0, "nonnegative edge-avoidance term")
    exact = (
        353088 * h * h
        + 472302 * h * w
        + 118656 * w * w
        + avoidance
        - 352116 * moment
        - 1170 * mixed
        + 18 * endpoint_square
    )
    lower = 6 * (-33 * h * h + 274 * h * w + 214 * w * w)
    endpoint_defect = w * w + 4 * h * w + 3 * h * h - 3 * moment
    need(endpoint_defect >= 0, "endpoint-energy cap")
    need(mixed <= (w + h) * h, "mixed endpoint cap")
    debt = (
        avoidance
        + 117372 * endpoint_defect
        + 1170 * ((w + h) * h - mixed)
        + 18 * endpoint_square
    )
    need(exact - lower == debt >= 0, "canonical lower-bound debt")
    return exact, lower, debt


def main() -> None:
    # Symbolic reconstruction of the exact coordinate change.
    h, w = sp.symbols("h w", nonnegative=True)
    v, t, m, p, tau = sp.symbols("v t m p tau", nonnegative=True)
    v0, v1 = (v - t) / 2, (v + t) / 2
    rooted_linear = sp.expand(1134 * v0 + 1152 * v1)
    need(rooted_linear == 1143 * v + 9 * t, "linear coordinate change")

    m00 = (m - 2 * p + tau) / 4
    m01 = (m - tau) / 4
    m11 = (m + 2 * p + tau) / 4
    rooted_quadratic = sp.expand(
        -341856 * m00 - 695052 * m01 - 353268 * m11
    )
    need(
        rooted_quadratic == -347544 * m - 5706 * p - 18 * tau,
        "quadratic coordinate change",
    )

    fan_rewrite = sp.expand(
        rooted_quadratic
        - 4 * (1143 * m - 1134 * p - 9 * tau)
    )
    need(
        fan_rewrite == -352116 * m - 1170 * p + 18 * tau,
        "fan-to-avoidance rewrite",
    )

    lower = 6 * (-33 * h**2 + 274 * h * w + 214 * w**2)
    u = sp.symbols("u", nonnegative=True)
    need(
        sp.expand(
            lower.subs(w, h + u)
            - 2730 * h**2
            - 12 * u * (351 * h + 107 * u)
        )
        == 0,
        "non-singleton marked-Hamilton floor",
    )

    one = base.tournament_data(base.transitive(1))
    p2 = base.tournament_data(base.transitive(2))
    cycle = base.tournament_data(base.parse("101", 3))
    q = {
        tail: cycle
        if tail == 0
        else base.ordinal_data(cycle, base.tournament_data(base.transitive(tail)))
        for tail in range(0, 8)
    }

    contexts = tuple(
        data for order in range(1, 6) for data in response.labelled(order)
    )
    need(len(contexts) == 1099, "labelled context count through order five")
    exact_rows = 0
    nonsingleton_rows = 0
    propagated_rows = 0
    minimum_gap: int | None = None
    digest = hashlib.sha256()
    for data in contexts:
        direct = base.remainder(q[5], data)
        exact, certified, debt = exact_decomposition(data)
        need(direct == exact, "exact Q5 unary decomposition")
        need(direct >= certified, "Q5 certified polynomial floor")
        if len(data.out) == 1:
            need(direct == -180, "singleton boundary")
        else:
            need(data.mass >= data.hamilton, "marked-Hamilton non-singleton bound")
            need(direct >= 2730 * data.hamilton**2 > 0,
                 "non-singleton positivity")
            gap = direct - 2730 * data.hamilton**2
            minimum_gap = gap if minimum_gap is None else min(minimum_gap, gap)
            nonsingleton_rows += 1
        for tail in (6, 7):
            later = base.remainder(q[tail], data)
            need(later >= 10764 * data.hamilton**2 > 0,
                 "later-tail universal positivity control")
            propagated_rows += 1
        digest.update(
            f"{base.label(data.out)}|{direct}|{certified}|{debt}\n".encode()
        )
        exact_rows += 1

    p2_hostiles = tuple(base.remainder(q[tail], p2) for tail in range(5))
    need(
        p2_hostiles == (-288, -684, -1368, -2232, -1512),
        "sharp non-singleton pre-threshold hostiles",
    )
    singleton_hostiles = tuple(base.remainder(q[tail], one) for tail in range(6))
    need(
        singleton_hostiles == (-72, -216, -468, -900, -1332, -180),
        "sharp all-context pre-threshold hostiles",
    )

    # A small direct control for the strengthened THM-4213 corollary: a
    # tail-five-separated two-cycle word now needs no final extra singleton.
    two_cycle = ordinal(q[5], q[5])
    no_sink_controls = tuple(
        data
        for order in range(3, 5)
        for data in response.labelled(order)
        if not base.has_sink(data.out)
    )
    need(len(no_sink_controls) == 34, "labelled no-sink controls")
    for data in no_sink_controls:
        need(base.remainder(two_cycle, data) > 0,
             "final-tail-five multi-cycle OS+ control")

    lines = [
        "theorem=THM-4216",
        "symbolic_coordinate_change=PASS",
        "symbolic_fan_to_avoidance_rewrite=PASS",
        "certified_floor=6*(-33*H^2+274*H*W+214*W^2)",
        "nonsingleton_uniform_floor=2730*H^2",
        f"labelled_contexts_orders_1_to_5={len(contexts)}",
        f"exact_decomposition_rows={exact_rows}",
        f"nonsingleton_positive_rows={nonsingleton_rows}",
        f"later_tail_rows={propagated_rows}",
        f"minimum_finite_gap_over_2730H2={minimum_gap}",
        "p2_hostiles_n0_to_n4=" + ",".join(map(str, p2_hostiles)),
        "p1_hostiles_n0_to_n5=" + ",".join(map(str, singleton_hostiles)),
        f"final_tail_five_no_sink_controls={len(no_sink_controls)}",
        f"semantic_sha256={digest.hexdigest()}",
    ]
    print("\n".join(lines))


if __name__ == "__main__":
    main()
