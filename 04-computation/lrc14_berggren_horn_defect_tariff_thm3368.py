#!/usr/bin/env python3
"""Exact companion for THM-3368's Berggren Horn/clock separation.

The general proof is algebraic: THM-3357's weighted Horn inequality and
``G >= -max(0,-G)`` give the defect tariff.  This script independently
checks the underlying polynomial identities, exercises the tariff on hostile
integer decks, and freezes the three exact examples used by THM-3368.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd, lcm

import sympy as sp


KAPPAS = (0, 1, 7, 91, 137)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def branches(m, n):
    return (n, 2 * n - m), (n, 2 * n + m), (m, 2 * m + n)


def pythagorean_weights(m, n):
    return n * n - m * m, 2 * m * n, n * n + m * m


def correction(m, n):
    return 2 * m * (n - m) * (3 * n * n + 4 * m * n - m * m)


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def delta(v, deck):
    return max(abs(det(v, column)) for column in deck)


def score(v, deck, kappa=91):
    return v[0] * v[0] + v[1] * v[1] - kappa * delta(v, deck)


def physical_row(v, deck):
    return tuple(v[0] * column[0] + v[1] * column[1] for column in deck)


def symbolic_control():
    m, n = sp.symbols("m n")
    x = (n, 2 * n - m)
    y = (n, 2 * n + m)
    z = (m, 2 * m + n)
    a = n**2 - m**2
    b = 2 * m * n
    c = n**2 + m**2
    K = 2 * m * (n - m) * (3 * n**2 + 4 * m * n - m**2)
    for coordinate in range(2):
        require(sp.expand(a * x[coordinate] + b * z[coordinate] - c * y[coordinate]) == 0,
                ("circuit", coordinate))
    norm = lambda v: v[0] ** 2 + v[1] ** 2
    require(sp.expand(c * norm(y) - a * norm(x) - b * norm(z) - K) == 0,
            "norm correction")


def hostile_sweep():
    decks = (
        ((1, 0), (0, 1)),
        tuple((i, 1 if i == 12 else 0) for i in range(1, 14)),
        tuple(((2 * i * i + i) % 9 - 4, (i * i + 3 * i) % 11 - 5)
              for i in range(1, 14)),
    )
    packets = 0
    score_checks = 0
    tariff_hits = 0
    for m in range(1, 25):
        for n in range(m + 1, 36):
            if gcd(m, n) != 1:
                continue
            x, y, z = branches(m, n)
            a, b, c = pythagorean_weights(m, n)
            K = correction(m, n)
            require((a * x[0] + b * z[0], a * x[1] + b * z[1]) ==
                    (c * y[0], c * y[1]), ("integer circuit", m, n))
            require(c * sum(t * t for t in y) -
                    a * sum(t * t for t in x) -
                    b * sum(t * t for t in z) == K, ("integer K", m, n))
            for deck in decks:
                for kappa in KAPPAS:
                    g_l, g_m, g_r = (score(v, deck, kappa) for v in (x, y, z))
                    require(c * g_m >= a * g_l + b * g_r + K,
                            ("Horn", m, n, kappa, deck))
                    A_l, A_r = max(0, -g_l), max(0, -g_r)
                    require(c * g_m >= K - a * A_l - b * A_r,
                            ("tariff", m, n, kappa, deck))
                    if a * A_l + b * A_r <= K:
                        require(g_m >= 0, ("tariff implication", m, n, kappa, deck))
                        tariff_hits += 1
                    score_checks += 1
            packets += 1
    return packets, score_checks, tariff_hits


def ap_controls():
    deck = tuple((i, 1 if i == 12 else 0) for i in range(1, 14))
    core = tuple(range(1, 12)) + (13,)

    # Positive packet: both outer gates fail, but their weighted debt is
    # smaller than K and the middle gate passes exactly at the tariff floor.
    m, n = 353, 356
    x, y, z = branches(m, n)
    a, b, c = pythagorean_weights(m, n)
    K = correction(m, n)
    scores = tuple(score(v, deck) for v in (x, y, z))
    require((x, y, z) == ((356, 359), (356, 1065), (353, 1062)), "AP children")
    require((a, b, c) == (2127, 251336, 251345), "AP weights")
    require(K == 1606017978, "AP K")
    require(scores == (-169080, 1066, -3893), ("AP scores", scores))
    bill = a * max(0, -scores[0]) + b * max(0, -scores[2])
    require(bill == 1338084208 and bill < K, ("AP bill", bill))
    require(c * scores[1] == K - bill == 267933770, "AP exact tariff")
    rows = tuple(physical_row(v, deck) for v in (x, y, z))
    expected_tails = (4631, 5337, 5298)
    for v, row, tail in zip((x, y, z), rows, expected_tails):
        scale = v[0]
        require(tuple(value for i, value in enumerate(row, start=1) if i != 12) ==
                tuple(scale * i for i in core), ("AP core", v))
        require(row[11] == tail == 12 * v[0] + v[1], ("AP tail", v))

    # Logical hostile: all three rows belong to the same proved-safe AP plane,
    # although all three determinant gates fail.
    root_children = branches(1, 2)
    root_scores = tuple(score(v, deck) for v in root_children)
    require(root_children == ((2, 3), (2, 5), (1, 4)), "root children")
    require(root_scores == (-3536, -5886, -4715), ("root scores", root_scores))
    root_tails = tuple(12 * v[0] + v[1] for v in root_children)
    require(root_tails == (27, 29, 16), root_tails)
    return scores, bill, K - bill, root_scores, root_tails


def complement_separation_control():
    body = (1, 2, 3, 4, 6, 12)
    L_scale = 14 * lcm(*body)
    numerators = (3, 5, 9, 11, 13, 15)
    denominators = (14,) * 6
    drifts = tuple(L_scale * numerator // denominator
                   for numerator, denominator in zip(numerators, denominators))
    row = body + (L_scale,) + drifts
    require(L_scale == 168, L_scale)
    require(drifts == (36, 60, 108, 132, 156, 180), drifts)
    require(len(set(row)) == 13 and all(value > 0 for value in row), row)
    for numerator in numerators:
        require(gcd(numerator, 14) == 1, ("reduced drift", numerator))
    quotient_speeds = tuple(14 * value // L_scale for value in drifts)
    require(quotient_speeds == numerators, quotient_speeds)

    d = branches(1, 2)[0]
    require(d == (2, 3), d)
    base = []
    for value in row:
        parity = value % 2
        base.append(((value - 3 * parity) // 2, parity))
    base = tuple(base)
    require(det(base[0], base[1]) == -1, "saturated deck witness")
    require(physical_row(d, base) == row, "base physical row")

    controls = (0, 1, 10, 100, 1000)
    observed = []
    for q in controls:
        deck = list(base)
        u, z = deck[-1]
        deck[-1] = (u + 3 * q, z - 2 * q)
        deck = tuple(deck)
        require(det(deck[0], deck[1]) == -1, ("saturation", q))
        require(physical_row(d, deck) == row, ("fixed physical row", q))
        D = delta(d, deck)
        G = score(d, deck)
        require(D == 270 + 13 * q, ("Delta family", q, D))
        require(G == -24557 - 1183 * q, ("score family", q, G))
        observed.append((q, D, G))
    return row, quotient_speeds, tuple(observed)


def main():
    symbolic_control()
    packets, score_checks, tariff_hits = hostile_sweep()
    ap_scores, ap_bill, ap_slack, root_scores, root_tails = ap_controls()
    fixed_row, quotient_speeds, unbounded_controls = complement_separation_control()

    semantic = sha256(repr((
        KAPPAS, packets, score_checks, tariff_hits,
        ap_scores, ap_bill, ap_slack, root_scores, root_tails,
        fixed_row, quotient_speeds, unbounded_controls,
    )).encode()).hexdigest()
    print("LRC14 BERGGREN HORN DEFECT TARIFF AND COMPLEMENT SEPARATION")
    print("SYMBOLIC circuit_and_norm_correction PASS")
    print(f"HOSTILE_SWEEP packets={packets} score_checks={score_checks} tariff_hits={tariff_hits} PASS")
    print(f"AP_POSITIVE scores={ap_scores} bill={ap_bill} strict_slack={ap_slack} PASS")
    print(f"AP_SAFE_GATE_HOSTILE scores={root_scores} tails={root_tails} PASS")
    print(f"THM3363_FIXED_ROW row={fixed_row} quotient_speeds={quotient_speeds}")
    print(f"UNBOUNDED_DECK_CONTROLS {unbounded_controls} PASS")
    print("FORMULAS Delta_q=270+13*q G_q=-24557-1183*q")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
