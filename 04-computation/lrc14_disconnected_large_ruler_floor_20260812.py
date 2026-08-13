#!/usr/bin/env python3
"""Exact large-ruler pair floor for the disconnected-low LRC14 branch.

Universe
--------
Every feasible upper-median context with ruler L>=4592 and every primitive
high channel P<Q<8P, at every common dilation g>=1.

For P>=4 the linked projective floor and the THM-3350 midpoint envelope give
one uniform analytic certificate already at g=1.  The finitely many P<4
channels use the contextwise midpoint bound for their tails and the exact
THM-3352 mass engine for all preceding heads.
"""
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TAIL = ROOT / "04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py"
MASS = ROOT / "04-computation/lrc_general_reflected_pair_mass_thm3352.py"
EXPECTED_TAIL = "78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9"
EXPECTED_MASS = "afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b"
TARGET = F(1, 294)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def filehash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gamma(k):
    # THM-3350 has the sharper value 1/2 for even k.  The displayed formula
    # is a valid upper bound for every parity and, crucially, decreases in k;
    # that monotonicity is what the all-Q envelope below uses.
    return F(k * k + 1, 2 * k * k)


def phase_floor(P, Q):
    return max(F(1, 105), F(1, 49) - F(12, 49 * P * Q))


def midpoint_error(P, Q, g, L, e, f):
    C = abs(Q * e - P * f)
    return (
        gamma(P) * F(e * P, g * L * P - e)
        + gamma(Q) * F(f * Q, g * L * Q - f)
        + F(C * (C // L + 1), 2 * g * g * L * P * Q)
    )


def envelope(P, L, e, f):
    """Worst-Q midpoint loss for Q in [P,8P], at the worst dilation g=1."""
    determinant = max((e + f) ** 2, F((8 * e + f) ** 2, 8))
    return (
        gamma(P) * F(e * P, L * P - e)
        + gamma(P) * F(f * P, L * P - f)
        + F(determinant, 2 * L * L)
        + F(e + f, 2 * L * P)
    )


def envelope_margin(P, lane):
    L, e, f = lane
    return F(1, 49) - F(12, 49 * P * (P + 1)) - envelope(P, L, e, f) - TARGET


def main():
    require(filehash(TAIL) == EXPECTED_TAIL, ("tail hash", filehash(TAIL)))
    require(filehash(MASS) == EXPECTED_MASS, ("mass hash", filehash(MASS)))
    T = load("large_floor_tail", TAIL)
    M = load("large_floor_mass", MASS)

    contexts = set()
    for body, L in T.SEL.MS.body_universe():
        cell, *_ = T.SEL.body_geometry(body, L)
        if L >= 4592:
            for e in body:
                for f in body:
                    if e != f:
                        contexts.add((L, cell, e, f))
    contexts = tuple(sorted(contexts))
    require(len(contexts) == 1514, ("contexts", len(contexts)))

    # Error decreases with L, so retain only the smallest realized large
    # ruler for each ordered endpoint lane.  For real x>=1 the derivative of
    # gamma_x * ex/(Lx-e) is negative; the determinant and phase losses also
    # decrease.  Thus the exact P=4 checks certify every real P>=4, and every
    # g>=1 because each midpoint-error term decreases in g.
    lane_minimum = {}
    for L, _, e, f in contexts:
        lane_minimum[e, f] = min(lane_minimum.get((e, f), 10**100), L)
    lanes = tuple(sorted((L, e, f) for (e, f), L in lane_minimum.items()))
    require(len(lanes) == 180, ("lanes", len(lanes)))
    envelope_weakest = min((envelope_margin(4, lane), lane) for lane in lanes)
    require(envelope_weakest[0] > 0, ("P>=4 envelope", envelope_weakest))

    primitive = tuple(
        (P, Q)
        for P in range(1, 4)
        for Q in range(P + 1, 8 * P)
        if gcd(P, Q) == 1 and P + Q >= 8
    )
    require(len(primitive) == 19 and (3, 5) in primitive, primitive)

    controls = 0
    tail_histogram = {}
    tail_weakest = None
    weakest = None
    semantic = sha256()
    for context in contexts:
        L, cell, e, f = context
        for P, Q in primitive:
            threshold = 1
            while phase_floor(P, Q) - midpoint_error(P, Q, threshold, L, e, f) <= TARGET:
                threshold += 1
                require(threshold <= 1000, ("tail ceiling", context, P, Q))
            tail_histogram[threshold] = tail_histogram.get(threshold, 0) + 1
            tail_row = (
                phase_floor(P, Q) - midpoint_error(P, Q, threshold, L, e, f) - TARGET,
                threshold, P, Q, L, cell, e, f,
            )
            if tail_weakest is None or tail_row < tail_weakest:
                tail_weakest = tail_row
            # Independent exact positive control at the worst dilation g=1.
            # The theorem for all later g uses the decreasing analytic error,
            # not a monotonicity claim about the physical mass itself.
            value = M.mass(L, cell, e, P, f, Q)
            row = (value - TARGET, 1, P, Q, L, cell, e, f, value)
            if weakest is None or row < weakest:
                weakest = row
            require(value >= TARGET, ("literal failure", row))
            semantic.update((repr((tail_row, row)) + "\n").encode())
            controls += 1

    require(controls == 28766, ("control count", controls))
    require(set(tail_histogram) == {1}, ("unexpected finite tail", tail_histogram))
    require(tail_weakest is not None and tail_weakest[0] > 0, ("tail weakest", tail_weakest))
    require(weakest is not None and weakest[0] > 0, ("weakest", weakest))
    print("LRC14 DISCONNECTED LARGE-RULER UNIVERSAL PAIR FLOOR")
    print("contexts", len(contexts), "rulers", len(set(x[0] for x in contexts)), "endpoint_lanes", len(lanes))
    print("P_ge_4_envelope_weakest", envelope_weakest)
    print("P_lt_4_primitive_channels", len(primitive), "literal_g1_controls", controls)
    print("tail_threshold_histogram", tuple(sorted(tail_histogram.items())))
    print("P_lt_4_tail_weakest", tail_weakest)
    print("weakest_literal_margin", weakest)
    print("semantic_sha256", semantic.hexdigest())
    print("conclusion=for every large-ruler context, primitive high P<Q<8P, and g>=1, physical overlap>=1/294")


if __name__ == "__main__":
    main()
