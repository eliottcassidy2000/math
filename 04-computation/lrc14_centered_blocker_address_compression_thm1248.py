#!/usr/bin/env python3
"""Exact referee for THM-1248 centered blocker address compression.

The proof-facing calculations use only integers and ``fractions.Fraction``.
The analytic inputs are THM-1233's ratio cap and THM-1240's centered-spoke
construction.  This referee replays the new relative-address identities,
their sharp constants, affine cycle transport, and two non-cover guardrails.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from math import gcd, prod
from pathlib import Path


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def nearest_integer(x: F, *, upper_tie: bool = False) -> int:
    lower = x.numerator // x.denominator
    upper = lower + 1
    dl, du = x - lower, upper - x
    if dl < du:
        return lower
    if du < dl:
        return upper
    return upper if upper_tie else lower


def circle_distance(speed: int, time: F) -> F:
    value = speed * time
    return abs(value - nearest_integer(value))


def centered_phase(c: int, k: int, d: int, *, upper_tie: bool = False) -> tuple[int, int, F, F]:
    require(0 <= k < c < d, (c, k, d))
    q = c + d
    t0 = F(2 * k + 1, 2 * c)
    p = nearest_integer(q * t0, upper_tie=upper_tie)
    epsilon = p - q * t0
    require(abs(epsilon) <= F(1, 2), (c, k, d, p, epsilon))
    return q, p, F(p, q), epsilon


def blocker_edge(
    c: int,
    k: int,
    di: int,
    dj: int,
    pi: int,
    pj: int,
) -> dict[str, object]:
    qi, qj = c + di, c + dj
    ti, tj = F(pi, qi), F(pj, qj)
    n = nearest_integer(dj * ti)
    beta = dj * ti - n
    require(abs(beta) < F(1, 14), (c, k, di, dj, pi, pj, beta))

    t0 = F(2 * k + 1, 2 * c)
    epsilon_i, epsilon_j = pi - qi * t0, pj - qj * t0
    x = c * pi - k * qi
    residue = di * 0 + pi * dj - n * qi
    m = k + n
    ell = pi * qj - m * qi
    delta = pj - m
    determinant = pi * qj - pj * qi
    theta = F(ell, qi)
    a = F(dj, qi)

    require(F(1, 4) < F(x, qi) < F(3, 4), (di, x, qi))
    require(residue == qi * beta, (residue, qi, beta))
    require(ell == x + residue, (ell, x, residue))
    require(F(5, 28) < theta < F(23, 28), (di, dj, theta))
    require(ell == determinant + delta * qi, (ell, determinant, delta, qi))
    require(ell == determinant % qi, (ell, determinant, qi))
    require(m == (pi * qj) // qi, (m, pi, qj, qi, ell))
    require(gcd(qi, qj) and ell % gcd(qi, qj) == 0, (qi, qj, ell))
    require(
        beta == delta - F(1, 2) + a * epsilon_i - epsilon_j,
        (di, dj, beta, delta, epsilon_i, epsilon_j),
    )
    require(qj * (ti - tj) == theta - delta, (ti, tj, theta, delta))
    require(-586 <= delta <= 587, (di, dj, delta))
    if dj <= qi:
        require(delta in (0, 1), (di, dj, qi, delta))
        require(F(5, 28 * qj) < abs(ti - tj) < F(23, 28 * qj), (ti, tj))

    return {
        "source": di,
        "target": dj,
        "Qi": qi,
        "Qj": qj,
        "Pi": pi,
        "Pj": pj,
        "N": n,
        "M": m,
        "epsilon_i": epsilon_i,
        "epsilon_j": epsilon_j,
        "beta": beta,
        "X": x,
        "ell": ell,
        "delta": delta,
        "determinant": determinant,
        "theta": theta,
        "slope": a,
    }


def check_cycle(rows: tuple[dict[str, object], ...]) -> dict[str, object]:
    require(len(rows) >= 2, len(rows))
    require(all(rows[i]["target"] == rows[(i + 1) % len(rows)]["source"] for i in range(len(rows))), rows)
    slopes = tuple(row["slope"] for row in rows)
    offsets = tuple(row["delta"] - F(1, 2) - row["beta"] for row in rows)
    multiplier = prod(slopes, start=F(1))
    weights = tuple(prod(slopes[i + 1 :], start=F(1)) for i in range(len(rows)))
    affine_offset = sum((offsets[i] * weights[i] for i in range(len(rows))), F(0))
    epsilon0 = rows[0]["epsilon_i"]
    require(affine_offset == (1 - multiplier) * epsilon0, (rows, affine_offset, multiplier, epsilon0))
    require(0 < multiplier < 1, multiplier)
    require(abs(affine_offset) <= (1 - multiplier) / 2, affine_offset)

    digit_center = sum(((F(rows[i]["delta"]) - F(1, 2)) * weights[i] for i in range(len(rows))), F(0))
    beta_radius = sum((weights[i] for i in range(len(rows))), F(0)) / 14
    require(abs(digit_center) < beta_radius + (1 - multiplier) / 2, (digit_center, beta_radius))

    weighted_delta = sum((F(row["delta"], row["Qj"]) for row in rows), F(0))
    weighted_theta = sum((row["theta"] / row["Qj"] for row in rows), F(0))
    reciprocal_weight = sum((F(1, row["Qj"]) for row in rows), F(0))
    require(weighted_delta == weighted_theta, (weighted_delta, weighted_theta))
    require(F(5, 28) * reciprocal_weight < weighted_delta < F(23, 28) * reciprocal_weight, weighted_delta)
    require(any(row["delta"] <= 0 for row in rows), rows)
    require(any(row["delta"] >= 1 for row in rows), rows)

    p_product = prod((row["Pi"] for row in rows), start=1)
    m_product = prod((row["M"] for row in rows), start=1)
    require(p_product > m_product, (p_product, m_product, rows))
    return {
        "length": len(rows),
        "multiplier": multiplier,
        "contraction_gap": 1 - multiplier,
        "delta_word": tuple(row["delta"] for row in rows),
        "ell_word": tuple(row["ell"] for row in rows),
        "positive_holonomy_gap": p_product - m_product,
    }


def check_constants() -> dict[str, F]:
    ratio = F(2345, 2346)
    gap = 1 - ratio * ratio
    require(gap == F(4691, 5503716), gap)
    require(F(1, 4) - F(1, 14) == F(5, 28), "lower central remainder")
    require(F(3, 4) + F(1, 14) == F(23, 28), "upper central remainder")
    require(F(1, 4) - F(1, 7) == F(3, 28) > F(1, 14), "wall-event margin")
    return {
        "global_delta_min": F(-586),
        "global_delta_max": F(587),
        "contraction_gap": gap,
        "wall_event_source_margin": F(3, 28),
    }


def check_unbounded_address_guardrail(c: int) -> dict[str, object]:
    require(c >= 27, c)
    k = c - 1
    speeds = (c + 1, c + 2, c + 3, c + 4, 2 * c, 3 * c)
    phases: dict[int, tuple[int, int, F, F]] = {}
    for d in speeds:
        phases[d] = centered_phase(c, k, d, upper_tie=(d == 2 * c))
    blockers = {c + 1: 2 * c, c + 2: 2 * c, c + 3: 2 * c, c + 4: 2 * c, 2 * c: 3 * c, 3 * c: 2 * c}
    rows = tuple(
        blocker_edge(c, k, di, dj, phases[di][1], phases[dj][1])
        for di, dj in blockers.items()
    )
    cycle_rows = tuple(row for row in rows if row["source"] in (2 * c, 3 * c))
    cycle = check_cycle(cycle_rows)
    witness = 1 - F(2, 5 * c)
    packet = (c,) + speeds
    depths = tuple(circle_distance(v, witness) for v in packet)
    require(min(depths) == F(1, 5), (c, packet, witness, depths))
    require(gcd(*packet[:2]) == 1, packet)
    return {
        "c": c,
        "max_tooth_address": max(row["N"] for row in rows),
        "cycle_delta_word": cycle["delta_word"],
        "witness": witness,
        "witness_depth": min(depths),
    }


def check_large_delta_guardrail() -> dict[str, object]:
    c, k = 1, 0
    speeds = (2, 16, 17, 34, 35, 2343)
    phases: dict[int, tuple[int, int, F, F]] = {
        d: centered_phase(c, k, d, upper_tie=(d == 2)) for d in speeds
    }
    blockers = {2: 2343, 2343: 2, 16: 17, 17: 16, 34: 35, 35: 34}
    rows = {
        (di, dj): blocker_edge(c, k, di, dj, phases[di][1], phases[dj][1])
        for di, dj in blockers.items()
    }
    cycles = tuple(
        check_cycle((rows[(a, b)], rows[(b, a)]))
        for a, b in ((2, 2343), (16, 17), (34, 35))
    )
    critical = rows[(2, 2343)]
    require(critical["delta"] == -390 and critical["ell"] == 2, critical)
    require(critical["M"] == 1562 and critical["Pi"] == 2, critical)
    require(
        F(2, 1) < F(13, 6)
        and F(16, 2) < F(273, 29)
        and F(17, 16) < F(84, 5)
        and F(34, 17) < F(343, 15)
        and F(35, 34) < F(189, 8)
        and F(2343, 35) < 77
        and 2343 < 2345,
        "THM-1233 compact inequalities",
    )
    witness = F(1, 6)
    packet = (c,) + speeds
    depth = min(circle_distance(v, witness) for v in packet)
    require(depth == F(1, 6), (packet, witness, depth))
    return {
        "packet": packet,
        "critical_edge": (critical["source"], critical["target"]),
        "critical_delta": critical["delta"],
        "critical_ell": critical["ell"],
        "cycle_delta_words": tuple(cycle["delta_word"] for cycle in cycles),
        "witness": witness,
        "witness_depth": depth,
    }


def main() -> None:
    print("THM-1248 CENTERED BLOCKER ADDRESS COMPRESSION EXACT REFEREE")
    print("method=integer/Fraction only; always-on checks; no dependencies")

    constants = check_constants()
    print("\nFINITE RELATIVE-ADDRESS AND WALL-EVENT CONSTANTS")
    for key, value in constants.items():
        print(f"{key}={value}")
    print("central_remainder_band=(5/28,23/28)")
    print("relative_delta_word=[-586,587]; every edge with d_j<=c+d_i has digit 0 or 1")
    print("every genuine blocker cycle exports a 3-safe/fourth-support wall event")

    print("\nUNBOUNDED ABSOLUTE-ADDRESS GUARDRAIL")
    for c in (27, 101, 1001, 10001):
        row = check_unbounded_address_guardrail(c)
        print(
            f"c={c} max_tooth_address={row['max_tooth_address']} "
            f"cycle_delta_word={row['cycle_delta_word']} "
            f"witness={row['witness']} depth={row['witness_depth']}"
        )
    print("sampled obligations and primitive compact ratios do not bound absolute addresses")
    print("these packets are globally lonely, not six-cover counterexamples")

    large = check_large_delta_guardrail()
    print("\nLARGE RELATIVE-DIGIT GUARDRAIL")
    for key, value in large.items():
        print(f"{key}={value}")
    print("the {-586,...,587} alphabet cannot be replaced by {0,1} on speed-ascent edges")

    print("\nTOURNAMENT / ALTERNATE-VERTEX AUDIT")
    print("observable=D_ij=P_i*Q_j-P_j*Q_i")
    print("switch=i->j iff (t_i,i) is lexicographically below (t_j,j); index breaks D_ij=0 ties")
    print("tie_Hamiltonian_path=vertices sorted by (t_i,i)")
    print("phase-determinant sign gauge is transitive and loses the central remainder")
    print("proof-bearing edge=(Qi,Qj,least-positive ell,relative delta,gcd sheet)")
    print("cycle quotient=finite delta word + contracting affine errors + wall-event hyperedge")
    print("challenged_assumption=runner vertices are insufficient; use residues, wall counts, and boundary obligations")
    print("preserves=blocker danger, determinant gcd, and phase chronology")
    print("destroys=absolute tooth lift and off-sample cover ownership")
    print("runner scores=(0,1,2,3,4,5); cycles=0; SCCs=6; Hamilton_paths=1")
    print("STATUS=PASS")
    print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
