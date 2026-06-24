#!/usr/bin/env python3
"""
pi/unital/flower carrier atlas (codex, 2026-06-23).

Prompt ingredients:
  - 22/7 and cbrt(31) as low-complexity pi approximants,
  - flower/petal families from turning by 1/pi,
  - geometric unitals 2-(q^3+1, q+1, 1),
  - algebraic/functional-analytic unital maps preserving identity.

Purpose:
  Keep units and binary relations labelled before scalarizing.  The "22
  families" statement is true for a turn fraction 1/pi ~= 7/22.  If the angle
  is literally 1/pi radians, the turn fraction is 1/(2*pi^2), whose first
  denominator is near 20, not 22.  This is the same guardrail as a unital map:
  preserve the chosen unit object.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction


PI = math.pi
SPECIAL_H = {3, 7, 13, 21, 31, 43, 73, 91}
LRC_SHELLS = {14, 21, 27, 31, 33, 41, 43, 84, 91}


def continued_fraction(x: float, depth: int = 14):
    out = []
    y = x
    for _ in range(depth):
        a = math.floor(y)
        out.append(a)
        frac = y - a
        if abs(frac) < 1e-15:
            break
        y = 1.0 / frac
    return out


def convergents_from_cf(cf):
    conv = []
    p_nm2, p_nm1 = 0, 1
    q_nm2, q_nm1 = 1, 0
    for a in cf:
        p = a * p_nm1 + p_nm2
        q = a * q_nm1 + q_nm2
        conv.append(Fraction(p, q))
        p_nm2, p_nm1 = p_nm1, p
        q_nm2, q_nm1 = q_nm1, q
    return conv


def best_convergents(x: float, max_den: int = 500):
    rows = []
    for f in convergents_from_cf(continued_fraction(x, 24)):
        if f.denominator <= max_den:
            rows.append(f)
    return rows


def factor(n: int):
    d = 2
    out = []
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            out.append((d, e))
        d += 1 if d == 2 else 2
    if n > 1:
        out.append((n, 1))
    return out


def factor_str(n: int) -> str:
    fs = factor(n)
    if not fs:
        return "1"
    parts = []
    for p, e in fs:
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(parts)


def approximation_rows():
    rows = [
        ("22/7", 22 / 7, "continued-fraction pi convergent"),
        ("cuberoot(31)", 31 ** (1 / 3), "nearest-integer pi^3 carrier"),
        ("355/113", 355 / 113, "higher rational comparator"),
    ]
    return [(name, value, value - PI, abs(value - PI), note) for name, value, note in rows]


def flower_rows():
    turn_fraction = 1 / PI
    radian_fraction = 1 / (2 * PI * PI)
    rows = []
    for label, alpha in [
        ("1/pi turns", turn_fraction),
        ("1/pi radians", radian_fraction),
    ]:
        conv = best_convergents(alpha, 400)
        chosen = []
        for f in conv:
            q = f.denominator
            drift_turns = abs(q * alpha - round(q * alpha))
            chosen.append((f, drift_turns, 360 * drift_turns))
        rows.append((label, alpha, chosen))
    return rows


def unital_parameter_rows(q_min=2, q_max=14):
    rows = []
    for q in range(q_min, q_max + 1):
        points = q**3 + 1
        block = q + 1
        h = q * q - q + 1
        replication = q * q
        blocks = q * q * h
        tags = []
        if h in SPECIAL_H:
            tags.append(f"h={h}")
        if h in LRC_SHELLS:
            tags.append("LRC-shell")
        if points == 28:
            tags.append("C(8,2)")
        if points + 7 == 72:
            tags.append("65+7=72")
        if block == 7:
            tags.append("7-block")
        if h == 31:
            tags.append("pi^3-near")
        if h == 91:
            tags.append("C(14,2)")
        rows.append(
            {
                "q": q,
                "points": points,
                "block": block,
                "h": h,
                "replication": replication,
                "blocks": blocks,
                "factor_points": factor_str(points),
                "factor_h": factor_str(h),
                "tags": ",".join(tags) if tags else "-",
            }
        )
    return rows


def route_tournament():
    vertices = [
        "unit_preserving_guardrail",
        "C27_shell_transfer",
        "q3_unital_pair_frame",
        "q6_formal_unital_31",
        "pi_turn_22_flower",
        "literal_radian_correction",
        "q4_unital_65_plus_7_code_hint",
        "raw_pi_numerology",
    ]
    # tuple = (unit fidelity, LRC relevance, incidence strength, novelty, anti-scalar guard)
    scores = {
        "unit_preserving_guardrail": (5, 4, 2, 3, 5),
        "C27_shell_transfer": (4, 5, 3, 4, 5),
        "q3_unital_pair_frame": (3, 4, 5, 3, 4),
        "q6_formal_unital_31": (3, 3, 3, 5, 3),
        "pi_turn_22_flower": (3, 3, 2, 4, 3),
        "literal_radian_correction": (5, 2, 1, 3, 5),
        "q4_unital_65_plus_7_code_hint": (3, 2, 3, 4, 4),
        "raw_pi_numerology": (1, 1, 1, 2, 1),
    }
    order = {v: i for i, v in enumerate(vertices)}
    edges = set()
    vertex_scores = Counter()
    for a, b in itertools.combinations(vertices, 2):
        key_a = (scores[a], -order[a])
        key_b = (scores[b], -order[b])
        winner, loser = (a, b) if key_a >= key_b else (b, a)
        edges.add((winner, loser))
        vertex_scores[winner] += 1
        vertex_scores.setdefault(loser, vertex_scores[loser])

    def beats(a, b):
        return (a, b) in edges

    cycles3 = 0
    for a, b, c in itertools.combinations(vertices, 3):
        if (beats(a, b) and beats(b, c) and beats(c, a)) or (
            beats(a, c) and beats(c, b) and beats(b, a)
        ):
            cycles3 += 1
    hamiltonian_paths = 0
    for perm in itertools.permutations(vertices):
        if all(beats(perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            hamiltonian_paths += 1
    return vertices, scores, dict(vertex_scores), cycles3, hamiltonian_paths


def main():
    print("=" * 88)
    print("pi / unital / flower carrier atlas")
    print("=" * 88)

    print("\n[1] Low-complexity pi approximants")
    print(f"    pi = {PI:.15f}")
    print(f"    pi^3 = {PI**3:.15f}; pi^3 - 31 = {PI**3 - 31:.12f}")
    for name, value, signed_error, abs_error, note in approximation_rows():
        print(
            f"    {name:13s} value={value:.15f} "
            f"signed_error={signed_error:+.12e} abs_error={abs_error:.12e} :: {note}"
        )

    print("\n[2] Flower turn families and the unit guardrail")
    for label, alpha, conv in flower_rows():
        print(f"    {label}: alpha={alpha:.15f} turns per petal")
        for f, drift_turns, drift_degrees in conv[:7]:
            print(
                f"        {f.numerator:4d}/{f.denominator:<4d} "
                f"den={f.denominator:<4d} drift={drift_turns:.9e} turns "
                f"= {drift_degrees:.6f} deg"
            )
    print("    Readout: the 22-family belongs to 1/pi as a turn fraction (7/22).")
    print("             Literal 1/pi radians means alpha=1/(2*pi^2), whose first")
    print("             family is near denominator 20, not 22.")

    print("\n[3] Formal unital parameter rows 2-(q^3+1, q+1, 1)")
    print("    h=q^2-q+1, points=(q+1)h, replication=q^2, blocks=q^2 h")
    for row in unital_parameter_rows():
        print(
            f"    q={row['q']:2d}: points={row['points']:4d}={row['factor_points']:<8s} "
            f"block={row['block']:2d} h={row['h']:3d}={row['factor_h']:<6s} "
            f"r={row['replication']:3d} b={row['blocks']:5d} tags={row['tags']}"
        )

    print("\n[4] Interpretation")
    print("    q=3 gives the known unital pair-frame: 28=C(8,2), block size 4.")
    print("    q=6 is only a formal parameter row in finite-geometry terms, but it")
    print("    is the clean numerological bridge: h=31, points=7*31, block size 7,")
    print("    and cbrt(31) approximates pi because pi^3 is just above 31.")
    print("    q=4 gives 65 points; 65+7=72 is a code-length hint, not a construction.")
    print("    q=10 gives h=91=C(14,2), echoing the LRC14 pair-slot shell.")

    print("\n[5] Proof-carrier Tournament Analysis")
    vertices, scores, vertex_scores, cycles3, hpaths = route_tournament()
    ordered = sorted(vertices, key=lambda v: (-vertex_scores[v], v))
    print(f"    vertices={vertices}")
    print(f"    role_scores={scores}")
    print(f"    vertex_scores={dict(sorted(vertex_scores.items(), key=lambda kv: (-kv[1], kv[0])))}")
    print(f"    Hamiltonian_order={ordered}")
    print(f"    score_hist={dict(sorted(Counter(vertex_scores.values()).items()))}")
    print(f"    directed_3cycles={cycles3}")
    print(f"    SCC_sizes={[1 for _ in vertices]}")
    print(f"    Hamiltonian_path_count={hpaths}")

    print("\nVerdict:")
    print("    The shared word 'unital' is not itself a theorem, but it is a guardrail:")
    print("    preserve the unit of measurement / identity / pair-incidence before")
    print("    quotienting.  The 22-flower, cbrt(31), q=6 formal unital row, q=3")
    print("    pair frame, and q=4 length-72 hint are useful only as labelled carriers.")


if __name__ == "__main__":
    main()
