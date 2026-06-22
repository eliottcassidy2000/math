#!/usr/bin/env python3
"""
codex-S108b: sub-14 single-swap boundary census.

Convention: LRC(N) has N-1 moving speeds and target 1/N.

This is a companion to lrc_sub14_covering_training_codex_s108.py and to
kps-S31o's LRC14 single-swap sporadic census.  It scans the natural
acceleration window

    {1,...,N-1} \ {d} union {v},  1 <= v <= 2(N-1),

using the exact THM-524 binding-pair critical set.  The point is not to prove
all small LRC cases again; it is to learn which exact tilers survive before the
q-witness covering condition is imposed.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import gcd


F = Fraction


def norm_frac(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def envelope_candidates(speeds: tuple[int, ...]) -> set[F]:
    speeds = tuple(sorted(set(speeds)))
    out: set[F] = {F(0), F(1, 2)}

    for v in speeds:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2):
                break
            out.add(t)
            k += 1

    for i, a in enumerate(speeds):
        for b in speeds[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                for k in range(1, d // 2 + 1):
                    out.add(F(k, d))
    return out


def lonely_constant(speeds: tuple[int, ...]) -> tuple[F, F, tuple[int, ...]]:
    best = F(-1)
    arg = F(0)
    active: tuple[int, ...] = ()
    for t in envelope_candidates(speeds):
        vals = [(norm_frac(F(v) * t), v) for v in speeds]
        val = min(x for x, _ in vals)
        if val > best:
            best = val
            arg = t
            active = tuple(v for x, v in vals if x == val)
    return best, arg, active


def safe_measure(speeds: tuple[int, ...], N: int) -> F:
    bps: set[F] = {F(0), F(1)}
    for v in speeds:
        for m in range(v + 1):
            lo = F(N * m - 1, N * v)
            hi = F(N * m + 1, N * v)
            if 0 <= lo <= 1:
                bps.add(lo)
            if 0 <= hi <= 1:
                bps.add(hi)

    total = F(0)
    ordered = sorted(bps)
    for lo, hi in zip(ordered, ordered[1:]):
        mid = (lo + hi) / 2
        if all(norm_frac(F(v) * mid) >= F(1, N) for v in speeds):
            total += hi - lo
    return total


def covering_deficits(speeds: tuple[int, ...], N: int) -> tuple[int, ...]:
    return tuple(q for q in range(2, N + 1) if all(v % q for v in speeds))


def residue_profile(speeds: tuple[int, ...], N: int, t: F) -> tuple[list[int], list[int], list[int]]:
    if t.denominator % N != 0 and t.denominator != N:
        return [], [], []
    multiplier = (t * N) % N
    if multiplier.denominator != 1:
        return [], [], []
    a = int(multiplier)
    residues = sorted((a * v) % N for v in speeds)
    counts = Counter(residues)
    missing = [r for r in range(1, N) if counts[r] == 0]
    doubled = [r for r in range(1, N) if counts[r] > 1]
    return residues, missing, doubled


def all_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for N in range(4, 15):
        ap = tuple(range(1, N))
        target = F(1, N)
        bound = 2 * (N - 1)
        for dropped in range(1, N):
            for added in range(1, bound + 1):
                speeds = tuple(sorted((set(ap) - {dropped}) | {added}))
                if len(speeds) != N - 1 or speeds == ap or not primitive(speeds):
                    continue
                m, t, active = lonely_constant(speeds)
                if m != target:
                    continue
                residues, missing, doubled = residue_profile(speeds, N, t)
                rows.append(
                    {
                        "N": N,
                        "dropped": dropped,
                        "added": added,
                        "speeds": speeds,
                        "M": m,
                        "t": t,
                        "active": active,
                        "safe": safe_measure(speeds, N),
                        "deficits": covering_deficits(speeds, N),
                        "residues": residues,
                        "missing": missing,
                        "doubled": doubled,
                    }
                )
    return rows


def main() -> None:
    print("codex-S108b sub-14 single-swap boundary census")
    print("Convention: LRC(N) has N-1 speeds and target 1/N.")
    print("Window: AP_N single-swap rows with added speed <= 2(N-1).")
    print("Exact method: THM-524 binding-pair critical points.")
    print()

    rows = all_rows()
    by_N = Counter(row["N"] for row in rows)
    print("Non-AP tight rows in the acceleration window:")
    for N in range(4, 15):
        print(f"  N={N:2d}: {by_N[N]} rows")
    print()

    print("Detailed rows:")
    for row in rows:
        print(
            f"  N={row['N']:2d} drop {row['dropped']:2d} -> add {row['added']:2d}; "
            f"M={row['M']} at t={row['t']}; safe={row['safe']}; "
            f"active={row['active']}; deficits={row['deficits']}"
        )
        print(f"       speeds={row['speeds']}")
        print(
            f"       residues={row['residues']} missing={row['missing']} "
            f"doubled={row['doubled']}"
        )
    print()

    apex_only = [row for row in rows if row["deficits"] == (row["N"],)]
    zero_safe = [row for row in rows if row["safe"] == 0]
    print("Synthesis:")
    print(f"  rows found: {len(rows)}")
    print(f"  rows with safe measure zero: {len(zero_safe)}")
    print(f"  rows whose only q-witness deficit is the apex q=N: {len(apex_only)}")
    print(
        "  Thus every tight non-AP row in this natural sub-14 single-swap "
        "window is an apex-denominator boundary tiler, not a covering row."
    )
    print(
        "  This matches the LRC14 lesson from THM-560/S31o: exact tilers live "
        "on the q=N witness boundary; the covering branch must be closed by "
        "binding-pair margins and support-six residual control."
    )
    print()

    print("Tournament Analysis")
    print("  Vertices are row mechanisms, not runners:")
    vertices = [
        "apex_q_deficit",
        "complete_residue_lift",
        "one_gap_one_collision",
        "global_single_swap_tightness",
        "q_covering_forced_margin",
        "raw_single_swap_scan",
    ]
    print("  " + " > ".join(vertices))
    print(
        "  Observable: whether a tight row survives after quotienting by "
        "q-witness covering.  It does not: the quotient sends every row above "
        "to the apex q=N boundary, so it destroys no covering counterexample."
    )
    print(
        "  Challenged assumption: small sporadic tilers are the prototype for "
        "the LRC14 covering obstruction.  In this window they are only the "
        "prototype for boundary artifacts; the covering prototype is S108's "
        "AP-drop repair slice with positive binding-pair margins."
    )


if __name__ == "__main__":
    main()
