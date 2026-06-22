#!/usr/bin/env python3
"""
codex-S108: exact sub-14 LRC training atlas.

Convention: LRC(N) has N-1 moving speeds and target 1/N.

This extends HYP-2649 with the most useful current LRC14 proof tools:

* THM-523 q-witness covering reduction.
* THM-524 binding-pair exact maximization.
* HYP-2893 Goddyn-Wong/Jacobsthal acceleration atoms.
* HYP-2649 support-floor ladder.

The script is intentionally a proof-pattern atlas, not a brute-force proof of
all below-14 cases.  It asks what the known tools do on the smaller cases and
what first becomes sharp at N=14.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from math import gcd, lcm


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
    """THM-524 binding-pair candidate set in [0,1/2]."""

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
    """Exact measure of {t: ||v t|| >= 1/N for all v}."""

    bps: set[F] = {F(0), F(1)}
    for v in speeds:
        for m in range(v + 1):
            lo = F(N * m - 1, N * v)
            hi = F(N * m + 1, N * v)
            if 0 <= lo <= 1:
                bps.add(lo)
            if 0 <= hi <= 1:
                bps.add(hi)

    ordered = sorted(bps)
    total = F(0)
    for lo, hi in zip(ordered, ordered[1:]):
        mid = (lo + hi) / 2
        if all(norm_frac(F(v) * mid) >= F(1, N) for v in speeds):
            total += hi - lo
    return total


def covering_deficits(speeds: tuple[int, ...], N: int) -> tuple[int, ...]:
    return tuple(q for q in range(2, N + 1) if all(v % q for v in speeds))


def covering_repair(N: int, dropped: int) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    """AP_N with one speed dropped and the minimal q-witness covering repair."""

    base = tuple(v for v in range(1, N) if v != dropped)
    deficits = covering_deficits(base, N)
    repair = 1
    for q in deficits:
        repair = lcm(repair, q)
    speeds = tuple(sorted(base + (repair,)))
    assert primitive(speeds)
    assert not covering_deficits(speeds, N)
    return repair, deficits, speeds


def gw_accelerable_velocities(n: int) -> tuple[int, ...]:
    """HYP-2893/Tao-Goddyn-Wong nonunit interval test."""

    good = []
    for v in range(n // 2 + 1, n + 1):
        lo = n - v + 1
        hi = 2 * n - 2 * v + 1
        if lo <= hi and all(gcd(v, j) > 1 for j in range(lo, hi + 1)):
            good.append(v)
    return tuple(good)


def accelerated_row(n: int, accelerations: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(v for v in range(1, n + 1) if v not in accelerations) + sorted(2 * v for v in accelerations))


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def print_ap_and_gw_rows() -> None:
    print("=" * 96)
    print("A. AP tight rows and Goddyn-Wong acceleration atoms for N<=14")
    print("=" * 96)
    print("  N  target     AP_M       GW acc velocities        GW_M / safe measure")
    for N in range(3, 15):
        target = F(1, N)
        ap = tuple(range(1, N))
        ap_m, ap_t, _ = lonely_constant(ap)
        n = N - 1
        acc = gw_accelerable_velocities(n)
        if acc:
            gw = accelerated_row(n, acc)
            gw_m, gw_t, gw_active = lonely_constant(gw)
            gw_safe = safe_measure(gw, N) if max(gw) <= 80 else None
            gw_text = f"M={gw_m}, t={gw_t}, active={gw_active}, safe={gw_safe}"
        else:
            gw_text = "-"
        print(f" {N:2d}  {str(target):>7}   {str(ap_m):>7}   {str(acc):<24} {gw_text}")
    print()


def print_covering_repairs() -> list[dict[str, object]]:
    print("=" * 96)
    print("B. q-witness covering repairs of AP_N: drop d, add minimal repair w")
    print("=" * 96)
    print("Rows cover every q=2..N, so THM-523's easy witness is disabled.")
    print("The comparison line N=14 is included to show the first-open contrast.")
    print("This is an AP-drop repair slice, not a global covering-set minimization.")
    print()
    records: list[dict[str, object]] = []
    for N in range(4, 15):
        target = F(1, N)
        rows = []
        for dropped in range(1, N):
            repair, deficits, speeds = covering_repair(N, dropped)
            m, t, active = lonely_constant(speeds)
            sm = safe_measure(speeds, N)
            rows.append((m, sm, dropped, repair, deficits, speeds, t, active))
        rows.sort()
        best = rows[0]
        m, sm, dropped, repair, deficits, speeds, t, active = best
        margin = m - target
        records.append(
            {
                "N": N,
                "target": target,
                "M": m,
                "safe": sm,
                "dropped": dropped,
                "repair": repair,
                "deficits": deficits,
                "speeds": speeds,
                "t": t,
                "active": active,
                "margin": margin,
            }
        )
        print(
            f"N={N:2d} target={str(target):>5s} best_drop={dropped:2d}->w={repair:<4d} "
            f"deficits={deficits!s:<16s} M={str(m):>8s} margin={str(margin):>10s} "
            f"safe={str(sm):>12s} t={str(t):>8s} active={active}"
        )
        print(f"     row={speeds}")
    print()
    return records


def print_margin_story(records: list[dict[str, object]]) -> None:
    print("=" * 96)
    print("C. What smaller N teaches")
    print("=" * 96)
    print("Covering-repair rows are the finite lab where the q-witness route is disabled.")
    print("  N   N*(M-1/N)       active switch type           repair reading")
    for rec in records:
        N = int(rec["N"])
        scaled = F(N) * rec["margin"]  # type: ignore[arg-type]
        active = rec["active"]
        repair = rec["repair"]
        dropped = rec["dropped"]
        if N < 14:
            reading = "loose training row"
        else:
            reading = "LRC14 AP-drop covering comparison"
        print(f" {N:2d}   {str(scaled):>14s}   {str(active):<26s} drop {dropped}->w={repair}: {reading}")
    print()
    print("Observation: AP and GW rows are exact tilers, but they are not covering rows.")
    print("Once q-witness covering is forced, every AP-drop repair in this atlas is loose.")
    print("The AP-drop N=14 line is not a new phenomenon in kind; it is the first line where the")
    print("support-six signed tail can coexist with a still-small covering margin.")
    print()


def print_support_floor() -> None:
    print("=" * 96)
    print("D. Even-denominator support floor")
    print("=" * 96)
    for N in range(4, 16, 2):
        q = N // 2
        floor = q - 1
        tag = "below support-six" if N < 14 else ("first support-six" if N == 14 else "beyond")
        print(f"N={N:2d}: q={q}, support floor={floor}, {tag}")
    print()


def print_tournament_analysis() -> None:
    print("=" * 96)
    print("Tournament Analysis")
    print("=" * 96)
    vertices = [
        "q_witness_covering_reduction",
        "binding_pair_covering_margin",
        "AP_tight_equal_spacing",
        "Goddyn_Wong_acceleration_atom",
        "support_floor_q_minus_1",
        "raw_speed_search",
    ]
    print("Vertices are proof carriers rather than runners.")
    print("Observable: how much LRC proof structure survives before scalarizing.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("Challenged assumption: below-14 proof should be learned from raw finite")
    print("enumeration.  The useful training objects are q-cover deficits, binding")
    print("switches, acceleration windows, and the support floor.")


def main() -> None:
    print("codex-S108 exact sub-14 LRC covering/tiler training atlas")
    print("Convention: LRC(N) has N-1 speeds and target 1/N.")
    print_ap_and_gw_rows()
    records = print_covering_repairs()
    print_margin_story(records)
    print_support_floor()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
