#!/usr/bin/env python3
"""HYP-2790/S76: bridge period-max pressure to Boolean/type slack.

This is a focused scout, not a replacement for the incoming S6/S7 period-max
scripts.  It overlays HYP-2791's low-depth missed-sector functional

    Phi_low = 21*T1 + 57*T2sep + 2*T2adj

on the bounded bases that are already known to matter:

* the highest-plateau bounded bases with 0 in B and max(B) <= 14;
* the S6 broad-scan worst row;
* the S21 skipped-row frontier with large endpoint periods.

This script is coordinate-only: exact period-max values are delegated to the
incoming S6/S7 scripts.  It reports plateau, margins, endpoint period, and
Boolean features, so the period-max rows can be cross-read against Boolean/type
coordinates without pretending to close their period maxima here.

Tournament vertices are proof lenses rather than runners: direct period scan,
AP/dilation filter, Boolean low-depth slack, q0 slack, skipped-period audit,
and canonical-vs-strict cap normalization.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd, sqrt

sys.stdout.reconfigure(line_buffering=True)


CAPS_CANONICAL = {
    8: F(2243, 5880),
    9: F(1979, 4004),
    10: F(55, 91),
}

# Older period-max scripts use the proved floor 4/7 for k=10; it is stricter
# than the canonical cap 55/91.  Report both ratios so cap mismatches cannot
# masquerade as proof progress.
CAPS_STRICT = dict(CAPS_CANONICAL)
CAPS_STRICT[10] = F(4, 7)

LOW_FEATURES = ((1, (1,)), (2, (1, 1)), (2, (2,)))
LOW_COEFFS = (F(21), F(57), F(2))

FORCED_ROWS = {
    8: {
        (0, 4, 6, 8, 10, 12, 14),  # S6 broad checked worst ratio row.
        (0, 8, 9, 10, 11, 12, 13),
        (0, 1, 2, 10, 11, 12, 13),
        (0, 1, 3, 5, 9, 11, 13),
    },
    9: {
        (0, 1, 3, 5, 7, 9, 11, 13),
        (0, 1, 3, 5, 6, 7, 9, 11),
        (0, 6, 7, 8, 9, 10, 11, 12),
        (0, 8, 9, 10, 11, 12, 13, 14),
    },
    10: {
        (0, 2, 4, 6, 8, 10, 11, 12, 14),
        (0, 1, 2, 3, 4, 5, 6, 7, 13),
        (0, 2, 4, 6, 8, 10, 12, 13, 14),
        (0, 1, 2, 3, 4, 5, 6, 7, 8),
        (0, 2, 4, 6, 7, 8, 10, 12, 14),
    },
}

TOP_PLATEAU_ROWS = {
    # From 05-knowledge/results/lrc_periodmax_general_macmini_0621s6.out.
    8: [
        (0, 2, 4, 6, 8, 10, 12),
        (0, 1, 2, 3, 4, 5, 6),
        (0, 2, 3, 4, 5, 6, 8),
        (0, 2, 4, 5, 6, 8, 10),
        (0, 1, 3, 4, 5, 6, 9),
        (0, 1, 2, 3, 4, 6, 12),
        (0, 2, 3, 4, 6, 8, 12),
        (0, 1, 2, 4, 5, 6, 10),
        (0, 1, 4, 5, 6, 9, 10),
        (0, 2, 3, 5, 6, 8, 11),
        (0, 3, 4, 5, 6, 8, 9),
        (0, 3, 5, 6, 8, 9, 11),
    ],
    # From 05-knowledge/results/lrc_periodmax_general_macmini_0621s6.out.
    9: [
        (0, 2, 4, 6, 8, 10, 12, 14),
        (0, 1, 2, 3, 4, 5, 6, 7),
        (0, 2, 3, 4, 5, 6, 7, 8),
        (0, 2, 4, 6, 7, 8, 10, 12),
        (0, 1, 3, 4, 5, 6, 7, 9),
        (0, 1, 3, 5, 7, 9, 11, 13),
        (0, 1, 2, 3, 4, 5, 6, 14),
        (0, 1, 2, 3, 4, 5, 6, 8),
        (0, 2, 4, 5, 6, 7, 8, 10),
        (0, 1, 2, 4, 5, 6, 7, 10),
        (0, 3, 4, 5, 6, 7, 8, 9),
        (0, 2, 4, 5, 6, 8, 10, 12),
    ],
    # From 05-knowledge/results/lrc_periodmax_skipped_audit_thread5.out.
    10: [
        (0, 1, 2, 3, 4, 5, 6, 7, 8),
        (0, 2, 4, 6, 7, 8, 10, 12, 14),
        (0, 1, 2, 3, 4, 5, 6, 7, 9),
        (0, 1, 2, 3, 4, 5, 6, 7, 10),
        (0, 2, 3, 4, 6, 8, 10, 12, 14),
    ],
}


def fmt(x: F | None) -> str:
    if x is None:
        return "NA"
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def ffloat(x: F | None) -> str:
    if x is None:
        return "NA"
    return f"{float(x):.6f}"


def lcm(a: int, b: int) -> int:
    return abs(a * b) // gcd(a, b)


def frac_part(x: F) -> F:
    return x - (x.numerator // x.denominator)


def sector_of(x: F) -> int:
    return int(7 * frac_part(x))


def mask_runs(mask: int) -> tuple[int, ...]:
    bits = [(mask >> i) & 1 for i in range(6)]
    if sum(bits) == 0:
        return ()
    if all(bits):
        return (6,)
    candidates = []
    for seq in (bits, list(reversed(bits))):
        for shift in range(6):
            row = seq[shift:] + seq[:shift]
            if row[-1] == 0 and row[0] == 1:
                lens = []
                i = 0
                while i < 6:
                    if row[i]:
                        j = i
                        while j < 6 and row[j]:
                            j += 1
                        lens.append(j - i)
                        i = j
                    else:
                        i += 1
                candidates.append(tuple(lens))
    return min(candidates)


def breakpoints(E: tuple[int, ...]) -> list[F]:
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(7 * e + 1):
            bps.add(F(m, 7 * e))
    return sorted(bps)


def base_data(E: tuple[int, ...]):
    bps = breakpoints(E)
    p0 = F(0)
    p1 = F(0)
    atoms: dict[int, F] = defaultdict(F)
    arcs = []
    cur = {}
    for a, b in zip(bps, bps[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        hit = {sector_of(e * mid) for e in E}
        missing = [s for s in range(1, 7) if s not in hit]
        if len(hit) == 7:
            p0 += b - a
        if len(missing) == 1:
            p1 += b - a
        mask = 0
        for s in range(1, 7):
            if s not in hit:
                mask |= 1 << (s - 1)
        atoms[mask] += b - a
        one_miss = missing[0] if len(missing) == 1 else None
        for j in range(1, 7):
            active = one_miss == j
            if active and j not in cur:
                cur[j] = a
            if (not active) and j in cur:
                arcs.append((j, cur.pop(j), a))
    for j, a in cur.items():
        arcs.append((j, a, F(1)))
    profile: dict[tuple[int, tuple[int, ...]], F] = defaultdict(F)
    for mask, val in atoms.items():
        profile[(mask.bit_count(), mask_runs(mask))] += val
    return p0 + p1 / 7, dict(profile), arcs


def plateau_only(E: tuple[int, ...]) -> F:
    bps = breakpoints(E)
    p0 = F(0)
    p1 = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        hit = {sector_of(e * mid) for e in E}
        missing = [s for s in range(1, 7) if s not in hit]
        if len(hit) == 7:
            p0 += b - a
        if len(missing) == 1:
            p1 += b - a
    return p0 + p1 / 7


def low_phi(profile: dict[tuple[int, tuple[int, ...]], F]) -> F:
    return sum(c * profile.get(f, F(0)) for c, f in zip(LOW_COEFFS, LOW_FEATURES))


def q0(profile: dict[tuple[int, tuple[int, ...]], F]) -> F:
    return profile.get((0, ()), F(0))


def Sj_raw(t: F, j: int) -> F:
    t = frac_part(t)
    a = F(j, 7)
    b = F(j + 1, 7)
    if t <= a:
        return -t / 7
    if t <= b:
        return -a / 7 + F(6, 7) * (t - a)
    return -a / 7 + F(6, 7) * (b - a) - F(1, 7) * (t - b)


def mean_sj(j: int) -> F:
    a = F(j, 7)
    b = F(j + 1, 7)
    pts = [F(0), a, b, F(1)]
    vals = [Sj_raw(x, j) for x in pts]
    return sum((vals[i] + vals[i + 1]) / 2 * (pts[i + 1] - pts[i]) for i in range(3))


MEAN = {j: mean_sj(j) for j in range(1, 7)}


def Sc(t: F, j: int) -> F:
    return Sj_raw(t, j) - MEAN[j]


def endpoint_period(B: tuple[int, ...]) -> int:
    out = 1
    for e in B:
        if e:
            out = lcm(out, e)
    return 7 * out


def period_max(arcs, P: int) -> tuple[F, int]:
    best = F(-10)
    best_w = -1
    for w in range(15, 15 + P):
        val = F(0)
        for j, a, b in arcs:
            val += Sc(w * b, j) - Sc(w * a, j)
        if val > best:
            best = val
            best_w = w
    return best, best_w


def corr(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def tournament(order_scores: dict[str, tuple[int, ...]]):
    names = list(order_scores)
    wins = {n: 0 for n in names}
    edge = {}
    for i, j in combinations(range(len(names)), 2):
        a, b = names[i], names[j]
        winner = a if order_scores[a] >= order_scores[b] else b
        edge[(a, b)] = winner
        wins[winner] += 1
    cycles = 0
    for a, b, c in combinations(names, 3):
        out = {a: 0, b: 0, c: 0}
        for u, v in ((a, b), (a, c), (b, c)):
            key = (u, v) if (u, v) in edge else (v, u)
            out[edge[key]] += 1
        if sorted(out.values()) == [1, 1, 1]:
            cycles += 1
    hp = 0
    for perm in permutations(names):
        ok = True
        for i in range(len(perm) - 1):
            key = (perm[i], perm[i + 1])
            if key in edge:
                ok &= edge[key] == perm[i]
            else:
                ok &= edge[(perm[i + 1], perm[i])] == perm[i]
        hp += int(ok)
    return wins, cycles, hp


def top_plateau_bases(k: int, topn: int) -> list[tuple[int, ...]]:
    if k in TOP_PLATEAU_ROWS:
        return TOP_PLATEAU_ROWS[k][:topn]
    rows = []
    for combo in combinations(range(1, 15), k - 2):
        B = (0,) + combo
        rows.append((plateau_only(B), B))
    rows.sort(reverse=True, key=lambda x: x[0])
    return [B for _, B in rows[:topn]]


def is_dilation_of_ap(B: tuple[int, ...]) -> bool:
    if len(B) <= 1:
        return True
    step = B[1]
    if step <= 0:
        return False
    return B == tuple(step * i for i in range(len(B)))


def scan(k: int, topn: int, period_limit: int):
    exact_limit = 0
    ap = tuple(range(k - 1))
    ap_plat, ap_profile, _ = base_data(ap)
    ap_phi = low_phi(ap_profile)
    ap_q0 = q0(ap_profile)

    candidates = set(top_plateau_bases(k, topn))
    candidates.add(ap)
    if 2 * (k - 2) <= 14:
        candidates.add(tuple(2 * i for i in range(k - 1)))
    candidates.update(FORCED_ROWS.get(k, set()))

    rows = []
    for B in sorted(candidates):
        plat, profile, arcs = base_data(B)
        margin_canonical = CAPS_CANONICAL[k] - plat
        margin_strict = CAPS_STRICT[k] - plat
        P = endpoint_period(B)
        phi_gap = low_phi(profile) - ap_phi
        q0_gap = ap_q0 - q0(profile)
        ap_type = profile == ap_profile
        dil_ap = is_dilation_of_ap(B)
        pm = None
        w = None
        ratio_canonical = None
        ratio_strict = None
        if P <= exact_limit:
            pm, w = period_max(arcs, P)
            if margin_canonical > 0:
                ratio_canonical = pm / margin_canonical
            if margin_strict > 0:
                ratio_strict = pm / margin_strict
        rows.append(
            {
                "B": B,
                "plat": plat,
                "margin_canonical": margin_canonical,
                "margin_strict": margin_strict,
                "P": P,
                "pm": pm,
                "w": w,
                "ratio_canonical": ratio_canonical,
                "ratio_strict": ratio_strict,
                "phi_gap": phi_gap,
                "q0_gap": q0_gap,
                "ap_type": ap_type,
                "dil_ap": dil_ap,
            }
        )
    rows.sort(
        key=lambda r: (
            r["ratio_strict"] is not None,
            float(r["ratio_strict"] or r["ratio_canonical"] or F(-1)),
            float(r["plat"]),
        ),
        reverse=True,
    )
    return ap_plat, ap_phi, ap_q0, rows


def row_line(row) -> str:
    ratio_c = row["ratio_canonical"]
    ratio_s = row["ratio_strict"]
    return (
        f"B={row['B']!s:<34} P={row['P']:<7} Plat={ffloat(row['plat'])} "
        f"mc={ffloat(row['margin_canonical'])} ms={ffloat(row['margin_strict'])} "
        f"pm={fmt(row['pm']):>10} w={str(row['w']):<7} "
        f"ratioC={ffloat(ratio_c):>9} ratioS={ffloat(ratio_s):>9} "
        f"phi_gap={fmt(row['phi_gap']):>12} q0_gap={fmt(row['q0_gap']):>10} "
        f"ap_type={str(row['ap_type']):<5} dil_ap={row['dil_ap']}"
    )


def main():
    topn = 12
    period_limit = 0
    print("HYP-2790/S76 period-max vs Boolean/type bridge")
    print(f"topn_by_plat={topn}; period_limit={period_limit} (coordinate-only overlay)")
    print("features: Phi_low=21*T1+57*T2sep+2*T2adj on bounded base B")
    print("caps: canonical k10=55/91; strict comparison also reports k10=4/7")
    print("exact period maxima are delegated to S6/S7 scripts")

    all_rows = []
    for k in (8, 9, 10):
        print("\n" + "=" * 96)
        print(f"scanning k={k}")
        ap_plat, ap_phi, ap_q0, rows = scan(k, topn, period_limit)
        checked = [r for r in rows if r["pm"] is not None]
        skipped = [r for r in rows if r["pm"] is None]
        print(
            f"k={k}: AP base={tuple(range(k-1))}; AP Plat={fmt(ap_plat)}; "
            f"AP Phi_low={fmt(ap_phi)}; AP q0={fmt(ap_q0)}"
        )
        print(f"  candidates={len(rows)} checked={len(checked)} skipped_large_period={len(skipped)}")
        print("  checked rows by strict/canonical pressure:")
        for row in checked[:14]:
            print("    " + row_line(row))
        if skipped:
            print("  skipped rows: Boolean/type coordinates without exact period-max")
            shown = sorted(skipped, key=lambda r: r["margin_strict"])
            if checked:
                shown = shown[:10]
            for row in shown:
                print("    " + row_line(row))

        non_ap_checked = [r for r in checked if not r["ap_type"]]
        xs = [float(r["ratio_strict"] or r["ratio_canonical"]) for r in non_ap_checked]
        phi = [float(r["phi_gap"]) for r in non_ap_checked]
        q0s = [float(r["q0_gap"]) for r in non_ap_checked]
        print(f"  corr(nonAP checked ratio, phi_gap)={corr(xs, phi)}")
        print(f"  corr(nonAP checked ratio, q0_gap)={corr(xs, q0s)}")

        best_non_ap = non_ap_checked[0] if non_ap_checked else None
        if best_non_ap is not None:
            print(
                "  best non-AP checked pressure: "
                f"ratioS={ffloat(best_non_ap['ratio_strict'])}, "
                f"ratioC={ffloat(best_non_ap['ratio_canonical'])}, "
                f"phi_gap={fmt(best_non_ap['phi_gap'])}, B={best_non_ap['B']}"
            )
        all_rows.extend((k, row) for row in rows)

    frontier_non_ap = [(k, r) for k, r in all_rows if not r["ap_type"]]
    print("\n" + "=" * 96)
    print("Synthesis")
    print("- This is a coordinate overlay.  Exact period maxima remain in the S6/S7")
    print("  period-max scripts and are not recomputed here.")
    print("- AP/dilation rows are zero-gap orbit labels and still explain the largest")
    print("  pressure rows that S6/S7 already close.")
    if frontier_non_ap:
        min_phi = min(r["phi_gap"] for _, r in frontier_non_ap)
        min_q0 = min(r["q0_gap"] for _, r in frontier_non_ap)
        print(f"- global non-AP frontier min Phi_low gap = {fmt(min_phi)}")
        print(f"- global non-AP frontier min q0 gap = {fmt(min_q0)}")
        for kk in (8, 9, 10):
            sub = [r for k, r in frontier_non_ap if k == kk]
            if sub:
                print(
                    f"- k={kk}: min Phi_low gap={fmt(min(r['phi_gap'] for r in sub))}; "
                    f"min q0 gap={fmt(min(r['q0_gap'] for r in sub))}"
                )
    print("- Therefore the naive transfer of HYP-2791's three-term Phi_low cut to the")
    print("  one-far base ledger is false at k=8.  The safer bridge is AP/dilation")
    print("  filtering plus q0/cover slack on bases, with Phi_low reserved for final-row")
    print("  Boolean laws or for a size-shifted k>=9 subledger.")
    print("- k=10 cap normalization must be kept explicit: 4/7 strict floor is harder")
    print("  than the canonical cap 55/91 used by the final LRC cap ledger.")

    scores = {
        "periodmax_direct_scan": (6, 6, 2, 6, 5),
        "ap_dilation_filter": (5, 6, 6, 5, 5),
        "skipped_period_audit": (4, 4, 5, 5, 4),
        "q0_slack": (4, 4, 6, 4, 4),
        "boolean_low_depth_slack": (3, 4, 5, 3, 5),
        "cap_normalization_guard": (3, 3, 6, 5, 5),
        "raw_plateau_margin": (2, 2, 6, 2, 3),
    }
    wins, cycles, hp = tournament(scores)
    print("\nTournament Analysis over proof lenses")
    print("observable=(preserves period predicate, preserves missed-shape, tractability,")
    print("            compatibility with THM-563/HYP-2788/S21, formalizability)")
    print(f"wins={wins}")
    print(f"directed_3cycles={cycles}; Hamiltonian_paths={hp}")
    print(f"pressure order={sorted(scores, key=lambda n: scores[n], reverse=True)}")


if __name__ == "__main__":
    main()
