#!/usr/bin/env python3
"""HYP-2790: period-max / Boolean-type bridge scout.

THM-563 reduces the single-far signed deviation for a bounded base B to a
finite endpoint-period maximum:

    periodmax(B) = max_w w * Delta_w(B).

The single-far branch closes whenever

    periodmax(B) < 15 * (cap_k - Plat(B)),

where Plat(B)=p0(B)+p1(B)/7 and k=|B|+1 is the final row size.

This script asks whether the dangerous bounded bases are visible in the
Boolean/type quotient from HYP-2751 or nearby bounded-row quotients.  It is an
explanatory scout, not a proof: ratios are exact, correlations are diagnostic.
"""

from __future__ import annotations

import functools
import itertools
import math
import statistics
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
print = functools.partial(print, flush=True)

P = 7
INNER = tuple(range(1, P))
LOW_FEATURES = ((1, (1,)), (2, (1, 1)), (2, (2,)))
LOW_COEFFS = (21, 57, 2)
KPS3_CERT = {
    (1,): F(584995, 9117927),
    (2,): F(9719, 77931),
    (1, 1, 1): F(2971571, 3039309),
}
CAP = {
    8: F(2243, 5880),
    9: F(1979, 4004),
    10: F(55, 91),
    11: F(66, 91),
    12: F(6, 7),
    13: F(1),
}


def fmt(q: F) -> str:
    if q.denominator == 1:
        return str(q.numerator)
    return f"{q.numerator}/{q.denominator}"


def fdec(q: F, digits: int = 9) -> str:
    return f"{float(q):.{digits}f}"


def lcm(a: int, b: int) -> int:
    return abs(a * b) // math.gcd(a, b)


def sector_of(x: F) -> int:
    return int((x % 1) * P)


def breakpoints(E: tuple[int, ...]) -> list[F]:
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(P * e + 1):
            bps.add(F(m, P * e))
    return sorted(bps)


def exact_mask_atoms(E: tuple[int, ...]) -> dict[int, F]:
    q: dict[int, F] = defaultdict(F)
    bps = breakpoints(E)
    for a, b in zip(bps, bps[1:]):
        if b <= a:
            continue
        mid = (a + b) / 2
        hit = {sector_of(e * mid) for e in E}
        mask = 0
        for s in INNER:
            if s not in hit:
                mask |= 1 << (s - 1)
        q[mask] += b - a
    return dict(q)


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


def type_profile(E: tuple[int, ...]) -> dict[tuple[int, tuple[int, ...]], F]:
    out: dict[tuple[int, tuple[int, ...]], F] = defaultdict(F)
    for mask, mass in exact_mask_atoms(E).items():
        out[(mask.bit_count(), mask_runs(mask))] += mass
    return dict(out)


def containment_type_profile(E: tuple[int, ...]) -> dict[tuple[int, ...], F]:
    q = exact_mask_atoms(E)
    out: dict[tuple[int, ...], F] = defaultdict(F)
    for A in range(1, 64):
        mass = sum(v for M, v in q.items() if (M & A) == A)
        if mass:
            out[mask_runs(A)] += mass
    return dict(out)


def phi_low(profile: dict[tuple[int, tuple[int, ...]], F]) -> F:
    return sum(c * profile.get(f, F(0)) for c, f in zip(LOW_COEFFS, LOW_FEATURES))


def kps3_cut(profile: dict[tuple[int, ...], F]) -> F:
    return sum(c * profile.get(t, F(0)) for t, c in KPS3_CERT.items())


def plat_and_arcs(E: tuple[int, ...]) -> tuple[F, list[tuple[int, F, F]], F, F]:
    bps = breakpoints(E)
    p0 = F(0)
    p1 = F(0)
    arcs: list[tuple[int, F, F]] = []
    cur: dict[int, F] = {}
    for a, b in zip(bps, bps[1:]):
        if b <= a:
            continue
        mid = (a + b) / 2
        hit = {sector_of(e * mid) for e in E}
        miss = [j for j in INNER if j not in hit]
        if len(hit) == P:
            p0 += b - a
        one_miss = miss[0] if len(miss) == 1 else None
        if one_miss is not None:
            p1 += b - a
        for j in INNER:
            active = one_miss == j
            if active and j not in cur:
                cur[j] = a
            if (not active) and j in cur:
                arcs.append((j, cur.pop(j), a))
    for j, a in cur.items():
        arcs.append((j, a, F(1)))
    return p0 + p1 / P, arcs, p0, p1


def Sj_raw(t: F, j: int) -> F:
    t = t - math.floor(t)
    a = F(j, P)
    b = F(j + 1, P)
    if t <= a:
        return -t / P
    if t <= b:
        return -a / P + F(P - 1, P) * (t - a)
    return -a / P + F(P - 1, P) * (b - a) - F(1, P) * (t - b)


def mean_sj(j: int) -> F:
    a = F(j, P)
    b = F(j + 1, P)
    pts = [F(0), a, b, F(1)]
    vals = [Sj_raw(x, j) for x in pts]
    total = F(0)
    for i in range(3):
        total += (vals[i] + vals[i + 1]) * (pts[i + 1] - pts[i]) / 2
    return total


MEAN = {j: mean_sj(j) for j in INNER}


def Sc(t: F, j: int) -> F:
    return Sj_raw(t, j) - MEAN[j]


def endpoint_period(E: tuple[int, ...]) -> int:
    out = 1
    for e in E:
        if e:
            out = lcm(out, e)
    return P * out


def raw_saw_num(grid_num: int, period: int, j: int) -> int:
    """Return Sj_raw(grid_num/period,j) with common denominator 7*period."""
    if 7 * grid_num <= j * period:
        return -grid_num
    if 7 * grid_num <= (j + 1) * period:
        return 6 * grid_num - j * period
    return period - grid_num


def period_max(arcs: list[tuple[int, F, F]], period: int, w0: int = 15) -> tuple[F, int]:
    # The centering constants cancel inside Sc(w*b,j)-Sc(w*a,j).  Since every
    # endpoint denominator divides period, we can scan exact integer numerators
    # over the period grid and divide once by 7*period.
    endpoint_nums: list[tuple[int, int, int]] = []
    for j, a, b in arcs:
        an = (a.numerator * (period // a.denominator)) % period
        bn = (b.numerator * (period // b.denominator)) % period
        endpoint_nums.append((j, an, bn))
    best_num = -10**18
    best_w = w0
    for w in range(w0, w0 + period):
        total_num = 0
        for j, an, bn in endpoint_nums:
            total_num += raw_saw_num((w * bn) % period, period, j)
            total_num -= raw_saw_num((w * an) % period, period, j)
        if total_num > best_num:
            best_num = total_num
            best_w = w
    return F(best_num, P * period), best_w


def Wvec(E: tuple[int, ...]) -> tuple[F, ...]:
    """Six cell cover widths W_a for cells centered at a/7."""
    vals = []
    for cell in INNER:
        lo = F(cell, P) - F(1, 2 * P)
        hi = F(cell, P) + F(1, 2 * P)
        bps = {lo, hi}
        for e in E:
            if e == 0:
                continue
            den = P * abs(e)
            j0 = math.floor(lo * den)
            j1 = math.ceil(hi * den)
            for j in range(j0 - 1, j1 + 2):
                x = F(j, den)
                if lo <= x <= hi:
                    bps.add(x)
        bps = sorted(bps)
        total = F(0)
        for a, b in zip(bps, bps[1:]):
            if b <= a:
                continue
            mid = (a + b) / 2
            if len({sector_of(e * mid) for e in E}) == P:
                total += b - a
        vals.append(total)
    return tuple(vals)


def sorted_prefixes(vec: tuple[F, ...]) -> tuple[F, ...]:
    total = F(0)
    out = []
    for x in sorted(vec, reverse=True):
        total += x
        out.append(total)
    return tuple(out)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    sx = statistics.pstdev(xs)
    sy = statistics.pstdev(ys)
    if sx == 0 or sy == 0:
        return None
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / len(xs) / sx / sy


def slope_line(rows: list[dict[str, object]], xkey: str, ykey: str) -> tuple[float, float] | None:
    xs = [float(r[xkey]) for r in rows]
    ys = [float(r[ykey]) for r in rows]
    if len(xs) < 2:
        return None
    vx = statistics.pvariance(xs)
    if vx == 0:
        return None
    mx = statistics.mean(xs)
    my = statistics.mean(ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / len(xs)
    slope = cov / vx
    return slope, my - slope * mx


def bounded_bases(k: int):
    for rest in itertools.combinations(range(1, 15), k - 2):
        yield (0,) + rest


def analyze_k(k: int, top_n: int, period_cap: int) -> list[dict[str, object]]:
    cap = CAP[k]
    consec = tuple(range(k - 1))
    consec_phi = phi_low(type_profile(consec))
    consec_kps3 = kps3_cut(containment_type_profile(consec))
    consec_w = sorted_prefixes(Wvec(consec))
    rows = []
    for B in bounded_bases(k):
        plat, arcs, p0, p1 = plat_and_arcs(B)
        margin = cap - plat
        period = endpoint_period(B)
        rows.append(
            {
                "B": B,
                "plat": plat,
                "margin": margin,
                "p0": p0,
                "p1": p1,
                "arcs": arcs,
                "arc_count": len(arcs),
                "period": period,
            }
        )
    rows.sort(key=lambda r: (r["plat"], -r["period"]), reverse=True)
    checked = []
    skipped = 0
    for r in rows[:top_n]:
        if r["period"] > period_cap:
            skipped += 1
            continue
        pm, best_w = period_max(r["arcs"], r["period"])
        prof = type_profile(r["B"])
        phi = phi_low(prof)
        cprof = containment_type_profile(r["B"])
        kps3 = kps3_cut(cprof)
        wprefix = sorted_prefixes(Wvec(r["B"]))
        r["phi"] = phi
        r["phi_slack"] = phi - consec_phi
        r["kps3_cut"] = kps3
        r["kps3_deficit"] = consec_kps3 - kps3
        r["w_leak_top1"] = wprefix[0] - consec_w[0]
        r["w_leak_top2"] = wprefix[1] - consec_w[1]
        r["w_leak_total"] = wprefix[-1] - consec_w[-1]
        r["periodmax"] = pm
        r["best_w"] = best_w
        r["ratio"] = pm / r["margin"] if r["margin"] > 0 else F(999)
        r["closes"] = pm < 15 * r["margin"]
        checked.append(r)
    print("\n" + "=" * 96)
    print(f"k={k}: top_n_by_Plat={top_n}, checked={len(checked)}, skipped_period>{period_cap}={skipped}")
    print(f"consec_base={consec}; Phi_low(consec_base)={fmt(consec_phi)}; KPS3cut(consec_base)={fmt(consec_kps3)}")
    failures = [r for r in checked if not r["closes"]]
    print(f"period inequality failures among checked: {len(failures)}")
    if checked:
        worst = max(checked, key=lambda r: r["ratio"])
        print(
            "worst ratio: "
            f"ratio={fdec(worst['ratio'])} "
            f"periodmax={fmt(worst['periodmax'])} margin={fmt(worst['margin'])} "
            f"P={worst['period']} best_w={worst['best_w']} B={worst['B']}"
        )
    print("\nTop checked rows by Plat:")
    for r in checked[:12]:
        print(
            f"  B={str(r['B']):<32} Plat={fdec(r['plat'])} margin={fdec(r['margin'])} "
            f"P={r['period']:<6} pm={fdec(r['periodmax'])} ratio={fdec(r['ratio'])} "
            f"PhiSlack={fdec(r['phi_slack'])} KPS3def={fdec(r['kps3_deficit'])} "
            f"Wtop2={fdec(r['w_leak_top2'])} "
            f"{'CLOSE' if r['closes'] else 'FAIL'}"
        )
    print("\nTop checked rows by period pressure ratio:")
    for r in sorted(checked, key=lambda x: x["ratio"], reverse=True)[:12]:
        print(
            f"  ratio={fdec(r['ratio'])} B={str(r['B']):<32} Plat={fdec(r['plat'])} "
            f"margin={fdec(r['margin'])} pm={fdec(r['periodmax'])} best_w={r['best_w']} "
            f"PhiSlack={fdec(r['phi_slack'])} KPS3def={fdec(r['kps3_deficit'])} W1={fdec(r['w_leak_top1'])} "
            f"W2={fdec(r['w_leak_top2'])} arcs={r['arc_count']}"
        )

    if len(checked) >= 3:
        print("\nDiagnostic Pearson correlations with ratio(B)=periodmax/margin:")
        for key in ("plat", "margin", "phi_slack", "kps3_deficit", "w_leak_top1", "w_leak_top2", "w_leak_total", "arc_count", "period"):
            corr = pearson([float(r["ratio"]) for r in checked], [float(r[key]) for r in checked])
            if corr is None:
                print(f"  {key:<13}: n/a")
            else:
                print(f"  {key:<13}: {corr:+.4f}")
        line = slope_line(checked, "phi_slack", "ratio")
        if line is not None:
            slope, intercept = line
            print(f"  least-squares ratio ~ {slope:+.6f}*PhiSlack + {intercept:+.6f}")

    # Tournament Analysis over explanatory lenses.
    if checked:
        pressure_keys = {
            "periodmax_ratio": max(float(r["ratio"]) for r in checked),
            "plat_margin": max(float(r["plat"]) for r in checked),
            "boolean_type_slack": abs(pearson([float(r["ratio"]) for r in checked], [float(r["phi_slack"]) for r in checked]) or 0),
            "kps_containment_cut": abs(pearson([float(r["ratio"]) for r in checked], [float(r["kps3_deficit"]) for r in checked]) or 0),
            "sorted_cell_leak": abs(pearson([float(r["ratio"]) for r in checked], [float(r["w_leak_top2"]) for r in checked]) or 0),
            "endpoint_arc_count": abs(pearson([float(r["ratio"]) for r in checked], [float(r["arc_count"]) for r in checked]) or 0),
            "raw_endpoint_period": abs(pearson([float(r["ratio"]) for r in checked], [float(r["period"]) for r in checked]) or 0),
        }
        ordered = sorted(pressure_keys, key=lambda x: (pressure_keys[x], x), reverse=True)
        print("\nTournament Analysis over explanatory lenses")
        print("  vertices: periodmax_ratio, plat_margin, boolean_type_slack, kps_containment_cut, sorted_cell_leak, endpoint_arc_count, raw_endpoint_period")
        print("  pairwise observable: stronger absolute correlation/pressure with ratio(B); tie by name")
        print(f"  order={ordered}")
        wins = {name: len(ordered) - 1 - i for i, name in enumerate(reversed(ordered))}
        print(f"  score_hist={dict(Counter(wins.values()))}; directed_3cycles=0; Hamiltonian_paths=1")
    return checked


def main() -> None:
    print("HYP-2790 period-max / Boolean-type bridge scout")
    print("Exact fractions are used for Plat, margin, periodmax, and ratio.")
    print("Diagnostic correlations are float summaries only.")
    all_checked = []
    settings = {8: (40, 30000), 9: (40, 30000), 10: (40, 30000), 11: (30, 30000), 12: (30, 30000)}
    for k, (top_n, period_cap) in settings.items():
        all_checked.extend(analyze_k(k, top_n, period_cap))

    print("\n" + "=" * 96)
    print("Cross-k synthesis")
    print(f"checked rows total={len(all_checked)}")
    fails = [r for r in all_checked if not r["closes"]]
    print(f"period inequality failures total={len(fails)}")
    if all_checked:
        global_worst = max(all_checked, key=lambda r: r["ratio"])
        print(
            f"global worst checked ratio={fdec(global_worst['ratio'])} "
            f"B={global_worst['B']} k={len(global_worst['B'])+1} "
            f"periodmax={fmt(global_worst['periodmax'])} margin={fmt(global_worst['margin'])}"
        )
    print("\nWorking verdict:")
    print("  If boolean_type_slack has strong negative correlation with pressure, HYP-2790")
    print("  gets a plausible finite-ledger route. If not, the endpoint-period numerator")
    print("  should be treated as its own coordinate rather than forced through HYP-2751.")


if __name__ == "__main__":
    main()
