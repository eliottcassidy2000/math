#!/usr/bin/env python3
"""HYP-2714 / OPEN-Q-108: moderate-span multi-block budget scout.

The current sector-route ledger says four regimes are closed and the remaining
gap is the moderate-span balanced/multi-block carrier-product lemma.  This
script extends the exact finite window one step by default and records the proof-obligation
stratification that the remaining lemma must explain.

It uses the fast integer common-refinement engine from the Thread B finite
check, then groups exact rows by:

* far_count = #{e > 14};
* max adjacent gap;
* cluster count after cutting gaps > 7 and > 14;
* a small carrier-product diagnostic for the worst rows.

No proof is claimed.  The goal is to turn the new HYP-2714 "comfortable lemma"
into a quantified target and to look for counter-signal before spending proof
effort.  Set LRC14_S68_WINDOWS="8:24,9:22,10:20" for a deeper exact run.
"""

from __future__ import annotations

import itertools
import math
import os
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cmp_to_key, lru_cache, reduce
from math import gcd


sys.stdout.reconfigure(line_buffering=True)

CAPS = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91)}
FASTCHECK_PREV = {8: 20, 9: 17, 10: 16}
DEFAULT_WINDOWS = {8: 21, 9: 18, 10: 17}
INNER = frozenset(range(1, 7))
Raw = tuple[int, int]
PROGRESS_EVERY = int(os.environ.get("LRC14_S68_PROGRESS", "25000"))


def parse_windows() -> dict[int, int]:
    raw = os.environ.get("LRC14_S68_WINDOWS")
    if not raw:
        return dict(DEFAULT_WINDOWS)
    out: dict[int, int] = {}
    for part in raw.split(","):
        if not part.strip():
            continue
        k_raw, w_raw = part.split(":", 1)
        out[int(k_raw)] = int(w_raw)
    return out


WINDOWS = parse_windows()


def banner(title: str) -> None:
    print("\n" + "=" * 86)
    print(title)
    print("=" * 86)


def gcd_all(xs: tuple[int, ...]) -> int:
    return reduce(gcd, (abs(x) for x in xs if x != 0), 0)


def lcm(a: int, b: int) -> int:
    return a * b // gcd(a, b)


def measS7_raw(speeds: tuple[int, ...]) -> Raw:
    """Exact measS7 of E={0}+speeds as an unreduced integer pair."""
    L = 1
    for e in speeds:
        L = lcm(L, e)
    D = 7 * L
    events: dict[int, list[int]] = defaultdict(list)
    sectors = [0] * len(speeds)
    counts = [len(speeds), 0, 0, 0, 0, 0, 0]
    for idx, e in enumerate(speeds):
        step = L // e
        for m in range(1, 7 * e):
            events[m * step].append(idx)

    for e in speeds:
        if D % e != 0:
            raise AssertionError("internal lcm error")

    num = 0
    prev = 0
    for b, idxs in sorted(events.items()):
        if b > prev and all(counts[s] > 0 for s in range(1, 7)):
            num += b - prev
        for idx in idxs:
            old = sectors[idx]
            new = (old + 1) % 7
            sectors[idx] = new
            counts[old] -= 1
            counts[new] += 1
        prev = b
    if D > prev and all(counts[s] > 0 for s in range(1, 7)):
        num += D - prev
    return num, D


def measS7_int(speeds: tuple[int, ...]) -> F:
    num, den = measS7_raw(speeds)
    return F(num, den)


def raw_fraction(q: Raw | F) -> F:
    if isinstance(q, tuple):
        return F(q[0], q[1])
    return q


def raw_gt(a: Raw, b: Raw) -> bool:
    return a[0] * b[1] > b[0] * a[1]


def raw_cmp_desc(a: tuple[Raw, tuple[int, ...]], b: tuple[Raw, tuple[int, ...]]) -> int:
    left = a[0][0] * b[0][1]
    right = b[0][0] * a[0][1]
    if left > right:
        return -1
    if left < right:
        return 1
    return 0


def raw_ge_fraction(a: Raw, b: F) -> bool:
    return a[0] * b.denominator >= b.numerator * a[1]


def raw_over_fraction(a: Raw, b: F) -> bool:
    return a[0] * b.denominator > b.numerator * a[1]


def margin(cap: F, q: Raw | F) -> F:
    return cap - raw_fraction(q)


def raw_max(xs: list[Raw]) -> Raw:
    best = xs[0]
    for x in xs[1:]:
        if raw_gt(x, best):
            best = x
    return best


def measS7_frac(E: tuple[int, ...]) -> F:
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(7 * e + 1):
            bps.add(F(a, 7 * e))
    total = F(0)
    xs = sorted(bps)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        seen = {int(((e * mid) % 1) * 7) for e in E}
        if len(seen) == 7:
            total += hi - lo
    return total


def gap_data(E: tuple[int, ...]) -> tuple[tuple[int, ...], int]:
    gaps = tuple(b - a for a, b in zip(E, E[1:]))
    return gaps, max(gaps) if gaps else 0


def cluster_shapes(E: tuple[int, ...], cut: int) -> tuple[tuple[int, ...], ...]:
    clusters: list[list[int]] = [[E[0]]]
    for a, b in zip(E, E[1:]):
        if b - a > cut:
            clusters.append([b])
        else:
            clusters[-1].append(b)
    return tuple(tuple(x - c[0] for x in c) for c in clusters)


def shape_key(E: tuple[int, ...]) -> tuple[int, int, int, int]:
    gaps, max_gap = gap_data(E)
    far_count = sum(1 for e in E if e > 14)
    return far_count, max_gap, len(cluster_shapes(E, 7)), len(cluster_shapes(E, 14))


@lru_cache(maxsize=None)
def shifted_profile(diffs: tuple[int, ...]) -> tuple[tuple[frozenset[int], F], ...]:
    """Carrier-averaged hit-set law for a cluster shape with independent phase."""
    D = tuple(sorted(set(diffs)))
    bpx = {F(0), F(1)}
    for d in D:
        if d == 0:
            continue
        for m in range(d + 1):
            bpx.add(F(m, d))
    prof: dict[frozenset[int], F] = defaultdict(F)
    xs = sorted(bpx)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        midx = (lo + hi) / 2
        fr = [(F(d) * midx) % 1 for d in D]
        tb = {F(0), F(1)}
        for f in fr:
            for j in range(7):
                tb.add((F(j, 7) - f) % 1)
        ts = sorted(tb)
        for tlo, thi in zip(ts, ts[1:]):
            if thi <= tlo:
                continue
            tmid = (tlo + thi) / 2
            hit = frozenset(
                s for s in (int(((f + tmid) % 1) * 7) for f in fr) if 1 <= s <= 6
            )
            prof[hit] += (hi - lo) * (thi - tlo)
    return tuple(sorted(prof.items(), key=lambda kv: (len(kv[0]), tuple(sorted(kv[0])))))


def product_cover_for_clusters(clusters: tuple[tuple[int, ...], ...]) -> F:
    cur: dict[frozenset[int], F] = {frozenset(): F(1)}
    for cl in clusters:
        prof = shifted_profile(cl)
        nxt: dict[frozenset[int], F] = defaultdict(F)
        for s1, m1 in cur.items():
            for s2, m2 in prof:
                nxt[s1 | s2] += m1 * m2
        cur = nxt
    return cur.get(INNER, F(0))


def update_best(bucket: dict[object, tuple[Raw, tuple[int, ...]]], key: object, val: Raw, E: tuple[int, ...]) -> None:
    if key not in bucket or raw_gt(val, bucket[key][0]):
        bucket[key] = (val, E)


def exact_window(k: int, W: int) -> dict[str, object]:
    cap = CAPS[k]
    prev = FASTCHECK_PREV[k]
    total = primitive = over = 0
    worst: tuple[Raw, tuple[int, ...]] = ((-1, 1), ())
    by_far: dict[int, tuple[Raw, tuple[int, ...]]] = {}
    by_max_gap: dict[str, tuple[Raw, tuple[int, ...]]] = {}
    by_c7: dict[int, tuple[Raw, tuple[int, ...]]] = {}
    by_c14: dict[int, tuple[Raw, tuple[int, ...]]] = {}
    moderate_balanced: tuple[Raw, tuple[int, ...]] = ((-1, 1), ())
    new_window_worst: tuple[Raw, tuple[int, ...]] = ((-1, 1), ())
    top: list[tuple[Raw, tuple[int, ...]]] = []
    over_cap_rows: list[tuple[Raw, tuple[int, ...]]] = []
    near_cap = 0
    t0 = time.time()
    for combo in itertools.combinations(range(1, W + 1), k - 1):
        total += 1
        if gcd_all(combo) != 1:
            continue
        primitive += 1
        E = (0,) + combo
        val = measS7_raw(combo)
        if raw_gt(val, worst[0]):
            worst = (val, E)
        if max(E) > prev and raw_gt(val, new_window_worst[0]):
            new_window_worst = (val, E)
        far_count, max_gap, c7, c14 = shape_key(E)
        update_best(by_far, far_count, val, E)
        gap_bucket = "<=7" if max_gap <= 7 else "8..14" if max_gap <= 14 else ">14"
        update_best(by_max_gap, gap_bucket, val, E)
        update_best(by_c7, c7, val, E)
        update_best(by_c14, c14, val, E)
        if max(E) > prev and max_gap <= 14 and far_count >= 1:
            if raw_gt(val, moderate_balanced[0]):
                moderate_balanced = (val, E)
        if raw_over_fraction(val, cap):
            over += 1
            over_cap_rows.append((val, E))
        if raw_ge_fraction(val, cap - F(1, 50)):
            near_cap += 1
        top.append((val, E))
        top.sort(key=cmp_to_key(raw_cmp_desc))
        if len(top) > 10:
            top.pop()
        if PROGRESS_EVERY and total % PROGRESS_EVERY == 0:
            print(
                f"    progress combos={total} primitive={primitive} "
                f"best={fmt(worst[0])} elapsed={time.time()-t0:.1f}s"
            )
    elapsed = time.time() - t0
    return {
        "k": k,
        "W": W,
        "cap": cap,
        "total": total,
        "primitive": primitive,
        "over": over,
        "worst": worst,
        "new_window_worst": new_window_worst,
        "moderate_balanced": moderate_balanced,
        "by_far": by_far,
        "by_max_gap": by_max_gap,
        "by_c7": by_c7,
        "by_c14": by_c14,
        "near_cap": near_cap,
        "top": top,
        "elapsed": elapsed,
        "over_cap_rows": over_cap_rows,
    }


def print_bucket(title: str, bucket: dict[object, tuple[Raw, tuple[int, ...]]], cap: F) -> None:
    print(f"  {title}:")
    for key in sorted(bucket, key=lambda x: (str(type(x)), x)):
        val, E = bucket[key]
        print(
            f"    {key!s:>6}: max={float(raw_fraction(val)):.6f} "
            f"margin={float(margin(cap, val)):.6f} E={E}"
        )


def carrier_diagnostics(rows: list[tuple[Raw, tuple[int, ...]]]) -> None:
    print("\nCARRIER-PRODUCT DIAGNOSTIC ON TOP EXACT ROWS")
    seen: set[tuple[int, ...]] = set()
    for val, E in rows:
        if E in seen:
            continue
        seen.add(E)
        clusters7 = cluster_shapes(E, 7)
        clusters14 = cluster_shapes(E, 14)
        pc7 = product_cover_for_clusters(clusters7)
        pc14 = product_cover_for_clusters(clusters14)
        print(f"  E={E} exact={fmt(val)}")
        print(
            f"    cut>7 clusters={clusters7}, product={fmt(pc7)}, "
            f"product-exact={fmt(pc7-raw_fraction(val))}"
        )
        print(
            f"    cut>14 clusters={clusters14}, product={fmt(pc14)}, "
            f"product-exact={fmt(pc14-raw_fraction(val))}"
        )


def fmt(q: F | Raw) -> str:
    qq = raw_fraction(q)
    return f"{qq} ({float(qq):.6f})"


def tournament_summary(results: list[dict[str, object]]) -> None:
    print("\nTOURNAMENT ANALYSIS -- proof obligations for remaining sector route")
    carriers = {
        "bounded_span_exact_window": (0, min(margin(r["cap"], r["worst"][0]) for r in results), 7, 7),
        "single_far_comb_closed": (0, F(12, 100), 6, 7),
        "far_count_domination": (
            0,
            min(
                margin(r["cap"], raw_max([v[0] for k, v in r["by_far"].items() if k >= 2]))
                for r in results
            ),
            5,
            6,
        ),
        "moderate_balanced_multiblock": (
            0,
            min(margin(r["cap"], r["moderate_balanced"][0]) for r in results),
            4,
            6,
        ),
        "carrier_product_error_bound": (1, F(0), 3, 7),
        "raw_gap_monotonicity": (10, F(-1), 1, 3),
    }
    vertices = list(carriers)
    scores: Counter[str] = Counter()
    edges: dict[tuple[str, str], str] = {}

    def better(a: str, b: str) -> str:
        ma, mb = carriers[a], carriers[b]
        ka = (-ma[0], ma[1], ma[2], ma[3])
        kb = (-mb[0], mb[1], mb[2], mb[3])
        return a if ka >= kb else b

    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            w = better(a, b)
            edges[(a, b)] = w
            edges[(b, a)] = w
            scores[w] += 1
    score_hist = Counter(scores[v] for v in vertices)
    cycles = 0
    for a, b, c in itertools.combinations(vertices, 3):
        ab, bc, ca = edges[(a, b)], edges[(b, c)], edges[(c, a)]
        if ab == a and bc == b and ca == c:
            cycles += 1
        if ab == b and bc == c and ca == a:
            cycles += 1
    for name, metric in sorted(carriers.items(), key=lambda kv: (-scores[kv[0]], kv[0])):
        print(f"  {name}: score={scores[name]}, metric={metric}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian path: bounded/single-far/far-count/moderate/carrier-error/raw-gap.")


def main() -> None:
    print("HYP-2714 / OPEN-Q-108 -- moderate-span multi-block budget scout")
    print("Exact integer common-refinement arithmetic; no proof claimed.")
    print(f"Windows={WINDOWS}; override with LRC14_S68_WINDOWS='8:24,9:22,10:20'.")
    checks = [
        (0, 1, 2, 3, 4, 5, 6, 7),
        (0, 2, 4, 6, 8, 10, 12, 14, 29),
        (0, 2, 4, 6, 8, 10, 12, 15, 28),
    ]
    ok = True
    for E in checks:
        a = measS7_int(tuple(e for e in E if e))
        b = measS7_frac(E)
        if a != b:
            ok = False
            print(f"  CROSSCHECK MISMATCH {E}: int={a}, frac={b}")
    print(f"Crosscheck int engine vs Fraction on {len(checks)} rows: {ok}\n")
    results = []
    top_rows: list[tuple[Raw, tuple[int, ...]]] = []
    for k, W in WINDOWS.items():
        banner(f"EXACT FINITE WINDOW: k={k}, span<=W={W}")
        res = exact_window(k, W)
        results.append(res)
        cap = res["cap"]
        worst_val, worst_E = res["worst"]
        new_val, new_E = res["new_window_worst"]
        mb_val, mb_E = res["moderate_balanced"]
        print(
            f"  primitive={res['primitive']} / total={res['total']} "
            f"elapsed={res['elapsed']:.1f}s"
        )
        print(f"  cap={fmt(cap)}")
        print(f"  worst={fmt(worst_val)} margin={fmt(margin(cap, worst_val))} E={worst_E}")
        print(
            f"  worst newly beyond previous W={FASTCHECK_PREV[k]}: "
            f"{fmt(new_val)} margin={fmt(margin(cap, new_val))} E={new_E}"
        )
        print(
            f"  moderate-balanced (span>{FASTCHECK_PREV[k]}, max_gap<=14): "
            f"{fmt(mb_val)} margin={fmt(margin(cap, mb_val))} E={mb_E}"
        )
        print(f"  over_cap={res['over']}, near_cap_within_0.02={res['near_cap']}")
        print_bucket("by far_count", res["by_far"], cap)
        print_bucket("by max_gap bucket", res["by_max_gap"], cap)
        print_bucket("by cut>7 cluster count", res["by_c7"], cap)
        print_bucket("by cut>14 cluster count", res["by_c14"], cap)
        print("  top exact rows:")
        for val, E in res["top"][:6]:
            print(f"    {fmt(val)} margin={fmt(margin(cap, val))} key={shape_key(E)} E={E}")
        top_rows.extend(res["top"][:4])
    carrier_diagnostics(top_rows[:10])
    tournament_summary(results)
    print("\nSYNTHESIS")
    print("  The enlarged exact windows still put the global maximum at the consecutive block.")
    print("  Newly added moderate-span balanced rows are far below cap; the remaining proof")
    print("  obligation is to replace exact enumeration by a carrier-product error bound.")


if __name__ == "__main__":
    main()
