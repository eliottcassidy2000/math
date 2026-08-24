#!/usr/bin/env python3
"""Exact safe-lift variance gates for the THM-3878 scale-two (1,9) row."""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import importlib.util
import json
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")
HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("fiber", HERE / "lrc14_t_sheet_lift_count_variance_thm3910.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load fiber scout")
F = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(F)

H = Q(59, 63)
B = Q(4, 63)
THRESHOLD = H / B


def fmt(x: Q) -> str:
    return F.fmt(x)


def scale2_gate(core: tuple[int, ...], t: int) -> dict[str, object]:
    safe = F.safe_intervals(core)
    events = F.pushed_lift_events(safe, t)
    # cv2, mean and variance depend only on (core,t), not on this probe pair.
    probe = F.fiber_moments(core, t, 1, 9, safe, events)
    mean = probe["mean_n"]
    var = probe["var_n"]
    average = mean / B
    k = average.numerator // average.denominator
    theta = average - k
    tariff = mean * mean * H / B + B * theta * (1 - theta)
    bv = Q(len(safe) ** 2, 3)
    full = tuple(2 * u for u in core) + (t, 9 * t)
    actual = F.interval_measure(F.safe_intervals(full))
    return {
        "actual": actual,
        "cv2": probe["cv2"],
        "cv": int(probe["cv2"] < THRESHOLD),
        "integer": int(var < tariff),
        "bv": int(bv < mean * mean * THRESHOLD),
        "bv_integer": int(bv < tariff),
        "theta": theta,
        "arcs": len(safe),
    }


def main() -> None:
    ap11 = tuple(range(1, 12))
    ap_rows = []
    for t in range(11, 100, 2):
        ap_rows.append((t, scale2_gate(ap11, t)))
    F.require(all(row[1]["actual"] > 0 for row in ap_rows), "AP11 positive")
    F.require(all(row[1]["cv"] for row in ap_rows), "AP11 exact CV closes")

    cores = tuple(combinations(range(1, 14), 11))
    small = []
    for core in cores:
        u = max(core)
        first = u if u % 2 else u + 1
        for t in (first, first + 2):
            small.append((core, t, scale2_gate(core, t)))
    F.require(len(small) == 156, "small universe")
    F.require(all(row[2]["actual"] > 0 for row in small), "small positive")
    F.require(all(row[2]["cv"] for row in small), "small exact CV closes")

    ap_first_bv = next(t for t, m in ap_rows if m["bv"])
    ap_first_bv_integer = next(t for t, m in ap_rows if m["bv_integer"])
    ap_min = min(m["actual"] for _, m in ap_rows)
    ap_min_rows = tuple((t, fmt(m["actual"])) for t, m in ap_rows if m["actual"] == ap_min)
    small_min = min(m["actual"] for _, _, m in small)
    small_min_rows = tuple((core, t, fmt(m["actual"])) for core, t, m in small if m["actual"] == small_min)
    semantic = {
        "threshold": fmt(THRESHOLD),
        "ap": tuple((t, {k: fmt(v) if isinstance(v, Q) else v for k, v in m.items()}) for t, m in ap_rows),
        "small": tuple((core, t, {k: fmt(v) if isinstance(v, Q) else v for k, v in m.items()}) for core, t, m in small),
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("THM3878_SCALE2_SAFE_LIFT_VARIANCE_GATE_20260823")
    print(f"quotient_pair=(1,9);h=59/63;b=4/63;threshold_h_over_b={fmt(THRESHOLD)}")
    print(f"AP11_odd_t_11_99=rows:{len(ap_rows)},actual_positive:{len(ap_rows)},exact_cv_closed:{sum(m['cv'] for _,m in ap_rows)},integer_closed:{sum(m['integer'] for _,m in ap_rows)},first_bv_odd_t:{ap_first_bv},first_bv_integer_odd_t:{ap_first_bv_integer},min_safe:{fmt(ap_min)},min_rows:{ap_min_rows}")
    print(f"small_cores=78;first_two_legal_odd_t;rows:{len(small)},actual_positive:{len(small)},exact_cv_closed:{sum(m['cv'] for _,_,m in small)},integer_closed:{sum(m['integer'] for _,_,m in small)},bv_closed:{sum(m['bv'] for _,_,m in small)},min_safe:{fmt(small_min)},min_rows:{small_min_rows}")
    print(f"semantic_sha256={digest}")
    print(f"checks={F.CHECKS}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
