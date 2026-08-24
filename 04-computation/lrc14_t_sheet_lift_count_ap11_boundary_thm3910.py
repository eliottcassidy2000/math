#!/usr/bin/env python3
"""Exact AP11 boundary table for sheet-CV, integer tariff, and BV gates."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
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


def fmt(x: Q) -> str:
    return F.fmt(x)


def main() -> None:
    core = tuple(range(1, 12))
    safe = F.safe_intervals(core)
    mu = F.interval_measure(safe)
    F.require(mu == Q(10931, 194040), "AP11 mass")
    F.require(len(safe) == 14, "AP11 positive arcs")

    by_t = []
    equality_rows = []
    theta_zero_rows = []
    records = []
    for t in range(12, 101):
        events = F.pushed_lift_events(safe, t)
        rows = []
        for p, q in F.SCALE1:
            m = F.fiber_moments(core, t, p, q, safe, events)
            rows.append((p, q, m))
            records.append((t, p, q, m))
            if m["cv2"] == m["threshold"] or m["cv2"] == m["integer_failure_cv2"] or m["cv2"] == m["bv_cv2_upper"]:
                equality_rows.append((t, p, q, fmt(m["cv2"]), fmt(m["threshold"]), fmt(m["integer_failure_cv2"]), fmt(m["bv_cv2_upper"])))
            if m["integer_theta"] == 0:
                theta_zero_rows.append((t, p, q, m["integer_k"]))
        by_t.append((
            t,
            sum(m["cs_closes"] for _, _, m in rows),
            sum(m["integer_closes"] for _, _, m in rows),
            sum(m["bv_closes"] for _, _, m in rows),
            sum(m["bv_integer_closes"] for _, _, m in rows),
            min(m["safe"] for _, _, m in rows),
        ))

    first_all = {}
    for label, index in (("exact_cv", 1), ("exact_integer", 2), ("bv", 3), ("bv_integer", 4)):
        first_all[label] = next((row[0] for row in by_t if row[index] == 57), None)
    eventual_all = {}
    for label, index in (("exact_cv", 1), ("exact_integer", 2), ("bv", 3), ("bv_integer", 4)):
        eventual_all[label] = next(
            (row[0] for j, row in enumerate(by_t) if all(later[index] == 57 for later in by_t[j:])),
            None,
        )

    exceptional = tuple(
        (t, cv, integer, bv, bv_integer, fmt(min_safe))
        for t, cv, integer, bv, bv_integer, min_safe in by_t
        if (cv, integer, bv, bv_integer) != (57, 57, 57, 57)
    )

    # Exact per-pair monotone BV threshold.  This is independent of the
    # realized N_t and witnesses the all-large-t theorem.
    per_pair_bv = []
    for p, q in F.SCALE1:
        probe = next(m for t, pp, qq, m in records if pp == p and qq == q)
        threshold = probe["threshold"]
        t0 = 1
        while Q(14 * 14, 3) / (t0 * t0 * mu * mu) >= threshold:
            t0 += 1
        per_pair_bv.append((p, q, threshold, t0))
    max_t0 = max(t0 for _, _, _, t0 in per_pair_bv)
    worst_pairs = tuple((p, q, fmt(threshold), t0) for p, q, threshold, t0 in per_pair_bv if t0 == max_t0)
    min_threshold = min(threshold for _, _, threshold, _ in per_pair_bv)
    min_threshold_pairs = tuple((p, q) for p, q, threshold, _ in per_pair_bv if threshold == min_threshold)

    integer_open_hist = tuple(sorted(Counter(t for t, _, _, m in records if not m["integer_closes"]).items()))
    semantic = {
        "core": core,
        "mu": fmt(mu),
        "by_t": tuple((t, cv, integer, bv, bv_integer, fmt(ms)) for t, cv, integer, bv, bv_integer, ms in by_t),
        "per_pair_bv": tuple((p, q, fmt(th), t0) for p, q, th, t0 in per_pair_bv),
        "equalities": tuple(equality_rows),
        "theta_zero": tuple(theta_zero_rows),
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("THM3878_AP11_MIXED_MOMENT_GATE_BOUNDARY_20260823")
    print(f"core=AP11;mu={fmt(mu)};positive_arcs={len(safe)};t=12..100;pairs=57")
    print(f"first_all_pairs={first_all};eventual_all_pairs={eventual_all}")
    print(f"integer_open_t_histogram={integer_open_hist}")
    print(f"exceptional_t_rows={exceptional}")
    print(f"bv_uniform_pair_threshold=min_h_over_b:{fmt(min_threshold)},pairs:{min_threshold_pairs},first_integer_t:{max_t0},worst:{worst_pairs}")
    print(f"gate_equalities={tuple(equality_rows)}")
    print(f"integer_theta_zero_rows={tuple(theta_zero_rows)}")
    print(f"semantic_sha256={digest}")
    print(f"checks={F.CHECKS}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
