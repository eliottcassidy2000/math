#!/usr/bin/env python3
r"""Exact THM-3910 scout for THM-3878 pack-conditioned two-danger moments.

For a core u and relative scale t, let G_u be the closed 1/14-safe set and

    N_t(w) = #{a mod t : (w+a)/t lies in G_u}.

Then every scale-one pair (p,q) has

    mu(G_u \ (D_{tp} union D_{tq}))
      = (1/t) integral N_t(w) 1_{S_p intersect S_q}(w) dw.

The code constructs N_t exactly by pushing the rational safe arcs of G_u
through y -> t*y, and audits M0-M1+M2, the t-sheet projection identity, and
the Cauchy--Schwarz/CV certificate on the 57 THM-3878 scale-one survivors.
It also scans AP11 through a finite range of legal t as a hostile control.
"""

from __future__ import annotations

from collections import defaultdict
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import ceil, floor, gcd
import sys


sys.stdout.reconfigure(newline="\n")

DELTA = Q(1, 14)
CHECKS = 0

SCALE1 = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (2, 15), (2, 21), (2, 23),
    (3, 7), (3, 8), (3, 14), (3, 17), (3, 19), (3, 20),
    (3, 22), (3, 26), (3, 31), (3, 38),
    (4, 7), (4, 13), (4, 19), (4, 21), (4, 25), (4, 37),
    (4, 43), (4, 49), (4, 51),
    (5, 6), (5, 12), (5, 17), (5, 18), (5, 24), (5, 29),
    (5, 36), (5, 39), (5, 41), (5, 42), (5, 48), (5, 53),
    (5, 54), (5, 63),
    (6, 11), (6, 17), (6, 19), (6, 23), (6, 41), (6, 47),
    (6, 53), (6, 65),
    (7, 10), (7, 13), (7, 15), (7, 22),
    (8, 9), (8, 21), (9, 11),
)


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def fmt(x: Q) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def dist(x: Q) -> Q:
    x %= 1
    return min(x, 1 - x)


def danger_intervals(speed: int) -> list[tuple[Q, Q]]:
    radius = Q(1, 14 * speed)
    ans = []
    for k in range(speed + 1):
        left = max(Q(0), Q(k, speed) - radius)
        right = min(Q(1), Q(k, speed) + radius)
        if left < right:
            ans.append((left, right))
    return ans


def merge(intervals: list[tuple[Q, Q]]) -> list[tuple[Q, Q]]:
    ans: list[list[Q]] = []
    for left, right in sorted(intervals):
        if not ans or left >= ans[-1][1]:
            ans.append([left, right])
        elif right > ans[-1][1]:
            ans[-1][1] = right
    return [(left, right) for left, right in ans]


def safe_intervals(speeds: tuple[int, ...]) -> list[tuple[Q, Q]]:
    danger = merge(sum((danger_intervals(v) for v in speeds), []))
    ans = []
    cursor = Q(0)
    for left, right in danger:
        if cursor < left:
            ans.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        ans.append((cursor, Q(1)))
    return ans


def interval_measure(intervals: list[tuple[Q, Q]]) -> Q:
    return sum((right - left for left, right in intervals), Q(0))


def pushed_lift_events(core_safe: list[tuple[Q, Q]], t: int) -> dict[Q, int]:
    """Events of N_t on [0,1]; endpoints are measure-zero and immaterial."""
    events: dict[Q, int] = defaultdict(int)
    for left, right in core_safe:
        lifted_left = t * left
        lifted_right = t * right
        first = floor(lifted_left)
        last = ceil(lifted_right) - 1
        for k in range(first, last + 1):
            a = max(lifted_left, Q(k)) - k
            b = min(lifted_right, Q(k + 1)) - k
            if a < b:
                events[a] += 1
                events[b] -= 1
    events[Q(0)] += 0
    events[Q(1)] += 0
    return dict(events)


def pair_flags(p: int, q: int, w: Q) -> tuple[bool, bool]:
    return dist(p * w) < DELTA, dist(q * w) < DELTA


def fiber_moments(
    core: tuple[int, ...],
    t: int,
    p: int,
    q: int,
    core_safe: list[tuple[Q, Q]] | None = None,
    events: dict[Q, int] | None = None,
) -> dict[str, Q | int]:
    if core_safe is None:
        core_safe = safe_intervals(core)
    m0 = interval_measure(core_safe)
    if events is None:
        events = pushed_lift_events(core_safe, t)
    pair_walls = {Q(0), Q(1)}
    for v in (p, q):
        for left, right in danger_intervals(v):
            pair_walls.add(left)
            pair_walls.add(right)
    walls = sorted(set(events) | pair_walls)
    n = 0
    int_n = Q(0)
    int_n2 = Q(0)
    int_safe = Q(0)
    int_dp = Q(0)
    int_dq = Q(0)
    int_both = Q(0)
    min_n = None
    max_n = None
    zero_safe_measure = Q(0)
    positive_safe_measure = Q(0)
    for left, right in zip(walls, walls[1:]):
        n += events.get(left, 0)
        if left == right:
            continue
        length = right - left
        mid = (left + right) / 2
        dp, dq = pair_flags(p, q, mid)
        pair_safe = not dp and not dq
        int_n += n * length
        int_n2 += n * n * length
        int_dp += n * length * dp
        int_dq += n * length * dq
        int_both += n * length * (dp and dq)
        if pair_safe:
            int_safe += n * length
            if n == 0:
                zero_safe_measure += length
            else:
                positive_safe_measure += length
        min_n = n if min_n is None else min(min_n, n)
        max_n = n if max_n is None else max(max_n, n)
    n += events.get(Q(1), 0)
    require(n == 0, ("event balance", core, t, p, q, n))
    require(int_n == t * m0, ("pushforward mass", core, t, p, q))

    m1 = (int_dp + int_dq) / t
    m2 = int_both / t
    safe = int_safe / t
    require(safe == m0 - m1 + m2, ("M0-M1+M2", core, t, p, q))

    pair_danger = merge(danger_intervals(p) + danger_intervals(q))
    obstruction = interval_measure(pair_danger)
    h = 1 - obstruction
    mean_n = int_n
    var_n = int_n2 - mean_n * mean_n
    cv2 = var_n / (mean_n * mean_n) if mean_n else Q(10**100)
    threshold = h / (1 - h)
    # Exact squared Cauchy--Schwarz test avoids introducing irrational roots.
    cs_closes = cv2 < threshold
    # If safe mass vanished, integer N would be zero on the h-set and live on
    # a set of measure b=1-h.  At fixed mean m, convexity on Z gives the sharp
    # two-level lower second moment: N=k or k+1 on that b-set.
    b = obstruction
    average_on_obstruction = mean_n / b
    k = floor(average_on_obstruction)
    theta = average_on_obstruction - k
    integer_failure_var = mean_n * mean_n * h / b + b * theta * (1 - theta)
    integer_failure_cv2 = integer_failure_var / (mean_n * mean_n)
    integer_closes = var_n < integer_failure_var

    # For f=1_G with r positive arcs, BV(f)=2r and
    # |fhat(kt)|<=r/(pi|k|t).  Parseval gives Var(N_t)<=r^2/3.
    arc_count = len(core_safe)
    bv_var_upper = Q(arc_count * arc_count, 3)
    require(var_n <= bv_var_upper, ("BV variance bound", core, t, p, q))
    bv_cv2_upper = bv_var_upper / (mean_n * mean_n)
    bv_closes = bv_cv2_upper < threshold
    bv_integer_closes = bv_var_upper < integer_failure_var
    return {
        "m0": m0,
        "m1": m1,
        "m2": m2,
        "safe": safe,
        "pair_safe": h,
        "mean_n": mean_n,
        "var_n": var_n,
        "cv2": cv2,
        "threshold": threshold,
        "cs_closes": int(cs_closes),
        "integer_failure_cv2": integer_failure_cv2,
        "integer_closes": int(integer_closes),
        "integer_k": k,
        "integer_theta": theta,
        "arc_count": arc_count,
        "bv_var_upper": bv_var_upper,
        "bv_cv2_upper": bv_cv2_upper,
        "bv_closes": int(bv_closes),
        "bv_integer_closes": int(bv_integer_closes),
        "min_n": -1 if min_n is None else min_n,
        "max_n": -1 if max_n is None else max_n,
        "zero_safe_measure": zero_safe_measure,
        "positive_safe_measure": positive_safe_measure,
    }


def record(core: tuple[int, ...], t: int, p: int, q: int, m: dict[str, Q | int]) -> dict[str, object]:
    return {
        "core": core,
        "t": t,
        "p": p,
        "q": q,
        **{key: (fmt(value) if isinstance(value, Q) else value) for key, value in m.items()},
    }


def main() -> None:
    require(len(SCALE1) == 57 and len(set(SCALE1)) == 57, "survivor universe")
    ap11 = tuple(range(1, 12))
    ap_m0 = interval_measure(safe_intervals(ap11))
    require(ap_m0 == Q(10931, 194040), "AP11 safe mass")

    rows = []
    for t in range(11, 85):
        for p, q in SCALE1:
            full = ap11 + (t * p, t * q)
            if len(set(full)) != 13:
                continue
            m = fiber_moments(ap11, t, p, q)
            require(m["safe"] > 0, ("AP11 strict positive", t, p, q))
            rows.append(record(ap11, t, p, q, m))

    require(rows, "nonempty AP11 scan")
    min_safe = min(Q(row["safe"]) for row in rows)
    min_safe_rows = tuple(
        (row["t"], row["p"], row["q"], row["safe"], row["m1"], row["m2"])
        for row in rows if Q(row["safe"]) == min_safe
    )
    max_cover_ratio = max((Q(row["m1"]) - Q(row["m2"])) / Q(row["m0"]) for row in rows)
    max_cover_rows = tuple(
        (row["t"], row["p"], row["q"], row["safe"])
        for row in rows
        if (Q(row["m1"]) - Q(row["m2"])) / Q(row["m0"]) == max_cover_ratio
    )
    cs_closed = [row for row in rows if row["cs_closes"]]
    cs_open = [row for row in rows if not row["cs_closes"]]
    integer_closed = [row for row in rows if row["integer_closes"]]
    integer_open = [row for row in rows if not row["integer_closes"]]
    bv_closed = [row for row in rows if row["bv_closes"]]
    bv_integer_closed = [row for row in rows if row["bv_integer_closes"]]
    zero_fiber_rows = [row for row in rows if Q(row["zero_safe_measure"]) > 0]
    cs_open_t_hist = tuple(sorted(Counter(row["t"] for row in cs_open).items()))
    worst_cs_ratio = max(Q(row["cv2"]) / Q(row["threshold"]) for row in rows)
    worst_cs_rows = tuple(
        (row["t"], row["p"], row["q"], row["cv2"], row["threshold"], row["safe"], row["min_n"], row["max_n"])
        for row in rows if Q(row["cv2"]) / Q(row["threshold"]) == worst_cs_ratio
    )
    first_cs_hostiles = tuple(
        (row["t"], row["p"], row["q"], row["cv2"], row["threshold"], row["safe"], row["zero_safe_measure"])
        for row in cs_open[:12]
    )

    # A compact boundary table at t=12, the first scale above AP11's U.
    boundary = []
    for p, q in SCALE1:
        full = ap11 + (12 * p, 12 * q)
        if len(set(full)) == 13:
            boundary.append(record(ap11, 12, p, q, fiber_moments(ap11, 12, p, q)))

    semantic = {
        "universe": "AP11; 57 THM3878 scale-one survivors; 11<=t<=84; distinct full rows",
        "rows": rows,
        "boundary_t12": boundary,
    }
    digest = sha256(json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()

    print("THM3878_PACK_CONDITIONED_MIXED_MOMENT_FIBER_SCOUT_20260823")
    print("scope=AP11_control_only;scale1_57_certificate_survivors;t=11..84;LRC14=OPEN")
    print(f"fiber_identity=mu_safe=(1/t)int_Nt*pair_safe=M0-M1+M2;AP11_M0={fmt(ap_m0)}")
    print(f"rows={len(rows)};strict_positive={len(rows)}")
    print(f"min_safe={fmt(min_safe)};rows={min_safe_rows}")
    print(f"max_conditioned_cover_ratio={fmt(max_cover_ratio)};rows={max_cover_rows}")
    print(f"cs_cv_certificate=closed:{len(cs_closed)},open:{len(cs_open)}")
    print(f"exact_integer_second_moment_certificate=closed:{len(integer_closed)},open:{len(integer_open)}")
    print(f"bv_certificate=plain_closed:{len(bv_closed)},integer_refined_closed:{len(bv_integer_closed)},rows:{len(rows)}")
    print(f"cs_open_t_histogram={cs_open_t_hist}")
    print(f"cs_worst_ratio={fmt(worst_cs_ratio)};rows={worst_cs_rows}")
    print(f"cs_first_hostiles={first_cs_hostiles}")
    print(f"pair_safe_cells_with_Nt_zero_rows={len(zero_fiber_rows)}")
    print("boundary_t12=" + repr(tuple(
        (row["p"], row["q"], row["safe"], row["m1"], row["m2"], row["cv2"], row["threshold"], row["cs_closes"], row["integer_failure_cv2"], row["integer_closes"], row["bv_cv2_upper"], row["bv_closes"], row["bv_integer_closes"], row["min_n"], row["max_n"])
        for row in boundary
    )))
    print(f"semantic_sha256={digest}")
    print(f"checks={CHECKS}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
