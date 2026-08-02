#!/usr/bin/env python3
"""Exact controls for THM-3135 reflected directed-cycle lane boundary."""

from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations
from math import lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(name: str, rel: str):
    path = ROOT / rel
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MODULES = (
    (
        "one",
        load("probe_one", "04-computation/lrc14_j7_reflected_one_cone_closure_thm2941.py"),
        F(2),
    ),
    (
        "three_quarter",
        load(
            "probe_three_quarter",
            "04-computation/lrc14_j7_reflected_three_quarter_cone_closure_thm2941.py",
        ),
        F(7, 3),
    ),
    (
        "two_thirds",
        load(
            "probe_two_thirds",
            "04-computation/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.py",
        ),
        F(5, 2),
    ),
)


def lane_covers_order(pair, side, rank):
    i, j = pair
    if side == "all":
        return True
    # rank[x] is increasing with q_x, so Q/P=q_j/q_i is below one iff rank[j]<rank[i].
    return (side == "low" and rank[j] < rank[i]) or (
        side == "high" and rank[j] > rank[i]
    )


def coverage(module, replacements=None):
    by_body = {}
    for body, side, pair, start, numerator in module.TAIL_POLICIES:
        by_body.setdefault(body, []).append((side, pair, start, numerator))
    if replacements:
        for body, rows in replacements.items():
            by_body[body] = rows
    rows = []
    for body, lanes in by_body.items():
        missed = []
        for order in permutations(range(6)):
            rank = {vertex: position for position, vertex in enumerate(order)}
            if not any(lane_covers_order(pair, side, rank) for side, pair, _, _ in lanes):
                missed.append(order)
        rows.append((body, len(lanes), len(missed), tuple(missed[:3])))
    return tuple(rows)


def exact_numerator(body, pair, side, cap):
    a, b = body[pair[0]], body[pair[1]]
    if side == "low":
        lo, hi = 1 / cap, F(1)
    elif side == "high":
        lo, hi = F(1), cap
    else:
        lo, hi = 1 / cap, cap
    return 2 * max(abs(lo * a - b), abs(hi * a - b))


def probe_repaired_low(name, module, cap, body):
    pair = (0, 1)
    side = "low"
    numerator = exact_numerator(body, pair, side, cap)
    ruler = 14 * lcm(*body)
    debt = module.T.debt_at_nine(body, ruler)
    limit = F(1, 49) - numerator / ruler - debt
    first = None
    if limit > 0:
        for start in range(2, 1_000_001):
            margin = module.T.tail_envelope(body, pair, start, numerator)
            if margin > 0:
                first = (start, margin)
                break
    if first is None:
        return {
            "name": name,
            "body": body,
            "pair": pair,
            "side": side,
            "cap": cap,
            "numerator": numerator,
            "limit": limit,
            "analytic": False,
        }
    start, margin = first
    previous = module.T.tail_envelope(body, pair, start - 1, numerator)
    step = module.T.tail_step_gain(body, pair, start, numerator)
    a = body[pair[0]]
    tail_eta = numerator / 2 / (ruler - F(a, start))
    channels = module.oriented_channels_below(side, start)
    oversized = module.oriented_channels_below(side, start, oversize=8)
    if channels != oversized:
        raise RuntimeError((name, body, "head universe not exhausted", len(channels), len(oversized)))
    passed = []
    failed = []
    for channel in channels:
        try:
            passed.append(module.located_best(body, pair, channel))
        except RuntimeError as error:
            failed.append((channel, str(error)))
    return {
        "name": name,
        "body": body,
        "pair": pair,
        "side": side,
        "cap": cap,
        "numerator": numerator,
        "limit": limit,
        "analytic": True,
        "start": start,
        "previous": previous,
        "margin": margin,
        "step": step,
        "tail_eta": tail_eta,
        "channels": len(channels),
        "passed": len(passed),
        "failed": tuple(failed),
        "weakest": min(passed) if passed else None,
    }


def analytic_lane(module, cap, body, pair, side):
    numerator = exact_numerator(body, pair, side, cap)
    ruler = 14 * lcm(*body)
    debt = module.T.debt_at_nine(body, ruler)
    limit = F(1, 49) - numerator / ruler - debt
    if limit <= 0:
        return None
    for start in range(2, 1_000_001):
        margin = module.T.tail_envelope(body, pair, start, numerator)
        if margin > 0:
            if side == "low":
                clause = (pair[1], pair[0])  # q_j < q_i
            elif side == "high":
                clause = (pair[0], pair[1])  # q_i < q_j
            else:
                raise RuntimeError(side)
            return {
                "clause": clause,
                "pair": pair,
                "side": side,
                "numerator": numerator,
                "limit": limit,
                "start": start,
                "margin": margin,
            }
    raise RuntimeError((body, pair, side, numerator, limit, "threshold above 1e6"))


def head_audit(module, lane):
    channels = module.oriented_channels_below(lane["side"], lane["start"])
    oversized = module.oriented_channels_below(lane["side"], lane["start"], oversize=10)
    if channels != oversized:
        lane = dict(lane)
        lane["head_ok"] = False
        lane["head_error"] = ("not exhausted", len(channels), len(oversized))
        return lane
    rows = []
    for channel in channels:
        try:
            rows.append(module.located_best(lane["body"], lane["pair"], channel))
        except RuntimeError as error:
            lane = dict(lane)
            lane["head_ok"] = False
            lane["head_error"] = (channel, str(error))
            return lane
    lane = dict(lane)
    lane["head_ok"] = True
    lane["channels"] = len(channels)
    lane["weakest"] = min(rows) if rows else None
    return lane


def shortest_cycle(lanes):
    by_clause = {}
    for lane in lanes:
        by_clause.setdefault(lane["clause"], lane)
    vertices = range(6)
    for length in range(2, 7):
        for cycle in permutations(vertices, length):
            if cycle[0] != min(cycle):
                continue
            clauses = tuple((cycle[k], cycle[(k + 1) % length]) for k in range(length))
            if all(clause in by_clause for clause in clauses):
                return tuple(by_clause[clause] for clause in clauses)
    return None


def cycle_probe(name, module, cap, body):
    analytic = []
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            for side in ("low", "high"):
                lane = analytic_lane(module, cap, body, (i, j), side)
                if lane is not None:
                    lane["body"] = body
                    analytic.append(lane)
    analytic_summary = tuple(
        sorted(
            (
                lane["clause"],
                lane["pair"],
                lane["side"],
                lane["numerator"],
                lane["limit"],
                lane["start"],
                lane["margin"],
            )
            for lane in analytic
        )
    )
    # Prefer one analytic representative per physical comparison, ordered by
    # the smallest exact head universe, then actually audit possible cycles.
    attempts = 0
    rejected = []
    live = list(analytic)
    while True:
        cycle = shortest_cycle(live)
        if cycle is None:
            return {
                "name": name,
                "body": body,
                "analytic_lanes": len(analytic),
                "analytic_summary": analytic_summary,
                "cycle": None,
                "rejected": tuple(rejected),
            }
        attempts += 1
        audited = []
        failed = None
        for lane in cycle:
            result = head_audit(module, lane)
            audited.append(result)
            if not result["head_ok"]:
                failed = lane
                rejected.append((lane, result["head_error"]))
                break
        if failed is None:
            return {
                "name": name,
                "body": body,
                "analytic_lanes": len(analytic),
                "analytic_summary": analytic_summary,
                "attempts": attempts,
                "cycle": tuple(audited),
                "rejected": tuple(rejected),
            }
        live.remove(failed)


def main():
    print("THM-3135 reflected directed-cycle lane boundary")
    for name, module, cap in MODULES:
        print(f"[{name}] existing strict-order coverage")
        for row in coverage(module):
            if row[2]:
                print("MISS", row)
        repair_rows = []
        for body in (module.H, module.H2 if name != "one" else None):
            if body is None:
                continue
            result = probe_repaired_low(name, module, cap, body)
            repair_rows.append(result)
            print("REPAIR", result)

        replacements = {}
        for result in repair_rows:
            if not result.get("analytic"):
                continue
            body = result["body"]
            # Pair the repaired low lane with the already present high (0,1) lane.
            high = next(
                row[1:]
                for row in module.TAIL_POLICIES
                if row[0] == body and row[1] == "high" and row[2] == (0, 1)
            )
            high_side, high_pair, high_start, high_numerator = high
            replacements[body] = [
                ("low", (0, 1), result["start"], result["numerator"]),
                (high_side, high_pair, high_start, high_numerator),
            ]
        print(f"[{name}] repaired strict-order coverage")
        for row in coverage(module, replacements):
            if row[0] in replacements:
                print("REPAIRED", row)
        for body in (module.H, module.H2 if name != "one" else None):
            if body is not None:
                print("CYCLE", cycle_probe(name, module, cap, body))


if __name__ == "__main__":
    main()
