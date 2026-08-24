#!/usr/bin/env python3
"""Exact safe-component response controls for THM-3910.

The symbolic theorem says that a closed body-safe arc of length lambda maps
to a compact arc of length t*lambda.  In a counterexample it must fit inside
one open pair-obstruction component, forcing the strict inequality
t*lambda<beta.  This companion checks the 57 beta_1 values, the scale-two
beta_2 value, and the canonical AP11 arc with its endpoint owners.
"""

from __future__ import annotations

from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")
ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "cyclic_probe",
    ROOT / "04-computation" / "lrc14_eleven_plus_two_cyclic_slack_thm3878.py",
)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load THM-3878 cyclic probe")
P = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(P)


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

DELTA = Q(1, 14)
CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def distance_to_integer(value: Q) -> Q:
    phase = value % 1
    return min(phase, 1 - phase)


def body_safe(value: Q) -> bool:
    return all(distance_to_integer(speed * value) >= DELTA for speed in range(1, 12))


def main() -> None:
    require(len(SCALE1) == 57 and len(set(SCALE1)) == 57, "57 pairs")

    left = Q(1, 14)
    right = Q(13, 154)
    length = right - left
    require(length == Q(1, 77), "AP11 arc length")
    require(body_safe(left) and body_safe(right), "closed AP11 endpoints")
    require(all(body_safe(left + (right - left) * Q(k, 100)) for k in range(101)), "AP11 arc samples")
    # The endpoint equalities identify the owners.  On the interior all
    # eleven inequalities are strict because every wall has been crossed.
    require(distance_to_integer(left) == DELTA, "left owner 1")
    require(distance_to_integer(11 * right) == DELTA, "right owner 11")
    middle = (left + right) / 2
    require(all(distance_to_integer(speed * middle) > DELTA for speed in range(1, 12)), "strict interior")

    beta_rows = tuple((p, q, P.widths(p, q)[0]) for p, q in SCALE1)
    beta_max = max(beta for _, _, beta in beta_rows)
    beta_max_pairs = tuple((p, q) for p, q, beta in beta_rows if beta == beta_max)
    require(beta_max == Q(1, 7), "maximum beta1")
    require(beta_max_pairs == ((1, 3), (1, 4), (1, 9), (1, 10)), "beta1 equality pairs")

    staircase_thresholds = (Q(1, 35), Q(1, 28), Q(1, 21), Q(1, 14), Q(1, 7))
    staircase = tuple(
        (threshold, sum(beta <= threshold for _, _, beta in beta_rows))
        for threshold in staircase_thresholds
    )
    require(tuple(count for _, count in staircase) == (15, 34, 46, 52, 57), "response staircase")

    u = 11
    u_lambda = u * length
    beta2 = Q(2, 63)
    require(u_lambda == Q(1, 7), "AP11 normalized arc")
    require(all(u_lambda >= beta for _, _, beta in beta_rows), "AP11 closes 57")
    require(u_lambda > beta2, "AP11 closes scale two")

    semantic = {
        "ap11_arc": [str(left), str(right), str(length), "owners=1,11"],
        "beta_rows": [[p, q, str(beta)] for p, q, beta in beta_rows],
        "staircase": [[str(threshold), count] for threshold, count in staircase],
        "scale2_beta": str(beta2),
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("LRC14_SAFE_COMPONENT_RESPONSE_THM3910")
    print("scope=THM3878_t>=U_shape_filter;AP11_control;LRC14=OPEN")
    print(f"AP11_closed_arc=[{left},{right}];length={length};owners=1,11;Ulambda={u_lambda}")
    print(f"beta1_max={beta_max};equality_pairs={beta_max_pairs}")
    print("response_staircase=" + repr(tuple((str(threshold), count) for threshold, count in staircase)))
    print(f"AP11_scale1_closed={len(beta_rows)};scale2_beta={beta2};AP11_scale2_closed=1")
    print("strict_boundary=closed_body_arc_cannot_embed_at_equal_length_in_open_obstruction")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
