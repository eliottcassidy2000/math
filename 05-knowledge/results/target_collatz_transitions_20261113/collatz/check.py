#!/usr/bin/env python3
"""Exact finite controls for the positive affine debt theorem, not orbit evidence."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def apply_action(action, x):
    return 3 * x + 1 if action == "U" else (x - 1) / 3


def trace(word, start):
    values = [F(start)]
    for action, k in word:
        values.append(apply_action(action, values[-1]))
        for _ in range(k):
            values.append(values[-1] / 2)
    return values


def affine(word):
    slope, offset = F(1), F(0)
    for action, k in word:
        if action == "U":
            slope, offset = 3 * slope, 3 * offset + 1
        else:
            slope, offset = slope / 3, (offset - 1) / 3
        slope, offset = slope / 2**k, offset / 2**k
    return slope, offset


def ordinary(x, y):
    return x > 0 and y > 0 and y == (x // 2 if x % 2 == 0 else 3 * x + 1)


def model_a(x, y):
    return ordinary(x, y) or ordinary(y, x)


def model_l(x, y):
    if x <= 0 or y <= 0:
        return False
    return (ordinary(x, y) or (x % 2 == 0 and y == 2 * x)
            or (x > 1 and x % 6 == 1 and y == (x - 1) // 3))


def main():
    alphabet = tuple(product("UV", (1, 2, 3)))
    starts = (F(1), F(7), F(100), F(10**8))
    totals = {"all_words": 0, "paid_words": 0,
              "positive_paid_traces": 0, "positive_formal_fixed_points": 0}
    for length in range(1, 6):
        for word in product(alphabet, repeat=length):
            totals["all_words"] += 1
            a = sum(action == "U" for action, _ in word)
            b = length - a
            slope, offset = affine(word)
            require(slope == F(3) ** (a - b) / 2 ** sum(k for _, k in word),
                    "incorrect affine slope")
            if b < a:
                continue
            totals["paid_words"] += 1
            require(0 < slope < 1, "paid word not contractive")
            for start in starts:
                values = trace(word, start)
                require(values[-1] == slope * start + offset, "trace/affine mismatch")
                if min(values) > 0:
                    totals["positive_paid_traces"] += 1
                    require(values[-1] < start, "positive debt descent failed")
            fixed = offset / (1 - slope)
            if fixed > 0:
                totals["positive_formal_fixed_points"] += 1
                values = trace(word, fixed)
                require(values[-1] == fixed and min(values) <= 0,
                        "all-positive paid cycle")

    hostile = trace((("V", 1), ("U", 1)), F(1, 3))
    require(hostile[-1] == F(1, 3) and hostile[2] == F(-1, 9),
            "formal fixed-point hostile changed")
    guards = {
        "7_to_2": {"A": model_a(7, 2), "L": model_l(7, 2)},
        "40_to_13": {"A": model_a(40, 13), "L": model_l(40, 13)},
        "5_to_10": {"A": model_a(5, 10), "L": model_l(5, 10)},
        "1_to_0": {"A": model_a(1, 0), "L": model_l(1, 0)},
    }
    require(guards == {
        "7_to_2": {"A": False, "L": True},
        "40_to_13": {"A": True, "L": False},
        "5_to_10": {"A": True, "L": False},
        "1_to_0": {"A": False, "L": False},
    }, "graph-type hostile failed")
    for path in ((1, 4, 2, 1), (2, 4, 2), (3, 10, 5, 16, 8, 4, 2, 1)):
        require(all(model_l(x, y) for x, y in zip(path, path[1:])),
                "named cycle/descent path invalid")

    family_controls = []
    for r in range(1, 6):
        for q in range(1, 31):
            n = 2 * 4**r * q + (4**r - 1) // 3
            path = [n, 3 * n + 1]
            for _ in range(2 * r):
                require(path[-1] % 2 == 0, "non-even halving source")
                path.append(path[-1] // 2)
            require(path[-1] == 6 * q + 1, "wrong extra-inverse source")
            path.extend((2 * q, q))
            require(len(path) - 1 == 2 * r + 3 and max(path) == 3 * n + 1,
                    "family length or peak mismatch")
            require(all(model_l(x, y) for x, y in zip(path, path[1:])),
                    "family has illegal literal edge")
            family_controls.append({"r": r, "q": q, "source": n,
                                    "edges": len(path) - 1, "peak": max(path)})
        source = (4**r - 1) // 3
        endpoint = 3 * source + 1
        for _ in range(2 * r):
            require(endpoint % 2 == 0, "q0 boundary halving failed")
            endpoint //= 2
        require(endpoint == 1 and not model_l(endpoint, 0), "q0 boundary failed")

    receipt = {
        "status": "FINITE-EXACT controls of independently proved all-word statements",
        "universe": {"actions": ["U", "V"], "word_lengths": [1, 5],
                     "halving_counts_per_action": [1, 2, 3],
                     "positive_test_starts": [str(x) for x in starts]},
        "counts": totals,
        "formal_endpoint_hostile": [str(x) for x in hostile],
        "graph_guard_controls": guards,
        "source_family": {"r": [1, 5], "q": [1, 30],
                          "total_controls": len(family_controls),
                          "first": family_controls[0], "last": family_controls[-1]},
        "scope": "No orbit census; no claim of universal classical or literal reachability.",
        "reproduce": "python3 -B " + str(Path(__file__).relative_to(Path.cwd())),
    }
    Path(__file__).with_name("receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
