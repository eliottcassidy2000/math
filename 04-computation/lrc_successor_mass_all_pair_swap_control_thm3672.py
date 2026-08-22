#!/usr/bin/env python3
"""Exact successor-mass ledger for the THM-3669 typed non-cover row."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys


sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
REFEREE = ROOT / (
    "04-computation/"
    "lrc14_169_twist_two_twist_referee_thm2334.py"
)

P = 13
UNITS = tuple(range(6))
EXPECTED_PARENT_SHA256 = (
    "0e4a9e181263647e13d2a6738b6996c"
    "45df901d9d2b37d4d589dfddfbdd91480"
)
EXPECTED_SUCCESSOR_SHA256 = (
    "a41ec0c3502221eb1e0e17bd1d941817"
    "fdba134f93c796050a4b454772fee668"
)
EXPECTED_DEFECT_SHA256 = (
    "9129d6dc7c7c03db0b99cd55f8ee504"
    "f65448b87ebdac008e0cc66356787014b"
)

EXPECTED_VALUES = (
    ("zero", (110219232915792, 60084076348296, 9948919780800, 188056)),
    ("a0", (99299969997228, 54135630512964, 8971291028700, 169431)),
    ("a1", (99299969997228, 54135630512964, 8971291028700, 169431)),
    ("a2", (99299969997228, 54135630512964, 8971291028700, 169431)),
    ("a3", (99050735597814, 54000012243177, 8949288888540, 169006)),
    ("a4", (96962693739912, 52861450868226, 8760207996540, 165443)),
    ("a5", (96754696164126, 52747826144013, 8740956123900, 165088)),
    ("b0", (113880033524816, 61887542465528, 9895051406240, 186674)),
    ("b1", (113880033524816, 61887542465528, 9895051406240, 186674)),
    ("b2", (113880033524816, 61887542465528, 9895051406240, 186674)),
    ("b3", (113525486798282, 61695142630901, 9864798463520, 186093)),
    ("b4", (111210186545380, 60436521934750, 9662857324120, 182297)),
    ("b5", (111203438664916, 60432849886708, 9662261108500, 182285)),
)


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


def digest(value):
    encoded = json.dumps(value, separators=(",", ":")).encode("ascii")
    return sha256(encoded).hexdigest()


def load_referee():
    spec = importlib.util.spec_from_file_location("thm2334_referee", REFEREE)
    require(spec is not None and spec.loader is not None, "cannot load referee")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def negative_target_dipole(target, graft):
    """Return -e_target+e_graft, the THM-2365 positive-coordinate shift."""
    vector = [0] * 9
    vector[target] = -1 % P
    vector[graft] = 1
    return tuple(vector)


def marked_intervals(ref, present, word):
    """Intersect E_(s,t) with {x:R*x mod 1 lies in Q} on the NN grid."""
    present = [(ref.R * left, ref.R * right) for left, right in present]
    output = []
    present_index = 0

    for q_left, q_right in ref.word_preimages(word):
        while (
            present_index < len(present)
            and present[present_index][1] <= q_left
        ):
            present_index += 1

        scan = present_index
        while scan < len(present) and present[scan][0] < q_right:
            e_left, e_right = present[scan]
            left = max(e_left, q_left)
            right = min(e_right, q_right)
            if left < right:
                if output and output[-1][1] == left:
                    output[-1] = (output[-1][0], right)
                else:
                    output.append((left, right))
            scan += 1

    return output


def danger_prefix(x, period, window):
    """Length of d(13*c*u)=1 in [0,x] on the NN numerator grid."""
    quotient, residue = divmod(x, period)
    return (
        quotient * (2 * window)
        + min(residue, window)
        + max(0, residue - (period - window))
    )


def main():
    parent_hash = sha256(REFEREE.read_bytes()).hexdigest()
    require(
        parent_hash == EXPECTED_PARENT_SHA256,
        ("THM-2334 referee hash drift", parent_hash),
    )

    ref = load_referee()
    require(
        ref.W == (1, 14, 27, 40, 53, 66, 13, 2197, 742586),
        ref.W,
    )
    require(
        ref.R == 169 and ref.NN == 50334435734703120,
        (ref.R, ref.NN),
    )

    zero = (0,) * 9
    word = ref.build_boolean_set(ref.PATTERN_QA, zero)

    successor_speed = P * ref.W[ref.TARGET_B]
    require(
        ref.NN % (14 * successor_speed) == 0,
        (ref.NN, successor_speed),
    )
    window = ref.NN // (14 * successor_speed)
    period = 14 * window
    require(
        (successor_speed, window, period)
        == (9653618, 372432060, 5214048840),
        (successor_speed, window, period),
    )

    def successor(marked):
        mass = sum(right - left for left, right in marked)
        danger = sum(
            danger_prefix(right, period, window)
            - danger_prefix(left, period, window)
            for left, right in marked
        )
        return 2 * mass - danger, mass, danger, len(marked)

    values = {
        "zero": successor(
            marked_intervals(
                ref,
                ref.build_boolean_set(ref.PATTERN_E, zero),
                word,
            )
        )
    }

    for graft in UNITS:
        values[f"a{graft}"] = successor(
            marked_intervals(
                ref,
                ref.build_boolean_set(
                    ref.PATTERN_E,
                    negative_target_dipole(ref.TARGET_A, graft),
                ),
                word,
            )
        )
        values[f"b{graft}"] = successor(
            marked_intervals(
                ref,
                ref.build_boolean_set(
                    ref.PATTERN_E,
                    negative_target_dipole(ref.TARGET_B, graft),
                ),
                word,
            )
        )

    labels = (
        ("zero",)
        + tuple(f"a{k}" for k in UNITS)
        + tuple(f"b{k}" for k in UNITS)
    )
    ledger = tuple((label, values[label]) for label in labels)
    require(ledger == EXPECTED_VALUES, ("successor ledger changed", ledger))
    require(
        digest(ledger) == EXPECTED_SUCCESSOR_SHA256,
        ("successor digest", digest(ledger)),
    )

    s0 = values["zero"][0]
    avec = tuple(values[f"a{k}"][0] for k in UNITS)
    bvec = tuple(values[f"b{k}"][0] for k in UNITS)
    rows = tuple(
        (ka, kb, s0 + avec[ka] - 2 * bvec[kb])
        for ka in UNITS
        for kb in UNITS
        if ka != kb
    )

    require(len(rows) == 30, len(rows))
    require(
        all(defect < 0 for _, _, defect in rows),
        ("all-negative check failed", rows),
    )
    require(
        digest(rows) == EXPECTED_DEFECT_SHA256,
        ("defect digest", digest(rows)),
    )

    canonical = next(row for row in rows if row[:2] == (1, 2))
    require(canonical == (1, 2, -18240864136612), canonical)

    print(f"denominator={ref.NN}")
    print(f"successor_speed={successor_speed}")
    print(f"window={window}")
    print(f"period={period}")
    for label, value in ledger:
        print(label, value)
    for ka, kb, defect in rows:
        print(f"defect[{ka},{kb}]={defect}/{ref.NN}")
    print(f"canonical_reduced={Fraction(canonical[2], ref.NN)}")
    print(
        "all_30_negative="
        f"{all(defect < 0 for _, _, defect in rows)}"
    )
    print(f"successor_ledger_sha256={digest(ledger)}")
    print(f"defect_ledger_sha256={digest(rows)}")
    print("PASS")


if __name__ == "__main__":
    main()
