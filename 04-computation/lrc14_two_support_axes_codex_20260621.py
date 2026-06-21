#!/usr/bin/env python3
"""HYP-2725: reconcile the two odd-support axes in the LRC14 Weyl error.

There are now two nearby uses of "odd support":

1. Factorial-support parity in HYP-2718/HYP-2720:
   Q0 = sum_even W_j - sum_odd W_j, where
   W_j = Delta E[binom(T,j)].

2. Relation-support parity in the Fourier/Weyl lattice:
   corr(E) = sum_{0 != n in Lambda(E)} K(n), stratified by |supp(n)|.

Incoming S9 work refutes relation-support odd dominance as a signed rule:
K(-n)=conj(K(n)), so reverse-pairs reinforce.  This script keeps the exact S69
factorial-support evidence and turns it into a proof-order ledger:

    relation-support filter first,
    factorial odd-L1 envelope second,
    finite even-led packets third,
    then evaluate Q0.

No proof is claimed.  The output is an exact rational triage table for the
HYP-2718 row bank.
"""

from __future__ import annotations

import importlib.util
from collections import Counter
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
S69_PATH = HERE / "lrc14_odd_support_weyl_error_codex_s69.py"
spec = importlib.util.spec_from_file_location("s69_odd_support", S69_PATH)
s69 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(s69)


BLOCK4 = (0, 1, 2, 3)
CASES = [
    ("two 4-blocks, ratio 2:1", ((14, BLOCK4), (28, BLOCK4))),
    ("two 4-blocks, moderate gap", ((14, BLOCK4), (30, BLOCK4))),
    ("two 4-blocks, wider gap", ((30, BLOCK4), (80, BLOCK4))),
    ("two 4-blocks, high relation phase", ((15, BLOCK4), (31, BLOCK4))),
    ("two 4-blocks, positive Q0", ((40, BLOCK4), (61, BLOCK4))),
    ("two 4-blocks, positive Q0 high", ((80, BLOCK4), (121, BLOCK4))),
    ("5+3 split", ((20, (0, 1, 2, 3, 4)), (55, (0, 1, 2)))),
    ("3+3+2 split", ((18, (0, 1, 2)), (45, (0, 1, 2)), (90, (0, 1)))),
    (
        "five 2-blocks",
        (
            (15, (0, 1)),
            (30, (0, 1)),
            (46, (0, 1)),
            (63, (0, 1)),
            (81, (0, 1)),
        ),
    ),
    (
        "seven singleton carriers",
        (
            (19, (0,)),
            (31, (0,)),
            (44, (0,)),
            (58, (0,)),
            (73, (0,)),
            (89, (0,)),
            (106, (0,)),
        ),
    ),
]


def fmt(q: F | None) -> str:
    if q is None:
        return "n/a"
    return f"{q} ({float(q):.9f})"


def ledger_class(row: dict[str, object]) -> str:
    odd_ge_even = row["odd_l1"] >= row["even_l1"]
    if odd_ge_even and row["signed_odd_dominates"]:
        return "odd_L1_signed_phase"
    if odd_ge_even:
        return "odd_L1_envelope"
    return "finite_even_L1_ledger"


def main() -> None:
    rows = [s69.row_report(name, blocks) for name, blocks in CASES]
    print("HYP-2725 two-support-axis ledger for the LRC14 Weyl/origin error")
    print("Exact S69 factorial-support arithmetic; relation-support lesson from S9.")
    print()
    print("AXIS WARNING")
    print("  relation-support parity is not the proof rule:")
    print("    K(-n)=conj(K(n)), so reverse relation-pairs reinforce as 2 Re K(n).")
    print("  keep support-size as cut/cycle structure, not as odd/even cancellation.")
    print()
    print("FACTORIAL-SUPPORT TRIAGE")
    print(
        "  class                      odd_share     |Q0|/(cap-product)     "
        "Q0 sign  row"
    )
    rows_sorted = sorted(
        rows,
        key=lambda r: (
            ledger_class(r) != "odd_L1_signed_phase",
            ledger_class(r) != "odd_L1_envelope",
            -(r["ratio"] or F(0)),
        ),
    )
    for row in rows_sorted:
        sign = "+" if row["q0"] > 0 else "-" if row["q0"] < 0 else "0"
        print(
            f"  {ledger_class(row):26s} "
            f"{float(row['odd_l1_share']):.9f}   "
            f"{fmt(row['ratio']):28s} {sign:>2s}     {row['name']}"
        )

    print()
    print("AGGREGATES")
    class_counts = Counter(ledger_class(row) for row in rows)
    for cls in sorted(class_counts):
        print(f"  {cls}: {class_counts[cls]}")

    for cls in ("odd_L1_signed_phase", "odd_L1_envelope", "finite_even_L1_ledger"):
        group = [row for row in rows if ledger_class(row) == cls]
        if not group:
            continue
        q0_l1 = sum(abs(row["q0"]) for row in group)
        odd_l1 = sum(row["odd_l1"] for row in group)
        even_l1 = sum(row["even_l1"] for row in group)
        max_ratio = max((row["ratio"] or F(0) for row in group), default=F(0))
        print(
            f"  {cls}: sum|Q0|={fmt(q0_l1)} "
            f"odd_L1={fmt(odd_l1)} even_L1={fmt(even_l1)} "
            f"max_cap_risk={fmt(max_ratio)}"
        )

    print()
    print("PROOF ORDER")
    print("  1. Relation-support axis: split support-2 cut terms, support-3 Schur")
    print("     packets, and higher support; do not expect reverse-pair cancellation.")
    print("  2. Factorial-support axis: on the high-height/nonresonant tail, prove")
    print("     an odd-L1 envelope for W_1+W_3+W_5 versus W_0+W_2+W_4+W_6.")
    print("  3. Route finite_even_L1_ledger rows and signed-even-led packets through")
    print("     the HYP-2717/HYP-2714 low-height ledger.")
    print("  4. Only after those filters evaluate Q0=even_support-odd_support.")
    print()
    print("SYNTHESIS")
    print("  Odd support remains useful only on the factorial axis, as an L1")
    print("  envelope.  On the relation-lattice axis, parity is a false friend;")
    print("  support size, especially support-3 Schur packets, is the live geometry.")


if __name__ == "__main__":
    main()
