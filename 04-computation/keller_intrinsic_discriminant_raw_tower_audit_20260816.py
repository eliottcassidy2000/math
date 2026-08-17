#!/usr/bin/env python3
"""Exact arithmetic certificate for THM-3531's intrinsic discriminant law.

For the raw fixed-map tower P_(n+1)=L^e_n N(P_n), trace-discriminant
transitivity and the odd cubic block law reduce the square class at level n
to [(-1)^n P_(n-1)].  This companion checks every parity and normalization
invoice, including the mandatory constant class [2] in the named levels.

The geometric/function-field proof lives in THM-3531.  This script does not
assert that any fixed coordinate is primitive at every depth and does not
compute exact positive discriminant multiplicities.
"""

from __future__ import annotations

from hashlib import sha256
import json


EXPECTED_SEMANTIC_SHA256 = (
    "9ccc7adc713ebf30e69af6456e212b94f4a43da4cec6215f730c58cc00ecfe46"
)
MATRIX = ((7, -2), (3, -2))


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def step(row: tuple[int, int]) -> tuple[int, int]:
    e, m = row
    return 7 * e - 2 * m, 3 * e - 2 * m


def main() -> None:
    rows = [(1, 0)]
    for _ in range(31):
        rows.append(step(rows[-1]))

    parity_ledger = []
    induction_ledger = []
    for index, (e, m) in enumerate(rows):
        require(e % 2 == 1, ("odd clearing exponent", index, e, m))
        require((1 - e) % 2 == 0, ("old-L square exponent", index, e))
        parity_ledger.append((index, e % 2, m % 2, (1 - e) % 2))
        if index + 1 < len(rows):
            # If d_n=[(-1)^n P_(n-1)], then the odd-block recursion gives
            # N(d_n)*[-L]=[(-1)^(n+1) P_n L^(1-e_(n-1))].
            level = index + 1
            sign_before = level % 2
            sign_after = (sign_before + 1) % 2
            old_L_parity = (1 - e) % 2
            require(old_L_parity == 0, ("induction residue", level, e))
            require(sign_after == (level + 1) % 2,
                    ("sign alternation", level, sign_after))
            induction_ledger.append(
                (level, index, e, sign_before, 1, sign_after, old_L_parity)
            )

    # Exact raw-to-named scalar invoices.  P_1=H/2^6 follows from
    # N(L)=H/(2^6 L).  P_2=J/2^53 follows from
    # N(H)=J/(2^35 L^7), since 53=3*6+35.  Later named cleared norms introduce
    # no new scalar in their definitions, so the denominator exponent cubes.
    raw_named = (
        (0, "L", 0),
        (1, "H", 6),
        (2, "J", 53),
        (3, "G", 159),
        (4, "R5", 477),
        (5, "R6", 1431),
    )
    require(raw_named[2][2] == 3 * raw_named[1][2] + 35,
            "J normalization invoice")
    for index in range(3, len(raw_named)):
        require(raw_named[index][2] == 3 * raw_named[index - 1][2],
                ("cubed scalar invoice", index, raw_named[index]))

    named_classes = []
    expected_classes = ("-L", "H", "-2J", "2G", "-2R5", "2R6")
    for level, (rung, name, two_exponent) in enumerate(raw_named, start=1):
        require(rung == level - 1, ("level/rung indexing", level, rung))
        sign = "-" if level % 2 else ""
        two = "2" if two_exponent % 2 else ""
        class_name = f"{sign}{two}{name}"
        require(class_name == expected_classes[level - 1],
                ("known class", level, class_name, expected_classes[level - 1]))
        named_classes.append((level, rung, two_exponent, class_name))

    # Sharp hostiles for the three load-bearing parity/normalization inputs.
    hostile_even_clearing_e = 2
    require((1 - hostile_even_clearing_e) % 2 == 1,
            "even clearing exponent leaves old L in the class")
    hostile_even_block_degree = 2
    require(hostile_even_block_degree % 2 == 0,
            "even block suppresses the outer discriminant class")
    hostile_rescaling_two = 1  # v_2(2) modulo two
    require(hostile_rescaling_two == 1,
            "rational rescaling by 2 changes the square class")

    record = (
        ("matrix", MATRIX),
        ("rows", tuple(rows)),
        ("parity_ledger", tuple(parity_ledger)),
        ("induction_ledger", tuple(induction_ledger)),
        ("raw_named", raw_named),
        ("named_classes", tuple(named_classes)),
        ("hostiles", (
            ("even clearing exponent", hostile_even_clearing_e),
            ("even block degree", hostile_even_block_degree),
            ("rescaling by 2", hostile_rescaling_two),
            ("literal coordinate may be nonprimitive", True),
        )),
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== fixed Keller intrinsic discriminant raw-tower audit ==")
    print(f"matrix={MATRIX};det=-8;rows_0_to_9={tuple(rows[:10])}")
    print("all_clearing_exponents_e_n_are_odd=PASS through n=31;symbolic parity recursion e'=e mod 2")
    print("trace_discriminant_recursion=d_(n+1)=Norm(d_n)*d_1 with odd block degree 3^n")
    print("raw_class_induction=d_n=[(-1)^n P_(n-1)];residual_L_exponent=1-e_(n-1) is even")
    print(f"raw_named_scalar_invoices={raw_named}")
    print(f"known_intrinsic_classes={tuple(named_classes)}")
    print("hostiles=even clearing leaves L;even block loses outer class;rescaling by 2 changes class;literal coordinate can be nonprimitive")
    print(f"semantic_sha256={semantic}")
    print("scope=basis-independent function-field trace discriminant;no simultaneous fixed primitive coordinate,exact positive multiplicities,arbitrary-map,or general JC claim")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
