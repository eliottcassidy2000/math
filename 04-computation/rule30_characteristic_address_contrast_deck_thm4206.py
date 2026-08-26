#!/usr/bin/env python3
"""Exact Rule-30 controls for THM-4206's characteristic contrast deck.

The universal theorem is proved symbolically in the theorem document.  This
script exhausts its smallest covariance-blind hostile and computes the
slope-one contrast correlations through depth twelve by a bit-parallel truth
table.  There are no random samples or floating-point values.
"""

from fractions import Fraction


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(label)


def variable_truth_table(variable: int, variables: int) -> int:
    """Return a packed truth table whose bit i is bit `variable` of i."""
    run = 1 << variable
    period = run << 1
    length = 1 << variables
    one_block = ((1 << run) - 1) << run
    repeat_mask = ((1 << length) - 1) // ((1 << period) - 1)
    return one_block * repeat_mask


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def hostile_census() -> None:
    # A=a_0(1), B=a_0(0), C=a_1(1), with addresses (1,0,0).
    triple = [0] * 8
    pair_ab = [0] * 4
    pair_ac = [0] * 4
    pair_bc = [0] * 4
    triple_walsh = 0
    conditional_supports = {0: set(), 1: set()}
    deck_rows = []

    for x0 in (0, 1):
        for x1 in (0, 1):
            for x2 in (0, 1):
                a = x1
                b = x0
                c = rule30(x0, x1, x2)
                contrast = b ^ c
                require(contrast == (x1 | x2), "contrast formula")
                code = (a << 2) | (b << 1) | c
                triple[code] += 1
                pair_ab[(a << 1) | b] += 1
                pair_ac[(a << 1) | c] += 1
                pair_bc[(b << 1) | c] += 1
                triple_walsh += 1 if (a ^ b ^ c) == 0 else -1
                conditional_supports[x2].add(code)
                deck_rows.append((x2, a, contrast))

    require(triple == [1, 1, 1, 1, 0, 2, 2, 0], "triple histogram")
    require(pair_ab == pair_ac == [2, 2, 2, 2], "cross-address pairs")
    require(pair_bc == [1, 3, 3, 1], "same-address pair")
    require(triple_walsh == 4, "third Walsh numerator")
    require(all(len(support) == 4 for support in conditional_supports.values()),
            "two conditional innovation bits")
    require(sorted(set(deck_rows)) == [(0, 0, 0), (0, 1, 1),
                                      (1, 0, 1), (1, 1, 1)],
            "directed contrast deck")

    print(
        "hostile=addresses:(1,0,0);"
        "triple:(1,1,1,1,0,2,2,0);"
        "AB:(2,2,2,2);AC:(2,2,2,2);BC:(1,3,3,1);"
        "cov_AB:0;cov_AC:0;cov_BC:-1/2;third_walsh:1/2"
    )
    print(
        "deck=R:(A,B);N:x2;D:C_xor_B=A_or_x2;"
        "conditional_supports:(4,4);H_given_N:2;H_joint:5/2;H_D_given_R:1/2"
    )


def slope_one_contrast_count(depth: int) -> tuple[int, int]:
    # Set the shared pivot x_0=0.  The final cell is then the contrast
    # D_t=a_t(t) xor x_0 as a function of x_1,...,x_(2t).
    variables = 2 * depth
    denominator = 1 << variables
    row = [0]
    row.extend(variable_truth_table(j, variables) for j in range(variables))
    for _ in range(depth):
        row = [rule30(row[j], row[j + 1], row[j + 2])
               for j in range(len(row) - 2)]
    require(len(row) == 1, "cone closes to one cell")
    return row[0].bit_count(), denominator


def slope_one_census() -> None:
    expected = (
        Fraction(-1, 2), Fraction(1, 4), Fraction(-1, 4), Fraction(5, 32),
        Fraction(-5, 64), Fraction(77, 1024), Fraction(-141, 2048),
        Fraction(39, 512), Fraction(-3273, 65536), Fraction(2785, 131072),
        Fraction(-21759, 1048576), Fraction(27905, 2097152),
    )
    observed = []
    counts = []
    for depth in range(1, 13):
        ones, total = slope_one_contrast_count(depth)
        correlation = Fraction(total - 2 * ones, total)
        observed.append(correlation)
        counts.append((ones, total))
    require(tuple(observed) == expected, "slope-one correlation bank")
    print("slope_one_ones=" + ";".join(f"{a}/{b}" for a, b in counts))
    print("slope_one_correlations=" + ",".join(str(value) for value in observed))
    print("finite_scope=depths:1..12;alternating_sign_observed_only:true")


def main() -> None:
    hostile_census()
    slope_one_census()
    print("THM4206_PRIMARY_PASS")


if __name__ == "__main__":
    main()
