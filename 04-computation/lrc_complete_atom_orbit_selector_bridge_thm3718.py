#!/usr/bin/env python3
"""Exact checks for the THM-2452 -> THM-3713 -> THM-2531 bridge.

This does not enumerate LRC rows.  It independently checks the new linear
identity on the anchored 23-ray path cone and the quantitative denominator
arithmetic inherited from the proved theorems.
"""

from fractions import Fraction
from random import Random


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def three_site(table, s, t):
    return table[s][t] + table[(s + 1) % P][t] - 2 * table[s][(t + 1) % P]


def s_edge(table, s, t):
    return table[(s + 1) % P][t] - table[s][t]


def t_edge(table, s, t):
    return table[s][(t + 1) % P] - table[s][t]


def main() -> None:
    rng = Random(3718)

    # alpha_j and beta_j are singleton and adjacent-pair masses in the
    # anchored path chart.  The unused slots beta_0=beta_12=0 are explicit.
    alpha = [
        [[rng.randrange(0, 1000) for _ in range(P)] for _ in range(P)]
        for _ in range(P)
    ]
    beta = [[[0 for _ in range(P)] for _ in range(P)] for _ in range(P)]
    for j in range(1, 12):
        for s in range(P):
            for t in range(P):
                beta[j][s][t] = rng.randrange(0, 1000)

    h = [[[0 for _ in range(P)] for _ in range(P)] for _ in range(P)]
    gamma_plus = [[[0 for _ in range(P)] for _ in range(P)] for _ in range(P)]
    gamma_minus = [[[0 for _ in range(P)] for _ in range(P)] for _ in range(P)]
    for j in range(1, P):
        for s in range(P):
            for t in range(P):
                h[j][s][t] = alpha[j][s][t] + beta[j - 1][s][t] + beta[j][s][t]
                gamma_plus[j][s][t] = alpha[j][s][t] + beta[j][s][t]
                gamma_minus[j][s][t] = alpha[j][s][t] + beta[j - 1][s][t]

                reconstructed_beta = sum(
                    gamma_plus[k][s][t] - gamma_minus[k][s][t]
                    for k in range(1, j)
                )
                require(reconstructed_beta == beta[j - 1][s][t], "beta telescope failed")
                reconstructed_h = gamma_plus[j][s][t] + reconstructed_beta
                require(reconstructed_h == h[j][s][t], "diagonal telescope failed")

    # Every target-cell linear mask commutes with the root-index telescope.
    operators = (s_edge, t_edge, three_site)
    for operator in operators:
        for j in range(1, P):
            for s in range(P):
                for t in range(P):
                    rhs = operator(gamma_plus[j], s, t)
                    rhs += sum(
                        operator(gamma_plus[k], s, t)
                        - operator(gamma_minus[k], s, t)
                        for k in range(1, j)
                    )
                    require(operator(h[j], s, t) == rhs, "masked telescope failed")

    # Denominators in the cover-derived quantitative consequence.
    require(2**14 == 16_384, "THM-2452 bare drift denominator")
    require(2**16 == 65_536, "THM-2452 word drift denominator")
    require(2**13 == 8_192, "bare edge denominator")
    require(2**10 == 1_024, "bare three-site denominator")
    require(2**15 == 32_768, "word edge denominator")
    require(2**12 == 4_096, "word three-site denominator")
    require(23**2 == 529, "bidirectional selector loss")
    require(8_192 * 529 == 4_333_568, "selector edge denominator")
    require(1_024 * 529 == 541_696, "selector three-site denominator")
    require(32_768 * 529 == 17_334_272, "word selector edge denominator")
    require(4_096 * 529 == 2_166_784, "word selector three-site denominator")

    # A rational bookkeeping control for the adaptive coefficient branch.
    require(576 * 26_208 == 15_095_808, "adaptive colour denominator")
    require(Fraction(16, 15_095_808) == Fraction(1, 943_488), "adaptive M tariff")

    print("bridge_identity=PASS")
    print("target_masks_checked=3")
    print("root_indices_checked=12")
    print("target_cells_checked=169")
    print("bare_drift_floor=D0/16384")
    print("bare_max_edge_sq_floor=sin(pi/13)^2*D0/8192")
    print("bare_max_three_site_sq_floor=sin(pi/13)^4*D0/1024")
    print("selector_edge_sq_floor=sin(pi/13)^2*D0/4333568")
    print("selector_three_site_sq_floor=sin(pi/13)^4*D0/541696")
    print("word_selector_edge_sq_floor=mu(Q)^2*sin(pi/13)^2*D0/17334272")
    print("word_selector_three_site_sq_floor=mu(Q)^2*sin(pi/13)^4*D0/2166784")


if __name__ == "__main__":
    main()
