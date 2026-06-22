#!/usr/bin/env python3
"""Exact arithmetic audit for the one-large-speed interval peeler.

Workspace convention: LRC(n) means n-1 speeds and threshold 1/n.

If a seed B with max speed m has a witness tau with
    ||b tau|| >= alpha > 1/n
for every b in B, then the seed stays threshold-1/n safe on a connected arc
of length 2*(alpha-1/n)/m.  A new speed v has unsafe comb teeth of length
2/(n*v).  If the seed-safe arc is longer than one tooth, it cannot be
contained in the unsafe set of v, so B union {v} is LRC(n)-safe.
"""

from fractions import Fraction


def least_integer_strictly_greater(q: Fraction) -> int:
    return q.numerator // q.denominator + 1


def threshold(n: int, alpha: Fraction, m: int) -> Fraction:
    """Strict threshold v > m / (n*(alpha - 1/n))."""
    margin = alpha - Fraction(1, n)
    assert margin > 0
    return Fraction(m, 1) / (n * margin)


def arc_length(n: int, alpha: Fraction, m: int) -> Fraction:
    return 2 * (alpha - Fraction(1, n)) / m


def tooth_length(n: int, v: int) -> Fraction:
    return Fraction(2, n * v)


def print_case(name: str, n: int, alpha: Fraction, m: int) -> None:
    th = threshold(n, alpha, m)
    v0 = least_integer_strictly_greater(th)
    print(f"{name}")
    print(f"  n={n}, alpha={alpha}, seed max m={m}")
    print(f"  seed-safe arc length >= {arc_length(n, alpha, m)}")
    print(f"  strict threshold v > {th}")
    print(f"  least integer v certified = {v0}")
    print(f"  tooth length at v0 = {tooth_length(n, v0)}")
    print(f"  arc > tooth? {arc_length(n, alpha, m) > tooth_length(n, v0)}")
    print()


def main() -> None:
    print("ONE-LARGE-SPEED INTERVAL PEELER (codex S117b)")
    print()
    print("General lemma:")
    print("  If seed max is m and seed witness margin is alpha>1/n,")
    print("  then adding speed v is safe whenever")
    print("      v > m / (n*(alpha-1/n)).")
    print("  The proof only compares one connected seed-safe arc with one")
    print("  connected tooth of the new speed's danger comb.")
    print()

    for n in range(4, 15):
        alpha = Fraction(1, n - 1)
        th_per_m = threshold(n, alpha, 1)
        assert th_per_m == n - 1
    print("LRC(n-1) corollary checked for n=4..14:")
    print("  alpha=1/(n-1) gives the sharp symbolic threshold v > (n-1)*m.")
    print()

    print_case(
        "LRC14 generic seed from LRC13",
        n=14,
        alpha=Fraction(1, 13),
        m=13,
    )
    print_case(
        "LRC14 AP-core seed {1,...,11,13} using tau=1/12",
        n=14,
        alpha=Fraction(1, 12),
        m=13,
    )
    print_case(
        "Pure dilation hard core {b,2b,...,12b,V} with b=1",
        n=14,
        alpha=Fraction(1, 13),
        m=12,
    )

    committed = [30030, 60060, 510510]
    ap_threshold = threshold(14, Fraction(1, 12), 13)
    print("Committed AP-core lcm speeds:")
    for v in committed:
        print(f"  v={v}: certified by local AP-core peeler? {v > ap_threshold}")
    print()

    print("Tournament Analysis carriers:")
    carriers = [
        "connected_seed_safe_arc",
        "single_comb_tooth_width",
        "smaller_LRC_margin",
        "scale_ratio_gate",
        "global_component_budget",
        "runner_count_only_induction",
        "raw_runner_vertices",
    ]
    for i, carrier in enumerate(carriers, start=1):
        print(f"  {i}. {carrier}")
    print("  Hamiltonian path is transitive: the local connected-arc carrier")
    print("  preserves the LRC predicate with less topology than the global")
    print("  component-budget route, while runner-count-only induction forgets")
    print("  the needed scale ratio.")


if __name__ == "__main__":
    main()
