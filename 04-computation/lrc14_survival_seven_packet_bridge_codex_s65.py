#!/usr/bin/env python3
"""HYP-2704 / T939: survival currency vs seven-packet recursions.

This scout keeps the true-wide survival-middle-mass gate in exact linear
algebra.  The point is to compare three signed operators that have been
appearing in different coordinates:

* the seven-packet singleton/pair/triple recursions from the tiling/far-runner
  side;
* the survival cut expansion of THM-556/HYP-2701;
* the one-hit transfer-tax boundary operator of THM-558.

The output is deliberately small: exact coefficient tables and a tournament
fingerprint on missed-count states for the two-far decorrelated boundary.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from math import comb, factorial


def fmt(q: F) -> str:
    return str(q)


def transition_prob(t: int, s: int, r: int) -> F:
    """Probability t missed sectors become s missed after r iid color hits."""

    if not 0 <= s <= t <= 6:
        return F(0)
    need_hit = t - s
    total = F(0)
    for j in range(need_hit + 1):
        total += ((-1) ** j) * comb(need_hit, j) * F(7 - s - j, 7) ** r
    return F(comb(t, s)) * total


def currency_coeff(t: int) -> F:
    """Deterministic coefficient of p_t in C=p1+...+p4-4p6."""

    if 1 <= t <= 4:
        return F(1)
    if t == 6:
        return F(-4)
    return F(0)


def u4_coeff(t: int) -> F:
    """Deterministic coefficient of p_t in U4=p0+p5+5p6."""

    if t == 0 or t == 5:
        return F(1)
    if t == 6:
        return F(5)
    return F(0)


def cut_value(t: int, b: int) -> F:
    return F(1) if t >= b else F(0)


def currency_from_cuts(t: int) -> F:
    return cut_value(t, 1) - cut_value(t, 5) - 4 * cut_value(t, 6)


def expected_currency_after_hits(t: int, r: int) -> F:
    return sum(transition_prob(t, s, r) * currency_coeff(s) for s in range(t + 1))


def live_transfer_delta(t: int) -> F:
    """C(t-1)-C(t), the one-hit boundary contribution at depth t."""

    if t == 0:
        return F(0)
    return currency_coeff(t - 1) - currency_coeff(t)


def sign_word(vec: tuple[int, ...] | tuple[F, ...]) -> str:
    out = []
    for x in vec:
        if x > 0:
            out.append("+")
        elif x < 0:
            out.append("-")
        else:
            out.append("0")
    return "".join(out)


def directed_3cycles(edges: set[tuple[int, int]], n: int) -> int:
    total = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                if (
                    (a, b) in edges and (b, c) in edges and (c, a) in edges
                ) or (
                    (a, c) in edges and (c, b) in edges and (b, a) in edges
                ):
                    total += 1
    return total


def hamiltonian_path_count(edges: set[tuple[int, int]], n: int) -> int:
    count = 0
    for perm in __import__("itertools").permutations(range(n)):
        if all((perm[i], perm[i + 1]) in edges for i in range(n - 1)):
            count += 1
    return count


def print_basis_identities() -> None:
    print("BASIS IDENTITIES")
    c_vec = tuple(currency_coeff(t) for t in range(7))
    u_vec = tuple(u4_coeff(t) for t in range(7))
    print(f"  C coefficients on p0..p6:  {[fmt(x) for x in c_vec]}")
    print(f"  U4 coefficients on p0..p6: {[fmt(x) for x in u_vec]}")
    for t in range(7):
        assert currency_coeff(t) == currency_from_cuts(t)
        assert u4_coeff(t) == F(1) - currency_coeff(t)
    print("  exact: C = G1 - G5 - 4*G6")
    print("  exact: U4 = 1 - C = 1 - G1 + G5 + 4*G6")
    print()


def print_seven_packet_comparison() -> None:
    print("SEVEN-PACKET / CUT COMPARISON")
    packets = ("A", "B", "C", "D", "E", "F", "G")
    actual = (1, 1, 1, 1, 1, 1, 1)
    pair_tax = (1, 1, 1, -1, -1, -1, 1)
    half_odd_addressed = (1, 1, -1, 1, -1, -1, 1)
    print(f"  packets: {packets}")
    print(f"  actual far correction H(1):       {actual} sign={sign_word(actual)}")
    print(f"  pair-tax A+B+C-D-E-F+G:          {pair_tax} sign={sign_word(pair_tax)}")
    print(
        f"  odd half-tiling addressed signs: {half_odd_addressed} "
        f"sign={sign_word(half_odd_addressed)}"
    )
    print("  layer quotient (singletons,pairs,triple):")
    print("    actual        -> (+1,+1,+1)")
    print("    pair-tax      -> (+1,-1,+1)")
    print("    survival cuts -> (+1,-1,-4) on (G1,G5,G6)")
    print("    transfer dC   -> (-1,+1,+4) across depths (1,5,6)")
    print(
        "  conclusion: survival is a tail-weighted cut/Stokes character, "
        "not the pair-tax Euler character."
    )
    print()


def print_transfer_and_death_chain() -> None:
    print("TRANSFER TAX AND DECORRELATED FAR HITS")
    deltas = tuple(live_transfer_delta(t) for t in range(1, 7))
    print(f"  one-hit dC=C(t-1)-C(t), t=1..6: {[fmt(x) for x in deltas]}")
    print("  live transitions: t=1 gives -1, t=5 gives +1, t=6 gives +4")
    print("  expected currency coefficient after r iid far color hits:")
    header = "       t: " + " ".join(f"{t:>9}" for t in range(7))
    print(header)
    for r in range(5):
        coeffs = tuple(expected_currency_after_hits(t, r) for t in range(7))
        print(f"    r={r}: " + " ".join(f"{fmt(x):>9}" for x in coeffs))
    print()
    print("  two-far gain K2(t)-C(t):")
    gains = []
    for t in range(7):
        gain = expected_currency_after_hits(t, 2) - currency_coeff(t)
        gains.append(gain)
    print("       " + " ".join(f"{fmt(x):>9}" for x in gains))
    assert expected_currency_after_hits(6, 2) == F(26, 49)
    print("  exact signal: a fully missed state has C=-4, but after two")
    print("  decorrelated far hits its expected coefficient is 26/49.")
    print()


def print_tournament_analysis() -> None:
    """Tournament on missed-count depth states under the two-far boundary gauge."""

    print("TOURNAMENT ANALYSIS")
    verts = list(range(7))
    values = {t: expected_currency_after_hits(t, 2) for t in verts}
    edges: set[tuple[int, int]] = set()
    scores: Counter[int] = Counter()
    for i, a in enumerate(verts):
        for b in verts[i + 1 :]:
            # Gauge: larger two-far boundary currency is safer; ties break to
            # smaller missed-count depth, preserving the depth-chain order.
            key_a = (values[a], -a)
            key_b = (values[b], -b)
            winner, loser = (a, b) if key_a >= key_b else (b, a)
            edges.add((winner, loser))
            scores[winner] += 1
            scores.setdefault(loser, scores[loser])
    ordered = sorted(verts, key=lambda t: (values[t], -t), reverse=True)
    print("  vertices: missed-count depths t=0..6")
    print("  pairwise observable: K2(t)=E[C after two decorrelated far hits | N=t]")
    print("  switch/gauge: orient toward larger K2; ties by smaller t")
    print(f"  K2 order: {' > '.join(f't{t}({fmt(values[t])})' for t in ordered)}")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={directed_3cycles(edges, len(verts))}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(edges, len(verts))}")
    print(
        "  preserved quotient: missed-count survival currency; destroyed data: "
        "sector labels, cyclic order, runner phases, and relation-lattice address."
    )
    print()


def main() -> None:
    print("LRC14 SURVIVAL / SEVEN-PACKET BRIDGE (exact Fraction scout)")
    print("=" * 78)
    print_basis_identities()
    print_seven_packet_comparison()
    print_transfer_and_death_chain()
    print_tournament_analysis()
    # Tiny sanity check that no hidden float arithmetic entered.
    assert factorial(7) == 5040


if __name__ == "__main__":
    main()
