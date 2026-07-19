#!/usr/bin/env python3
"""Exact referee for THM-1213's peel-first/three-comb hybrid cone.

The first four-comb killer leaves a complete safe gap; THM-1097 is then
applied to the remaining three.  This script checks the resulting rational
gate, its positive-dispersion identity, and its interaction with the exact
THM-1148 primitive residual predicate through height 64.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd


HEIGHT = 64
IDENTITY_HEIGHT = 120


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def transfer_phi_pair(x_num: int, x_den: int) -> tuple[int, int] | None:
    if x_num < x_den:
        return None
    if 7 * x_num >= 13 * x_den:
        return 6, 7
    numerator = 7 * x_num - x_den
    denominator = 14 * x_den
    common = gcd(numerator, denominator)
    return numerator // common, denominator // common


def exact_transfer_gate(shape: tuple[int, int, int, int]) -> bool:
    numerator, denominator = 6, 7
    for old, new in zip(shape, shape[1:]):
        output = transfer_phi_pair(new * numerator, old * denominator)
        if output is None:
            return False
        numerator, denominator = output
    return 7 * numerator > denominator


def primitive(shape: tuple[int, int, int, int]) -> bool:
    a, b, c, d = shape
    return gcd(gcd(a, b), gcd(c, d)) == 1


def thm1148_residual(shape: tuple[int, int, int, int]) -> bool:
    a, b, c, d = shape
    return not (8 * a > 7 * d or 2 * d > a + b + c or exact_transfer_gate(shape))


def hybrid_q(shape: tuple[int, int, int, int]) -> F:
    a, b, c, d = shape
    return F(6 * (3 * d - b - c), a) - F(6 * d, b) - F(6 * d, c) - 37


def peel_two_then_two_comb(shape: tuple[int, int, int, int]) -> bool:
    """Natural next recursion: two exact peels, then THM-1094's weak gate."""
    a, b, c, d = shape
    c2 = transfer_phi_pair(6 * b, 7 * a)
    if c2 is None:
        return False
    numerator, denominator = c2
    return numerator * c * (28 * d - 7 * c) > denominator * b * (6 * d + 29 * c)


def peel_three_then_one_comb(shape: tuple[int, int, int, int]) -> bool:
    """Three exact peels, then the sharp one-comb component threshold x>3/5."""
    a, b, c, d = shape
    output = transfer_phi_pair(6 * b, 7 * a)
    if output is None:
        return False
    numerator, denominator = output
    output = transfer_phi_pair(c * numerator, b * denominator)
    if output is None:
        return False
    numerator, denominator = output
    return 5 * d * numerator > 3 * c * denominator


def dispersion_q(shape: tuple[int, int, int, int]) -> F:
    a, b, c, d = shape
    return (
        F(6 * d, a)
        - 49
        + F(6 * (b - a) * (d - b), a * b)
        + F(6 * (c - a) * (d - c), a * c)
    )


def alternate_peel_q(shape: tuple[int, int, int, int], peel: str) -> F:
    """The same THM-1097 consumer after peeling b, c, or d first."""
    a, b, c, d = shape
    if peel == "b":
        return F(6 * (3 * d - a - c), b) - F(6 * d, a) - F(6 * d, c) - 37
    if peel == "c":
        return F(6 * (3 * d - a - b), c) - F(6 * d, a) - F(6 * d, b) - 37
    if peel == "d":
        return F(6 * (3 * c - a - b), d) - F(6 * c, a) - F(6 * c, b) - 37
    raise RuntimeError("unknown peel")


def core_entry_audit() -> tuple[F, F]:
    """Check the universal saturated first transfer at every least legal scale."""
    minimum = None
    for maximum in range(8, 13):
        ell_floor = F(72, 35 * (13 * maximum + 1))
        for a in range(1, 301):
            m = (13 * maximum) // a + 1
            normalized = m * a * ell_floor
            require(normalized >= F(72, 35), "integer legal floor weakened")
            minimum = normalized if minimum is None else min(minimum, normalized)
    require(F(72, 35) > F(13, 7), "first transfer is not saturated")
    require(minimum == F(72, 35), "finite audit missed the sharp legal floor")
    return minimum, minimum - F(13, 7)


def identity_audit(limit: int = IDENTITY_HEIGHT) -> int:
    rows = 0
    for d in range(4, limit + 1):
        for a, b, c in combinations(range(1, d), 3):
            shape = a, b, c, d
            left = (
                6 * (3 * d - b - c) * b * c
                - 6 * d * a * c
                - 6 * d * a * b
                - 37 * a * b * c
            )
            right = (
                (6 * d - 49 * a) * b * c
                + 6 * (b - a) * (d - b) * c
                + 6 * (c - a) * (d - c) * b
            )
            require(left == right, "cleared Q identity failed")
            # Keep the 8.2-million-row identity audit denominator-cleared.
            # All common denominators below are positive, so these integer
            # comparisons are exactly equivalent to the Fraction identities.
            q_b_cleared = (
                6 * (3 * d - a - c) * a * c
                - 6 * d * b * c
                - 6 * d * a * b
                - 37 * a * b * c
            )
            q_c_cleared = (
                6 * (3 * d - a - b) * a * b
                - 6 * d * b * c
                - 6 * d * a * c
                - 37 * a * b * c
            )
            q_d_cleared = (
                6 * (3 * c - a - b) * a * b
                - 6 * c * b * d
                - 6 * c * a * d
                - 37 * a * b * d
            )
            require(
                left - q_b_cleared == 6 * (b - a) * (4 * d - a - b - c) * c > 0,
                "peel-a versus peel-b dominance failed",
            )
            require(
                left - q_c_cleared == 6 * (c - a) * (4 * d - a - b - c) * b > 0,
                "peel-a versus peel-c dominance failed",
            )
            require(q_d_cleared < -19 * a * b * d < 0, "peel-d negativity failed")
            if 6 * d >= 49 * a:
                require(right > 0, "clean 49/6 cone lost strictness")
            rows += 1
    return rows


def residual_census() -> dict[str, object]:
    residual = 0
    hybrid = []
    clean = []
    exact_only = []
    recursive_two = []
    recursive_one = []
    for d in range(4, HEIGHT + 1):
        for a, b, c in combinations(range(1, d), 3):
            shape = a, b, c, d
            if not primitive(shape) or not thm1148_residual(shape):
                continue
            residual += 1
            q = hybrid_q(shape)
            require(q == dispersion_q(shape), "residual Q identity failed")
            if q > 0:
                hybrid.append((q, shape))
            if 6 * d >= 49 * a:
                clean.append(shape)
                require(q > 0, "residual clean cone did not hybrid-close")
            elif q > 0:
                exact_only.append((q, shape))
            if peel_two_then_two_comb(shape):
                recursive_two.append(shape)
            if peel_three_then_one_comb(shape):
                recursive_one.append(shape)

    require(residual == 95336, "height-64 residual census changed")
    require(len(hybrid) == 484, "hybrid success count changed")
    require(len(clean) == 351, "clean cone count changed")
    require(len(exact_only) == 133, "exact-only hybrid count changed")
    require(not recursive_two, "two-peel/two-comb recursion gained a residual row")
    require(not recursive_one, "three-peel/one-comb recursion gained a residual row")
    require(min(hybrid) == (F(51, 559), (6, 39, 43, 44)), "tightest hybrid moved")
    return {
        "residual": residual,
        "hybrid": len(hybrid),
        "clean": len(clean),
        "exact_only": len(exact_only),
        "recursive_two": len(recursive_two),
        "recursive_one": len(recursive_one),
        "tightest": min(hybrid),
        "largest": max(hybrid),
    }


def near_top_rays() -> tuple[tuple[int, int, F], ...]:
    answer = []
    for a in range(3, 13):
        for d in range(a + 3, 2000):
            shape = a, d - 2, d - 1, d
            if hybrid_q(shape) > 0:
                require(thm1148_residual(shape), "named near-top ray left residual scope")
                answer.append((a, d, hybrid_q(shape)))
                break
    expected_d = (22, 31, 39, 47, 55, 63, 71, 80, 88, 96)
    require(tuple(d for _, d, _ in answer) == expected_d, "near-top thresholds moved")
    return tuple(answer)


def tournament_audit() -> tuple[tuple[int, ...], int, int]:
    # Shape vertices a<b<c<d, oriented by numerical position.  The tournament
    # is an index only; the two interiority weights are the proof data.
    scores = 3, 2, 1, 0
    cycles = 0
    reverse_flips = 6
    return scores, cycles, reverse_flips


def main() -> None:
    entry, entry_slack = core_entry_audit()
    identity_rows = identity_audit()
    census = residual_census()
    rays = near_top_rays()
    scores, cycles, flips = tournament_audit()

    print("THM-1213 peel-first sharp-three-comb hybrid referee")
    print("arithmetic=integers and fractions.Fraction")
    print("optimized_mode_guard=require() only; assert_statements=0")
    print(f"first-transfer normalized floor={entry}; slack above 13/7={entry_slack}")
    print("Q_hyb=6(3d-b-c)/a-6d/b-6d/c-37")
    print("dispersion identity=6d/a-49+6(b-a)(d-b)/(ab)+6(c-a)(d-c)/(ac)")
    print("one-peel-Q optimality=Q_a>Q_b,Q_c; Q_d<18-37<0")
    print(f"identity audit rows through d<={IDENTITY_HEIGHT}: {identity_rows}")
    print(f"height-64 THM-1148 residual census={census['residual']}")
    print(f"hybrid successes={census['hybrid']}; clean 6d>=49a={census['clean']}; exact-only={census['exact_only']}")
    print(f"deeper recursive gates on THM-1148 residual: two-peel/two-comb={census['recursive_two']}; three-peel/one-comb={census['recursive_one']}")
    print(f"tightest hybrid={census['tightest']}; largest={census['largest']}")
    print(f"near-top rays first closed (a,d,Q)={rays}")
    print("Tournament Analysis:")
    print(f"  runner/shape order scores={scores}; cycles={cycles}; SCCs=(1,1,1,1); HP=1")
    print(f"  reverse-order gauge flips={flips}")
    print("  pairwise observable=numerical shape order; switch=reverse order; tie path=a->b->c->d")
    print("  faithful object=first-gap proof state plus two weighted interiority energies")
    print("  destroyed by tournament=metric gap|Q margin|endpoint dispersion|core entry width")
    print("  challenged vertices=runners|shape coordinates|gaps|walls|interval states|proof obligations")
    print("VERDICT: Q_hyb>0 closes; in particular 6d>=49a closes uniformly at every legal scale")


if __name__ == "__main__":
    main()
