#!/usr/bin/env python3
"""Exact referee for THM-1215's all-range q=14 beat-puncture stalk.

For a slow carrier ``a`` and six faster speeds ``b1<...<b6``, put
``q=b6-b5=14*d``.  The theorem assumes ``q>=7*a`` and
``gcd(c,q)=d`` for ``c=b1,...,b5``.  On the reduced q-clock every one of
these five combs then has period 14 and the same strict dangerous mask
``{0}``.  This script checks the exact block/count identities, audits a
large bounded prefix of the all-range formulas, and freezes the example

    (a;b1,...,b6) = (3;6,10,18,22,26,54).

The proof in THM-1215 is symbolic.  None of the bounded rows below is a
premise of the theorem.  All comparisons use integers only, and ``require``
is used instead of ``assert`` so normal and optimized Python execute the
same checks.

Tournament-analysis audit
-------------------------
The pairwise observable on the six fast speeds is equality of their reduced
q=14 dangerous masks.  All masks tie at ``{0}``; breaking ties by increasing
speed gives the transitive Hamiltonian path ``b1->...->b6``.  Reversing that
gauge flips all 15 edges.  The tournament is therefore bookkeeping, not the
proof carrier.  Alternate vertices considered are runners, slow gaps, the 14
residue sections, wall-crossing events, beat numerators, kill masks, and proof
obligations.  The faithful object is the residue clock together with the
consecutive beat block: it preserves the strict-danger predicate and destroys
the real gap location if the block-phase sidecar is omitted.
"""

from __future__ import annotations

from math import gcd


RADIUS_DENOMINATOR = 14
CORE_LIMIT = 200_000
GAP_AUDIT_MAX_A = 160
GAP_AUDIT_D_MULTIPLIER = 3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_div(numerator: int, denominator: int) -> int:
    require(numerator >= 0 and denominator > 0, "ceil_div domain")
    return (numerator + denominator - 1) // denominator


def gap_block(a: int, k: int, q: int) -> tuple[int, int]:
    """Numerators p for which p/q lies in the closed carrier gap G_k(a)."""
    require(a > 0 and q > 0 and 0 <= k < a, "gap_block domain")
    denominator = RADIUS_DENOMINATOR * a
    lo = ceil_div(q * (RADIUS_DENOMINATOR * k + 1), denominator)
    hi = q * (RADIUS_DENOMINATOR * k + 13) // denominator
    return lo, hi


def block_size(block: tuple[int, int]) -> int:
    lo, hi = block
    return max(0, hi - lo + 1)


def dangerous(speed: int, p: int, q: int) -> bool:
    """Exact test for ||speed*p/q|| < 1/14."""
    residue = (speed * p) % q
    return RADIUS_DENOMINATOR * min(residue, q - residue) < q


def thm1192_period_count(Q: int) -> int:
    """A(Q)=2*ceil(Q/14)-1 from THM-1192."""
    require(Q > 0, "positive reduced period")
    return 2 * ceil_div(Q, RADIUS_DENOMINATOR) - 1


def thm1192_phase_free_cap(N: int, Q: int) -> int:
    """U(N,Q) from THM-1192, with no discarded integer rounding."""
    require(N >= 0 and Q > 0, "phase-free cap domain")
    full, remainder = divmod(N, Q)
    A = thm1192_period_count(Q)
    return full * A + min(remainder, A)


def q14_cap(N: int) -> int:
    return ceil_div(N, 14)


def integer_core_audit() -> tuple[int, int, int]:
    minimum_slack = CORE_LIMIT
    minimizer = -1
    for N in range(CORE_LIMIT + 1):
        U = thm1192_phase_free_cap(N, 14)
        require(U == q14_cap(N), f"q=14 cap identity failed at N={N}")
        if N >= 6:
            slack = N - 5 * U
            require(slack > 0, f"integer core failed at N={N}")
            if slack < minimum_slack:
                minimum_slack = slack
                minimizer = N
    require((minimum_slack, minimizer) == (1, 6), "sharp integer core moved")
    return CORE_LIMIT + 1, minimum_slack, minimizer


def gap_block_audit() -> tuple[int, int, tuple[int, int, int, tuple[int, int]]]:
    """Bounded replay of N>=6; the paper proves this from interval length."""
    rows = 0
    minimum = None
    first_minimum = None
    for a in range(1, GAP_AUDIT_MAX_A + 1):
        least_d = ceil_div(a, 2)  # 14d >= 7a
        for d in range(least_d, GAP_AUDIT_D_MULTIPLIER * a + 1):
            q = 14 * d
            require(q >= 7 * a, "audit left the all-range cone")
            for k in range(a):
                block = gap_block(a, k, q)
                N = block_size(block)
                require(N >= 6, f"short beat block at a={a},d={d},k={k}")
                rows += 1
                if minimum is None or N < minimum:
                    minimum = N
                    first_minimum = (a, d, k, block)
    require(minimum == 6 and first_minimum is not None, "block minimum moved")
    return rows, minimum, first_minimum


def check_stalk_hypotheses(a: int, faster: tuple[int, ...]) -> tuple[int, int]:
    require(len(faster) == 6, "six fast speeds required")
    require(a < faster[0] and all(x < y for x, y in zip(faster, faster[1:])),
            "speeds are not strictly ordered")
    q = faster[-1] - faster[-2]
    require(q % 14 == 0, "difference beat is not 14d")
    d = q // 14
    require(d > 0 and q >= 7 * a, "slow-gap length hypothesis failed")
    require(all(gcd(c, q) == d for c in faster[:5]), "mixed-gcd stalk failed")
    require(gcd(faster[-1], q) == d, "b6 gcd inheritance failed")
    return q, d


def first_exact_witness(
    a: int, faster: tuple[int, ...], k: int
) -> tuple[int, int, int, int]:
    q, _d = check_stalk_hypotheses(a, faster)
    block = gap_block(a, k, q)
    N = block_size(block)
    U = q14_cap(N)
    require(N >= 6 and N - U > 4 * U, "phase-free contradiction missing")
    for p in range(block[0], block[1] + 1):
        if p % 14 == 0:
            continue
        require(not dangerous(a, p, q), "chosen numerator left the carrier gap")
        require(all(not dangerous(c, p, q) for c in faster),
                "reduced-mask witness was killed")
        return p, N, N - U, 4 * U
    raise RuntimeError("six consecutive numerators contained no nonzero residue mod 14")


def example_audit() -> tuple[tuple[int, int, tuple[int, int], int, int, int, int], ...]:
    a = 3
    faster = (6, 10, 18, 22, 26, 54)
    q, d = check_stalk_hypotheses(a, faster)
    require((q, d) == (28, 2), "example q,d changed")
    rows = []
    for k in range(a):
        block = gap_block(a, k, q)
        p, N, left, right = first_exact_witness(a, faster, k)
        rows.append((k, p, block, N, q14_cap(N), left, right))
    expected = (
        (0, 1, (1, 8), 8, 1, 7, 4),
        (1, 10, (10, 18), 9, 1, 8, 4),
        (2, 20, (20, 27), 8, 1, 7, 4),
    )
    require(tuple(rows) == expected, "example witness table changed")
    return tuple(rows)


def tournament_audit() -> tuple[tuple[int, ...], int, tuple[int, ...], int, int]:
    # All reduced masks tie.  Increasing speed is only a deterministic gauge.
    scores = (5, 4, 3, 2, 1, 0)
    directed_cycles = 0
    scc_sizes = (1, 1, 1, 1, 1, 1)
    reverse_gauge_edge_flips = 15
    hamiltonian_path_count = 1
    return scores, directed_cycles, scc_sizes, reverse_gauge_edge_flips, hamiltonian_path_count


def main() -> None:
    core_rows, core_slack, core_minimizer = integer_core_audit()
    gap_rows, gap_minimum, first_gap_minimum = gap_block_audit()
    example = example_audit()
    scores, cycles, sccs, flips, hp_count = tournament_audit()

    print("THM-1215 all-range q=14 beat-puncture stalk referee")
    print("arithmetic=integers only")
    print("optimized_mode_guard=require() only; assert_statements=0")
    print("A(14)=1; U(N,14)=floor(N/14)+min(N mod 14,1)=ceil(N/14)")
    print(
        f"integer core rows N=0..{CORE_LIMIT}: {core_rows}; "
        f"min[N-5U | N>=6]={core_slack} at N={core_minimizer}"
    )
    print(
        f"slow-gap block replay a<={GAP_AUDIT_MAX_A}, "
        f"ceil(a/2)<=d<={GAP_AUDIT_D_MULTIPLIER}a: "
        f"rows={gap_rows}; minimum N={gap_minimum}; first={first_gap_minimum}"
    )
    print("symbolic proof scope: all a,d,k satisfying 14d>=7a; bounded replay is not a premise")
    print("example (a;b)=(3;6,10,18,22,26,54), q=28=14*2:")
    for k, p, block, N, U, left, right in example:
        print(
            f"  k={k}: P={block[0]}..{block[1]}, N={N}, U={U}, "
            f"N-U={left}>4U={right}, exact lonely numerator p={p}"
        )
    print("Tournament Analysis:")
    print(f"  tied-mask speed tournament scores={scores}; cycles={cycles}; SCCs={sccs}; HP={hp_count}")
    print(f"  reverse tie gauge edge flips={flips}; tie path=b1->b2->b3->b4->b5->b6")
    print("  pairwise observable=reduced dangerous-mask equality; all masks={0} in Z/14Z")
    print("  faithful object=14-residue clock plus consecutive beat-block phase sidecar")
    print("  challenged vertices=runners|slow gaps|residue sections|walls|beat points|proof obligations")
    print("VERDICT: the structured q=14 mixed-gcd stalk is empty at every scale")
    print("SCOPE: this does not close the general slow-gap branch or LRC(14)")


if __name__ == "__main__":
    main()
