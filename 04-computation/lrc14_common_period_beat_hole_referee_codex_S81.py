#!/usr/bin/env python3
"""Exact referee for THM-1216's common-period beat-hole theorem.

The defining difference beat is q=b6-b5=dQ.  If every fast speed c has
gcd(c,q)=d, its strict-danger predicate on beat numerators has reduced
period Q and exactly A(Q)=2*ceil(Q/14)-1 residues.  The defining pair has
one common mask, so only five masks matter.  All five contain residue zero;
their union therefore has size at most 5*A(Q)-4.  THM-1216 uses the next
integer B(Q)=5*A(Q)-3 as a uniform consecutive-block escape threshold.

Every check below uses integer arithmetic.  The exhaustive mask audit is
diagnostic: the theorem itself is the cardinality argument above.

Tournament-analysis audit
-------------------------
The pairwise observable is common-zero mask intersection.  It is symmetric,
so increasing speed is only a tie-breaking gauge.  The induced tournament is
transitive with score histogram 0,1,2,3,4,5, no directed cycles, six singleton
SCCs, one Hamiltonian path, and 15 edge flips under the reverse gauge.  This
orientation proves nothing: the faithful object is the five-mask nerve with
its universal residue-zero vertex, together with the consecutive block phase.
"""

from __future__ import annotations

from itertools import combinations
from math import gcd


def require(condition: bool, message: object) -> None:
    """Always-on certificate check; unlike ``assert``, survives ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def ceil_div(n: int, d: int) -> int:
    return -((-n) // d)


def window_count(Q: int) -> int:
    """A(Q)=#{x mod Q: 14*min(x,-x)<Q}."""
    return 2 * ceil_div(Q, 14) - 1


def block_threshold(Q: int, mask_classes: int = 5) -> int:
    """B_r(Q)=2+r(A(Q)-1), one beyond the common-zero union cap."""
    require(1 <= mask_classes <= 5, (Q, mask_classes))
    return 2 + mask_classes * (window_count(Q) - 1)


def danger_mask(Q: int, unit: int) -> frozenset[int]:
    require(Q >= 1 and gcd(unit, Q) == 1, (Q, unit))
    return frozenset(
        p for p in range(Q)
        if 14 * min((unit * p) % Q, (-(unit * p)) % Q) < Q
    )


def euler_phi(Q: int) -> int:
    return sum(gcd(unit, Q) == 1 for unit in range(1, Q + 1))


def uniform_mask_class_cap(Q: int) -> int:
    """Exact total mask-class count, capped at the five relevant masks."""
    require(Q >= 2, Q)
    return 1 if Q <= 14 else min(5, euler_phi(Q) // 2)


def phase_free_cap(N: int, Q: int) -> int:
    A = window_count(Q)
    full, remainder = divmod(N, Q)
    return full * A + min(remainder, A)


def phase_free_tail(Q: int) -> int | None:
    """First N after the last phase-free passing size; None if no tail."""
    A = window_count(Q)
    beta = Q - 5 * A
    if beta <= 0:
        return None
    m = 4 * A // beta
    last = m * Q + min(Q - 1, 5 * A - m * beta)
    return last + 1


def longest_cyclic_run(values: set[int], Q: int) -> int:
    if len(values) == Q:
        return Q
    best = current = 0
    for index in range(2 * Q):
        if index % Q in values:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return min(best, Q)


def gap_block(a: int, k: int, q: int) -> tuple[int, int]:
    denominator = 14 * a
    lo = ceil_div(q * (14 * k + 1), denominator)
    hi = q * (14 * k + 13) // denominator
    return lo, hi


def check_window_and_budget(limit: int = 10_000) -> None:
    for Q in range(1, limit + 1):
        if Q <= 500:
            brute = len(danger_mask(Q, 1))
            require(brute == window_count(Q), (Q, brute, window_count(Q)))
        if Q >= 2:
            A = window_count(Q)
            B = block_threshold(Q)
            Bphi = block_threshold(Q, uniform_mask_class_cap(Q))
            require(Bphi <= B <= Q, (Q, Bphi, B))
            require(Q - (5 * A - 4) >= 1, (Q, A))


def check_mask_unions(limit: int = 35) -> dict[int, tuple[int, int]]:
    """Exhaust all choices of at most five distinct masks for small Q."""
    thresholds: dict[int, tuple[int, int]] = {}
    for Q in range(2, limit + 1):
        masks: list[frozenset[int]] = []
        for unit in range(1, Q):
            if gcd(unit, Q) != 1:
                continue
            mask = danger_mask(Q, unit)
            if mask not in masks:
                masks.append(mask)
        A = window_count(Q)
        B = block_threshold(Q)
        expected_classes = 1 if Q <= 14 else euler_phi(Q) // 2
        require(len(masks) == expected_classes,
                (Q, len(masks), expected_classes))
        require(all(len(mask) == A and 0 in mask for mask in masks), Q)
        best_run = 0
        for size in range(1, min(5, len(masks)) + 1):
            for chosen in combinations(masks, size):
                union = set().union(*chosen)
                threshold = block_threshold(Q, size)
                run = longest_cyclic_run(union, Q)
                require(len(union) <= 1 + size * (A - 1),
                        (Q, size, len(union)))
                require(run < threshold <= B,
                        (Q, size, run, threshold, B))
                best_run = max(best_run, run)
        class_bound = block_threshold(Q, min(5, len(masks)))
        thresholds[Q] = best_run + 1, class_bound
    return thresholds


def check_mask_stabilizers(limit: int = 100) -> None:
    """Classify when two normalized unit masks are equal."""
    for Q in range(2, limit + 1):
        units = [unit for unit in range(1, Q) if gcd(unit, Q) == 1]
        stabilizer = [
            unit for unit in units
            if danger_mask(Q, unit) == danger_mask(Q, 1)
        ]
        if Q <= 14:
            require(stabilizer == units, (Q, stabilizer, units))
            require(len({danger_mask(Q, unit) for unit in units}) == 1, Q)
        else:
            require(stabilizer == [1, Q - 1], (Q, stabilizer))
            require(len({danger_mask(Q, unit) for unit in units}) == len(units) // 2,
                    Q)


def check_phase_free_classification(limit: int = 300) -> None:
    exceptional = {1, 2, 3, 4, 5, 15}
    for Q in range(1, limit + 1):
        tail = phase_free_tail(Q)
        require((tail is None) == (Q in exceptional), (Q, tail))
        if tail is None:
            for multiplier in (1, 2, 20):
                N = multiplier * Q
                require(N <= 5 * phase_free_cap(N, Q), (Q, N))
            continue
        require(tail > 1, (Q, tail))
        require(tail - 1 <= 5 * phase_free_cap(tail - 1, Q), (Q, tail))
        for N in range(tail, tail + 10 * Q + 1):
            require(N > 5 * phase_free_cap(N, Q), (Q, N, tail))


def check_q14_sharp_corollary() -> None:
    # New equality-scale example: q/a=14/6=7/3, versus THM-1204's q/a>=7.
    a = 6
    q = 14
    speeds = (9, 11, 13, 15, 17, 31)
    require(speeds[-1] - speeds[-2] == q, speeds)
    require(all(gcd(speed, q) == 1 for speed in speeds), speeds)
    blocks = []
    for k in range(a):
        lo, hi = gap_block(a, k, q)
        block = list(range(lo, hi + 1))
        require(len(block) >= 2, (k, block))
        witness = next(p for p in block if p % 14 != 0)
        require(all(
            14 * min((speed * witness) % q, (-speed * witness) % q) >= q
            for speed in speeds
        ), (k, witness))
        blocks.append((k, lo, hi, witness))
    require(blocks == [
        (0, 1, 2, 1),
        (1, 3, 4, 3),
        (2, 5, 6, 5),
        (3, 8, 9, 8),
        (4, 10, 11, 10),
        (5, 12, 13, 12),
    ], blocks)

    # Sharpness of the two-point threshold as a phase-uniform claim.  For
    # odd a and even d, the middle slow gap is centred at t=1/2.  This row
    # has a singleton beat block p=14, which is the common dangerous residue.
    for d in range(2, 402, 2):
        a, k, q = 6 * d + 1, 3 * d, 14 * d
        sharp_speeds = tuple(d * unit for unit in (9, 11, 13, 15, 17, 31))
        require(sharp_speeds[-1] - sharp_speeds[-2] == q, (d, sharp_speeds))
        require(all(gcd(speed, q) == d for speed in sharp_speeds), d)
        require(gap_block(a, k, q) == (7 * d, 7 * d),
                (a, d, k, gap_block(a, k, q)))
        require((7 * d) % 14 == 0, d)


def print_selected_table() -> None:
    selected = [1, 2, 5, 6, 10, 14, 15, 16, 17, 21, 28, 29, 42, 56]
    print("selected common-period thresholds:")
    print("  Q   A   B_5  B_phi  aligned_N  phase_tail")
    for Q in selected:
        A = window_count(Q)
        B = block_threshold(Q)
        Bphi = block_threshold(Q, uniform_mask_class_cap(Q)) if Q >= 2 else None
        aligned = A + 1 if Q >= 2 else None
        phase = phase_free_tail(Q)
        print(
            f"  {Q:2d}  {A:2d}   {B:3d}    {str(Bphi):>4s}        "
            f"{str(aligned):>4s}       {str(phase):>4s}"
        )


def print_run_table(run_thresholds: dict[int, tuple[int, int]]) -> None:
    print("exact worst-run vs distinct-mask cardinal threshold, Q=15..30:")
    print("  Q   B_run  B_phi")
    for Q in range(15, 31):
        run, cardinal = run_thresholds[Q]
        print(f"  {Q:2d}     {run:2d}      {cardinal:2d}")


def main() -> None:
    check_window_and_budget()
    run_thresholds = check_mask_unions()
    check_mask_stabilizers()
    check_phase_free_classification()
    check_q14_sharp_corollary()

    print("THM-1216 COMMON-PERIOD BEAT-HOLE REFEREE")
    print("A(Q) brute Q=1..500; formula/budget Q=1..10000: PASS")
    print("Q>=2 => B(Q)=5*A(Q)-3 <= Q, Q=2..10000: PASS")
    print("forced common-zero hole Q-(5*A-4)>=1, Q=2..10000: PASS")
    print("class-sensitive bound B_r=2+r*(A-1), r=1..5: PASS")
    print("all unions of <=5 distinct unit masks, Q=2..35: PASS")
    print("mask classes: one for Q<=14; phi(Q)/2 for Q=15..100: PASS")
    print("phase-free relaxed-cap tail formula, Q=1..300: PASS")
    print("phase-free no-tail periods: 1,2,3,4,5,15")
    print("exact common-period no-go period: 1 only")
    print("q=14 sharpened corollary: N>=2 and q/a>=7/3: PASS")
    print("q=14 equality-scale row: (a;b)=(6;9,11,13,15,17,31): PASS")
    print("q=14 sharpness family: even d<=400, (a,k,P)=(6d+1,3d,{7d}): PASS")
    print_selected_table()
    print_run_table(run_thresholds)
    print("tournament gauge: increasing-speed tie path b1->...->b6")
    print("fingerprint: scores=0,1,2,3,4,5 cycles=0 SCCs=6 paths=1 flips=15")
    print("faithful object: five-mask nerve + universal residue 0 + block phase")


if __name__ == "__main__":
    main()
