#!/usr/bin/env python3
"""Independent audit replay for THM-1097's three-comb certificate.

The primary C++ certificate represents a safe region by successively
subtracting danger teeth.  This standard-library replay uses a distinct
representation: it builds the complete rational endpoint arrangement for all
speeds at once and classifies every open cell by an exact midpoint.  It also
re-derives every guarded-bank count and checks strict positivity at the first
omitted k1, k2, and k3 boundary.

The exact carrier is the endpoint word with coordinate and owner sidecars.
Endpoint order supplies a transitive tournament (score sequence 0,...,m-1,
no directed cycles, singleton SCCs, one sorted Hamiltonian path), but that
tournament quotient forgets the metric gap lengths needed by the LRC
predicate.  This challenges runners-as-vertices: proof obligations and
wall-crossing endpoints are the faithful vertices here, and even they require
metric sidecars.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PRIMARY = (
    ROOT
    / "05-knowledge/results/"
    / "lrc14_r4_three_comb_component_exact_codex_S73.out"
)


def danger(x: F, speed: int) -> bool:
    """Whether x lies in the open radius-1/14 danger comb."""
    y = x * speed
    integer = y.numerator // y.denominator
    residue = y - integer
    return min(residue, 1 - residue) < F(1, 14)


def safe_components(speeds: tuple[int, ...]) -> list[tuple[F, F]]:
    """Full endpoint-arrangement reconstruction, with no tooth subtraction."""
    breakpoints = {F(0), F(1)}
    for speed in speeds:
        for centre in range(speed + 1):
            for sign in (-1, 1):
                point = F(14 * centre + sign, 14 * speed)
                if 0 <= point <= 1:
                    breakpoints.add(point)

    cells: list[tuple[F, F]] = []
    ordered = sorted(breakpoints)
    for left, right in zip(ordered, ordered[1:]):
        if left == right:
            continue
        midpoint = (left + right) / 2
        if all(not danger(midpoint, speed) for speed in speeds):
            if cells and cells[-1][1] == left:
                cells[-1] = (cells[-1][0], right)
            else:
                cells.append((left, right))
    return cells


def longest_components(
    speeds: tuple[int, ...],
) -> tuple[F, list[tuple[F, F]]]:
    cells = safe_components(speeds)
    length = max((right - left for left, right in cells), default=F(0))
    return length, [cell for cell in cells if cell[1] - cell[0] == length]


def floor(x: F) -> int:
    return x.numerator // x.denominator


def tail_margin(ell: F, k1: int, k2: int, k3: int) -> F:
    """The exact sufficient-tail expression F from THM-1097."""
    return (
        21 * ell * k3
        - 7 * ell * (k1 + k2)
        - 6 * (F(k3, k1) + F(k3, k2))
        - 37
    )


def parse_primary(text: str):
    row_re = re.compile(
        r"^core=\[([^]]+)\] ell=([^ ]+) finite_triples=(\d+) "
        r"failures=(\d+) (.+)$",
        re.MULTILINE,
    )
    hardest_re = re.compile(
        r"hardest_triple=\((\d+),(\d+),(\d+)\) min_7k3L=(\S+)"
    )
    rows = {}
    for match in row_re.finditer(text):
        core = tuple(map(int, match.group(1).split(",")))
        tail = match.group(5)
        hardest = None
        hardest_match = hardest_re.search(tail)
        if hardest_match:
            hardest = (
                tuple(map(int, hardest_match.group(1, 2, 3))),
                F(hardest_match.group(4)),
            )
        else:
            assert tail == "hardest_triple=none min_7k3L=analytic"
        rows[core] = {
            "ell": F(match.group(2)),
            "count": int(match.group(3)),
            "failures": int(match.group(4)),
            "hardest": hardest,
        }
    return rows


def main() -> None:
    primary = Path(sys.argv[1]) if len(sys.argv) == 2 else DEFAULT_PRIMARY
    if len(sys.argv) > 2:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} [primary-output]")
    primary_bytes = primary.read_bytes()
    primary_text = primary_bytes.decode("utf-8")
    primary_rows = parse_primary(primary_text)

    cores = list(combinations(range(1, 13), 9))
    atlas = {}
    minimum_ell = None
    minimum_cores = []
    total = 0
    per_core_counts = {}

    first_omitted_k1_checks = 0
    first_omitted_k2_checks = 0
    first_omitted_k3_checks = 0
    minimum_omitted_k3_margin = None
    minimum_b_denominator = None
    minimum_c_denominator = None

    maximum_formal_cutoffs = [0, 0, 0]
    maximum_counted_indices = [0, 0, 0]

    for core in cores:
        ell, intervals = longest_components(core)
        atlas[core] = (ell, intervals)
        if minimum_ell is None or ell < minimum_ell:
            minimum_ell = ell
            minimum_cores = [core]
        elif ell == minimum_ell:
            minimum_cores.append(core)

        lo = 13 * max(core) + 1
        k1_max = floor(F(7) / ell) + 1
        maximum_formal_cutoffs[0] = max(maximum_formal_cutoffs[0], k1_max)

        # The first omitted integer is strictly beyond ell*k1=7.  Equality
        # would already suffice because 1<k2/k1<k3/k1 supplies strictness.
        k1_omitted = k1_max + 1
        assert ell * k1_omitted > 7
        first_omitted_k1_checks += 1

        core_count = 0
        for k1 in range(lo, k1_max + 1):
            x = ell * k1
            b_denominator = 14 * x - 6
            assert b_denominator > 0
            minimum_b_denominator = (
                b_denominator
                if minimum_b_denominator is None
                else min(minimum_b_denominator, b_denominator)
            )
            k2_cut = F(k1) * (7 * x + 43) / b_denominator
            k2_max = floor(k2_cut) + 1
            maximum_formal_cutoffs[1] = max(maximum_formal_cutoffs[1], k2_max)

            # At k3=k2 this is precisely the solved k2 cutoff.  The k3
            # coefficient is positive, hence every ordered triple beyond it
            # is analytic as well.
            k2_omitted = k2_max + 1
            assert k2_omitted > k2_cut
            assert tail_margin(ell, k1, k2_omitted, k2_omitted) > 0
            assert 21 * ell - F(6, k1) - F(6, k2_omitted) > 0
            first_omitted_k2_checks += 1

            for k2 in range(k1 + 1, k2_max + 1):
                c_denominator = 21 * ell - F(6, k1) - F(6, k2)
                assert c_denominator > 0
                minimum_c_denominator = (
                    c_denominator
                    if minimum_c_denominator is None
                    else min(minimum_c_denominator, c_denominator)
                )
                k3_cut = (7 * ell * (k1 + k2) + 37) / c_denominator
                k3_max = floor(k3_cut) + 1
                maximum_formal_cutoffs[2] = max(
                    maximum_formal_cutoffs[2], k3_max
                )

                k3_omitted = k3_max + 1
                assert k3_omitted > k3_cut
                margin = tail_margin(ell, k1, k2, k3_omitted)
                assert margin > 0
                minimum_omitted_k3_margin = (
                    margin
                    if minimum_omitted_k3_margin is None
                    else min(minimum_omitted_k3_margin, margin)
                )
                first_omitted_k3_checks += 1

                count = max(0, k3_max - k2)
                core_count += count
                total += count
                if count:
                    maximum_counted_indices[0] = max(
                        maximum_counted_indices[0], k1
                    )
                    maximum_counted_indices[1] = max(
                        maximum_counted_indices[1], k2
                    )
                    maximum_counted_indices[2] = max(
                        maximum_counted_indices[2], k3_max
                    )
        per_core_counts[core] = core_count

    # Match all 220 primary atlas/count rows, then independently reconstruct
    # every nonempty core's reported hardest row by the combined arrangement.
    assert set(primary_rows) == set(cores)
    replayed_hardest = 0
    analytic_only_rows = 0
    replay_global = None
    for core in cores:
        primary = primary_rows[core]
        ell = atlas[core][0]
        assert primary["ell"] == ell
        assert primary["count"] == per_core_counts[core]
        assert primary["failures"] == 0
        if primary["hardest"] is None:
            assert per_core_counts[core] == 0
            analytic_only_rows += 1
            continue
        triple, reported_metric = primary["hardest"]
        length, intervals = longest_components(core + triple)
        metric = 7 * triple[2] * length
        assert metric == reported_metric
        assert metric > 1
        replayed_hardest += 1
        row = (metric, core, triple, length, intervals)
        if replay_global is None or row[:3] < replay_global[:3]:
            replay_global = row

    assert len(cores) == 220
    assert total == 39_778_595
    assert minimum_ell == F(11, 1008)
    assert minimum_cores == [(1, 2, 3, 7, 8, 9, 10, 11, 12)]
    assert maximum_formal_cutoffs == [642, 642, 642]
    assert maximum_counted_indices == [639, 640, 641]
    assert first_omitted_k1_checks == 220
    assert first_omitted_k2_checks == 23_589
    assert first_omitted_k3_checks == 1_238_741
    assert minimum_omitted_k3_margin == F(28_429, 182_160)
    assert replayed_hardest == 200
    assert analytic_only_rows == 20
    assert replay_global is not None

    metric, core, triple, length, intervals = replay_global
    assert core == (1, 2, 3, 6, 7, 8, 10, 11, 12)
    assert triple == (162, 168, 174)
    assert length == F(1, 522)
    assert metric == F(7, 3)
    assert intervals == [
        (F(53, 252), F(517, 2436)),
        (F(1919, 2436), F(199, 252)),
    ]

    print("THM-1097 independent guard/endpoint replay")
    print("method=combined breakpoint arrangement plus exact midpoint classification")
    print(f"primary_output_sha256={sha256(primary_bytes).hexdigest()}")
    print(f"cores={len(cores)}")
    print(f"min_core_ell={minimum_ell}")
    for minimum_core in minimum_cores:
        print(f"min_ell_core={list(minimum_core)}")
    print(f"finite_triples_from_guards={total}")
    print(f"maximum_formal_cutoffs={tuple(maximum_formal_cutoffs)}")
    print(f"maximum_counted_indices={tuple(maximum_counted_indices)}")
    print(
        "first_omitted_checks="
        f"(k1:{first_omitted_k1_checks},"
        f"k2:{first_omitted_k2_checks},"
        f"k3:{first_omitted_k3_checks})"
    )
    print(f"minimum_b_denominator={minimum_b_denominator}")
    print(f"minimum_c_denominator={minimum_c_denominator}")
    print(f"minimum_first_omitted_k3_F={minimum_omitted_k3_margin}")
    print(f"replayed_nonempty_core_hardest_rows={replayed_hardest}")
    print(f"analytic_only_core_rows={analytic_only_rows}")
    print(f"finite_bank_extremal_core={list(core)}")
    print(f"finite_bank_extremal_triple={triple}")
    print(f"finite_bank_extremal_longest_component={length}")
    print(f"finite_bank_extremal_7k3L={metric}")
    print(f"finite_bank_extremal_R={1 / metric}")
    for left, right in intervals:
        print(f"finite_bank_extremal_interval=[{left},{right}]")
    print("tournament_carrier=endpoint order with exact coordinate/owner sidecars")
    print("tournament_fingerprint=transitive; no cycles; singleton SCCs; one sorted Hamiltonian path")
    print("challenged_assumption=runners alone are not faithful vertices for component length")
    print("certificate=PASS")


if __name__ == "__main__":
    main()
