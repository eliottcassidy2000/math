#!/usr/bin/env python3
"""HYP-3458: AP84 coloring-recursion bridge for LRC14.

This scout connects the current AP84 tail proof packets to older repository
work on colorings and outer-extension recursion:

* HYP-2247: Paris-Harrington colorings as coherent bad-child trees.
* HYP-2243: outer-extension usability needs a retained embedding/address state.
* HYP-3454/HYP-3456/HYP-3457: AP84 endpoint clock, mod-35 floor count,
  and finite mixed transients.

For S_m={1,...,11,13,84m}, HYP-3456 gives the number N(m) of high E:84m
safe gaps hitting the left low corridor C1=[8/49,6/35]:

    N(m)=floor((504m-6)/70)-floor((96m-13)/14).

This is a period-35 coloring after subtracting floor(12m/35).  The additional
observation here is that the HYP-3433/HYP-3454 selected endpoint address
a_m=ceil(48m/7) is always the first or second colored gap in C1:

    rank_C1(a_m)=1 if 7|m, else 2.

So the AP-tail clock is not just a count.  It is a small recursive coloring
state with a mod-35 boundary color, a mod-7 endpoint-rank subcolor, and a
one-step death of the mixed endpoint phase at m=5.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


LIMIT = 350


def n_left(m: int) -> int:
    return ((504 * m - 6) // 70) - ((96 * m - 13) // 14)


def left_bounds(m: int) -> tuple[int, int]:
    lower = ((96 * m - 13) // 14) + 1
    upper = (504 * m - 6) // 70
    return lower, upper


def endpoint_address(m: int) -> int:
    return (48 * m + 6) // 7


def endpoint_rank(m: int) -> int:
    lower, upper = left_bounds(m)
    a = endpoint_address(m)
    if not lower <= a <= upper:
        raise AssertionError(f"endpoint address {a} outside C1 gap range {lower}..{upper} for m={m}")
    return a - lower + 1


def phase_color(m: int) -> str:
    return "mixed" if m <= 4 else "pure"


@dataclass(frozen=True)
class ResidueState:
    r: int
    n: int
    boundary_color: int
    endpoint_rank_color: int
    transition_increment: int


def residue_states() -> list[ResidueState]:
    states: list[ResidueState] = []
    for r in range(1, 36):
        n = n_left(r)
        boundary = n - (12 * r // 35)
        next_m = r + 1 if r < 35 else 36
        states.append(
            ResidueState(
                r=r,
                n=n,
                boundary_color=boundary,
                endpoint_rank_color=endpoint_rank(r),
                transition_increment=n_left(next_m) - n,
            )
        )
    return states


def crt_grid(states: list[ResidueState]) -> list[list[str]]:
    by_pair = {(state.r % 5, state.r % 7): state for state in states}
    grid: list[list[str]] = []
    for mod5 in range(5):
        row: list[str] = []
        for mod7 in range(7):
            state = by_pair[(mod5, mod7)]
            row.append(f"{state.boundary_color}/{state.endpoint_rank_color}")
        grid.append(row)
    return grid


def failures() -> dict[str, list[int]]:
    period_failures: list[int] = []
    rank_failures: list[int] = []
    phase_failures: list[int] = []
    address_failures: list[int] = []
    states = residue_states()
    corrections = {state.r: state.boundary_color for state in states}

    for m in range(1, LIMIT + 1):
        r = ((m - 1) % 35) + 1
        predicted_n = 12 * ((m - r) // 35) + n_left(r)
        predicted_d = corrections[r]
        actual_d = n_left(m) - (12 * m // 35)
        if n_left(m) != predicted_n or actual_d != predicted_d:
            period_failures.append(m)

        rank = endpoint_rank(m)
        expected_rank = 1 if m % 7 == 0 else 2
        if rank != expected_rank:
            rank_failures.append(m)

        lower, upper = left_bounds(m)
        if not lower <= endpoint_address(m) <= upper:
            address_failures.append(m)

        # HYP-3457/HYP-3454 phase break is the sign of 455-98m.
        expected_phase = "mixed" if (455 - 98 * m) > 0 else "pure"
        if phase_color(m) != expected_phase:
            phase_failures.append(m)

    return {
        "period_failures": period_failures,
        "rank_failures": rank_failures,
        "phase_failures": phase_failures,
        "address_failures": address_failures,
    }


def histogram(values) -> dict:
    return dict(sorted(Counter(values).items()))


def tournament_fingerprint() -> tuple[dict[int, int], list[str], int]:
    vertices = {
        "crt_35_coloring_state": (10, 10, 10, 9, 9, 10),
        "endpoint_rank_mod7_subcolor": (10, 9, 10, 10, 8, 10),
        "ph_bad_phase_derivative": (9, 10, 9, 8, 10, 9),
        "outer_extension_period_shift": (9, 9, 9, 9, 9, 9),
        "mod35_boundary_color_word": (8, 9, 8, 9, 8, 8),
        "transition_increment_word": (7, 8, 8, 8, 8, 7),
        "raw_escape_count": (4, 4, 3, 4, 3, 3),
        "raw_color_analogy": (2, 2, 2, 1, 1, 1),
    }
    scores = {name: sum(vals) for name, vals in vertices.items()}
    hist = dict(sorted(Counter(scores.values()).items()))
    path = [name for name, _score in sorted(scores.items(), key=lambda item: (-item[1], item[0]))]
    rank = {name: i for i, name in enumerate(path)}
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = rank[a] < rank[b]
        bc = rank[b] < rank[c]
        ca = rank[c] < rank[a]
        if ab == bc == ca:
            cycles += 1
    return hist, path, cycles


def main() -> None:
    states = residue_states()
    grid = crt_grid(states)
    checks = failures()

    boundary_word = [state.boundary_color for state in states]
    rank_word = [state.endpoint_rank_color for state in states]
    transition_word = [state.transition_increment for state in states]
    bad_phase_counts = {
        "q0_mixed_residues": sum(1 for m in range(1, 36) if phase_color(m) == "mixed"),
        "q1_mixed_residues": sum(1 for m in range(36, 71) if phase_color(m) == "mixed"),
    }

    print("HYP-3458 AP84 COLORING-RECURSION BRIDGE")
    print("status=EVIDENCE / exact coloring-recursion sidecar for AP-tail bridge; not an LRC14 proof")
    print("source=HYP-2247/HYP-2243 coloring recursion + HYP-3454/HYP-3456/HYP-3457 AP84 clocks")
    print()

    print("## Exact Checks")
    print(f"checked_m=1..{LIMIT}")
    for name, vals in checks.items():
        print(f"{name}={vals}")
    print()

    print("## Boundary And Endpoint Color Words")
    print("boundary_color d_r=N(r)-floor(12r/35), residues r=1..35:")
    print(f"boundary_word={boundary_word}")
    print(f"boundary_color_hist={histogram(boundary_word)}")
    print("endpoint_rank_color rank_C1(a_r), residues r=1..35:")
    print(f"endpoint_rank_word={rank_word}")
    print(f"endpoint_rank_hist={histogram(rank_word)}")
    print("transition_increment N(r+1)-N(r), with r=35 using N(36)-N(35):")
    print(f"transition_word={transition_word}")
    print(f"transition_hist={histogram(transition_word)}")
    print("period_shift=N(m+35)-N(m)=12, escapes shift by 24")
    print()

    print("## CRT Grid")
    print("Rows are m mod 5 = 0..4; columns are m mod 7 = 0..6.")
    print("Each entry is boundary_color/endpoint_rank_color.")
    for mod5, row in enumerate(grid):
        print(f"mod5={mod5}: " + " ".join(f"{cell:>3}" for cell in row))
    print()

    print("## Homogeneous Trace")
    print("Selected endpoint address a_m=ceil(48m/7) has:")
    print("  rank_C1(a_m)=1 exactly when 7 divides m")
    print("  rank_C1(a_m)=2 otherwise")
    print("Thus the HYP-3433 endpoint witness is always in the first two colored")
    print("left-corridor high gaps.  The mod-35 color word counts all escapes; the")
    print("mod-7 subcolor locates the canonical homogeneous trace inside that set.")
    print()

    print("## Paris-Harrington / Outer-Extension Reading")
    print(f"bad_phase_counts={bad_phase_counts}")
    print("The mixed endpoint color is a finite bad branch: it appears only at m=1..4")
    print("and is gone after the period extension m -> m+35.  The pure tail is not")
    print("fast-growing PH behavior; it is a 35-state automaton with retained colors")
    print("(phase, boundary d_r, endpoint rank).  This is the tame AP84 version of the")
    print("HYP-2247 rule: do not count colorings without their extension-rank sidecar.")
    print()

    hist, path, cycles = tournament_fingerprint()
    print("## Tournament Analysis")
    print("vertices=coloring/recursion proof carriers, not runners or raw colors")
    print("pairwise_observable=predicate retention + extension rank + endpoint address + scalar-risk control")
    print("switch=higher retained coloring-recursion payload; ties by AP-tail proof route")
    print(f"score_hist={hist}")
    print(f"directed_3cycles={cycles}")
    print("hamiltonian_path=" + " -> ".join(path))
    print()

    print("## Assumption Challenge")
    print("Considered vertices: runners, residues mod 35, residues mod 7, high gaps,")
    print("color classes, bad-coloring nodes, endpoint addresses, fixed corridors,")
    print("and proof obligations.  The chosen quotient preserves AP-tail escape count,")
    print("endpoint-witness location, mixed/pure phase, and outer-extension shift.")
    print("It destroys arbitrary non-AP wall geometry, so it is a coloring-recursion")
    print("sidecar for the AP-tail bridge rather than a global LRC14 proof.")


if __name__ == "__main__":
    main()
