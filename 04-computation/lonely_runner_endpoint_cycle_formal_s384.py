#!/usr/bin/env python3
"""
lonely_runner_endpoint_cycle_formal_s384.py

codex-2026-05-31 S384

Formal/computational probe for the endpoint-cycle version of the Lonely
Runner endpoint program.

The point of this pass is deliberately small and structural:

* a nonempty protection core forces a directed endpoint-protection cycle;
* therefore any LRC counterexample must realize such a cycle with the integer
  labels from THM-357;
* but bare circular-arc topology has all-protected cores very easily, so the
  next proof/search object has to be an arithmetic labelled cycle, not just an
  abstract endpoint graph.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import ceil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S362 = SourceFileLoader(
    "lonely_runner_bohr_descent_s362",
    str(ROOT / "04-computation" / "lonely_runner_bohr_descent_s362.py"),
).load_module()


@dataclass(frozen=True, order=True)
class Arc:
    q: int
    start: int
    end: int

    @property
    def length(self) -> int:
        return (self.end - self.start) % self.q

    @property
    def boundary(self) -> frozenset[int]:
        return frozenset((self.start, self.end))

    @property
    def cells(self) -> frozenset[int]:
        return frozenset((self.start + j) % self.q for j in range(self.length))

    @property
    def interior_points(self) -> frozenset[int]:
        return frozenset(
            (self.start + j) % self.q for j in range(1, self.length)
        )

    def label(self) -> str:
        return f"({self.start}->{self.end};L{self.length})"


@dataclass(frozen=True)
class AbstractCover:
    q: int
    max_len: int | None
    arcs: tuple[Arc, ...]
    total_length: int
    endpoint_count: int
    core_endpoint_count: int
    core_arc_count: int
    cycle: tuple[int, ...]


@dataclass(frozen=True)
class LrcCycleRow:
    label: str
    n: int
    classification: str
    gap_ratio: Fraction
    endpoint_count: int
    unprotected: int
    peel_depth: int
    core_endpoint_count: int
    core_interval_count: int
    terminal_cycle: tuple[str, ...]
    first_removed_layer: str


def fmt_frac(x: Fraction | None) -> str:
    return S356.fmt_frac(x)


def fmt_float(x: Fraction) -> str:
    return f"{float(x):.6f}"


def all_arcs(q: int, max_len: int | None = None) -> tuple[Arc, ...]:
    arcs: list[Arc] = []
    for start in range(q):
        for length in range(1, q):
            if max_len is not None and length > max_len:
                continue
            arcs.append(Arc(q, start, (start + length) % q))
    return tuple(arcs)


def endpoints_of(arcs: tuple[Arc, ...]) -> set[int]:
    endpoints: set[int] = set()
    for arc in arcs:
        endpoints.update(arc.boundary)
    return endpoints


def full_cell_cover(q: int, arcs: tuple[Arc, ...]) -> bool:
    covered: set[int] = set()
    for arc in arcs:
        covered.update(arc.cells)
    return covered == set(range(q))


def protectors(arcs: tuple[Arc, ...], endpoints: set[int]) -> dict[int, set[Arc]]:
    return {
        endpoint: {arc for arc in arcs if endpoint in arc.interior_points}
        for endpoint in endpoints
    }


def peel_core(
    arcs: tuple[Arc, ...],
) -> tuple[set[int], set[Arc]]:
    remaining_endpoints = endpoints_of(arcs)
    remaining_arcs = set(arcs)

    while True:
        current_protectors = protectors(tuple(remaining_arcs), remaining_endpoints)
        dead = {
            endpoint
            for endpoint in remaining_endpoints
            if not current_protectors[endpoint]
        }
        if not dead:
            return remaining_endpoints, remaining_arcs

        remaining_endpoints -= dead
        remaining_arcs = {
            arc for arc in remaining_arcs if not (arc.boundary & dead)
        }


def endpoint_graph_from_arcs(
    arcs: tuple[Arc, ...],
    endpoints: set[int] | None = None,
) -> dict[int, set[int]]:
    if endpoints is None:
        endpoints = endpoints_of(arcs)
    graph: dict[int, set[int]] = {endpoint: set() for endpoint in endpoints}
    for endpoint in endpoints:
        for arc in arcs:
            if endpoint in arc.interior_points:
                graph[endpoint].update(arc.boundary & endpoints)
    return graph


def find_directed_cycle(graph: dict[int, set[int]]) -> tuple[int, ...]:
    seen: set[int] = set()
    stack: list[int] = []
    on_stack: set[int] = set()

    def visit(node: int) -> tuple[int, ...] | None:
        seen.add(node)
        stack.append(node)
        on_stack.add(node)
        for nxt in sorted(graph.get(node, ())):
            if nxt not in seen:
                found = visit(nxt)
                if found:
                    return found
            elif nxt in on_stack:
                i = stack.index(nxt)
                return tuple(stack[i:] + [nxt])
        stack.pop()
        on_stack.remove(node)
        return None

    for node in sorted(graph):
        if node not in seen:
            found = visit(node)
            if found:
                return found
    return tuple()


def is_all_protected_cover(q: int, arcs: tuple[Arc, ...]) -> bool:
    endpoints = endpoints_of(arcs)
    if not endpoints or not full_cell_cover(q, arcs):
        return False
    return all(protectors(arcs, endpoints)[endpoint] for endpoint in endpoints)


def find_minimal_abstract_cover(q: int, max_len: int | None) -> AbstractCover | None:
    arcs = all_arcs(q, max_len)
    for size in range(1, min(6, len(arcs)) + 1):
        best: tuple[Arc, ...] | None = None
        best_key: tuple[int, tuple[tuple[int, int], ...]] | None = None
        for combo in combinations(arcs, size):
            if not is_all_protected_cover(q, combo):
                continue
            key = (
                sum(arc.length for arc in combo),
                tuple((arc.start, arc.end) for arc in combo),
            )
            if best is None or key < best_key:
                best = combo
                best_key = key
        if best is not None:
            core_endpoints, core_arcs = peel_core(best)
            graph = endpoint_graph_from_arcs(tuple(core_arcs), core_endpoints)
            return AbstractCover(
                q=q,
                max_len=max_len,
                arcs=tuple(best),
                total_length=sum(arc.length for arc in best),
                endpoint_count=len(endpoints_of(best)),
                core_endpoint_count=len(core_endpoints),
                core_arc_count=len(core_arcs),
                cycle=find_directed_cycle(graph),
            )
    return None


def interval_name(interval) -> str:
    return f"{interval.speed}:{interval.center}"


def endpoint_name(endpoint: Fraction) -> str:
    return fmt_frac(endpoint)


def lrc_endpoint_graph(
    endpoints: set[Fraction],
    intervals: set[object],
    protectors_by_endpoint: dict[Fraction, set[object]],
    boundary_by_interval: dict[object, set[Fraction]],
) -> dict[Fraction, set[Fraction]]:
    graph: dict[Fraction, set[Fraction]] = {endpoint: set() for endpoint in endpoints}
    for endpoint in endpoints:
        for interval in protectors_by_endpoint[endpoint] & intervals:
            graph[endpoint].update(boundary_by_interval[interval] & endpoints)
    return graph


def find_fraction_cycle(graph: dict[Fraction, set[Fraction]]) -> tuple[Fraction, ...]:
    int_graph: dict[int, set[int]] = {}
    values = sorted(graph)
    index = {value: i for i, value in enumerate(values)}
    for value, targets in graph.items():
        int_graph[index[value]] = {index[target] for target in targets}
    cycle = find_directed_cycle(int_graph)
    if not cycle:
        return tuple()
    return tuple(values[i] for i in cycle)


def lrc_cycle_row(label: str, speeds: tuple[int, ...]) -> LrcCycleRow:
    report = S356.report(label, list(speeds))
    descent = S362.summarize(list(report.speeds))
    endpoints, intervals, _owners, protectors_by_endpoint, boundary = (
        S362.build_endpoint_system(report.speeds)
    )
    layers, core_endpoints, core_intervals = S362.peel_protection_core(
        descent.q, endpoints, intervals, protectors_by_endpoint, boundary
    )
    graph = lrc_endpoint_graph(
        core_endpoints, core_intervals, protectors_by_endpoint, boundary
    )
    cycle = find_fraction_cycle(graph)
    first_layer = "none"
    if layers:
        first = layers[0]
        first_layer = (
            f"removed={first.removed_endpoints},"
            f"subgroup={first.removed_subgroup_modulus}"
        )
    return LrcCycleRow(
        label=label,
        n=len(report.speeds) + 1,
        classification=descent.classification,
        gap_ratio=report.max_gap / report.threshold if report.threshold else Fraction(0),
        endpoint_count=descent.endpoint_count,
        unprotected=descent.unprotected_count,
        peel_depth=len(layers),
        core_endpoint_count=len(core_endpoints),
        core_interval_count=len(core_intervals),
        terminal_cycle=tuple(endpoint_name(value) for value in cycle),
        first_removed_layer=first_layer,
    )


def sample_speed_sets() -> tuple[tuple[str, tuple[int, ...]], ...]:
    return (
        ("initial n=8", tuple(range(1, 8))),
        ("sporadic tight n=8A", (1, 4, 5, 6, 7, 11, 13)),
        ("initial n=14", tuple(range(1, 14))),
        (
            "n14 seven-ladder",
            (1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91),
        ),
        ("n14 single-gate", tuple(list(range(1, 13)) + [14])),
        ("initial n=15", tuple(range(1, 15))),
        (
            "n15 3x5 ladder",
            (3, 5, 6, 9, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50),
        ),
    )


def print_formal_statement() -> None:
    print("FORMAL ENDPOINT-CYCLE REDUCTION")
    print("=" * 88)
    print(
        "Given a finite endpoint/interval system (E,I,B,P), draw an arrow "
        "e -> f when some interval protecting e has f as a boundary endpoint."
    )
    print(
        "Every nonempty protection core has at least one outgoing arrow from "
        "each remaining endpoint, hence contains a directed endpoint cycle."
    )
    print(
        "By THM-357 and THM-359, an LRC counterexample therefore must realize "
        "a directed cycle in the labelled integer endpoint-protection graph."
    )
    print()


def print_abstract_mirages() -> None:
    print("ABSTRACT CIRCULAR-ARC MIRAGES")
    print("=" * 88)
    print(
        "Bare topology is too weak: all-protected circular-arc covers appear "
        "immediately once arithmetic labels are forgotten."
    )
    print()
    print("q max_len min_arcs total_len endpoints coreE coreI cycle arcs")
    print("-" * 88)
    for q in range(3, 10):
        for max_len in (None, ceil(q / 2)):
            cover = find_minimal_abstract_cover(q, max_len)
            if cover is None:
                print(f"{q:>1} {str(max_len):>7} none")
                continue
            cycle = "->".join(str(x) for x in cover.cycle) if cover.cycle else "none"
            arcs = " ".join(arc.label() for arc in cover.arcs)
            max_label = "any" if max_len is None else str(max_len)
            print(
                f"{q:>1} {max_label:>7} {len(cover.arcs):>8} "
                f"{cover.total_length:>9} {cover.endpoint_count:>9} "
                f"{cover.core_endpoint_count:>5} {cover.core_arc_count:>5} "
                f"{cycle:<11} {arcs}"
            )
    print()


def print_lrc_cycle_audit() -> None:
    print("LRC TERMINAL-CORE CYCLE AUDIT")
    print("=" * 88)
    print(
        "The sampled LRC systems all peel to empty terminal core.  This is "
        "stronger than having no visible gap: even tiny-gap ladders fail to "
        "hold an endpoint cycle after protection peeling."
    )
    print()
    print(
        "label                  n class          gap/th endpoints unprot "
        "peel coreE coreI cycle first_layer"
    )
    print("-" * 112)
    for label, speeds in sample_speed_sets():
        row = lrc_cycle_row(label, speeds)
        cycle = "yes" if row.terminal_cycle else "no"
        print(
            f"{row.label:<22} {row.n:>2} {row.classification:<14} "
            f"{fmt_float(row.gap_ratio):>7} {row.endpoint_count:>9} "
            f"{row.unprotected:>7} {row.peel_depth:>4} "
            f"{row.core_endpoint_count:>5} {row.core_interval_count:>5} "
            f"{cycle:>5} {row.first_removed_layer}"
        )
    print()


def print_cycle_label_grammar() -> None:
    print("ARITHMETIC LABEL GRAMMAR FOR THE NEXT SEARCH")
    print("=" * 88)
    print(
        "A protected endpoint cycle is not just e0 -> e1 -> ... -> e0.  Each "
        "arrow must carry an owner speed u, a protector speed p, a center m, "
        "and a sign eps satisfying the strict integer inequality"
    )
    print()
    print("    | p*(n*m + eps) - a*n*u | < u")
    print()
    print(
        "for some integer a.  The cycle-first search should enumerate these "
        "labelled inequalities, then ask whether their boundary intervals can "
        "also produce full measure.  That separates real LRC arithmetic from "
        "abstract circular-arc mirages."
    )
    print()


def print_new_ideas() -> None:
    print("NEW IDEAS GENERATED")
    print("=" * 88)
    ideas = [
        "Endpoint-cycle rank: the length and slack profile of the first surviving protection cycle after partial peeling.",
        "Arithmetic mirage filter: generate abstract all-protected cores first, then reject them by speed-orbit and integer-inequality labels.",
        "Cycle slack potential: sum the strict margins u-|p*(n*m+eps)-a*n*u| around a cycle; prove primitive cycles have positive leak.",
        "Unit-layer trap: any cycle touching unit endpoints must pass through speeds divisible by n, so quotient descent should start there.",
        "Torsion shortcut test: product-sum/torsion coordinates may be exactly the short labelled cycles that almost close but expose higher endpoints.",
        "LRC-TDA feature: terminal-core size is too coarse; record the largest preterminal directed cycle before its first endpoint leaf is removed.",
    ]
    for i, idea in enumerate(ideas, 1):
        print(f"{i}. {idea}")
    print()


def main() -> None:
    print("Lonely Runner endpoint-cycle formal session (codex-2026-05-31 S384)")
    print()
    print_formal_statement()
    print_abstract_mirages()
    print_lrc_cycle_audit()
    print_cycle_label_grammar()
    print_new_ideas()


if __name__ == "__main__":
    main()
