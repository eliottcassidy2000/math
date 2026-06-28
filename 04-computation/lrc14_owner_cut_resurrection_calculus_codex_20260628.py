#!/usr/bin/env python3
"""HYP-3414 owner-cut resurrection calculus for LRC14 mixed fibers.

This script sits directly on HYP-3410's exact HYP-3406 mixed-fiber rows.
It translates endpoint-owner separation into a finite clause/hitting-set
calculus:

* rows with different theorem exits generate separation clauses;
* an owner cut is a transversal of those clauses;
* the induced cut code must be theorem-exit pure.

The goal is not to prove LRC14 here.  The goal is to make the strongest current
proof-facing obligation precise enough that future scans can either certify it
or name the first residual.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
import importlib.util
from pathlib import Path
import sys


def load_hyp3410():
    path = Path(__file__).with_name(
        "lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628.py"
    )
    spec = importlib.util.spec_from_file_location("hyp3410_charal", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import HYP-3410 script from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


H3410 = load_hyp3410()
FIBERS = H3410.FIBERS
Row = H3410.Row


def owner_labels(rows: tuple[Row, ...]) -> tuple[str, ...]:
    return tuple(sorted({label for row in rows for label in row.owner_support}))


def exit_hist(rows: tuple[Row, ...]) -> dict[str, int]:
    return dict(sorted(Counter(row.exit for row in rows).items()))


def code(row: Row, cut: tuple[str, ...]) -> str:
    return "".join("1" if label in row.owner_support else "0" for label in cut)


def code_partition(rows: tuple[Row, ...], cut: tuple[str, ...]) -> dict[str, list[Row]]:
    part: dict[str, list[Row]] = defaultdict(list)
    for row in rows:
        part[code(row, cut)].append(row)
    return dict(sorted(part.items()))


def code_exit_pure(rows: tuple[Row, ...], cut: tuple[str, ...]) -> bool:
    return all(
        len({row.exit for row in bucket}) == 1
        for bucket in code_partition(rows, cut).values()
    )


def cross_exit_clauses(rows: tuple[Row, ...]) -> list[tuple[tuple[str, str], tuple[str, ...]]]:
    clauses: list[tuple[tuple[str, str], tuple[str, ...]]] = []
    for left, right in combinations(rows, 2):
        if left.exit == right.exit:
            continue
        delta = tuple(sorted(set(left.owner_support) ^ set(right.owner_support)))
        clauses.append(((left.name, right.name), delta))
    return clauses


def min_hitting_sets(clauses: list[tuple[tuple[str, str], tuple[str, ...]]]) -> list[tuple[str, ...]]:
    labels = sorted({label for _, clause in clauses for label in clause})
    if not clauses:
        return [()]
    for size in range(1, len(labels) + 1):
        winners = []
        for subset in combinations(labels, size):
            test = set(subset)
            if all(test.intersection(clause) for _, clause in clauses):
                winners.append(subset)
        if winners:
            return winners
    return []


def singleton_conflicts(rows: tuple[Row, ...]) -> list[tuple[str, str, tuple[str, ...]]]:
    out = []
    for label in owner_labels(rows):
        partition = code_partition(rows, (label,))
        for pattern, bucket in partition.items():
            exits = tuple(sorted({row.exit for row in bucket}))
            if len(exits) > 1:
                out.append((label, pattern, exits))
                break
    return out


def cut_core_and_union(cuts: list[tuple[str, ...]]) -> tuple[tuple[str, ...], tuple[str, ...], Counter]:
    if not cuts:
        return (), (), Counter()
    core = set(cuts[0])
    union = set()
    freq: Counter[str] = Counter()
    for cut in cuts:
        labels = set(cut)
        core &= labels
        union |= labels
        freq.update(cut)
    return tuple(sorted(core)), tuple(sorted(union)), freq


def bdh_top(rows: tuple[Row, ...], limit: int = 5) -> list[tuple[str, Fraction]]:
    return [(label, score) for score, label in H3410.channel_variance(rows)[:limit]]


def print_owner_cut_certificates() -> None:
    print("## Owner-cut clause certificates")
    for fiber_name, rows in FIBERS.items():
        clauses = cross_exit_clauses(rows)
        hitting_sets = min_hitting_sets(clauses)
        h3410_size, h3410_cuts = H3410.min_owner_cut(rows)
        core, union, freq = cut_core_and_union(hitting_sets)
        first_cut = hitting_sets[0]

        print(f"\n### {fiber_name}")
        print(f"rows={len(rows)} exit_hist={exit_hist(rows)}")
        print(f"cross_exit_clause_count={len(clauses)}")
        print(f"minimum_transversal_size={len(first_cut)}")
        print(f"HYP3410_min_owner_cut_size={h3410_size}")
        print(f"minimum_cut_count={len(hitting_sets)}")
        print(f"minimum_cut_core={core}")
        print(f"minimum_cut_union={union}")
        print(f"minimum_cut_participation={dict(freq.most_common())}")
        print(f"first_minimum_cut={first_cut}")
        print(f"first_cut_exit_pure={code_exit_pure(rows, first_cut)}")
        print("first_cut_code_table=")
        for pattern, bucket in code_partition(rows, first_cut).items():
            exits = sorted({row.exit for row in bucket})
            names = [row.name for row in bucket]
            print(f"  {pattern}: exits={exits} rows={names}")
        print("top_BDH_style_owner_variance=")
        for label, score in bdh_top(rows):
            print(f"  {label}: {score}")

        conflicts = singleton_conflicts(rows)
        if conflicts:
            print("singleton_cut_failures=")
            for label, pattern, exits in conflicts[:8]:
                print(f"  {label} has mixed code {pattern} with exits={exits}")
        else:
            print("singleton_cut_failures=none")

        if sorted(hitting_sets) != sorted(h3410_cuts):
            print("WARNING: clause transversals differ from HYP-3410 min cuts")


def print_resurrection_api() -> None:
    print("\n## Mixed-fiber resurrection API")
    stages = [
        (
            "Q0 coarse packet",
            "preserves theorem-row substrate and known exits",
            "destroys residue/height/owner refinements; mixed fibers are expected",
            "if mixed, emit first destroyed coordinate as sidecar debt",
        ),
        (
            "Q1 residue/C3/Qsqrt(-7) skeleton",
            "preserves the 7-adic binding layer and unit-pair symmetry",
            "destroys 2-adic magnitude and endpoint-owner support",
            "repairs curated-bank residue leaks but not expanded-bank leaks",
        ),
        (
            "Q2 residue plus v2/height",
            "preserves the 2-adic hinge and tropical height wall",
            "destroys owner currents and SC accessory parameters",
            "kills the first 12-family leak but not height-persistent owner leaks",
        ),
        (
            "Q3 residue plus owner_support",
            "preserves endpoint-owner cut coordinates",
            "destroys raw row identity and scalar analogies",
            "exact on HYP-3406 through (72,20); next task is first-failure search",
        ),
        (
            "Q4 owner-cut dual certificate",
            "preserves only the owner labels needed by a cut code",
            "destroys irrelevant owner labels and row names",
            "legal when every cut-code bucket has one theorem exit",
        ),
    ]
    for name, preserves, destroys, repair in stages:
        print(f"- {name}:")
        print(f"    preserves={preserves}")
        print(f"    destroys={destroys}")
        print(f"    repair_rule={repair}")


def print_graph_geometry_synthesis() -> None:
    print("\n## Geometry/topology/graph synthesis")
    print(
        "incidence_complex=rows are 0-cells, cross-exit pairs are obstruction "
        "1-cells, owner labels are cut hyperplanes, and a legal quotient is a "
        "collapse whose fibers do not identify two theorem exits."
    )
    print(
        "menger_to_farkas=each cross-exit pair contributes a clause requiring at "
        "least one distinguishing owner label; a minimum owner cut is the finite "
        "dual current hitting all obstruction clauses."
    )
    print(
        "schwarz_christoffel=turn words give the polygonal skeleton, but the "
        "accessory parameter is exactly the owner label set needed by the cut."
    )
    print(
        "bdh=variance is a search heuristic for energetic owner channels; it is "
        "not a certificate unless the induced cut code is exit-pure."
    )
    print(
        "bring=the five theorem exits are a branch alphabet; proof progress is "
        "making the branch label single-valued on the chosen quotient."
    )
    print(
        "recursive_charal=+14 ladders should be tested by whether child decks "
        "preserve the clause transversal number and the terminal exit router."
    )


WEIGHTS = {
    "finite_exactness": 14,
    "exit_purity": 13,
    "dual_certificate": 12,
    "owner_cut": 12,
    "mixed_fiber_api": 11,
    "globalization": 10,
    "terminal_router": 10,
    "recursive_charal": 9,
    "collar_to_bank": 9,
    "height_wall": 8,
    "bdh_prefilter": 5,
    "sc_accessory": 5,
    "bring_branch": 4,
    "raw_analogy": 1,
}


@dataclass(frozen=True)
class Obligation:
    code: str
    name: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    evidence: int
    priority: int

    def score(self) -> int:
        keep = sum(WEIGHTS[x] for x in self.preserves)
        loss = sum(max(1, WEIGHTS[x] // 3) for x in self.destroys)
        return keep + self.evidence + self.priority - loss


def obligations() -> list[Obligation]:
    return [
        Obligation(
            "O00",
            "owner_cut_dual_certificate",
            frozenset({"finite_exactness", "exit_purity", "dual_certificate", "owner_cut"}),
            frozenset({"raw_analogy"}),
            evidence=18,
            priority=18,
        ),
        Obligation(
            "O01",
            "mixed_fiber_resurrection_api",
            frozenset({"finite_exactness", "exit_purity", "mixed_fiber_api", "terminal_router"}),
            frozenset({"raw_analogy"}),
            evidence=15,
            priority=17,
        ),
        Obligation(
            "O02",
            "bounded_owner_cut_theorem",
            frozenset({"owner_cut", "dual_certificate", "globalization", "terminal_router"}),
            frozenset({"height_wall"}),
            evidence=12,
            priority=16,
        ),
        Obligation(
            "O03",
            "owner_support_exactness_extension",
            frozenset({"finite_exactness", "owner_cut", "mixed_fiber_api", "globalization"}),
            frozenset({"recursive_charal"}),
            evidence=10,
            priority=15,
        ),
        Obligation(
            "O04",
            "terminal_chamber_router",
            frozenset({"terminal_router", "exit_purity", "globalization"}),
            frozenset({"raw_analogy"}),
            evidence=8,
            priority=14,
        ),
        Obligation(
            "O05",
            "collar_to_bank_owner_lift",
            frozenset({"collar_to_bank", "height_wall", "owner_cut", "finite_exactness"}),
            frozenset({"globalization"}),
            evidence=7,
            priority=13,
        ),
        Obligation(
            "O06",
            "charal_child_deck_no_new_kernel",
            frozenset({"recursive_charal", "owner_cut", "terminal_router"}),
            frozenset({"finite_exactness"}),
            evidence=5,
            priority=12,
        ),
        Obligation(
            "O07",
            "bdh_exception_ledger",
            frozenset({"bdh_prefilter", "globalization", "terminal_router"}),
            frozenset({"dual_certificate", "finite_exactness"}),
            evidence=3,
            priority=11,
        ),
        Obligation(
            "O08",
            "schwarz_christoffel_accessory_lemma",
            frozenset({"sc_accessory", "owner_cut", "recursive_charal"}),
            frozenset({"exit_purity"}),
            evidence=2,
            priority=10,
        ),
        Obligation(
            "O09",
            "bring_branch_normal_form",
            frozenset({"bring_branch", "terminal_router", "exit_purity"}),
            frozenset({"owner_cut", "dual_certificate"}),
            evidence=1,
            priority=9,
        ),
        Obligation(
            "O10",
            "raw_analogy_rejection",
            frozenset({"raw_analogy"}),
            frozenset(
                {
                    "finite_exactness",
                    "exit_purity",
                    "dual_certificate",
                    "owner_cut",
                    "mixed_fiber_api",
                    "terminal_router",
                    "globalization",
                }
            ),
            evidence=0,
            priority=1,
        ),
    ]


def orient(a: Obligation, b: Obligation) -> tuple[str, str]:
    ka = (a.score(), -len(a.destroys), a.priority, a.code)
    kb = (b.score(), -len(b.destroys), b.priority, b.code)
    return (a.code, b.code) if ka >= kb else (b.code, a.code)


def adjacency(items: list[Obligation]) -> dict[str, set[str]]:
    adj = {item.code: set() for item in items}
    for a, b in combinations(items, 2):
        u, v = orient(a, b)
        adj[u].add(v)
    return adj


def directed_3cycles(adj: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    cycles = []
    for a, b, c in combinations(sorted(adj), 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            cycles.append((a, b, c))
        if c in adj[a] and b in adj[c] and a in adj[b]:
            cycles.append((a, c, b))
    return cycles


def hamiltonian_path_count(adj: dict[str, set[str]]) -> int:
    names = tuple(sorted(adj))
    index = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            last_name = names[last]
            for nxt_name in adj[last_name]:
                nxt = index[nxt_name]
                if mask & (1 << nxt):
                    continue
                dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def priority_path(items: list[Obligation]) -> list[Obligation]:
    return sorted(items, key=lambda item: (item.score(), item.priority, item.code), reverse=True)


def print_tournament() -> None:
    items = obligations()
    adj = adjacency(items)
    print("\n## Tournament Analysis")
    print("vertices=proof-facing obligations, not runners, arcs, or owner labels")
    print("alternate_vertices_considered=runners,gaps,sections,boundary events,residues,cover arcs,Fourier modes,matroid circuits,owner labels,proof obligations")
    print("chosen_vertices=proof obligations generated by mixed-fiber resurrection")
    print("pairwise_observable=weighted retained proof payload minus destroyed sidecars plus exact evidence")
    print("switch_gauge=higher score; ties by fewer destroyed coordinates and priority")
    print("preserved_predicate=theorem-exit purity on legal quotient fibers")
    print("destroyed_information=row order, raw runner identity, irrelevant owner labels, scalar analogies")
    print(f"vertex_count={len(items)}")
    print(f"score_hist={dict(sorted(Counter(item.score() for item in items).items()))}")
    print(f"directed_3cycles={len(directed_3cycles(adj))}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print("priority_hamiltonian_path=")
    for item in priority_path(items):
        print(
            f"  {item.code} {item.name}: score={item.score()} "
            f"keeps={sorted(item.preserves)} destroys={sorted(item.destroys)}"
        )


def print_theorem_targets() -> None:
    print("\n## Strongest proof-facing theorem targets")
    targets = [
        (
            "owner-cut resurrection lemma",
            "For every legal mixed theorem-exit fiber after the current q/AP/GW "
            "pre-routing, the cross-exit clause hypergraph has a bounded owner "
            "transversal whose cut code is theorem-exit pure, or the fiber emits "
            "a named residual debt.",
        ),
        (
            "owner-support exactness extension",
            "Push HYP-3406 beyond (72,20).  If residue+owner_support first fails, "
            "save the failing mixed fiber and run this clause calculus on it.",
        ),
        (
            "terminal chamber router",
            "Every cut-pure fiber must route to AP/GW boundary, strict/positive "
            "open mass, q-witness, state-lift/H7, exact-period exception, or a "
            "new finite residual; owner cuts alone are not terminal exits.",
        ),
        (
            "charal child-deck stability",
            "Endpoint deletion, mirror swap, and +14 child decks should preserve "
            "transversal number and cut-pure terminal codes unless they cross a "
            "recorded owner/accessory wall.",
        ),
    ]
    for name, text in targets:
        print(f"- {name}: {text}")


def main() -> None:
    print("HYP-3414 LRC14 owner-cut resurrection calculus")
    print("=" * 78)
    print("substrate=HYP-3410 exact HYP-3406 mixed fibers")
    print("calculus=cross-exit clauses -> minimum owner transversals -> exit-pure cut codes")
    print_owner_cut_certificates()
    print_resurrection_api()
    print_graph_geometry_synthesis()
    print_theorem_targets()
    print_tournament()
    print("\n## Assumption challenge")
    print(
        "The challenged assumption is that the next LRC14 recursion should be "
        "over runners, residues, or special-function names.  The exact data says "
        "to recurse over legal quotients and the finite owner clauses created "
        "when those quotients identify distinct theorem exits."
    )


if __name__ == "__main__":
    main()
