#!/usr/bin/env python3
"""HYP-3419: charal owner-cut recursion prototype for LRC14.

This is a second pass after HYP-3410.  It keeps the same exact mixed fibers but
pushes the Menger/charal signal one step closer to a finite lemma API:

    charal quotient -> mixed theorem-exit fiber -> owner cut sidecar
    -> binary cut recursion -> terminal theorem exit or named debt.

Famous outside objects remain guardrails.  Bring supplies the branch alphabet;
Schwarz-Christoffel supplies turn/accessory language; BDH supplies finite
channel variance; Menger supplies the cut theorem shape; Krasner,
Sophie-Germain, HLW, Soldner, and Meissel-Mertens are admitted only after they
become exact packet fields or terminal normalizers.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import log2

import lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628 as h3410


Row = h3410.Row


def exit_set(rows: tuple[Row, ...]) -> frozenset[str]:
    return frozenset(row.exit for row in rows)


def is_pure(rows: tuple[Row, ...]) -> bool:
    return len(exit_set(rows)) <= 1


def exit_name(rows: tuple[Row, ...]) -> str:
    exits = sorted(exit_set(rows))
    return exits[0] if len(exits) == 1 else "MIXED:" + ",".join(exits)


def all_labels(rows: tuple[Row, ...]) -> tuple[str, ...]:
    return tuple(sorted({label for row in rows for label in row.owner_support}))


def key_without_owner_labels(row: Row) -> tuple[object, ...]:
    return tuple(item for item in h3410.charal_signature(row) if item[0] != "owner_labels")


def mixed_fiber_count(rows: tuple[Row, ...], key_fn) -> tuple[int, int, int]:
    fibers: dict[object, list[Row]] = defaultdict(list)
    for row in rows:
        fibers[key_fn(row)].append(row)
    mixed = sum(1 for fiber_rows in fibers.values() if not is_pure(tuple(fiber_rows)))
    largest = max((len(fiber_rows) for fiber_rows in fibers.values()), default=0)
    return mixed, len(fibers), largest


def all_min_cut_sets(rows: tuple[Row, ...]) -> tuple[int, list[tuple[str, ...]]]:
    labels = all_labels(rows)
    for size in range(1, len(labels) + 1):
        winners: list[tuple[str, ...]] = []
        for subset in combinations(labels, size):
            mixed, _, _ = mixed_fiber_count(
                rows,
                lambda row, s=subset: tuple(label for label in s if label in row.owner_support),
            )
            if mixed == 0:
                winners.append(subset)
        if winners:
            return size, winners
    return 0, []


def cut_core(cuts: list[tuple[str, ...]]) -> tuple[str, ...]:
    if not cuts:
        return ()
    core = set(cuts[0])
    for cut in cuts[1:]:
        core &= set(cut)
    return tuple(sorted(core))


def cut_frequency(cuts: list[tuple[str, ...]]) -> list[tuple[int, str]]:
    counts = Counter(label for cut in cuts for label in cut)
    return sorted(((count, label) for label, count in counts.items()), reverse=True)


def channel_variance_map(rows: tuple[Row, ...]) -> dict[str, Fraction]:
    return {label: score for score, label in h3410.channel_variance(rows)}


def exit_entropy(rows: tuple[Row, ...]) -> float:
    n = len(rows)
    counts = Counter(row.exit for row in rows)
    return -sum((count / n) * log2(count / n) for count in counts.values())


def impurity_cost(groups: tuple[tuple[Row, ...], ...]) -> tuple[int, int, float]:
    mixed_rows = sum(len(group) for group in groups if not is_pure(group))
    mixed_nodes = sum(1 for group in groups if not is_pure(group))
    entropy = sum(len(group) * exit_entropy(group) for group in groups)
    return mixed_rows, mixed_nodes, entropy


@dataclass(frozen=True)
class Tree:
    label: str | None
    rows: tuple[Row, ...]
    present: "Tree | None" = None
    absent: "Tree | None" = None
    residual: bool = False

    @property
    def depth(self) -> int:
        if self.label is None:
            return 0
        return 1 + max(self.present.depth if self.present else 0, self.absent.depth if self.absent else 0)

    @property
    def leaf_count(self) -> int:
        if self.label is None:
            return 1
        return (self.present.leaf_count if self.present else 0) + (self.absent.leaf_count if self.absent else 0)

    @property
    def residual_count(self) -> int:
        if self.label is None:
            return int(self.residual)
        return (self.present.residual_count if self.present else 0) + (
            self.absent.residual_count if self.absent else 0
        )


def tree_lines(tree: Tree, indent: str = "") -> list[str]:
    names = ", ".join(row.name for row in tree.rows)
    if tree.label is None:
        marker = "RESIDUAL" if tree.residual else "EXIT"
        return [f"{indent}{marker} {exit_name(tree.rows)} rows={len(tree.rows)} names=[{names}]"]
    lines = [f"{indent}test owner_label {tree.label} rows={len(tree.rows)} exits={sorted(exit_set(tree.rows))}"]
    lines.extend(tree_lines(tree.present, indent + "  yes -> ") if tree.present else [])
    lines.extend(tree_lines(tree.absent, indent + "  no  -> ") if tree.absent else [])
    return lines


def optimal_tree(rows: tuple[Row, ...], labels: tuple[str, ...]) -> Tree:
    row_names = tuple(row.name for row in rows)
    row_lookup = {row.name: row for row in rows}
    variance = channel_variance_map(rows)

    @lru_cache(maxsize=None)
    def build(names: tuple[str, ...], available: tuple[str, ...]) -> tuple[tuple[int, int, int, Fraction, str], Tree]:
        node_rows = tuple(row_lookup[name] for name in names)
        if is_pure(node_rows):
            tree = Tree(label=None, rows=node_rows)
            return (0, 1, 0, Fraction(0), ""), tree
        if not available:
            tree = Tree(label=None, rows=node_rows, residual=True)
            return (999, 1, len(node_rows), Fraction(0), ""), tree

        best_key: tuple[int, int, int, Fraction, str] | None = None
        best_tree: Tree | None = None
        for label in available:
            present_names = tuple(name for name in names if label in row_lookup[name].owner_support)
            absent_names = tuple(name for name in names if label not in row_lookup[name].owner_support)
            if not present_names or not absent_names:
                continue
            next_available = tuple(x for x in available if x != label)
            present_key, present_tree = build(present_names, next_available)
            absent_key, absent_tree = build(absent_names, next_available)
            depth = 1 + max(present_tree.depth, absent_tree.depth)
            leaves = present_tree.leaf_count + absent_tree.leaf_count
            residuals = present_tree.residual_count + absent_tree.residual_count
            # Prefer shallow exact trees; use finite BDH variance only as a tie-break.
            key = (depth, leaves, residuals, -variance.get(label, Fraction(0)), label)
            if best_key is None or key < best_key:
                best_key = key
                best_tree = Tree(label=label, rows=node_rows, present=present_tree, absent=absent_tree)
        if best_tree is None:
            tree = Tree(label=None, rows=node_rows, residual=True)
            return (999, 1, len(node_rows), Fraction(0), ""), tree
        return best_key, best_tree

    return build(row_names, labels)[1]


def greedy_tree(rows: tuple[Row, ...], labels: tuple[str, ...]) -> Tree:
    if is_pure(rows):
        return Tree(label=None, rows=rows)
    if not labels:
        return Tree(label=None, rows=rows, residual=True)
    variance = channel_variance_map(rows)
    candidates = []
    for label in labels:
        present = tuple(row for row in rows if label in row.owner_support)
        absent = tuple(row for row in rows if label not in row.owner_support)
        if not present or not absent:
            continue
        cost = impurity_cost((present, absent))
        candidates.append((cost, -variance.get(label, Fraction(0)), label, present, absent))
    if not candidates:
        return Tree(label=None, rows=rows, residual=True)
    _, _, label, present, absent = min(candidates)
    next_labels = tuple(x for x in labels if x != label)
    return Tree(
        label=label,
        rows=rows,
        present=greedy_tree(present, next_labels),
        absent=greedy_tree(absent, next_labels),
    )


def purity_table(rows: tuple[Row, ...], cut: tuple[str, ...]) -> list[tuple[str, int, int, int]]:
    keys = [
        ("raw_route", lambda row: row.route),
        ("turn_word", h3410.sc_turn_word),
        ("owner_group_counts", h3410.bdh_channel_vector),
        ("turn+groups+apex", lambda row: (h3410.sc_turn_word(row), h3410.bdh_channel_vector(row), "7:g7" in row.owner_support)),
        ("charal_without_owner_labels", key_without_owner_labels),
        ("first_min_owner_cut", lambda row, s=cut: tuple(label for label in s if label in row.owner_support)),
        ("full_owner_labels", lambda row: tuple(sorted(row.owner_support))),
        ("full_charal_signature", h3410.charal_signature),
    ]
    return [(name, *mixed_fiber_count(rows, key_fn)) for name, key_fn in keys]


@dataclass(frozen=True)
class Module:
    code: str
    name: str
    keeps: tuple[str, ...]
    destroys: tuple[str, ...]
    proof_use: str
    scores: dict[str, int]

    @property
    def total(self) -> int:
        weights = {
            "finite_exactness": 5,
            "cut_power": 5,
            "recursion": 4,
            "sidecar_hygiene": 3,
            "analytic_guardrail": 2,
            "risk": -4,
        }
        return sum(weights[key] * value for key, value in self.scores.items())


def modules() -> list[Module]:
    return [
        Module(
            "M00",
            "bounded owner-cut theorem",
            ("exit_status", "endpoint_owner", "finite Menger cut", "terminal debt"),
            ("raw scalar analogy",),
            "Prove every mixed charal fiber has a bounded owner cut or a named residual.",
            {
                "finite_exactness": 5,
                "cut_power": 5,
                "recursion": 4,
                "sidecar_hygiene": 4,
                "analytic_guardrail": 1,
                "risk": 0,
            },
        ),
        Module(
            "M01",
            "charal decision-tree API",
            ("charal signature", "binary owner tests", "terminal exits"),
            ("global theorem without recursion",),
            "Turn mixed fibers into explicit yes/no owner-label decision trees.",
            {
                "finite_exactness": 5,
                "cut_power": 4,
                "recursion": 5,
                "sidecar_hygiene": 4,
                "analytic_guardrail": 1,
                "risk": 0,
            },
        ),
        Module(
            "M02",
            "finite BDH label variance",
            ("owner channel discrepancy", "cut-label priority"),
            ("prime-distribution black box",),
            "Use mean-square label imbalance only to prioritize exact cut sidecars.",
            {
                "finite_exactness": 4,
                "cut_power": 4,
                "recursion": 3,
                "sidecar_hygiene": 4,
                "analytic_guardrail": 3,
                "risk": 0,
            },
        ),
        Module(
            "M03",
            "Krasner owner-stability gate",
            ("local contact/root stability", "endpoint_owner", "same-residue lifts"),
            ("raw p-adic closeness",),
            "Admit local stability only when the owner/contact packet is unchanged.",
            {
                "finite_exactness": 3,
                "cut_power": 3,
                "recursion": 4,
                "sidecar_hygiene": 4,
                "analytic_guardrail": 3,
                "risk": 1,
            },
        ),
        Module(
            "M04",
            "Schwarz-Christoffel accessory reconstruction",
            ("turn word", "accessory owner debt", "contact polygon"),
            ("turn-only proof",),
            "Treat owner labels as accessory parameters hidden by a polygon turn word.",
            {
                "finite_exactness": 3,
                "cut_power": 3,
                "recursion": 3,
                "sidecar_hygiene": 5,
                "analytic_guardrail": 2,
                "risk": 1,
            },
        ),
        Module(
            "M05",
            "Sophie-Germain height-factor channel",
            ("quartic height/flex debt", "two quadratic channels"),
            ("factor identity as terminal proof",),
            "Split a live height/flex debt after the cut recursion names it.",
            {
                "finite_exactness": 3,
                "cut_power": 2,
                "recursion": 3,
                "sidecar_hygiene": 4,
                "analytic_guardrail": 2,
                "risk": 1,
            },
        ),
        Module(
            "M06",
            "Bring branch alphabet",
            ("five theorem exits", "branch normalization"),
            ("quintic formula chase", "owner labels"),
            "Use Bring only to normalize the branch problem after sidecars are known.",
            {
                "finite_exactness": 2,
                "cut_power": 1,
                "recursion": 2,
                "sidecar_hygiene": 3,
                "analytic_guardrail": 2,
                "risk": 2,
            },
        ),
        Module(
            "M07",
            "Soldner-HLW-Mertens scalar firewall",
            ("zero-level hygiene", "no-scalar-shadow", "tail entropy labels"),
            ("finite packet certificate", "endpoint owner"),
            "Allow scalar constants only as normalization or tail labels after finite cuts.",
            {
                "finite_exactness": 1,
                "cut_power": 1,
                "recursion": 1,
                "sidecar_hygiene": 5,
                "analytic_guardrail": 5,
                "risk": 2,
            },
        ),
    ]


def ordered_modules() -> list[Module]:
    return sorted(modules(), key=lambda module: (-module.total, module.code))


def directed_3cycles(ms: list[Module]) -> int:
    index = {module.code: idx for idx, module in enumerate(ms)}
    count = 0
    for a, b, c in combinations(ms, 3):
        ab = index[a.code] < index[b.code]
        bc = index[b.code] < index[c.code]
        ca = index[c.code] < index[a.code]
        ba = index[b.code] < index[a.code]
        cb = index[c.code] < index[b.code]
        ac = index[a.code] < index[c.code]
        if (ab and bc and ca) or (ba and cb and ac):
            count += 1
    return count


def print_cut_report() -> None:
    print("## Owner-Cut Recursion")
    for fiber_name, rows in h3410.FIBERS.items():
        cut_size, cuts = all_min_cut_sets(rows)
        core = cut_core(cuts)
        frequencies = cut_frequency(cuts)
        first_cut = cuts[0] if cuts else ()
        greedy = greedy_tree(rows, all_labels(rows))
        optimal = optimal_tree(rows, all_labels(rows))
        print()
        print(f"### {fiber_name}")
        print(f"rows={len(rows)} exit_hist={dict(Counter(row.exit for row in rows))}")
        print(f"min_owner_cut_size={cut_size}")
        print(f"min_owner_cut_count={len(cuts)}")
        print(f"min_owner_cut_core={core}")
        print(f"top_labels_across_min_cuts={frequencies[:8]}")
        print(f"first_min_cut={first_cut}")
        print(f"greedy_tree_depth={greedy.depth} leaves={greedy.leaf_count} residuals={greedy.residual_count}")
        print(f"optimal_tree_depth={optimal.depth} leaves={optimal.leaf_count} residuals={optimal.residual_count}")
        print("optimal_tree=")
        for line in tree_lines(optimal):
            print("  " + line)


def print_purity_report() -> None:
    print("\n## Charal Signature Purity Levels")
    for fiber_name, rows in h3410.FIBERS.items():
        cut_size, cuts = all_min_cut_sets(rows)
        cut = cuts[0] if cuts else ()
        print()
        print(f"### {fiber_name}")
        for name, mixed, fibers, largest in purity_table(rows, cut):
            print(f"{name:32s} mixed={mixed} fibers={fibers} largest_fiber={largest}")


def print_guardrail_report() -> None:
    print("\n## Guardrail-Gated Proof Roles")
    roles = [
        (
            "Bring radical",
            "five-exit branch alphabet",
            "not a quintic formula; branch ambiguity must be killed by exact sidecars",
        ),
        (
            "Schwarz-Christoffel",
            "turn word plus accessory owner debt",
            "turns alone are demoted whenever owner labels separate exits",
        ),
        (
            "Barban-Davenport-Halberstam",
            "finite owner-channel variance",
            "variance ranks labels but never replaces the cut certificate",
        ),
        (
            "Menger cuts",
            "bounded endpoint-owner separator",
            "this is the theorem-shaped core after HYP-3410",
        ),
        (
            "Krasner's lemma",
            "owner/contact stability gate",
            "same-residue or p-adic-near moves are legal only with stable contact packets",
        ),
        (
            "Sophie Germain identity",
            "quartic height/flex split into two quadratic channels",
            "used only after a live height sidecar is named by the recursion",
        ),
        (
            "Hermite-Lindemann-Weierstrass",
            "no-scalar-shadow firewall",
            "transcendental shadows cannot certify finite rational packets",
        ),
        (
            "Ramanujan-Soldner constant",
            "zero-level normalization hygiene",
            "renormalizes branch potentials but does not split a mixed fiber",
        ),
        (
            "Meissel-Mertens constant",
            "tail-entropy label after finite exceptions",
            "global averages start only after owner-cut residuals are named",
        ),
    ]
    for name, admitted, blocked in roles:
        print(f"- {name}: admitted={admitted}; guardrail={blocked}.")


def print_theorem_target() -> None:
    print("\n## Candidate Finite Lemma Shape")
    print("For every expanded-bank mixed charal fiber after residue/height sidecars:")
    print("  1. either a bounded endpoint-owner cut makes theorem_exit a function;")
    print("  2. or the fiber admits a dual owner-current / Farkas certificate;")
    print("  3. or a stable charal ladder preserves positive-Haar-open exit;")
    print("  4. or it terminates at AP/GW, strict-open mass, q-witness, H7/state lift,")
    print("     Schwarz-Christoffel accessory debt, exact-period/BDH exception,")
    print("     height-factor debt, or a newly named finite residual.")
    print()
    print("Known finite evidence from HYP-3410/HYP-3419:")
    print("  height_leak_12_family cut_size=1")
    print("  persistent_owner_leak_26_40_54_family cut_size=1")
    print("  height_persistent_owner_leak_10_20_drop_add_family cut_size=3")
    print("The new warning is that one-label owner theorems are already too optimistic;")
    print("the plausible next theorem is bounded owner-cut recursion.")


def print_tournament() -> None:
    print("\n## Tournament Analysis")
    ordered = ordered_modules()
    print("vertices=proof modules and guardrail gates, not runners/arcs/constants")
    print("pairwise_observable=finite exactness + cut power + recursion - scalar risk")
    print("switch=higher weighted score; ties by module code")
    print(f"vertex_count={len(ordered)}")
    print(f"score_hist={dict(sorted(Counter(module.total for module in ordered).items()))}")
    print(f"directed_3cycles={directed_3cycles(ordered)}")
    print("hamiltonian_path=" + " -> ".join(module.code for module in ordered))
    for module in ordered:
        print(
            f"  {module.code} score={module.total:2d} {module.name}; "
            f"keeps={module.keeps}; destroys={module.destroys}"
        )


def main() -> None:
    print("HYP-3419 CHARAL OWNER-CUT RECURSION PROTOTYPE")
    print("status=SYNTHESIS / finite cut-recursion API prototype; not an LRC14 proof")
    print("source=HYP-3410 exact mixed fibers + HYP-3408 guardrails + HYP-3409 recursion")
    print()
    print_cut_report()
    print_purity_report()
    print_guardrail_report()
    print_theorem_target()
    print_tournament()
    print("\n## Assumption Challenge")
    print("Considered vertices: runners, owner labels, cut sets, charal turns, BDH labels,")
    print("Bring branches, Schwarz-Christoffel prevertices, algebraic constants, and proof modules.")
    print("Chosen vertices: proof modules and guardrail gates.  The preserved LRC predicate is")
    print("theorem-exit purity after controlled forgetting; destroyed raw data is allowed only")
    print("when a later mixed fiber demands it as a named sidecar.")


if __name__ == "__main__":
    main()
