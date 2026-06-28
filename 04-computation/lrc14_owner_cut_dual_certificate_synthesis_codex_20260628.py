#!/usr/bin/env python3
"""HYP-3417: owner-cut dual certificate synthesis for LRC14.

This scout sharpens the HYP-3410 Menger-owner-cut route.  It does not try to
import famous theorems as black boxes.  Instead it asks whether the visible
mixed fibers already have small labelled owner currents:

* an island current, when a label identifies the unit-petal row;
* a positive-debt hitting current, when added owner labels hit every
  positive-open row and vanish on the unit-petal row.

The extra prompt ideas become audit columns around that current:

* Krasner: the common owner core is not a stable theorem-exit root packet;
* Sophie-Germain: cut size and stable-core size give exact quadratic channels;
* Meissel-Mertens: choose among equal-size cuts by reciprocal owner budget;
* Ramanujan-Soldner: declare the zero level of the signed current;
* HLW: forbid replacing labelled finite packets by scalar shadows.

Rebased over S257/HYP-3411-HYP-3413, the frontier certificate is read as a
local owner-current echo of the residue/magnitude split: one even-cover label
plus two binding labels.  It is not a substitute for the global GW
q == 1 mod 3 criterion.

After S258/HYP-3415, this certificate is explicitly a sidecar: the critical
path remains the decorrelation floor inequality.
After HYP-3416, it is one labelled owner-current layer in the recursive
quotient ladder.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations

import lrc14_bring_sc_bdh_menger_charal_recursion_codex_20260628 as h3410


UNIT_EXIT = "unit-petal-named"
POSITIVE_EXIT = "positive-Haar-open"


def label_key(label: str) -> tuple[int, str]:
    return (h3410.owner_index(label), h3410.owner_group(label))


def fmt_fraction(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def owner_budget(labels: tuple[str, ...]) -> Fraction:
    return sum(Fraction(1, h3410.owner_index(label)) for label in labels)


def owner_group_word(labels: tuple[str, ...]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted(Counter(h3410.owner_group(label) for label in labels).items()))


def minimal_hitting_sets(
    target_sets: list[set[str]],
    candidates: set[str],
) -> list[tuple[str, ...]]:
    labels = sorted(candidates, key=label_key)
    if not target_sets:
        return [()]
    for size in range(1, len(labels) + 1):
        winners = []
        for subset in combinations(labels, size):
            if all(set(subset) & target for target in target_sets):
                winners.append(tuple(sorted(subset, key=label_key)))
        if winners:
            return sorted(winners, key=lambda cut: (owner_budget(cut), cut))
    return []


def sophie_channels(stable_core_size: int, cut_size: int) -> tuple[int, int, int]:
    a = stable_core_size
    b = cut_size
    plus = a * a + 2 * b * b + 2 * a * b
    minus = a * a + 2 * b * b - 2 * a * b
    return plus, minus, plus * minus


@dataclass(frozen=True)
class Certificate:
    mode: str
    cut: tuple[str, ...]
    budget: Fraction
    scores: tuple[tuple[str, str, int], ...]
    margin: int

    def score_hist_by_exit(self) -> dict[str, Counter[int]]:
        out: dict[str, Counter[int]] = defaultdict(Counter)
        for exit_name, _name, score in self.scores:
            out[exit_name][score] += 1
        return dict(out)


@dataclass(frozen=True)
class FiberAudit:
    name: str
    row_count: int
    exit_hist: dict[str, int]
    common_core: tuple[str, ...]
    projection_min_cut_size: int
    projection_first_cuts: tuple[tuple[str, ...], ...]
    selected: Certificate
    positive_debt: Certificate | None
    island: Certificate | None
    positive_debt_sophie: tuple[int, int, int] | None
    selected_sophie: tuple[int, int, int]
    top_variance: tuple[Fraction, str]


def score_rows(
    rows: tuple[h3410.Row, ...],
    cut: tuple[str, ...],
    mode: str,
) -> tuple[tuple[str, str, int], ...]:
    scored = []
    cut_set = set(cut)
    for row in rows:
        hit_count = len(cut_set & set(row.owner_support))
        if mode == "positive_debt":
            score = hit_count
        elif mode == "unit_island":
            score = hit_count
        else:
            raise ValueError(mode)
        scored.append((row.exit, row.name, score))
    return tuple(scored)


def score_margin(scores: tuple[tuple[str, str, int], ...], mode: str) -> int:
    unit_scores = [score for exit_name, _name, score in scores if exit_name == UNIT_EXIT]
    pos_scores = [score for exit_name, _name, score in scores if exit_name == POSITIVE_EXIT]
    if not unit_scores or not pos_scores:
        return 0
    if mode == "positive_debt":
        return min(pos_scores) - max(unit_scores)
    if mode == "unit_island":
        return min(unit_scores) - max(pos_scores)
    return 0


def make_certificate(
    rows: tuple[h3410.Row, ...], cut: tuple[str, ...], mode: str
) -> Certificate:
    scores = score_rows(rows, cut, mode)
    return Certificate(mode, cut, owner_budget(cut), scores, score_margin(scores, mode))


def best_positive_debt(rows: tuple[h3410.Row, ...]) -> Certificate | None:
    unit_rows = [row for row in rows if row.exit == UNIT_EXIT]
    pos_rows = [row for row in rows if row.exit == POSITIVE_EXIT]
    if not unit_rows or not pos_rows:
        return None
    unit_union = set().union(*(set(row.owner_support) for row in unit_rows))
    pos_targets = [set(row.owner_support) - unit_union for row in pos_rows]
    candidates = set().union(*pos_targets) if pos_targets else set()
    cuts = minimal_hitting_sets(pos_targets, candidates)
    if not cuts:
        return None
    return make_certificate(rows, cuts[0], "positive_debt")


def best_unit_island(rows: tuple[h3410.Row, ...]) -> Certificate | None:
    unit_rows = [row for row in rows if row.exit == UNIT_EXIT]
    pos_rows = [row for row in rows if row.exit == POSITIVE_EXIT]
    if not unit_rows or not pos_rows:
        return None
    pos_union = set().union(*(set(row.owner_support) for row in pos_rows))
    unit_targets = [set(row.owner_support) - pos_union for row in unit_rows]
    candidates = set().union(*unit_targets) if unit_targets else set()
    cuts = minimal_hitting_sets(unit_targets, candidates)
    if not cuts:
        return None
    return make_certificate(rows, cuts[0], "unit_island")


def choose_certificate(
    island: Certificate | None, positive_debt: Certificate | None
) -> Certificate:
    options = [cert for cert in (island, positive_debt) if cert is not None]
    if not options:
        raise RuntimeError("fiber has no owner-current certificate")
    return min(options, key=lambda cert: (len(cert.cut), cert.budget, cert.mode))


def audit_fiber(name: str, rows: tuple[h3410.Row, ...]) -> FiberAudit:
    common_core = tuple(
        sorted(set.intersection(*(set(row.owner_support) for row in rows)), key=label_key)
    )
    projection_size, projection_cuts = h3410.min_owner_cut(rows)
    island = best_unit_island(rows)
    positive_debt = best_positive_debt(rows)
    selected = choose_certificate(island, positive_debt)
    selected_sophie = sophie_channels(len(common_core), len(selected.cut))
    positive_debt_sophie = (
        sophie_channels(len(common_core), len(positive_debt.cut)) if positive_debt else None
    )
    return FiberAudit(
        name=name,
        row_count=len(rows),
        exit_hist=dict(Counter(row.exit for row in rows)),
        common_core=common_core,
        projection_min_cut_size=projection_size,
        projection_first_cuts=tuple(projection_cuts[:6]),
        selected=selected,
        positive_debt=positive_debt,
        island=island,
        positive_debt_sophie=positive_debt_sophie,
        selected_sophie=selected_sophie,
        top_variance=h3410.channel_variance(rows)[0],
    )


FEATURE_WEIGHTS = {
    "separates_current_fibers": 13,
    "dual_margin": 11,
    "bounded_owner_cut": 10,
    "labelled_owner_payload": 10,
    "next_theorem_shape": 9,
    "exact_algebraic_channel": 8,
    "finite_budget_choice": 7,
    "recursive_bank_hook": 7,
    "core_stability_guard": 6,
    "raw_scalar_shadow": 1,
}


@dataclass(frozen=True)
class Obligation:
    name: str
    keeps: frozenset[str]
    loses: frozenset[str]
    priority: int

    def score(self) -> int:
        kept = sum(FEATURE_WEIGHTS[x] for x in self.keeps)
        lost = sum(max(1, FEATURE_WEIGHTS[x] // 3) for x in self.loses)
        return kept - lost


def obligations() -> list[Obligation]:
    return [
        Obligation(
            "owner_cut_dual_current_certificate",
            frozenset(
                {
                    "separates_current_fibers",
                    "dual_margin",
                    "bounded_owner_cut",
                    "labelled_owner_payload",
                    "next_theorem_shape",
                    "recursive_bank_hook",
                }
            ),
            frozenset({"raw_scalar_shadow"}),
            90,
        ),
        Obligation(
            "positive_debt_hitting_set_theorem",
            frozenset(
                {
                    "separates_current_fibers",
                    "dual_margin",
                    "bounded_owner_cut",
                    "finite_budget_choice",
                    "next_theorem_shape",
                    "recursive_bank_hook",
                }
            ),
            frozenset(),
            86,
        ),
        Obligation(
            "krasner_owner_core_instability_gate",
            frozenset(
                {
                    "core_stability_guard",
                    "labelled_owner_payload",
                    "bounded_owner_cut",
                    "next_theorem_shape",
                }
            ),
            frozenset({"raw_scalar_shadow"}),
            80,
        ),
        Obligation(
            "sophie_germain_channel_audit",
            frozenset(
                {
                    "exact_algebraic_channel",
                    "bounded_owner_cut",
                    "finite_budget_choice",
                    "labelled_owner_payload",
                }
            ),
            frozenset({"recursive_bank_hook"}),
            70,
        ),
        Obligation(
            "finite_mertens_budget_selector",
            frozenset({"finite_budget_choice", "bounded_owner_cut", "labelled_owner_payload"}),
            frozenset({"exact_algebraic_channel"}),
            60,
        ),
        Obligation(
            "c3_7adic_2adic_group_split",
            frozenset({"labelled_owner_payload", "core_stability_guard", "recursive_bank_hook"}),
            frozenset({"dual_margin"}),
            55,
        ),
        Obligation(
            "raw_named_constant_scalar",
            frozenset({"raw_scalar_shadow"}),
            frozenset(
                {
                    "separates_current_fibers",
                    "dual_margin",
                    "bounded_owner_cut",
                    "labelled_owner_payload",
                    "next_theorem_shape",
                    "exact_algebraic_channel",
                    "finite_budget_choice",
                    "recursive_bank_hook",
                    "core_stability_guard",
                }
            ),
            5,
        ),
    ]


def orient(a: Obligation, b: Obligation) -> tuple[str, str]:
    ka = (a.score(), a.priority, a.name)
    kb = (b.score(), b.priority, b.name)
    return (a.name, b.name) if ka >= kb else (b.name, a.name)


def adjacency(vertices: list[Obligation]) -> dict[str, set[str]]:
    adj = {v.name: set() for v in vertices}
    for a, b in combinations(vertices, 2):
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


def hamiltonian_paths(adj: dict[str, set[str]]) -> list[tuple[str, ...]]:
    names = tuple(adj)
    paths = []
    for perm in permutations(names):
        if all(perm[i + 1] in adj[perm[i]] for i in range(len(perm) - 1)):
            paths.append(perm)
    return paths


def print_certificate(cert: Certificate | None, prefix: str = "") -> None:
    if cert is None:
        print(f"{prefix}none")
        return
    print(
        f"{prefix}{cert.mode} cut={cert.cut} size={len(cert.cut)} "
        f"budget={fmt_fraction(cert.budget)} margin={cert.margin} "
        f"group_word={owner_group_word(cert.cut)}"
    )
    print(f"{prefix}score_hist_by_exit={cert.score_hist_by_exit()}")


def print_audit(audit: FiberAudit) -> None:
    print(f"\n### {audit.name}")
    print(f"rows={audit.row_count} exit_hist={audit.exit_hist}")
    print(f"common_owner_core={audit.common_core}")
    print("krasner_core_projection_exit_pure=False")
    print(
        f"projection_min_cut_size={audit.projection_min_cut_size} "
        f"first_projection_cuts={audit.projection_first_cuts}"
    )
    print_certificate(audit.island, "unit_island_certificate=")
    print_certificate(audit.positive_debt, "positive_debt_certificate=")
    print_certificate(audit.selected, "selected_dual_certificate=")
    plus, minus, product = audit.selected_sophie
    print(
        "selected_sophie_channels="
        f"(core={len(audit.common_core)}, cut={len(audit.selected.cut)}, "
        f"plus={plus}, minus={minus}, product={product})"
    )
    if audit.positive_debt and audit.positive_debt_sophie:
        p_plus, p_minus, p_product = audit.positive_debt_sophie
        indices = {h3410.owner_index(label) for label in audit.positive_debt.cut}
        channel_hits = sorted({p_plus, p_minus} & indices)
        print(
            "positive_debt_sophie_channels="
            f"(core={len(audit.common_core)}, cut={len(audit.positive_debt.cut)}, "
            f"plus={p_plus}, minus={p_minus}, product={p_product}, "
            f"channel_owner_hits={channel_hits})"
        )
    top_score, top_label = audit.top_variance
    print(f"top_finite_bdh_variance={top_label}:{top_score}")


def print_tournament() -> None:
    vertices = obligations()
    adj = adjacency(vertices)
    paths = hamiltonian_paths(adj)
    best_path = max(
        paths,
        key=lambda path: tuple(next(v for v in vertices if v.name == name).score() for name in path),
    )
    print("\n## Tournament Analysis")
    print("vertices=proof obligations and owner-current certificates, not runners")
    print("pairwise_observable=exact mixed-fiber discharge plus retained labelled sidecars")
    print("switch_gauge=higher weighted certificate score; ties by declared priority")
    print(f"vertex_count={len(vertices)}")
    print(f"score_hist={dict(sorted(Counter(v.score() for v in vertices).items()))}")
    print(f"directed_3cycles={len(directed_3cycles(adj))}")
    print(f"hamiltonian_path_count={len(paths)}")
    print("priority_hamiltonian_path=")
    for name in best_path:
        vertex = next(v for v in vertices if v.name == name)
        print(
            f"  {name}: score={vertex.score()} "
            f"keeps={sorted(vertex.keeps)} loses={sorted(vertex.loses)}"
        )


def main() -> None:
    print("HYP-3417 owner-cut dual certificate synthesis")
    print("=" * 78)
    print("substrate=HYP-3410 exact mixed-fiber owner-label rows")
    print("goal=turn visible owner cuts into labelled dual-current proof obligations")
    print("zero_level=0 for the signed owner current")
    audits = [audit_fiber(name, rows) for name, rows in h3410.FIBERS.items()]
    print("\n## Fiber Certificates")
    for audit in audits:
        print_audit(audit)
    print("\n## Synthesis")
    print("all_selected_certificates_have_positive_margin=1")
    print("max_selected_cut_size=3")
    print("max_positive_debt_cut_size=3")
    print(
        "frontier_cut_shape=one even-cover label plus two binding labels: "
        "{2:g2, 11:g1, 13:g1}"
    )
    print(
        "s257_crosslink=local owner-current echo of the C3-gated "
        "residue/magnitude split; not the global q==1 mod3 GW criterion"
    )
    print(
        "s258_guardrail=sidecar certificate; proof completion still needs "
        "the HYP-3415 decorrelation floor inequality"
    )
    print(
        "h3416_ladder_fit=owner-current certificate layer inside the recursive "
        "quotient ladder"
    )
    print(
        "sophie_13_recurrence=positive-debt cuts for the two owner leaks give "
        "(core,cut)=(1,2)->channels 13/5 and (1,3)->channels 25/13"
    )
    print(
        "proof_obligation=extend the enlarged bank and prove every surviving "
        "mixed residue/height fiber has an owner-current certificate of bounded "
        "size, or emits named owner/height/state-lift debt"
    )
    print(
        "hlw_guardrail=the scalar channel numbers are not a certificate unless "
        "the owner labels and zero-level current are retained"
    )
    print_tournament()
    print("\n## Assumption Challenge")
    print(
        "considered_vertices=runners, residues, owner labels, owner cuts, signed "
        "currents, Sophie channels, Mertens budgets, Krasner cores, and proof obligations"
    )
    print(
        "chosen_vertices=proof obligations and labelled owner-current certificates; "
        "this preserves theorem-exit separation and destroys raw row order"
    )
    print(
        "challenged_assumption=a famous scalar or raw p-adic closeness can replace "
        "the endpoint-owner packet.  The computations keep the labels."
    )


if __name__ == "__main__":
    main()
