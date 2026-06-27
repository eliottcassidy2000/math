#!/usr/bin/env python3
"""Hidden statement ledger for the current LRC14 proof stack.

This is a synthesis pass, not a full-bank recomputation.  It records small
statements that became visible only after several recent carriers were placed
beside older work: path/deletion topology, Haar product zeta, Ramanujan period
decks, analytic sieve blindness, automatic-word guardrails, and residual
capacitor cuts.

Tournament Analysis declaration:
  vertices: micro-statements / proof obligations, not runners;
  pairwise observable: predicate sharpness, localizer strength, noncircularity,
    past-work connection count, compression, and theorem-actionability;
  switch/gauge: majority comparison of observable vectors;
  tie Hamiltonian path: boundary H1 -> first tooth -> primitive deck ->
    capacitor cut -> zeta square -> topology-scale tooth -> stalk descent ->
    analytic blindness -> automaton shadow -> decoy generator.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class MicroStatement:
    ident: str
    handle: str
    source: str
    role: str
    statement: str
    evidence: str
    keeps: tuple[str, ...]
    destroys: tuple[str, ...]
    past_threads: tuple[str, ...]
    vector: tuple[int, int, int, int, int, int]


STATEMENTS: tuple[MicroStatement, ...] = (
    MicroStatement(
        "MS01",
        "owner_essential_boundary_h1",
        "HYP-3034/HYP-3030/HYP-3025",
        "front_gate",
        "Zero-open is not just beta1=1; on the tested collision surface it is an owner-essential closed danger-arc H1 representative.",
        "AP and GW are the only closed-H1 rows; each has a 58-edge representative and deleting any owner speed kills H1.",
        ("status", "topology", "owner", "boundary", "deletion"),
        ("route", "exact_magnitude", "runner_order"),
        ("path_homology", "deletion_persistence", "arc_cech"),
        (5, 5, 5, 5, 4, 5),
    ),
    MicroStatement(
        "MS02",
        "coarse_status_gate_before_routes",
        "HYP-3024/HYP-3028/HYP-3030/HYP-3034",
        "ordering_rule",
        "The coarse ET+Henselian gate is a status theorem before it is a route theorem.",
        "Full-bank coarse ET+unit fibers have 0 mixed boundary/open fibers; their route-mixed survivors are all strict-open.",
        ("status", "et_hensel", "residual", "ordering"),
        ("route", "exact_magnitude", "full_packet_identity"),
        ("fiber_zipper", "henselian_units", "erdos_turan"),
        (5, 3, 5, 4, 5, 4),
    ),
    MicroStatement(
        "MS03",
        "residual_first_tooth_owner_strip",
        "HYP-3035/HYP-3033/HYP-3029",
        "first_tooth_manifest",
        "The 15 coarse residual route fibers are a finite owner-strip atlas, not a residual count.",
        "Arc topology is first tooth for 13 fibers, coarse safe stalk for 2, and every first repair is owner_strip.",
        ("status", "route", "owner", "topology", "stalk", "residual"),
        ("route_labels_until_late", "exact_magnitude_until_backup"),
        ("safe_component_stalk", "barcode", "normal_fan", "haar_repair"),
        (5, 5, 4, 5, 4, 5),
    ),
    MicroStatement(
        "MS04",
        "primitive_q13_boundary",
        "HYP-3036/HYP-2979/HYP-2978",
        "period_scheduler",
        "The direct-witness layer is q<=13 primitive safe mass; q=14 belongs to the boundary/covering layer.",
        "On the 38 S194 residual packets, every Q-WITNESS row has positive primitive safe mass for q<=13 and every covering row has zero.",
        ("status", "route", "period", "ramanujan", "q_witness"),
        ("exact_magnitude", "safe_interval_length", "raw_trace_without_inequality"),
        ("ramanujan_projector", "exact_period", "divisor_packets"),
        (5, 4, 5, 5, 4, 5),
    ),
    MicroStatement(
        "MS05",
        "analytic_capacity_blindness",
        "HYP-3032/HYP-2982/HYP-2983",
        "analytic_guardrail",
        "Analytic sieve clocks are capacity and blindness meters, not LRC quotients.",
        "mu^2/phi kills C27 prime powers and fibbinary q=25, while squarefree q=23 still needs packet geometry.",
        ("analytic_clock", "capacity", "blindness", "residual"),
        ("prime_power_packets", "owner", "topology", "route"),
        ("large_sieve", "mobius", "kaczynski", "smoothing"),
        (3, 4, 5, 5, 5, 4),
    ),
    MicroStatement(
        "MS06",
        "residual_capacitor_min_cut",
        "HYP-3037/HYP-3027/HYP-3031",
        "cut_lattice",
        "After status purity, a mixed open-route pair is a capacitor whose proof content is the first nonzero cut.",
        "The petal/covering capacitor is cut by boundary topology; the K33/covering capacitor is cut by exact M+q.",
        ("status", "route", "cut", "topology", "magnitude", "residual"),
        ("universal_scalar_order", "runner_identity"),
        ("max_flow_min_cut", "repair_ladder", "haar_zeta"),
        (4, 5, 5, 5, 3, 5),
    ),
    MicroStatement(
        "MS07",
        "q23_diagonal_zeta_owner_strip",
        "HYP-3038/HYP-3032/HYP-3031",
        "local_square",
        "The squarefree q=23 analytic residual is a drop/add Haar square with a diagonal zeta and endpoint-owner strip.",
        "Diagonal corners keep M=2/23, off-diagonal corners open as q=10 and q=8 witnesses, and endpoint-owner strips split petal from covering.",
        ("route", "zeta", "owner", "drop_add_square", "endpoint", "residual"),
        ("raw_analytic_q23", "row_column_shadow", "endpoint_scalar_B18Z6"),
        ("haar_product", "fixed_margin", "drop_add_grid"),
        (4, 5, 5, 5, 3, 5),
    ),
    MicroStatement(
        "MS08",
        "topology_bucket_plus_unit_scale",
        "HYP-3033/HYP-3030",
        "compressed_scheduler",
        "A tiny topology bucket plus unit-scale tooth schedules all stored coarse residual Q/covering route debt.",
        "Topology alone leaves 3 mixed route classes, unit-scale alone leaves 1, but the joined tooth has 21 fibers and 0 route mixing.",
        ("status", "route", "topology", "unit_scale", "compression"),
        ("exact_magnitude", "coarse_et_address", "row_identity"),
        ("cech_topology", "henselian_unit_rule", "fiber_zipper"),
        (5, 4, 4, 4, 5, 4),
    ),
    MicroStatement(
        "MS09",
        "safe_stalk_as_geometry_not_magnitude",
        "HYP-3029/HYP-3035/HYP-3018",
        "geometric_descent",
        "The safe-component stalk is a geometric descent of magnitude, not a numerical fallback.",
        "Exact stalk splits the target automatic fiber like exact M, while coarse stalk supplies the two missing residual first teeth.",
        ("status", "route", "stalk", "owner", "safe_component"),
        ("full_barcode", "global_magnitude_scalar"),
        ("normal_fan", "barcode", "safe_component_stalk"),
        (4, 4, 4, 5, 4, 4),
    ),
    MicroStatement(
        "MS10",
        "automatic_words_are_row_column_shadows",
        "HYP-3023/HYP-3031/HYP-3008/HYP-3010",
        "shadow_guardrail",
        "Moser/fibbinary/automatic words are useful telemetry but unsafe quotients unless a mixed coordinate or packet label returns.",
        "The AP/GW word mixes routes badly; q=23 drop/add row-column shadows forget zeta; exact packet labels or zeta/owner strips repair the loss.",
        ("automaton", "zeta", "packet_label", "shadow", "sequence"),
        ("route", "boundary", "mixed_coordinate"),
        ("moser_de_bruijn", "fibbinary", "ostrowski_hadamard", "finite_automata"),
        (3, 3, 5, 5, 5, 3),
    ),
    MicroStatement(
        "MS11",
        "pair_good_decoys_are_generator_teeth",
        "HYP-3021/HYP-3022/HYP-3019",
        "decoy_normal_form",
        "Pair-good decoys are modular blocker teeth, not a raw abundance problem.",
        "The blocker rule 14*min(c*p mod q,q-c*p mod q)<q collapses thousands of decoys to bounded lane/blocker generators tied to barcode and normal-fan supports.",
        ("decoy", "blocker", "barcode", "normal_fan", "owner"),
        ("raw_decoy_count", "pair_good_boolean"),
        ("binding_pair_switch", "barcode", "normal_fan"),
        (3, 4, 5, 4, 5, 4),
    ),
)


TIE_PATH = {
    handle: index
    for index, handle in enumerate(
        (
            "owner_essential_boundary_h1",
            "residual_first_tooth_owner_strip",
            "primitive_q13_boundary",
            "residual_capacitor_min_cut",
            "q23_diagonal_zeta_owner_strip",
            "topology_bucket_plus_unit_scale",
            "safe_stalk_as_geometry_not_magnitude",
            "analytic_capacity_blindness",
            "automatic_words_are_row_column_shadows",
            "pair_good_decoys_are_generator_teeth",
            "coarse_status_gate_before_routes",
        )
    )
}


def compare(a: MicroStatement, b: MicroStatement) -> int:
    awins = sum(1 for x, y in zip(a.vector, b.vector) if x > y)
    bwins = sum(1 for x, y in zip(a.vector, b.vector) if y > x)
    if awins > bwins:
        return 1
    if bwins > awins:
        return -1
    return 1 if TIE_PATH[a.handle] < TIE_PATH[b.handle] else -1


def tournament_fingerprint(statements: tuple[MicroStatement, ...]) -> dict[str, object]:
    n = len(statements)
    edges: dict[tuple[int, int], bool] = {}
    scores = [0] * n
    for i, j in combinations(range(n), 2):
        if compare(statements[i], statements[j]) > 0:
            edges[(i, j)] = True
            scores[i] += 1
        else:
            edges[(j, i)] = True
            scores[j] += 1

    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (
            edges.get((i, j), False)
            and edges.get((j, k), False)
            and edges.get((k, i), False)
        ) or (
            edges.get((i, k), False)
            and edges.get((k, j), False)
            and edges.get((j, i), False)
        ):
            c3 += 1

    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if edges.get((last, nxt), False):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count

    full = (1 << n) - 1
    hp = sum(dp.get((full, last), 0) for last in range(n))
    order = sorted(range(n), key=lambda i: (-scores[i], TIE_PATH[statements[i].handle]))
    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3cycles": c3,
        "hamiltonian_path_count": hp,
        "score_order": [statements[i].handle for i in order],
        "scores": dict((statements[i].handle, scores[i]) for i in range(n)),
    }


def connection_weight(a: MicroStatement, b: MicroStatement) -> tuple[int, tuple[str, ...], tuple[str, ...]]:
    shared_keeps = tuple(sorted(set(a.keeps) & set(b.keeps)))
    shared_past = tuple(sorted(set(a.past_threads) & set(b.past_threads)))
    weight = 2 * len(shared_keeps) + len(shared_past)
    return weight, shared_keeps, shared_past


def main() -> None:
    print("=== LRC14 hidden statement ledger S203 ===")
    print("statements=", len(STATEMENTS), sep="")
    print()

    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, speeds, gaps, route labels, automatic words, exact M values,")
    print("    cover arcs, proof carriers, micro-statements, and proof obligations.")
    print("  chosen vertices:")
    print("    micro-statements / proof obligations extracted from recent and past")
    print("    LRC14 work, because the user asked for hidden details and connections.")
    print("  preserved LRC predicate:")
    print("    boundary/open status, post-status route schedulability, and the")
    print("    named repair route by which a quotient recovers lost information.")
    print("  destroyed information:")
    print("    raw row enumeration, scalar ranking, and route labels until a legal")
    print("    topology, period, owner, stalk, zeta, or cut carrier reattaches them.")
    print()

    print("[1] Micro-statement ledger")
    for st in STATEMENTS:
        print(f"  {st.ident} {st.handle} ({st.source})")
        print(f"    role={st.role}")
        print(f"    statement={st.statement}")
        print(f"    evidence={st.evidence}")
        print(f"    keeps={','.join(st.keeps)}")
        print(f"    destroys={','.join(st.destroys)}")
        print(f"    past_threads={','.join(st.past_threads)}")
    print()

    print("[2] Hidden connection edges")
    edges = []
    for a, b in combinations(STATEMENTS, 2):
        weight, keeps, past = connection_weight(a, b)
        if weight:
            edges.append((weight, a, b, keeps, past))
    for weight, a, b, keeps, past in sorted(edges, key=lambda item: (-item[0], item[1].ident, item[2].ident))[:18]:
        print(
            f"  {a.ident}-{b.ident} weight={weight} "
            f"shared_keeps={keeps or '-'} shared_past={past or '-'}"
        )
    print()

    print("[3] Field histogram")
    field_counts = Counter(field for st in STATEMENTS for field in st.keeps)
    for field, count in field_counts.most_common():
        print(f"  {field:<18} {count}")
    print()

    print("[4] Tournament Analysis over hidden statements")
    fp = tournament_fingerprint(STATEMENTS)
    print("  vertices_are=micro-statements / proof obligations, not runners")
    print("  score_hist=" + str(fp["score_hist"]))
    print("  directed_3cycles=" + str(fp["directed_3cycles"]))
    print("  hamiltonian_path_count=" + str(fp["hamiltonian_path_count"]))
    print("  score_order=" + " > ".join(fp["score_order"]))
    print("  scores=" + str(fp["scores"]))
    print()

    print("[5] Fundamental readout")
    print("  1. LRC14 is no longer looking like one master scalar.  It is a")
    print("     layered obstruction calculus: owner-essential boundary H1 first,")
    print("     status-preserving coarse gates second, route teeth third, and")
    print("     residual capacitors/zeta squares only at the end.")
    print("  2. The AP/GW equality atom is topological and deletion-critical,")
    print("     not merely a tight maximin value.")
    print("  3. The sharp q-boundary is q<=13 versus q=14: q<=13 primitive safe")
    print("     mass is direct-witness currency; q=14 is boundary/covering currency.")
    print("  4. Owner identity is repeatedly the missing coordinate: first teeth are")
    print("     owner strips, endpoint-current words can be lossy, and deleting any")
    print("     AP/GW owner kills the closed-H1 representative.")
    print("  5. Past scalar shadows become useful once labelled by their blindness:")
    print("     automatic words are row/column shadows, mu^2/phi is a squarefree")
    print("     capacity meter, and pair-good decoys are blocker-generator teeth.")
    print()

    print("[6] Sidecar theorem target")
    print("  Add a hidden-statement sidecar vocabulary to the HYP-2963 packet ledger:")
    print("    boundary_h1_owner_support, first_tooth, primitive_safe_deck_2_13,")
    print("    residual_capacitor_id, first_cut_stage, exact_M_zeta,")
    print("    endpoint_owner_strip, analytic_blindness_report, automaton_shadow_class.")
    print("  Candidate principle:")
    print("    Every primitive packet either carries the AP/GW owner-essential")
    print("    boundary H1 atom, has a q<=13 primitive witness, descends through")
    print("    an owner/stalk/topology tooth, or exposes a named capacitor/zeta")
    print("    cut leading to covering, petal, K33/THM-572, Haar/Fejer/Ramanujan")
    print("    annihilation, or F7 debt.")


if __name__ == "__main__":
    main()
