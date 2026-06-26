#!/usr/bin/env python3
"""S237: cycle-class and observability scout for two remaining LRC14 angles.

This is deliberately a proof-interface scout, not a proof of LRC14.  It takes
two live angles from HYP-3066 and HYP-3063/HYP-3065:

1. the exact HYP-2963 residual-fiber observability matrix already exposed by
   S199/S200, and
2. a small rational cycle-class matrix whose columns are named certificate
   generators and whose only intentionally ungenerated atom is F7/THM-572
   phantom debt.

Tournament Analysis vertices are proof carriers / certificate columns, not
runners.
"""

from __future__ import annotations

from ast import literal_eval
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from pathlib import Path
import re


REPO = Path(__file__).resolve().parents[1]
S199_OUT = REPO / "05-knowledge" / "results" / "lrc14_residual_tooth_atlas_codex_s199.out"
S200_OUT = REPO / "05-knowledge" / "results" / "lrc14_ramanujan_route_scheduler_codex_s200.out"


@dataclass(frozen=True)
class ResidualFiber:
    index: int
    size: int
    routes: dict[str, int]
    status: dict[str, int]
    q0_range: str
    first_tooth: str
    repair_class: str
    topology_pure: bool
    coarse_stalk_pure: bool
    exact_stalk_pure: bool


@dataclass(frozen=True)
class Generator:
    name: str
    vector: tuple[Fraction, ...]
    cost: int
    role: str


def parse_residual_fibers() -> list[ResidualFiber]:
    pattern = re.compile(
        r"^\s+fiber\[(?P<index>\d+)\]\s+"
        r"size=(?P<size>\d+)\s+"
        r"routes=(?P<routes>\{.*?\})\s+"
        r"status=(?P<status>\{.*?\})\s+"
        r"q0_range=(?P<q0_range>\S+)\s+"
        r"first=(?P<first>[^/\s]+)/(?P<repair>\S+)\s+"
        r"topology_pure=(?P<topology>True|False)\s+"
        r"coarse_stalk_pure=(?P<coarse>True|False)\s+"
        r"exact_stalk_pure=(?P<exact>True|False)"
    )
    fibers: list[ResidualFiber] = []
    for line in S199_OUT.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if not match:
            continue
        groups = match.groupdict()
        fibers.append(
            ResidualFiber(
                index=int(groups["index"]),
                size=int(groups["size"]),
                routes=literal_eval(groups["routes"]),
                status=literal_eval(groups["status"]),
                q0_range=groups["q0_range"],
                first_tooth=groups["first"],
                repair_class=groups["repair"],
                topology_pure=groups["topology"] == "True",
                coarse_stalk_pure=groups["coarse"] == "True",
                exact_stalk_pure=groups["exact"] == "True",
            )
        )
    if len(fibers) != 15 or sum(f.size for f in fibers) != 38:
        raise RuntimeError(f"unexpected S199 shape: fibers={len(fibers)} rows={sum(f.size for f in fibers)}")
    return fibers


def parse_s200_mixed_route(label: str) -> int:
    pattern = re.compile(rf"^\s+{re.escape(label)}\s+.*mixed_route=(\d+)")
    for line in S200_OUT.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line)
        if match:
            return int(match.group(1))
    raise RuntimeError(f"missing S200 line for {label}")


def observability_columns(fibers: list[ResidualFiber]) -> dict[str, list[int]]:
    primitive_pure = parse_s200_mixed_route("coarse_plus_primitive_deck_2_13") == 0
    first_q_pure = parse_s200_mixed_route("coarse_plus_first_q_2_13") == 0
    return {
        "arc_topology_compact": [int(f.topology_pure) for f in fibers],
        "coarse_safe_stalk": [int(f.coarse_stalk_pure) for f in fibers],
        "exact_safe_stalk": [int(f.exact_stalk_pure) for f in fibers],
        "magnitude_cocycle": [1 for _ in fibers],
        "first_primitive_safe_q_2_13": [int(first_q_pure) for _ in fibers],
        "primitive_safe_deck_2_13": [int(primitive_pure) for _ in fibers],
    }


def minimal_covers(columns: dict[str, list[int]]) -> list[tuple[str, ...]]:
    names = list(columns)
    target = [1] * len(next(iter(columns.values())))
    covers: list[tuple[str, ...]] = []
    for r in range(1, len(names) + 1):
        for combo in combinations(names, r):
            seen = [0] * len(target)
            for name in combo:
                seen = [int(a or b) for a, b in zip(seen, columns[name])]
            if seen == target:
                covers.append(combo)
        if covers:
            return covers
    return covers


def unit_vector(basis: list[str], *names: str) -> tuple[Fraction, ...]:
    selected = set(names)
    return tuple(Fraction(1 if name in selected else 0) for name in basis)


def add_vectors(*vectors: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(sum(values, Fraction(0)) for values in zip(*vectors))


def rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    mat = [row[:] for row in matrix]
    if not mat:
        return mat, []
    rows, cols = len(mat), len(mat[0])
    pivots: list[int] = []
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if mat[i][c]), None)
        if pivot is None:
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        scale = mat[r][c]
        mat[r] = [x / scale for x in mat[r]]
        for i in range(rows):
            if i != r and mat[i][c]:
                factor = mat[i][c]
                mat[i] = [x - factor * y for x, y in zip(mat[i], mat[r])]
        pivots.append(c)
        r += 1
        if r == rows:
            break
    return mat, pivots


def matrix_rank(columns: list[tuple[Fraction, ...]]) -> int:
    if not columns:
        return 0
    rows = [[col[i] for col in columns] for i in range(len(columns[0]))]
    _, pivots = rref(rows)
    return len(pivots)


def solve_columns(
    columns: list[tuple[Fraction, ...]], target: tuple[Fraction, ...]
) -> tuple[bool, list[Fraction]]:
    """Solve A x = target and return one exact solution with free vars zero."""
    if not columns:
        return (not any(target), [])
    m = len(target)
    n = len(columns)
    augmented = [[columns[j][i] for j in range(n)] + [target[i]] for i in range(m)]
    reduced, pivots = rref(augmented)
    for row in reduced:
        if not any(row[:n]) and row[n]:
            return False, []
    solution = [Fraction(0) for _ in range(n)]
    for i, pivot_col in enumerate(pivots):
        if pivot_col < n:
            solution[pivot_col] = reduced[i][n]
    return True, solution


def cycle_class_system() -> tuple[list[str], list[Generator], dict[str, tuple[Fraction, ...]]]:
    basis = [
        "ap_gw_closed_boundary_h1",
        "endpoint_owner_current",
        "primitive_period_mass",
        "haar_zeta_square",
        "observer_cut_payload",
        "rectangle_hourglass_residue",
        "partial_cube_theta",
        "simplex_face_boundary",
        "low_height_wall",
        "octahedral_face_curl",
        "toeplitz_scale_gate",
        "k33_state_lift",
        "phantom_f7_class",
    ]
    e = lambda *names: unit_vector(basis, *names)
    generators = [
        Generator("ap_gw_boundary_cycle", e("ap_gw_closed_boundary_h1"), 1, "boundary"),
        Generator("endpoint_owner_boundary", e("endpoint_owner_current"), 1, "owner"),
        Generator("ramanujan_period_character", e("primitive_period_mass"), 2, "period"),
        Generator("haar_zipper_square", e("haar_zeta_square"), 2, "haar"),
        Generator("observer_cut_boundary", e("observer_cut_payload"), 2, "observer"),
        Generator("rectangle_hourglass_cycle", e("rectangle_hourglass_residue"), 2, "cycle"),
        Generator("partial_cube_theta_class", e("partial_cube_theta"), 2, "partial_cube"),
        Generator("simplex_face_boundary", e("simplex_face_boundary"), 2, "simplex"),
        Generator("roth_minkowski_low_height_wall", e("low_height_wall"), 3, "diophantine"),
        Generator("octahedral_face_curl", e("octahedral_face_curl"), 2, "hodge"),
        Generator("toeplitz_square_scale_gate", e("toeplitz_scale_gate"), 2, "noncollapse"),
        Generator("k33_state_lift_incidence", e("k33_state_lift"), 3, "state_lift"),
    ]
    targets = {
        "ap_gw_boundary_atom": add_vectors(e("ap_gw_closed_boundary_h1"), e("endpoint_owner_current")),
        "residual_owner_tooth": e("endpoint_owner_current"),
        "primitive_period_scheduler": e("primitive_period_mass"),
        "q23_drop_add_square": add_vectors(e("haar_zeta_square"), e("endpoint_owner_current")),
        "observer_duodecimal_defect": add_vectors(e("observer_cut_payload"), e("rectangle_hourglass_residue")),
        "partial_cube_bridge_debt": add_vectors(
            e("partial_cube_theta"), e("simplex_face_boundary"), e("rectangle_hourglass_residue")
        ),
        "support_six_octahedral_debt": add_vectors(e("low_height_wall"), e("octahedral_face_curl")),
        "toeplitz_four_witness_gate": e("toeplitz_scale_gate"),
        "k33_state_lift_exit": e("k33_state_lift"),
        "f7_phantom_probe": add_vectors(e("k33_state_lift"), e("phantom_f7_class")),
    }
    return basis, generators, targets


def count_directed_3cycles(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        if ((a, b) in edges and (b, c) in edges and (c, a) in edges) or (
            (a, c) in edges and (c, b) in edges and (b, a) in edges
        ):
            total += 1
    return total


def scc_sizes(vertices: list[str], edges: set[tuple[str, str]]) -> list[int]:
    adjacency = defaultdict(list)
    reverse = defaultdict(list)
    for a, b in edges:
        adjacency[a].append(b)
        reverse[b].append(a)

    def reach(start: str, graph: dict[str, list[str]]) -> set[str]:
        seen = {start}
        todo = deque([start])
        while todo:
            node = todo.popleft()
            for nxt in graph[node]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        return seen

    remaining = set(vertices)
    sizes = []
    while remaining:
        start = next(iter(remaining))
        comp = reach(start, adjacency) & reach(start, reverse)
        sizes.append(len(comp))
        remaining -= comp
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(vertices: list[str], edges: set[tuple[str, str]]) -> int:
    n = len(vertices)
    index = {v: i for i, v in enumerate(vertices)}
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if (vertices[last], vertices[nxt]) in edges:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[-1])


def tournament_fingerprint(
    observability: dict[str, list[int]], generators: list[Generator], basis: list[str]
) -> tuple[list[str], set[tuple[str, str]], dict[str, tuple[int, int, int]]]:
    generator_by_name = {g.name: g for g in generators}
    vertices = [
        "coarse_safe_stalk",
        "arc_topology_compact",
        "primitive_safe_deck_2_13",
        "endpoint_owner_boundary",
        "haar_zipper_square",
        "observer_cut_boundary",
        "rectangle_hourglass_cycle",
        "partial_cube_theta_class",
        "simplex_face_boundary",
        "octahedral_face_curl",
        "roth_minkowski_low_height_wall",
        "toeplitz_square_scale_gate",
        "k33_state_lift_incidence",
        "phantom_f7_marker",
    ]
    support_scores: dict[str, tuple[int, int, int]] = {}
    for v in vertices:
        sep = sum(observability.get(v, []))
        if v in generator_by_name:
            atom_count = sum(1 for x in generator_by_name[v].vector if x)
            cost_bonus = 5 - generator_by_name[v].cost
        elif v == "coarse_safe_stalk":
            atom_count = 2
            cost_bonus = 3
        elif v == "arc_topology_compact":
            atom_count = 1
            cost_bonus = 4
        elif v == "primitive_safe_deck_2_13":
            atom_count = 1
            cost_bonus = 3
        else:
            atom_count = 0
            cost_bonus = 0
        support_scores[v] = (sep, atom_count, cost_bonus)
    priority = {v: i for i, v in enumerate(vertices)}
    edges: set[tuple[str, str]] = set()
    for a, b in combinations(vertices, 2):
        key_a = (support_scores[a], -priority[a])
        key_b = (support_scores[b], -priority[b])
        edges.add((a, b) if key_a > key_b else (b, a))
    return vertices, edges, support_scores


def fmt_fraction(fr: Fraction) -> str:
    if fr == 0:
        return "0"
    if fr == 1:
        return "1"
    if fr == -1:
        return "-1"
    return f"{fr.numerator}/{fr.denominator}" if fr.denominator != 1 else str(fr.numerator)


def main() -> None:
    fibers = parse_residual_fibers()
    columns = observability_columns(fibers)

    print("=== LRC14 cycle-class / observability scout S237 ===")
    print("source=S199 residual tooth atlas + S200 primitive-period scheduler")
    print(f"residual_fibers={len(fibers)} residual_packets={sum(f.size for f in fibers)}")
    print()
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, endpoint owners, residual fibers, automaton states,")
    print("    partial-cube Theta classes, rectangle/hourglass cycles, primitive")
    print("    period decks, Hodge/cohomology classes, state-lift sectors, and proof obligations.")
    print("  chosen vertices:")
    print("    proof obligations, separating sidecar columns, and certificate-cycle generators.")
    print("  preserved LRC predicate:")
    print("    boundary/open status and Q-WITNESS vs COVERING-MOMENT route separability")
    print("    inside the 15 strict-open HYP-2963 residual coarse ET+unit fibers.")
    print("  destroyed information:")
    print("    exact row identity and exact M unless later sidecars reattach them.")
    print("  challenged assumption:")
    print("    the remaining proof problem is not another scalar count; it is an")
    print("    observability / cycle-generation problem over labelled packet fibers.")

    print()
    print("[1] HYP-2963 residual-fiber observability matrix")
    header = "fiber  size  first_tooth                 " + "  ".join(f"{name:>28s}" for name in columns)
    print(header)
    for f in fibers:
        row = "  ".join(f"{columns[name][f.index]:28d}" for name in columns)
        print(f"{f.index:02d}     {f.size:1d}  {f.first_tooth + '/' + f.repair_class:28s}  {row}")
    print()
    for name, values in columns.items():
        print(f"  column {name:28s} separates {sum(values):2d}/{len(values)} residual fibers")
    print(f"  first_tooth_counts={dict(sorted(Counter(f.first_tooth for f in fibers).items()))}")
    print(f"  repair_class_counts={dict(sorted(Counter(f.repair_class for f in fibers).items()))}")
    print(f"  minimal_single_or_small_covers={minimal_covers(columns)[:8]}")
    print("  reading: arc topology is the first tooth in 13 fibers; the two topology")
    print("  collisions are repaired by the coarse safe-component stalk.  Primitive")
    print("  q<=13 decks and magnitude cocycles also separate all 15 fibers, but they")
    print("  are later/higher-payload exits rather than the first local proof tooth.")

    print()
    print("[2] Exact rational cycle-class matrix template")
    basis, generators, targets = cycle_class_system()
    generator_columns = [g.vector for g in generators]
    known_rank = matrix_rank(generator_columns)
    print(f"  basis_dimension={len(basis)} known_generator_count={len(generators)} known_rank={known_rank}")
    print(f"  unresolved_basis_atoms={[basis[i] for i in range(len(basis)) if all(g.vector[i] == 0 for g in generators)]}")
    print("  generators:")
    for g in generators:
        support = [basis[i] for i, x in enumerate(g.vector) if x]
        print(f"    {g.name:34s} cost={g.cost} role={g.role:13s} support={support}")
    print("  target rows:")
    generated = 0
    for target_name, target in targets.items():
        ok, solution = solve_columns(generator_columns, target)
        if ok:
            generated += 1
            terms = [
                f"{fmt_fraction(coeff)}*{generators[i].name}"
                for i, coeff in enumerate(solution)
                if coeff
            ]
            print(f"    {target_name:31s} status=generated_by_known_cycles  decomposition={terms}")
        else:
            support = [basis[i] for i, x in enumerate(target) if x]
            print(f"    {target_name:31s} status=phantom_unresolved        support={support}")
    print(f"  generated_targets={generated}/{len(targets)}")
    print("  reading: the matrix deliberately spans every named certificate atom except")
    print("  `phantom_f7_class`.  Therefore F7 should mean a concrete residual basis")
    print("  vector outside the current certificate image, not an anonymous hard row.")

    print()
    print("[3] Two proof angles extracted")
    print("  angle A: observability-first residual proof")
    print("    Work fiberwise.  The route leak after coarse ET+unit is killed by the")
    print("    first separating sidecar: arc topology for 13 fibers, safe-component")
    print("    stalk for the remaining 2.  Later primitive decks and magnitude data")
    print("    are fallback certificates, not replacements for the local tooth order.")
    print("  angle B: cycle-generation / Hodge-style proof")
    print("    Treat every residual cochain as a row in Cert(P)->H_res(P).  Partial")
    print("    cubes, simplex faces, rectangle residues, endpoint owners, Ramanujan")
    print("    periods, Haar squares, and state-lifts are columns.  The proof closes")
    print("    only when each row is generated, dual-annihilated, descended, boundary-")
    print("    stopped, or emitted as an explicit F7/THM-572 class.")

    print()
    print("[4] Tournament Analysis")
    vertices, edges, scores = tournament_fingerprint(columns, generators, basis)
    outdegree = Counter()
    for a, b in edges:
        outdegree[a] += 1
        outdegree.setdefault(b, outdegree[b])
    score_hist = Counter(outdegree[v] for v in vertices)
    hpath = sorted(vertices, key=lambda v: (scores[v], -vertices.index(v)), reverse=True)
    print("  vertices_are=proof carriers / sidecar columns / certificate generators, not runners")
    print("  pairwise_observable=(residual_fiber_separation_count, cycle_atom_support, inverse_payload_cost)")
    print("  gauge=orient toward the carrier retaining more route/status/cycle payload")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={count_directed_3cycles(vertices, edges)}")
    print(f"  scc_sizes={scc_sizes(vertices, edges)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(vertices, edges)}")
    print("  tie_hamiltonian_path=" + " > ".join(hpath))

    print()
    print("[5] Next exact job")
    print("  Replace the template target rows with actual HYP-2963 packet cochains:")
    print("  for each residual fiber, emit the cochain basis coordinates for topology,")
    print("  owner current, primitive deck, zeta square, observer-cut payload,")
    print("  partial-cube Theta/simplex sidecar, median owner/root fields, and")
    print("  state-lift target.  Then rerun the exact row-reduction and promote")
    print("  only rows with a legal image status.")


if __name__ == "__main__":
    main()
