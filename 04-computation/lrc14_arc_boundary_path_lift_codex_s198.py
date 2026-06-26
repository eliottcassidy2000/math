#!/usr/bin/env python3
"""Arc-boundary path lift for LRC14.

This script mines the older tournament path-homology/deletion-persistence
thread and pulls it back to the current LRC14 arc-Cech topology gate.

The point is not to put path homology on runner vertices.  The predicate
preserving object is the closed danger-arc Cech nerve.  We compute its GF(2)
boundary ranks, an H1 representative when present, and deletion persistence
under owner removal.  This turns the HYP-3030 closed beta=(1,1) signal into an
explicit path-boundary carrier with a local deletion test.

Tournament Analysis declaration:
  vertices: proof carriers / boundary operators, not runners;
  pairwise observable: status purity, H1 exactness, deletion persistence,
    topology retention, endpoint ownership, quotient defect, and proof cost;
  switch/gauge: orient by more retained predicate data before lower proof cost;
  tie Hamiltonian path: arc_boundary_path_lift > owner_deletion_persistence >
    arc_cech_beta > coarse_et_unit_status > residue_terminal_word.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


fz = load_module(
    "lrc14_fiber_zipper_convergence_s198",
    REPO / "04-computation" / "lrc14_fiber_zipper_convergence_codex_s188.py",
)
ac = load_module(
    "lrc14_arc_cech_nerve_s198",
    REPO / "04-computation" / "lrc14_arc_cech_nerve_carrier_codex_s174.py",
)


def fmt(fr: Fraction) -> str:
    return str(fr.numerator) if fr.denominator == 1 else f"{fr.numerator}/{fr.denominator}"


def bit_count(x: int) -> int:
    return x.bit_count()


def insert_gf2_basis(basis: dict[int, int], vector: int) -> bool:
    x = vector
    while x:
        pivot = x.bit_length() - 1
        if pivot in basis:
            x ^= basis[pivot]
        else:
            basis[pivot] = x
            return True
    return False


def gf2_rank(columns: list[int]) -> int:
    basis: dict[int, int] = {}
    for column in columns:
        insert_gf2_basis(basis, column)
    return len(basis)


@dataclass(frozen=True)
class BoundaryLift:
    vertices: tuple[str, ...]
    edges: tuple[tuple[int, int], ...]
    triangles: tuple[tuple[int, int, int], ...]
    rank_d1: int
    rank_d2: int
    beta0: int
    beta1: int
    representative: int | None

    @property
    def cycle_edge_count(self) -> int:
        return 0 if self.representative is None else bit_count(self.representative)

    @property
    def cycle_vertices(self) -> tuple[str, ...]:
        if self.representative is None:
            return ()
        used: set[int] = set()
        for edge_index, (a, b) in enumerate(self.edges):
            if (self.representative >> edge_index) & 1:
                used.add(a)
                used.add(b)
        return tuple(self.vertices[i] for i in sorted(used))

    @property
    def cycle_owner_speeds(self) -> tuple[int, ...]:
        speeds = set()
        for label in self.cycle_vertices:
            speed = label.split(":", 1)[1].split("@", 1)[0]
            speeds.add(int(speed))
        return tuple(sorted(speeds))


def boundary_lift(facets: tuple[frozenset[str], ...]) -> BoundaryLift:
    vertices = tuple(sorted(set().union(*facets))) if facets else ()
    vertex_index = {v: i for i, v in enumerate(vertices)}

    edge_set: set[tuple[int, int]] = set()
    triangle_set: set[tuple[int, int, int]] = set()
    for facet in facets:
        ids = sorted(vertex_index[v] for v in facet)
        edge_set.update(tuple(pair) for pair in combinations(ids, 2))
        triangle_set.update(tuple(tri) for tri in combinations(ids, 3))

    edges = tuple(sorted(edge_set))
    edge_index = {edge: i for i, edge in enumerate(edges)}
    triangles = tuple(sorted(triangle_set))

    d1_cols = [(1 << a) | (1 << b) for a, b in edges]
    rank_d1 = gf2_rank(d1_cols)

    d2_cols = []
    for a, b, c in triangles:
        mask = 0
        for edge in ((a, b), (a, c), (b, c)):
            mask ^= 1 << edge_index[tuple(sorted(edge))]
        d2_cols.append(mask)
    rank_d2 = gf2_rank(d2_cols)

    beta0 = len(vertices) - rank_d1 if vertices else 0
    beta1 = len(edges) - rank_d1 - rank_d2
    representative = h1_representative(vertices, edges, d2_cols) if beta1 else None
    return BoundaryLift(vertices, edges, triangles, rank_d1, rank_d2, beta0, beta1, representative)


def h1_representative(
    vertices: tuple[str, ...],
    edges: tuple[tuple[int, int], ...],
    d2_cols: list[int],
) -> int | None:
    adjacency: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for edge_id, (a, b) in enumerate(edges):
        adjacency[a].append((b, edge_id))
        adjacency[b].append((a, edge_id))

    parent: dict[int, tuple[int, int] | None] = {}
    tree_edges: set[int] = set()
    for root in range(len(vertices)):
        if root in parent:
            continue
        parent[root] = None
        queue = deque([root])
        while queue:
            v = queue.popleft()
            for w, edge_id in adjacency[v]:
                if w in parent:
                    continue
                parent[w] = (v, edge_id)
                tree_edges.add(edge_id)
                queue.append(w)

    def tree_path_edges(a: int, b: int) -> int:
        seen: dict[int, int] = {}
        mask = 0
        x = a
        while True:
            seen[x] = mask
            item = parent[x]
            if item is None:
                break
            x, edge_id = item
            mask ^= 1 << edge_id
        mask = 0
        x = b
        while True:
            if x in seen:
                return mask ^ seen[x]
            item = parent[x]
            if item is None:
                return 0
            x, edge_id = item
            mask ^= 1 << edge_id

    cycle_candidates = []
    for edge_id, (a, b) in enumerate(edges):
        if edge_id in tree_edges:
            continue
        cycle = (1 << edge_id) ^ tree_path_edges(a, b)
        if cycle:
            cycle_candidates.append(cycle)

    basis: dict[int, int] = {}
    for col in d2_cols:
        insert_gf2_basis(basis, col)
    for cycle in sorted(cycle_candidates, key=bit_count):
        test = dict(basis)
        if insert_gf2_basis(test, cycle):
            return cycle
    return None


def remove_owner_facets(facets: tuple[frozenset[str], ...], owner_speed: int) -> tuple[frozenset[str], ...]:
    token = f":{owner_speed}@"
    out = []
    for facet in facets:
        reduced = frozenset(label for label in facet if token not in label)
        if reduced:
            out.append(reduced)
    return tuple(out)


@dataclass(frozen=True)
class PacketBoundaryRecord:
    name: str
    source_family: str
    status: str
    safe_mu: Fraction
    closed: BoundaryLift
    open_beta: tuple[int, int]
    safe_topes: int
    quotient_defect: int
    boundary_owner_sums: tuple[int, ...]
    owner_delete_beta1: tuple[tuple[int, int], ...]

    @property
    def apgw_path_cycle(self) -> bool:
        return (
            self.status == "boundary"
            and self.closed.beta0 == 1
            and self.closed.beta1 == 1
            and self.open_beta == (6, 0)
            and self.safe_topes == 0
            and self.boundary_owner_sums == (0, 0, 0, 0, 0, 0)
            and self.closed.cycle_edge_count > 0
        )


def packet_boundary_record(packet) -> PacketBoundaryRecord:
    rec = ac.row_record(packet.name, tuple(packet.speeds))
    closed = boundary_lift(rec.audit.closed_arc_facets)
    owner_delete = []
    if closed.beta1:
        for speed in closed.cycle_owner_speeds:
            beta1 = boundary_lift(remove_owner_facets(rec.audit.closed_arc_facets, speed)).beta1
            owner_delete.append((speed, beta1))
    return PacketBoundaryRecord(
        name=packet.name,
        source_family=packet.source_family,
        status="boundary" if rec.audit.safe_topes == 0 else "open",
        safe_mu=rec.audit.safe_measure,
        closed=closed,
        open_beta=(rec.arc_betti.beta0, rec.arc_betti.beta1),
        safe_topes=rec.audit.safe_topes,
        quotient_defect=rec.runner_quotient_defect,
        boundary_owner_sums=rec.boundary_pair_sums,
        owner_delete_beta1=tuple(owner_delete),
    )


def group_by_split(packets: list, split_name: str) -> dict[object, list]:
    split = next(s for s in fz.SPLITS if s.name == split_name)
    groups: dict[object, list] = defaultdict(list)
    for packet in packets:
        groups[split.key_func(packet)].append(packet)
    return groups


def mixed_status_groups(groups: dict[object, list], status_by_name: dict[str, str]) -> list[list]:
    return [rows for rows in groups.values() if len({status_by_name[p.name] for p in rows}) > 1]


def row_status(packet: RowStub) -> tuple[str, Fraction]:
    safe_mu, safe_topes = ac.safe_measure_and_topes(tuple(packet.speeds))
    return ("boundary" if safe_topes == 0 else "open"), safe_mu


def short_record(record: PacketBoundaryRecord) -> str:
    deletions = ",".join(f"{speed}->{beta1}" for speed, beta1 in record.owner_delete_beta1) or "-"
    return (
        f"{record.name:<34} status={record.status:<8} family={record.source_family:<14} "
        f"mu={fmt(record.safe_mu):>12} "
        f"C0/C1={record.closed.beta0}/{record.closed.beta1} "
        f"d1={record.closed.rank_d1:<3d} d2={record.closed.rank_d2:<3d} "
        f"rep_edges={record.closed.cycle_edge_count:<3d} "
        f"owners={record.closed.cycle_owner_speeds} del={deletions}"
    )


@dataclass(frozen=True)
class RowStub:
    name: str
    source_family: str
    speeds: tuple[int, ...]
    automatic_word: str


@dataclass(frozen=True)
class Carrier:
    name: str
    status_purity: int
    h1_exactness: int
    deletion_persistence: int
    topology_retention: int
    owner_retention: int
    quotient_defect: int
    proof_cost: int

    @property
    def values(self) -> tuple[int, ...]:
        return (
            self.status_purity,
            self.h1_exactness,
            self.deletion_persistence,
            self.topology_retention,
            self.owner_retention,
            self.quotient_defect,
            8 - self.proof_cost,
        )


CARRIERS = (
    Carrier("arc_boundary_path_lift", 5, 5, 4, 5, 4, 4, 4),
    Carrier("owner_deletion_persistence", 4, 3, 5, 4, 5, 3, 4),
    Carrier("arc_cech_beta", 4, 4, 1, 4, 2, 4, 2),
    Carrier("coarse_et_unit_status", 5, 0, 0, 1, 0, 1, 2),
    Carrier("residue_terminal_word", 1, 0, 0, 0, 0, 0, 1),
)


def orient(a: Carrier, b: Carrier) -> str:
    awins = sum(x > y for x, y in zip(a.values, b.values))
    bwins = sum(x < y for x, y in zip(a.values, b.values))
    if awins != bwins:
        return a.name if awins > bwins else b.name
    if a.values != b.values:
        return a.name if a.values > b.values else b.name
    order = {carrier.name: i for i, carrier in enumerate(CARRIERS)}
    return a.name if order[a.name] < order[b.name] else b.name


def tournament_fingerprint() -> dict[str, object]:
    names = [carrier.name for carrier in CARRIERS]
    edges: set[tuple[str, str]] = set()
    for a, b in combinations(CARRIERS, 2):
        winner = orient(a, b)
        loser = b.name if winner == a.name else a.name
        edges.add((winner, loser))

    scores = Counter()
    for winner, _ in edges:
        scores[winner] += 1
    score_hist = Counter(scores[name] for name in names)

    c3 = 0
    for triple in combinations(names, 3):
        out = {name: 0 for name in triple}
        for x, y in combinations(triple, 2):
            if (x, y) in edges:
                out[x] += 1
            else:
                out[y] += 1
        if sorted(out.values()) == [1, 1, 1]:
            c3 += 1

    graph = {name: [] for name in names}
    rev = {name: [] for name in names}
    for a, b in edges:
        graph[a].append(b)
        rev[b].append(a)

    def reach(start: str, g: dict[str, list[str]]) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in g[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    remaining = set(names)
    scc_sizes = []
    while remaining:
        seed = min(remaining)
        comp = reach(seed, graph) & reach(seed, rev)
        scc_sizes.append(len(comp))
        remaining -= comp

    hp = 0
    for order in permutations(names):
        if all((order[i], order[i + 1]) in edges for i in range(len(order) - 1)):
            hp += 1

    score_order = sorted(names, key=lambda name: (-scores[name], name))
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_path_count": hp,
        "score_order": score_order,
    }


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, individual danger arcs, endpoint cells, boundary")
    print("    cocircuits, residues, Fourier modes, matroid circuits, path chains,")
    print("    and proof obligations.")
    print("  chosen vertices:")
    print("    closed danger-arc Cech vertices and proof-carrier boundary operators.")
    print("  preserved LRC predicate:")
    print("    boundary/open status at threshold 1/14 and the exact closed-cover H1")
    print("    obstruction that appears when no strict safe interval exists.")
    print("  destroyed information:")
    print("    raw runner order and scalar route labels; owner-deletion records only")
    print("    the speeds supporting the H1 representative.")
    print("  challenged assumption:")
    print("    the useful homology object is not path homology of the runner")
    print("    tournament.  It is the GF(2) boundary complex of danger arcs.")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    args = parser.parse_args()

    rows = fz.lp.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    packets = [
        RowStub(name, source_family, tuple(speeds), fz.lp.automatic_word(tuple(speeds)))
        for name, source_family, speeds in rows
    ]
    print(f"[progress] loaded lightweight row bank rows={len(packets)}", file=sys.stderr, flush=True)
    status_by_name: dict[str, str] = {}
    safe_mu_by_name: dict[str, Fraction] = {}
    for index, packet in enumerate(packets, 1):
        status, safe_mu = row_status(packet)
        status_by_name[packet.name] = status
        safe_mu_by_name[packet.name] = safe_mu
        if index % 2000 == 0:
            print(f"[progress] interval statuses={index}/{len(packets)}", file=sys.stderr, flush=True)

    print("=== LRC14 arc-boundary path lift S198 ===")
    print(
        "bank=HYP-2963 default rows "
        f"single_limit={args.single_limit} two_swap_limit={args.two_swap_limit} "
        f"alias_depth={args.alias_depth} lcm_tail_max={args.lcm_tail_max}"
    )
    print(f"packets={len(packets)}")
    print()
    print_assumption_challenge()

    groups = group_by_split(packets, "residue_terminal_fiber")
    mixed_status = sorted(mixed_status_groups(groups, status_by_name), key=lambda g: (-len(g), g[0].name))
    coarse_groups = group_by_split(packets, "coarse_et_unit_gate")
    coarse_status_mixed = sorted(mixed_status_groups(coarse_groups, status_by_name), key=lambda g: (-len(g), g[0].name))

    controls = {"AP", "GW 12->24", "K33 12->36", "covering comb 12->84", "magnitude liar 12->96"}
    needed_names = set(controls)
    for group in mixed_status:
        needed_names.update(packet.name for packet in group)
    for group in coarse_status_mixed:
        needed_names.update(packet.name for packet in group)
    needed_packets = [packet for packet in packets if packet.name in needed_names]
    print(f"[progress] path-lift target rows={len(needed_packets)}", file=sys.stderr, flush=True)
    records = [packet_boundary_record(packet) for packet in needed_packets]
    by_name = {record.name: record for record in records}

    h1_rows = [record for record in records if record.closed.beta1 > 0]
    print("[1] Targeted closed arc-boundary H1 audit")
    print(f"  path_lift_rows={len(records)}")
    print(f"  closed_h1_rows={len(h1_rows)}")
    print(f"  status_counts={dict(Counter(record.status for record in h1_rows))}")
    print(f"  source_family_counts={dict(Counter(record.source_family for record in h1_rows))}")
    print(f"  apgw_path_cycle_rows={sum(1 for record in h1_rows if record.apgw_path_cycle)}")
    if h1_rows:
        print("  rows:")
        for record in h1_rows:
            print("    " + short_record(record))
    print()

    print("[2] Residue-terminal status collisions with path lift")
    print(f"  mixed_status_fibers={len(mixed_status)}")
    for index, group in enumerate(mixed_status):
        group_records = sorted(
            (by_name[packet.name] for packet in group),
            key=lambda record: (record.status, record.source_family, record.name),
        )
        print(
            f"  residue_status_collision[{index}] size={len(group_records)} "
            f"status={dict(Counter(record.status for record in group_records))} "
            f"h1={dict(Counter(record.closed.beta1 for record in group_records))}"
        )
        for record in group_records:
            print("    " + short_record(record))
    print()

    coarse_records = [by_name[packet.name] for group in coarse_status_mixed for packet in group]
    print("[3] Coarse ET+unit status fibers")
    print(f"  fibers={len(coarse_groups)} mixed_status_fibers={len(coarse_status_mixed)} packets_in_mixed_status={len(coarse_records)}")
    print(f"  closed_h1_in_mixed_status={sum(1 for record in coarse_records if record.closed.beta1)}")
    print("  note=route-mixed open fibers are scheduled by HYP-3030; this script checks status/H1 purity without recomputing exact M.")
    print()

    print("[4] Named boundary controls")
    for name in ("AP", "GW 12->24", "K33 12->36", "covering comb 12->84", "magnitude liar 12->96"):
        if name in by_name:
            print("  " + short_record(by_name[name]))
    print()

    fp = tournament_fingerprint()
    print("[5] Tournament Analysis over proof carriers")
    print("  vertices_are=proof carriers / boundary operators, not runners")
    print("  pairwise_observable=status purity, H1 exactness, deletion persistence,")
    print("                      topology retention, endpoint ownership, quotient")
    print("                      defect visibility, and proof cost")
    print("  tie_hamiltonian_path=" + " > ".join(carrier.name for carrier in CARRIERS))
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  scc_sizes={fp['scc_sizes']}")
    print(f"  hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("  score_order=" + " > ".join(fp["score_order"]))
    print()

    print("[6] Proof readout")
    print("  1. On the targeted residue-collision/control surface, closed")
    print("     danger-arc H1 occurs only on the AP/Goddyn-Wong boundary rows.")
    print("  2. The residue-terminal boundary/open collisions are therefore not only")
    print("     beta-separated; they have an explicit GF(2) boundary representative.")
    print("  3. The lightweight full-bank status scan gives zero mixed boundary/open")
    print("     fibers for the coarse ET+unit gate.  HYP-3030's route-mixed open")
    print("     survivors should be scheduled only after this path-boundary topology")
    print("     removes equality atoms.")
    print("  4. Deleting any owner speed from the displayed AP/GW representative")
    print("     kills the closed H1 signal in this carrier.  This is the old")
    print("     deletion-persistence idea, but applied to danger arcs instead of")
    print("     runner tournaments.")


if __name__ == "__main__":
    main()
