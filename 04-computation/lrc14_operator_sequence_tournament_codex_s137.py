#!/usr/bin/env python3
"""S137: operator sequences and carrier tournaments for the LRC14 Farey route.

This is a synthesis layer over S130--S136.  It returns to the prompt's four
Farey mutations

    p + q, p*q, q^p, p^q

and treats them first as integer sequences along the LRC14 unit-excess chain

    p/q = p/(14p - 1),

then as binary-relation carriers.  The goal is not to prove LRC14 in this
file.  The goal is to preserve the information that the current LRC14 proof
route appears to need before any scalar quotient is applied.

Topics deliberately kept in one comparison:
  1. Farey unit-excess sequence fingerprints,
  2. the additive n+2 lane versus the multiplicative n*2/product lane,
  3. the C=27 shell-transfer quotient,
  4. the K_{p,q}/K33 incidence wall,
  5. the octahedral L(K4) support-six current carrier,
  6. the Clebsch / halved-cube covariance carrier,
  7. Paley-Zygmund as a second-moment gateway.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd, log
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
N = 14
THR = F(1, N)
C27 = 2 * N - 1


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s124 = load_module(
    "s124_apgw_for_s137",
    REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py",
)
s132 = load_module(
    "s132_graphs_for_s137",
    REPO / "04-computation" / "lrc14_farey_graph_pz_carriers_codex_s132.py",
)


def edge(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        raise ValueError("loop")
    if i < j:
        bit = i * n - i * (i + 1) // 2 + (j - i - 1)
        return bool(mask & (1 << bit))
    return not edge(mask, n, j, i)


def hamiltonian_path_count(mask: int, n: int) -> int:
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for state in range(1 << n):
        for last in range(n):
            val = dp[state][last]
            if not val:
                continue
            for nxt in range(n):
                if state & (1 << nxt):
                    continue
                if edge(mask, n, last, nxt):
                    dp[state | (1 << nxt)][nxt] += val
    return sum(dp[full])


def tournament_fingerprint(mask: int, n: int) -> dict[str, object]:
    scores = tuple(sum(1 for j in range(n) if i != j and edge(mask, n, i, j)) for i in range(n))
    c3 = 0
    for tri in combinations(range(n), 3):
        local = [sum(1 for j in tri if i != j and edge(mask, n, i, j)) for i in tri]
        if sorted(local) == [1, 1, 1]:
            c3 += 1

    adj = {i: {j for j in range(n) if i != j and edge(mask, n, i, j)} for i in range(n)}
    radj = {i: set() for i in range(n)}
    for i, outs in adj.items():
        for j in outs:
            radj[j].add(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    comps: list[int] = []

    def rdfs(v: int, comp: list[int]) -> None:
        seen.add(v)
        comp.append(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[int] = []
            rdfs(v, comp)
            comps.append(len(comp))

    return {
        "score_hist": tuple(sorted(Counter(scores).items())),
        "c3": c3,
        "scc": tuple(sorted(comps, reverse=True)),
        "hp": hamiltonian_path_count(mask, n),
    }


def finite_differences(values: list[int]) -> list[list[int]]:
    layers = [values]
    cur = values
    while len(cur) > 1:
        cur = [b - a for a, b in zip(cur, cur[1:])]
        layers.append(cur)
    return layers


def polynomial_degree(values: list[int], max_degree: int = 5) -> str:
    layers = finite_differences(values)
    for degree, layer in enumerate(layers[: max_degree + 1]):
        if len(layer) > 1 and len(set(layer)) == 1:
            return str(degree)
    return f">{max_degree}"


def unit_excess_chain(limit: int) -> list[dict[str, object]]:
    rows = []
    for p in range(1, limit + 1):
        q = N * p - 1
        rows.append(
            {
                "p": p,
                "value": F(p, q),
                "gap": F(1, N * q),
                "q": q,
                "sum": p + q,
                "product": p * q,
                "log_qpow": p * log(q),
                "log_ppow": q * log(p) if p > 1 else 0.0,
                "k_rank": "star" if p == 1 else "two-block" if p == 2 else "K33-wall",
            }
        )
    return rows


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]
    family: str


@dataclass(frozen=True)
class Shell:
    a: int
    gcd27: int
    lower_count: int
    upper_count: int

    @property
    def total(self) -> int:
        return self.lower_count + self.upper_count

    @property
    def detail(self) -> str:
        return f"{self.a}:g{self.gcd27}"


@dataclass(frozen=True)
class RowItem:
    label: str
    speeds: tuple[int, ...]
    M: F
    branch: str
    q_threshold: int
    shells: tuple[Shell, ...]
    zero_count: int

    @property
    def holes(self) -> tuple[Shell, ...]:
        return tuple(shell for shell in self.shells if shell.total == 0)

    @property
    def doubles(self) -> tuple[Shell, ...]:
        return tuple(shell for shell in self.shells if shell.total >= 2)


def candidate_rows(max_replacement: int) -> list[Row]:
    ap = tuple(range(1, 14))
    rows = [
        Row("AP", ap, "known tight"),
        Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "known tight"),
        Row("residue-liar 12->26", tuple(list(range(1, 12)) + [13, 26]), "q-threshold loose"),
        Row("near-miss 12->36", tuple(list(range(1, 12)) + [13, 36]), "Farey child loose"),
    ]
    for r in range(1, 14):
        for v in range(14, max_replacement + 1):
            row = tuple(sorted((set(ap) - {r}) | {v}))
            if len(row) == 13:
                rows.append(Row(f"swap {r}->{v}", row, "single AP replacement"))
    for m in range(2, 13):
        row = tuple(sorted(set(list(range(1, 12)) + [13, 12 * m])))
        if len(row) == 13:
            rows.append(Row(f"12m family m={m}", row, "12m tail"))

    seen: set[tuple[int, ...]] = set()
    out: list[Row] = []
    for row in rows:
        if row.speeds not in seen:
            seen.add(row.speeds)
            out.append(row)
    return out


def exact_M(speeds: tuple[int, ...]) -> F:
    M, _pts = s124.M_exact(speeds)
    return M


def farey_branch(M: F) -> str:
    p, q = M.numerator, M.denominator
    e = N * p - q
    if e == 0:
        return "tight-floor"
    if e == 1 and p == 1:
        return "q-parent-star"
    if e == 1 and p == 2 and q == C27:
        return "C27-petal-two-block"
    if e == 1 and p >= 3:
        return "K33-unit-excess"
    if e > 1:
        return "nonunit-excess"
    return "below-or-other"


def shell_profile(speeds: tuple[int, ...]) -> tuple[tuple[Shell, ...], int]:
    residues = Counter(v % C27 for v in speeds)
    shells = []
    for a in range(1, (C27 + 1) // 2):
        shells.append(
            Shell(
                a=a,
                gcd27=gcd(a, C27),
                lower_count=residues[a],
                upper_count=residues[C27 - a],
            )
        )
    return tuple(shells), residues[0]


def transfer_label(item: RowItem, detailed: bool = False) -> str:
    if not item.holes and not item.doubles and item.zero_count == 0:
        return "perfect-transversal"
    holes = ",".join(shell.detail for shell in item.holes) or "-"
    doubles = ",".join(
        shell.detail if not detailed else f"{shell.a}:g{shell.gcd27}:{shell.lower_count}{shell.upper_count}"
        for shell in item.doubles
    ) or "-"
    zero = str(item.zero_count) if item.zero_count else "-"
    return f"H[{holes}] D[{doubles}] Z[{zero}]"


def coarse_transfer_label(item: RowItem) -> str:
    if not item.holes and not item.doubles and item.zero_count == 0:
        return "perfect"
    h = tuple(shell.gcd27 for shell in item.holes)
    d = tuple(shell.gcd27 for shell in item.doubles)
    return f"Hg{h or '-'}->Dg{d or '-'};z={item.zero_count}"


def build_items(max_replacement: int) -> list[RowItem]:
    items = []
    for row in candidate_rows(max_replacement):
        M = exact_M(row.speeds)
        shells, zero_count = shell_profile(row.speeds)
        items.append(
            RowItem(
                label=row.label,
                speeds=row.speeds,
                M=M,
                branch=farey_branch(M),
                q_threshold=s124.q_threshold(row.speeds),
                shells=shells,
                zero_count=zero_count,
            )
        )
    return items


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, gaps, Farey fractions, fraction payloads, C=27 shells,")
    print("    Kpq incidence owners, octahedral packet currents, Clebsch cuts,")
    print("    second-moment events, exact-period packets, and proof obligations.")
    print("  chosen vertices:")
    print("    operator/carrier roles, plus a C=27 shell quotient when the p=2")
    print("    branch is isolated.")
    print("  quotient preserves:")
    print("    exact Farey excess e=14p-q, binding denominator q, additive and")
    print("    product side channels, low-gap shell transfer, and state-lift fit.")
    print("  quotient destroys:")
    print("    exact time geometry and row-specific witnesses; those remain attached")
    print("    through M(S), q-threshold, and the C27 transfer label.")
    print("  challenged assumption:")
    print("    a single runner tournament or scalar payload is enough.  The live")
    print("    predicate needs a labelled operator packet before scalarization.")


def print_sequence_table(limit: int) -> None:
    rows = unit_excess_chain(limit)
    print()
    print("[1] Unit-excess Farey sequence fingerprints")
    print("  chain: p/q = p/(14p-1), so M-1/14 = 1/(14q)")
    print(
        f"  {'p/q':>8s} {'q':>5s} {'p+q':>5s} {'p*q':>6s} "
        f"{'gap':>11s} {'Kpq':>9s} {'q mod27':>7s} {'p*q mod27':>9s}"
    )
    for row in rows:
        print(
            f"  {str(row['value']):>8s} {row['q']:5d} {row['sum']:5d} "
            f"{row['product']:6d} {str(row['gap']):>11s} {row['k_rank']:>9s} "
            f"{row['q'] % C27:7d} {row['product'] % C27:9d}"
        )
    print()
    for key in ("q", "sum", "product"):
        vals = [int(row[key]) for row in rows]
        layers = finite_differences(vals)
        first = layers[1][:6]
        second = layers[2][:6] if len(layers) > 2 else []
        print(
            f"  {key:7s}: detected polynomial degree {polynomial_degree(vals)}; "
            f"Delta={first}; Delta2={second}"
        )
    print("  q^p and p^q are intentionally not polynomial lanes; they are")
    print("  magnitude-stress tests for quotients that forget scale.")
    print()
    print("  Readout:")
    print("    q      is the theorem-level binding scale: +14 on the chain.")
    print("    p+q    is the additive/Stern-Brocot lane: +15 on the chain.")
    print("    p*q    is the product/coimage area lane: quadratic, Delta2=28.")
    print("    q^p,p^q amplify magnitude and should not be used as denominators.")


def print_row_bank(max_replacement: int) -> None:
    items = build_items(max_replacement)
    branch_counts = Counter(item.branch for item in items)
    low_341 = [item for item in items if item.M <= F(3, 41)]
    low_227 = [item for item in items if item.M <= F(2, 27)]
    transfer_counts = Counter(coarse_transfer_label(item) for item in items)
    print()
    print("[2] Row-bank bridge: exact Farey branch plus C=27 transfer")
    print(f"  S130 row bank through replacement {max_replacement}: {len(items)} rows")
    print("  branch counts:")
    for branch, count in branch_counts.most_common():
        print(f"    {branch:23s} {count:5d}")
    print()
    print("  low frontier:")
    print(f"    M <= 3/41 : {len(low_341)} rows -> {[item.label for item in low_341]}")
    print(f"    M <= 2/27 : {len(low_227)} rows -> {[item.label for item in low_227]}")
    print()
    print("  named transfer signatures:")
    for label in ("AP", "GW 12->24", "swap 10->20", "swap 13->26", "near-miss 12->36"):
        item = next((x for x in items if x.label == label), None)
        if item is None:
            continue
        print(
            f"    {item.label:18s} M={str(item.M):>5s} "
            f"branch={item.branch:23s} transfer={transfer_label(item, detailed=False)}"
        )
    print()
    print("  most common transfer labels, including loose repeats:")
    for label, count in transfer_counts.most_common(8):
        print(f"    {label:28s} {count:5d}")
    print("  Guardrail: transfer labels recur, so exact M/Farey branch remains the")
    print("  first coordinate.  The C27 quotient is a router, not a classifier.")


def print_graph_and_pz() -> None:
    print()
    print("[3] Graph and Paley-Zygmund carrier fingerprints")
    graphs = [
        ("octahedron L(K4)", s132.octahedron_line_k4()),
        ("Clebsch folded Q5", s132.folded_five_cube_clebsch()),
        ("halved 5-cube", s132.halved_cube(5)),
    ]
    for name, graph in graphs:
        stats = s132.graph_stats(graph)
        print(
            f"  {name:19s}: v={stats['vertices']:2d} e={stats['edges']:2d} "
            f"tri={stats['triangles']:3d} cycle_rank={stats['cycle_rank']:2d} "
            f"deg={stats['degree_hist']}"
        )
    print(
        "  halved 5-cube == complement(Clebsch) under even representatives: "
        f"{s132.halved_cube_is_complement_clebsch()}"
    )
    print()
    print("  PZ toy loss on six independent empty-sector indicators:")
    for k in (8, 10, 12):
        row = s132.binomial_miss_model(k)
        loss = row["exact_union"] - row["PZ"]
        print(
            f"    k={k}: exact={float(row['exact_union']):.6f} "
            f"PZ={float(row['PZ']):.6f} loss={float(loss):.6f}"
        )
    print("  Readout: octahedron carries support-six curl (cycle rank 7),")
    print("  Clebsch/half-cube carry folded covariance/cut data, and PZ is an")
    print("  existence gateway unless upgraded to labelled degree-4 moments.")


@dataclass(frozen=True)
class Carrier:
    key: str
    name: str
    role: str


CARRIERS = (
    Carrier("q", "q binding scale", "exact LRC gap denominator"),
    Carrier("sum", "p+q additive lane", "Stern-Brocot / n+2 recursion"),
    Carrier("prod", "p*q product lane", "Kpq area / n*2 coimage growth"),
    Carrier("C27", "C=27 shell", "p=2 second-gap transfer router"),
    Carrier("K33", "Kpq/K33 incidence", "p>=3 three-owner packet"),
    Carrier("oct", "octahedron L(K4)", "support-six current/curl"),
    Carrier("cleb", "Clebsch/half-cube", "folded covariance/cut carrier"),
    Carrier("PZ", "Paley-Zygmund", "second-moment existence gateway"),
    Carrier("power", "power payloads", "magnitude-leak stress tests"),
)


CRITERIA: dict[str, list[str]] = {
    "theorem_scale": ["q", "C27", "sum", "K33", "prod", "oct", "cleb", "PZ", "power"],
    "sequence_simplicity": ["q", "sum", "prod", "C27", "K33", "oct", "cleb", "PZ", "power"],
    "branch_separation": ["q", "C27", "K33", "sum", "prod", "oct", "cleb", "PZ", "power"],
    "state_lift_fit": ["K33", "oct", "cleb", "C27", "q", "prod", "PZ", "sum", "power"],
    "graph_packet_fit": ["oct", "cleb", "K33", "prod", "C27", "q", "PZ", "sum", "power"],
    "anti_scalar_guard": ["q", "C27", "K33", "oct", "cleb", "sum", "prod", "PZ", "power"],
    "magnitude_stress": ["power", "prod", "K33", "cleb", "PZ", "sum", "C27", "oct", "q"],
}


def majority_tournament() -> tuple[int, dict[str, set[str]], dict[tuple[str, str], tuple[str, int, int]]]:
    keys = [c.key for c in CARRIERS]
    ranks = {crit: {key: i for i, key in enumerate(order)} for crit, order in CRITERIA.items()}
    tie_path = {key: i for i, key in enumerate(CRITERIA["theorem_scale"])}
    wins: dict[str, set[str]] = {key: set() for key in keys}
    details: dict[tuple[str, str], tuple[str, int, int]] = {}
    mask = 0
    bit = 0
    for i, j in combinations(range(len(keys)), 2):
        a, b = keys[i], keys[j]
        av = bv = 0
        for crit in CRITERIA:
            if ranks[crit][a] < ranks[crit][b]:
                av += 1
            else:
                bv += 1
        if av > bv or (av == bv and tie_path[a] < tie_path[b]):
            winner, loser = a, b
        else:
            winner, loser = b, a
        wins[winner].add(loser)
        details[(a, b)] = (winner, av, bv)
        if winner == a:
            mask |= 1 << bit
        bit += 1
    return mask, wins, details


def conservative_tournament() -> int:
    keys = [c.key for c in CARRIERS]
    order = {key: i for i, key in enumerate(CRITERIA["theorem_scale"])}
    mask = 0
    bit = 0
    for i, j in combinations(range(len(keys)), 2):
        if order[keys[i]] < order[keys[j]]:
            mask |= 1 << bit
        bit += 1
    return mask


def directed_cycles(wins: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    cycles = []
    keys = [c.key for c in CARRIERS]
    for a, b, c in combinations(keys, 3):
        if b in wins[a] and c in wins[b] and a in wins[c]:
            cycles.append((a, b, c))
        if c in wins[a] and b in wins[c] and a in wins[b]:
            cycles.append((a, c, b))
    return cycles


def sccs(wins: dict[str, set[str]]) -> list[set[str]]:
    keys = [c.key for c in CARRIERS]

    def reach(start: str) -> set[str]:
        seen = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nxt in wins[cur]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen

    unseen = set(keys)
    comps = []
    for key in keys:
        if key not in unseen:
            continue
        fwd = reach(key)
        comp = {other for other in fwd if key in reach(other)}
        comps.append(comp)
        unseen -= comp
    return comps


def hamiltonian_paths(wins: dict[str, set[str]]) -> tuple[int, tuple[str, ...] | None]:
    keys = [c.key for c in CARRIERS]
    count = 0
    first = None
    for perm in permutations(keys):
        if all(perm[i + 1] in wins[perm[i]] for i in range(len(perm) - 1)):
            count += 1
            if first is None:
                first = perm
    return count, first


def name_map() -> dict[str, str]:
    return {c.key: c.name for c in CARRIERS}


def print_tournament_analysis() -> None:
    print()
    print("[4] Tournament Analysis on operator/carrier vertices")
    names = name_map()
    conservative = conservative_tournament()
    majority, wins, details = majority_tournament()
    cfp = tournament_fingerprint(conservative, len(CARRIERS))
    mfp = tournament_fingerprint(majority, len(CARRIERS))
    cycles = directed_cycles(wins)
    hp_count, first_hp = hamiltonian_paths(wins)
    score_hist = Counter(len(v) for v in wins.values())
    print("  vertices:")
    for carrier in CARRIERS:
        print(f"    {carrier.key:5s} {carrier.name:21s} -- {carrier.role}")
    print()
    print("  conservative gauge:")
    print("    observable: theorem scale / proof-router priority")
    print(
        f"    fingerprint: score_hist={cfp['score_hist']} c3={cfp['c3']} "
        f"scc={cfp['scc']} hp={cfp['hp']}"
    )
    print("    path: " + " > ".join(names[key] for key in CRITERIA["theorem_scale"]))
    print()
    print("  majority-of-criteria gauge:")
    print("    criteria: " + ", ".join(CRITERIA))
    print(
        f"    fingerprint: score_hist={dict(sorted(score_hist.items()))} "
        f"c3={mfp['c3']} scc={tuple(sorted((len(c) for c in sccs(wins)), reverse=True))} "
        f"hp={hp_count}"
    )
    if first_hp:
        print("    first Hamiltonian path: " + " > ".join(names[key] for key in first_hp))
    print(f"    directed 3-cycles: {len(cycles)}")
    for cyc in cycles[:6]:
        print("      " + " > ".join(names[key] for key in cyc) + f" > {names[cyc[0]]}")
    print("    nontrivial SCCs:")
    for comp in sccs(wins):
        if len(comp) > 1:
            print("      " + ", ".join(names[key] for key in sorted(comp)))
    print()
    print("  Edge-flip lesson:")
    print("    The conservative proof order is clean, but the majority gauge creates")
    print("    a real SCC among p+q, p*q, octahedron, and Clebsch/half-cube.")
    print("    That is useful: it marks the packet layer where no single scalar")
    print("    carrier should be trusted without its side-channel labels.")


def print_proof_readout() -> None:
    print()
    print("[5] Proof readout")
    print("  S137 keeps the current best interface and adds one warning.")
    print("  Interface:")
    print("    exact M/Farey branch first; C=27 shell transfer for p=2; Kpq/K33")
    print("    incidence for p>=3; octahedral/Clebsch packet data if a support-six")
    print("    or folded-mask residual remains; PZ only as an existence gateway.")
    print("  Warning:")
    print("    the additive lane and product/packet lanes form a majority-cycle once")
    print("    sequence simplicity, graph fit, state-lift fit, and magnitude stress")
    print("    are all allowed to vote.  This is exactly where scalarization is risky.")
    print("  Candidate lemma shape:")
    print("    after the standard LRC14 finite-core reductions, every low-gap")
    print("    non-AP/GW atom has either a unit-visible C27 hole discharged by")
    print("    petal/two-block rigidity, or a sign-visible K33/octahedral/Clebsch")
    print("    packet whose connected OCF state-lift would have activity-two value 7.")
    print("  This is still not an LRC14 proof.  It narrows the object that a proof")
    print("  must construct before invoking the HYP-2908 / THM-572 forbidden-H end.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain-limit", type=int, default=9)
    parser.add_argument("--max-replacement", type=int, default=70)
    args = parser.parse_args()

    print("S137 LRC14 OPERATOR SEQUENCE TOURNAMENT")
    print("=" * 78)
    print_assumption_challenge()
    print_sequence_table(args.chain_limit)
    print_row_bank(args.max_replacement)
    print_graph_and_pz()
    print_tournament_analysis()
    print_proof_readout()


if __name__ == "__main__":
    main()
