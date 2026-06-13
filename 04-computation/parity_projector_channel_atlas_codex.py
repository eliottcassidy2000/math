"""Parity projector atlas for midpoint and tournament reversal gates.

This script turns the prompt's slogan into exact finite checks:

* scalar midpoint anti-symmetrization keeps odd offset channels;
* tournament reversal/complement symmetrization keeps even Walsh channels;
* marked observers are transported by reversal, so their sum/difference
  split into even/odd channels.

The output is intentionally a research note, not a benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb
from typing import Dict, Iterable, List, Sequence, Tuple


def midpoint_support(max_degree: int = 8) -> List[Tuple[int, List[int], List[int]]]:
    """Return z-exponents surviving pair-sum and pair-difference gates."""
    rows = []
    for d in range(max_degree + 1):
        pair_sum = [r for r in range(d + 1) if r % 2 == 0 and comb(d, r)]
        pair_diff = [r for r in range(d + 1) if r % 2 == 1 and comb(d, r)]
        rows.append((d, pair_sum, pair_diff))
    return rows


def edge_pairs(n: int) -> List[Tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def adjacency_from_bits(n: int, bits: int, pairs: Sequence[Tuple[int, int]]) -> List[List[bool]]:
    adj = [[False] * n for _ in range(n)]
    for idx, (i, j) in enumerate(pairs):
        if (bits >> idx) & 1:
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def ham_counts(n: int, bits: int, pairs: Sequence[Tuple[int, int]]) -> Tuple[int, int, int]:
    """Return (all Hamiltonian paths, paths starting at 0, paths ending at 0)."""
    adj = adjacency_from_bits(n, bits, pairs)
    size = 1 << n

    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            row = adj[last]
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and row[nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    full = size - 1
    total = sum(dp[full])
    end0 = dp[full][0]

    dp0 = [[0] * n for _ in range(size)]
    dp0[1][0] = 1
    for mask in range(size):
        for last in range(n):
            val = dp0[mask][last]
            if not val:
                continue
            row = adj[last]
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and row[nxt]:
                    dp0[mask | (1 << nxt)][nxt] += val
    start0 = sum(dp0[full])
    return total, start0, end0


def c3_count(n: int, bits: int, pairs: Sequence[Tuple[int, int]]) -> int:
    adj = adjacency_from_bits(n, bits, pairs)
    total = 0
    for tri in combinations(range(n), 3):
        out = []
        for a in tri:
            out.append(sum(1 for b in tri if b != a and adj[a][b]))
        if sorted(out) == [1, 1, 1]:
            total += 1
    return total


def writhe(bits: int, m: int) -> int:
    return sum(1 if (bits >> idx) & 1 else -1 for idx in range(m))


def fwht(values: Sequence[int]) -> List[int]:
    out = list(values)
    h = 1
    n = len(out)
    while h < n:
        for i in range(0, n, 2 * h):
            for j in range(i, i + h):
                x = out[j]
                y = out[j + h]
                out[j] = x + y
                out[j + h] = x - y
        h *= 2
    return out


def walsh_degrees(values: Sequence[int], m: int) -> Dict[int, Dict[str, int]]:
    coeffs = fwht(values)
    summary: Dict[int, Dict[str, int]] = {}
    for mask, coeff in enumerate(coeffs):
        if coeff == 0:
            continue
        degree = mask.bit_count()
        bucket = summary.setdefault(degree, {"count": 0, "max_abs": 0, "energy": 0})
        bucket["count"] += 1
        bucket["max_abs"] = max(bucket["max_abs"], abs(coeff))
        bucket["energy"] += coeff * coeff
    return summary


def degree_list(summary: Dict[int, Dict[str, int]]) -> List[int]:
    return sorted(summary)


def parity_label(degrees: Iterable[int]) -> str:
    degs = list(degrees)
    has_even = any(d % 2 == 0 for d in degs)
    has_odd = any(d % 2 == 1 for d in degs)
    if has_even and has_odd:
        return "mixed"
    if has_even:
        return "even"
    if has_odd:
        return "odd"
    return "zero"


def tournament_walsh_audit(n_values: Sequence[int] = (3, 4, 5)) -> List[Dict[str, object]]:
    audits = []
    for n in n_values:
        pairs = edge_pairs(n)
        m = len(pairs)
        total_states = 1 << m
        mask_all = total_states - 1
        arrays = {
            "H": [],
            "c3": [],
            "writhe": [],
            "start0": [],
            "end0": [],
        }
        for bits in range(total_states):
            h_all, h_start0, h_end0 = ham_counts(n, bits, pairs)
            arrays["H"].append(h_all)
            arrays["c3"].append(c3_count(n, bits, pairs))
            arrays["writhe"].append(writhe(bits, m))
            arrays["start0"].append(h_start0)
            arrays["end0"].append(h_end0)

        arrays["start0_plus_end0"] = [
            arrays["start0"][i] + arrays["end0"][i] for i in range(total_states)
        ]
        arrays["start0_minus_end0"] = [
            arrays["start0"][i] - arrays["end0"][i] for i in range(total_states)
        ]
        arrays["edge0_flip_delta_H"] = [
            arrays["H"][i ^ 1] - arrays["H"][i] for i in range(total_states)
        ]
        arrays["edge0_oriented_grad_H"] = [
            (1 if (i & 1) else -1) * (arrays["H"][i ^ 1] - arrays["H"][i])
            for i in range(total_states)
        ]

        checks = {
            "H_complement_invariant": 0,
            "c3_complement_invariant": 0,
            "writhe_complement_anti": 0,
            "start0_to_end0_transport": 0,
            "edge0_delta_invariant": 0,
            "edge0_grad_anti": 0,
        }
        for bits in range(total_states):
            comp = mask_all ^ bits
            if arrays["H"][bits] != arrays["H"][comp]:
                checks["H_complement_invariant"] += 1
            if arrays["c3"][bits] != arrays["c3"][comp]:
                checks["c3_complement_invariant"] += 1
            if arrays["writhe"][bits] != -arrays["writhe"][comp]:
                checks["writhe_complement_anti"] += 1
            if arrays["start0"][bits] != arrays["end0"][comp]:
                checks["start0_to_end0_transport"] += 1
            if arrays["edge0_flip_delta_H"][bits] != arrays["edge0_flip_delta_H"][comp]:
                checks["edge0_delta_invariant"] += 1
            if arrays["edge0_oriented_grad_H"][bits] != -arrays["edge0_oriented_grad_H"][comp]:
                checks["edge0_grad_anti"] += 1

        walsh = {}
        for name, values in arrays.items():
            summary = walsh_degrees(values, m)
            walsh[name] = {
                "degrees": degree_list(summary),
                "parity": parity_label(summary),
                "nonzero_coeffs": sum(bucket["count"] for bucket in summary.values()),
                "max_abs_coeff": max((bucket["max_abs"] for bucket in summary.values()), default=0),
            }

        audits.append(
            {
                "n": n,
                "edges": m,
                "states": total_states,
                "checks": checks,
                "walsh": walsh,
            }
        )
    return audits


@dataclass(frozen=True)
class Carrier:
    name: str
    exactness: int
    parity_clarity: int
    lrc_leverage: int
    lift_readiness: int
    algorithmic_value: int
    cross_domain: int
    note: str

    def score_tuple(self) -> Tuple[int, ...]:
        return (
            self.lrc_leverage,
            self.lift_readiness,
            self.parity_clarity,
            self.exactness,
            self.algorithmic_value,
            self.cross_domain,
        )


def carrier_tournament(carriers: Sequence[Carrier]) -> Dict[str, object]:
    n = len(carriers)
    adj = [[False] * n for _ in range(n)]
    flips_vs_list_order = 0
    for i in range(n):
        for j in range(i + 1, n):
            left = carriers[i].score_tuple()
            right = carriers[j].score_tuple()
            if left >= right:
                adj[i][j] = True
            else:
                adj[j][i] = True
                flips_vs_list_order += 1

    scores = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]
    score_hist: Dict[int, int] = {}
    for s in scores:
        score_hist[s] = score_hist.get(s, 0) + 1

    c3 = 0
    for a, b, c in combinations(range(n), 3):
        edges = int(adj[a][b]) + int(adj[b][c]) + int(adj[c][a])
        if edges in (1, 2):
            # Cyclic iff each vertex has one out-edge in the triple.
            out = [
                int(adj[a][b]) + int(adj[a][c]),
                int(adj[b][a]) + int(adj[b][c]),
                int(adj[c][a]) + int(adj[c][b]),
            ]
            if out == [1, 1, 1]:
                c3 += 1

    # Hamiltonian paths in the carrier tournament.
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1 == 0 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    h_paths = sum(dp[size - 1])

    # Kosaraju SCCs.
    seen = [False] * n
    order: List[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in range(n):
            if adj[v][w] and not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    scc_sizes: List[int] = []

    def rdfs(v: int) -> int:
        seen[v] = True
        total = 1
        for w in range(n):
            if radj[v][w] and not seen[w]:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if not seen[v]:
            scc_sizes.append(rdfs(v))

    leaders = sorted(range(n), key=lambda idx: (-scores[idx], idx))
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted(scc_sizes, reverse=True),
        "hamiltonian_paths": h_paths,
        "flips_vs_list_order": flips_vs_list_order,
        "leaders": [(carriers[i].name, scores[i]) for i in leaders[:4]],
    }


def print_section(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def main() -> None:
    print("Parity projector channel atlas")
    print("Core rule: for an involution J, (f+fJ)/2 keeps invariant channels and (f-fJ)/2 keeps anti-invariant channels.")
    print("On a Boolean edge cube, global reversal x -> -x makes invariant functions even-Walsh and anti-invariant functions odd-Walsh.")

    print_section("1. Scalar midpoint gate")
    print("For a paired offset z around midpoint c:")
    for degree, pair_sum, pair_diff in midpoint_support(8):
        print(
            f"degree {degree}: pair-sum z-powers {pair_sum}; "
            f"pair-difference z-powers {pair_diff}"
        )
    print("Faulhaber anchor warning: the prompt's interval balance is paired offsets plus one fixed center term c^p.")
    print("Thus the surviving scalar data are odd moment channels plus a fixed-point scalar atom.")

    print_section("2. Exact tournament Walsh audit")
    audits = tournament_walsh_audit()
    names = [
        "H",
        "c3",
        "writhe",
        "start0",
        "start0_plus_end0",
        "start0_minus_end0",
        "edge0_flip_delta_H",
        "edge0_oriented_grad_H",
    ]
    for audit in audits:
        print(f"n={audit['n']} edges={audit['edges']} states={audit['states']}")
        print("  complement/transport mismatches:", audit["checks"])
        for name in names:
            row = audit["walsh"][name]
            print(
                f"  {name:24s} parity={row['parity']:5s} "
                f"degrees={row['degrees']} nonzero={row['nonzero_coeffs']} "
                f"max|W|={row['max_abs_coeff']}"
            )

    print_section("3. Cross-domain parity gate atlas")
    atlas_rows = [
        (
            "Faulhaber midpoint",
            "offset reversal z -> -z",
            "anti gate plus fixed center",
            "odd moments S_1,S_3,... and c^p",
            "build LRC wall ledgers as odd resources plus fixed atoms",
        ),
        (
            "Hamiltonian path count H",
            "tournament converse T -> T^op",
            "invariant gate",
            "even Walsh levels only",
            "scalar H cannot see odd owner/writhe channels",
        ),
        (
            "Rooted tournament perspective",
            "start observer <-> end observer",
            "transport gate",
            "start+end even, start-end odd",
            "observer-coupled LRC data should be paired before quotienting",
        ),
        (
            "H-gradient / flip sensitivity",
            "global converse plus marked edge",
            "derivative toggles parity",
            "raw flip delta even, oriented gradient odd",
            "odd side channels are where local obstruction pressure lives",
        ),
        (
            "Signed LRC cut gauge",
            "speed sign reversal",
            "observer-blind, pair-visible",
            "distance scalar even; sum/difference clocks are marked cuts",
            "ignore signs for one-runner distance, keep them for pair ledgers",
        ),
        (
            "Unit distance threshold",
            "distance equality vs oriented tiling",
            "symmetric metric gate plus directed support",
            "unit graph even; flipped/nonunit orientation is side data",
            "unit spine questions need support orientation, not only counts",
        ),
        (
            "Type II / code72 scalar gate",
            "Gleason invariant reversal",
            "even weight enumerator gate",
            "scalar enumerator positive; support design remains",
            "72 obstruction is support/matroid, not parity scalar",
        ),
        (
            "Polynomial sign cube",
            "global sign/reciprocal choices",
            "scalar chamber gate",
            "irreducibility needs residue and convolution lifts",
            "sign parity is a quotient; factor grids are the hidden carrier",
        ),
        (
            "OCF odd-cycle formula",
            "cycle-packet compatibility",
            "odd atoms with product packets",
            "alpha_k in I(Omega,2)",
            "model odd Faulhaber/LRC resources as compatible packets",
        ),
    ]
    for row in atlas_rows:
        print(" | ".join(row))

    print_section("4. Proof-carrier Tournament Analysis")
    carriers = [
        Carrier("marked_section_parity_toggle", 5, 5, 5, 4, 5, 5, "start/end and plus/minus split"),
        Carrier("lrc_sign_cut_side_channel", 4, 5, 5, 4, 5, 4, "observer-blind but pair-visible"),
        Carrier("faulhaber_odd_fixed_atom_gate", 5, 5, 4, 4, 4, 5, "odd moments plus center atom"),
        Carrier("h_gradient_derivative_gate", 5, 5, 4, 4, 5, 4, "derivative toggles even to odd"),
        Carrier("ocf_compatibility_packets", 5, 4, 4, 5, 4, 5, "odd atoms need alpha packets"),
        Carrier("code72_support_design_lift", 4, 4, 3, 5, 4, 5, "even scalar passes, support open"),
        Carrier("polynomial_convolution_lift", 5, 3, 3, 5, 5, 5, "scalar sign cube needs factor grid"),
        Carrier("unit_distance_spine_lift", 4, 4, 3, 4, 4, 4, "metric equality plus oriented support"),
        Carrier("tournament_even_walsh_gate", 5, 5, 3, 3, 5, 4, "H and c3 live on even Walsh levels"),
        Carrier("lrc14_q27_owner_carry_ledger", 4, 4, 5, 5, 4, 4, "separate scalar clocks from owner/carry channels"),
    ]
    fp = carrier_tournament(carriers)
    print("score_hist:", fp["score_hist"])
    print("directed_3cycles:", fp["directed_3cycles"])
    print("scc_sizes:", fp["scc_sizes"])
    print("hamiltonian_paths:", fp["hamiltonian_paths"])
    print("flips_vs_list_order:", fp["flips_vs_list_order"])
    print("leaders:", fp["leaders"])

    print_section("5. Working moves")
    moves = [
        "Do not quotient marked observer data until it is split into start+end and start-end analogues.",
        "For LRC14, tag every clock as invariant scalar, anti/owner channel, or transporter between perspectives.",
        "Use even scalar gates to shrink search; use odd derivative/owner channels to prove pressure or descent.",
        "Treat HYP-2458 odd atoms as incomplete until an OCF-style compatibility packet ledger is attached.",
        "When a scalar gate passes in code72/unit-distance/polynomials, immediately ask which support lift was forgotten.",
    ]
    for idx, move in enumerate(moves, 1):
        print(f"{idx}. {move}")


if __name__ == "__main__":
    main()
