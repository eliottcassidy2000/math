#!/usr/bin/env python3
"""
tournament_gap_transfer_exponent_s615.py

S615 scout for the user's "H gaps propagate" thesis.

The purpose is not to prove a new forbidden Hamiltonian-path value.  The repo
canon is stricter: only H=7 and H=21 are currently proved permanent gaps.  This
script formalizes a transfer lens around those gaps, the unit-distance n^1.014
construction, LRC residue/carry coimages, and Collatz-style correlated residue
walls.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class Gap:
    h: int
    theorem: str
    obstruction: str
    transfer_reading: str


GAPS = (
    Gap(
        7,
        "THM-029",
        "exactly three mutually conflicting odd-cycle packets cannot remain "
        "closed; tournament completion forces extra odd-cycle structure",
        "small forbidden evaluation becomes a certificate that a quotient has "
        "forgotten a forced completion channel",
    ),
    Gap(
        21,
        "THM-079",
        "all OCF decompositions of I(Omega,2)=21 are blocked by tournament "
        "forcing and base-case discharge",
        "larger forbidden evaluation becomes a multi-packet obstruction, not "
        "a Mersenne-pattern continuation",
    ),
)


@dataclass(frozen=True)
class Domain:
    name: str
    carrier: str
    observable: str
    hard_wall: str
    retained_side_channel: str


DOMAINS = (
    Domain(
        "tournament H-spectrum",
        "odd-cycle conflict graph Omega(T)",
        "H(T)=I(Omega(T),2)",
        "H=7 and H=21 are permanent forbidden evaluations",
        "which compatible odd-cycle packets completion forces",
    ),
    Domain(
        "LRC floor problem",
        "depth nerve / relation lattice / Res_(2n-1) coimage",
        "p0 and maximin radius M(S)",
        "p0=0 collapse is all-orders cancellation; finite moments are blind",
        "owner labels, carry cocycles, pinch witnesses, endpoint certificates",
    ),
    Domain(
        "unit-distance grids",
        "geometric graph plus Moser/CM arithmetic carrier",
        "number of unit-distance pairs",
        "graph-only or visible-grid quotients miss embeddability and tower data",
        "coordinates, unfaithful subgraphs, dense deletion cores, split primes",
    ),
    Domain(
        "Collatz/two-block residues",
        "carry/residue automaton over the log2-log3 near-lattice",
        "cycle or empty-residual certificate",
        "density of residue classes is blind to correlated carry constraints",
        "carry depth, determinant walls, CRT compatibility, valuation history",
    ),
)


@dataclass(frozen=True)
class Mechanism:
    name: str
    proof_status: int
    quotient_preservation: int
    transfer_bandwidth: int
    exponent_contact: int
    speculation_risk: int
    note: str

    @property
    def score(self) -> int:
        return (
            3 * self.proof_status
            + 2 * self.quotient_preservation
            + self.transfer_bandwidth
            + self.exponent_contact
            - 2 * self.speculation_risk
        )


MECHANISMS = (
    Mechanism(
        "H-gap completion forcing",
        5,
        4,
        5,
        2,
        1,
        "proved tournament side; explains why a tiny arithmetic evaluation can "
        "be impossible",
    ),
    Mechanism(
        "coimage side-channel retention",
        4,
        5,
        5,
        4,
        1,
        "shared with LRC Res_27, unit-distance n=22, and CM projection carriers",
    ),
    Mechanism(
        "two-block determinant walls",
        3,
        4,
        4,
        2,
        2,
        "local CRT incompatibility version of an unavailable global state",
    ),
    Mechanism(
        "exponent entropy normalization",
        2,
        4,
        5,
        5,
        3,
        "the right place to test the user's 1.014 tournament/unit-distance match",
    ),
    Mechanism(
        "Collatz carry-residue transfer",
        2,
        3,
        4,
        2,
        3,
        "promising analogy, but needs a sharper conserved observable",
    ),
    Mechanism(
        "raw scalar H matching",
        1,
        1,
        2,
        1,
        5,
        "dangerous: the repo already warns against promoting sampled H gaps",
    ),
)


@dataclass(frozen=True)
class ExponentNormalization:
    source: str
    value: str
    variable: str
    carrier: str
    status: str


EXPONENTS = (
    ExponentNormalization(
        "Sawin unit-distance lower bound",
        "1.014",
        "number n of planar points",
        "CM/class-field tower projected to the plane",
        "primary-source theorem statement; unit side is normalized",
    ),
    ExponentNormalization(
        "raw A000568 tournament growth",
        "not 1.014",
        "number of tournament vertices",
        "unlabelled tournament Burnside quotient",
        "wrong normalization for the user's claim; identity term dominates",
    ),
    ExponentNormalization(
        "tournament forbidden-carrier entropy",
        "unmeasured",
        "growth of feasible OCF/obligation carriers after forbidden packets",
        "Omega(T), proof obligations, residue/gap carriers",
        "candidate place where a 1.014-style exponent could genuinely live",
    ),
)


def ocf_arithmetic_decompositions(h: int, max_degree: int = 4) -> list[tuple[int, ...]]:
    """Return nonnegative alpha-tuples with 1+sum 2^k alpha_k = h."""
    target = h - 1
    weights = [2**k for k in range(1, max_degree + 1)]
    out: list[tuple[int, ...]] = []

    def rec(i: int, remaining: int, prefix: list[int]) -> None:
        if i == len(weights):
            if remaining == 0:
                out.append(tuple(prefix))
            return
        weight = weights[i]
        for count in range(remaining // weight + 1):
            prefix.append(count)
            rec(i + 1, remaining - count * weight, prefix)
            prefix.pop()

    rec(0, target, [])
    return sorted(out, key=lambda row: (sum(row), row))


def orient_tournament(items: tuple[Mechanism, ...]) -> list[list[int]]:
    n = len(items)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a, b = items[i], items[j]
        key_a = (a.score, -a.speculation_risk, a.proof_status, a.name)
        key_b = (b.score, -b.speculation_risk, b.proof_status, b.name)
        if key_a >= key_b:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        cyc = (
            adj[a][b] and adj[b][c] and adj[c][a]
        ) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        )
        total += int(cyc)
    return total


def strongly_connected_components(adj: list[list[int]]) -> list[list[int]]:
    n = len(adj)

    def reach_from(start: int, graph: list[list[int]]) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for u, edge in enumerate(graph[v]):
                if edge and u not in seen:
                    seen.add(u)
                    stack.append(u)
        return seen

    radj = [[adj[j][i] for j in range(n)] for i in range(n)]
    remaining = set(range(n))
    comps: list[list[int]] = []
    while remaining:
        v = min(remaining)
        comp = reach_from(v, adj) & reach_from(v, radj)
        comps.append(sorted(comp))
        remaining -= comp
    return comps


def hamiltonian_path_count(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1 << n) - 1])


def print_section(title: str) -> None:
    print()
    print("=" * 72)
    print(title)
    print("=" * 72)


def main() -> None:
    print("S615: tournament forbidden-H transfer and exponent normalization")
    print("Assumption challenge: vertices are transfer mechanisms/proof obligations,")
    print("not runners, arcs, unit-distance points, or Collatz residues by default.")

    print_section("Known permanent H gaps")
    for gap in GAPS:
        print(f"H={gap.h}: {gap.theorem}")
        print(f"  obstruction: {gap.obstruction}")
        print(f"  transfer:    {gap.transfer_reading}")
        decomps = ocf_arithmetic_decompositions(gap.h)
        print(
            f"  arithmetic alpha-tuples up to degree 4: {len(decomps)} "
            f"(first {min(8, len(decomps))}: {decomps[:8]})"
        )

    print()
    print("Guardrail: H=63 is known to be achieved at n=8; larger Mersenne-looking")
    print("values are not permanent gaps merely because 7 and 21 are forbidden.")

    print_section("Cross-domain carrier ledger")
    for domain in DOMAINS:
        print(f"{domain.name}:")
        print(f"  carrier:  {domain.carrier}")
        print(f"  observed: {domain.observable}")
        print(f"  wall:     {domain.hard_wall}")
        print(f"  retain:   {domain.retained_side_channel}")

    print_section("Tournament Analysis over mechanisms")
    adj = orient_tournament(MECHANISMS)
    outdegrees = [sum(row) for row in adj]
    score_hist: dict[int, int] = {}
    for d in outdegrees:
        score_hist[d] = score_hist.get(d, 0) + 1
    ranking = sorted(
        enumerate(MECHANISMS),
        key=lambda x: (outdegrees[x[0]], x[1].score, -x[1].speculation_risk),
        reverse=True,
    )
    for rank, (idx, mechanism) in enumerate(ranking, 1):
        print(
            f"{rank}. {mechanism.name}: out={outdegrees[idx]}, "
            f"score={mechanism.score}, risk={mechanism.speculation_risk}"
        )
        print(f"   {mechanism.note}")
    print(f"score histogram: {dict(sorted(score_hist.items()))}")
    print(f"directed 3-cycles: {directed_triangles(adj)}")
    print(f"SCCs: {strongly_connected_components(adj)}")
    print(f"Hamiltonian paths in mechanism tournament: {hamiltonian_path_count(adj)}")

    print_section("1.014 normalization audit")
    for item in EXPONENTS:
        print(f"{item.source}:")
        print(f"  value:    {item.value}")
        print(f"  variable: {item.variable}")
        print(f"  carrier:  {item.carrier}")
        print(f"  status:   {item.status}")

    print()
    print("Working target:")
    print("  If a tournament-side 1.014 exists, it should be measured as an")
    print("  entropy dividend/deficit of the feasible OCF or proof-obligation carrier,")
    print("  after forbidden H packets and side-channel constraints are retained.")
    print("  Raw scalar H counts and raw A000568 growth are the wrong coordinates.")

    print_section("Next proof probes")
    print("1. Complete-core H=21: prove the r=10 single-core signature gap or find")
    print("   a non-core route that explains why THM-079 needs its full decomposition.")
    print("2. Unit-distance transfer: map each killed graph-only state to the")
    print("   corresponding missing side-channel, not merely to an edge count.")
    print("3. LRC transfer: treat Res_(2n-1) floor atoms as forbidden-evaluation")
    print("   tests; a lift must preserve owner/carry/pinch probes Yoneda-style.")
    print("4. Exponent test: define one common carrier-size variable and measure the")
    print("   tournament entropy dividend against Sawin's n^1.014 unit-distance side.")


if __name__ == "__main__":
    main()
