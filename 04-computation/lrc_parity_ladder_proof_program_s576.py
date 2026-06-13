#!/usr/bin/env python3
"""
lrc_parity_ladder_proof_program_s576.py

codex-2026-06-03 S576

Use HYP-2091 as a proof-routing table for LRC.

HYP-2091 separates the runner tournament geometry into two parity ladders:
  * even LRC n: m=n-1 is odd, so the regular polygon outside is a clean
    rotational tournament;
  * odd LRC n: m=n-1 is even, so the regular polygon outside is an antipodal
    wall and the simplex mesh must resolve tie diagonals.

This script overlays that parity split with the clocks that HYP-2083/2088/2090
say matter: D_q divisibility clocks, U_a unit-shell clocks modulo 2n-1, N_j
floor clocks, converse-merged round nodes, nonunit summand strata, and endpoint
owner labels.

Tournament Analysis / assumption challenge:
  Vertices are LRC n-ladder rows, not runners.
  Pairwise observable is the proof-burden vector
      (wall tie pairs, nonunit shell count, D/U/N obligations,
       converse-merged round nodes, extended dihedral size).
  Switch/gauge orients toward the larger remaining burden; ties follow
  increasing n as the Hamiltonian path.
  This quotient preserves which proof lemma family is needed.  It destroys
  exact speed ownership and exact active pair-sum geometry, which must re-enter
  through D/U/N private labels and THM-397 endpoint blockers.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations
from math import factorial, gcd, log10


@dataclass(frozen=True)
class Row:
    n: int
    m: int
    c: int
    factors: str
    lane: str
    tie_pairs: int
    tie_choices: int
    d_obligations: int
    u_obligations: int
    n_obligations: int
    unit_shells: int
    nonunit_shells: int
    shell_strata: tuple[tuple[int, int], ...]
    round_classes: int
    self_converse_round: int
    converse_merged_round: int
    d_ext: int

    @property
    def total_obligations(self) -> int:
        return self.d_obligations + self.u_obligations + self.n_obligations

    @property
    def burden(self) -> tuple[int, int, int, int, int]:
        return (
            self.tie_pairs,
            self.nonunit_shells,
            self.total_obligations,
            self.converse_merged_round,
            self.d_ext,
        )


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def totient(n: int) -> int:
    return sum(1 for a in range(1, n + 1) if gcd(a, n) == 1)


def factorization(n: int) -> str:
    parts: list[str] = []
    d = 2
    x = n
    while d * d <= x:
        if x % d:
            d += 1
            continue
        e = 0
        while x % d == 0:
            x //= d
            e += 1
        parts.append(str(d) if e == 1 else f"{d}^{e}")
        d += 1
    if x > 1:
        parts.append(str(x))
    return "*".join(parts)


def round_classes(m: int) -> int:
    # A000016.  The valid labelled round words have a binary necklace quotient.
    return sum(totient(d) * 2 ** (m // d) for d in divisors(m) if d % 2 == 1) // (2 * m)


def self_converse_round_count(m: int) -> int:
    # HYP-2089/S575 audit: 2^floor((m-1)/2) through m=13.
    return 2 ** ((m - 1) // 2)


def shell_strata(c: int) -> tuple[tuple[int, int], ...]:
    # Antipodal shells {a,c-a}, 1 <= a <= (c-1)/2, grouped by gcd(a,c).
    k = (c - 1) // 2
    counts = Counter(gcd(a, c) for a in range(1, k + 1))
    return tuple(sorted(counts.items()))


def row_for_n(n: int) -> Row:
    m = n - 1
    c = 2 * n - 1
    strata = shell_strata(c)
    unit_shells = dict(strata).get(1, 0)
    total_shells = n - 1
    sc = self_converse_round_count(m)
    round_count = round_classes(m)
    return Row(
        n=n,
        m=m,
        c=c,
        factors=factorization(c),
        lane="clean_polygon" if m % 2 else "wall_mesh",
        tie_pairs=0 if m % 2 else m // 2,
        tie_choices=1 if m % 2 else 2 ** (m // 2),
        d_obligations=n - 2,
        u_obligations=totient(c) // 2,
        n_obligations=n - 1,
        unit_shells=unit_shells,
        nonunit_shells=total_shells - unit_shells,
        shell_strata=strata,
        round_classes=round_count,
        self_converse_round=sc,
        converse_merged_round=(round_count + sc) // 2,
        d_ext=2 * m,
    )


def burden_route(row: Row) -> str:
    if row.tie_pairs == 0 and row.nonunit_shells == 0:
        return "clean unit lane: round/converse seam + DUN exchange"
    if row.tie_pairs == 0:
        return "clean composite lane: add gcd descent to DUN exchange"
    if row.nonunit_shells == 0:
        return "wall lane: tie-discharge before unit-shell exchange"
    return "wall composite lane: tie-discharge plus gcd descent"


def fmt_strata(row: Row) -> str:
    return ",".join(f"g{g}:{count}" for g, count in row.shell_strata)


def log_or_int(x: int) -> str:
    if x < 100000:
        return str(x)
    return f"10^{log10(x):.2f}"


def tournament_fingerprint(rows: list[Row]) -> dict[str, object]:
    n = len(rows)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(rows):
        for j, b in enumerate(rows):
            if i == j:
                continue
            adj[i][j] = (a.burden, a.n) > (b.burden, b.n)

    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        ):
            c3 += 1

    def reach(start: int, graph: list[list[bool]]) -> set[int]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for w, edge in enumerate(graph[v]):
                if edge and w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        start = next(iter(remaining))
        forward = reach(start, adj)
        comp = {v for v in remaining if v in forward and start in reach(v, adj)}
        sccs.append(len(comp))
        remaining -= comp

    return {
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_path": [f"n={row.n}" for row in sorted(rows, key=lambda r: (r.burden, r.n), reverse=True)],
    }


def main() -> None:
    rows = [row_for_n(n) for n in range(4, 19)]

    print("S576 LRC parity-ladder proof program")
    print("=" * 78)
    print("HYP-2091 supplies the geometry split; HYP-2083/2088/2090 supply")
    print("the clock obligations that survive that split.")
    print()

    print(
        "n  lane           C       factors  ties choices  U/nonU  DUN  merged_round  route"
    )
    for row in rows:
        print(
            f"{row.n:2d} {row.lane:14s} {row.c:5d} {row.factors:8s} "
            f"{row.tie_pairs:5d} {log_or_int(row.tie_choices):>7s} "
            f"{row.unit_shells:2d}/{row.nonunit_shells:<2d} "
            f"{row.total_obligations:3d} {row.converse_merged_round:12d}  "
            f"{burden_route(row)}"
        )

    print("\nShell strata by the 2n-1 summand/unit clock")
    for row in rows:
        print(f"  n={row.n:2d} C={row.c:<3d} {row.factors:<8s} strata={fmt_strata(row)}")

    focus = row_for_n(14)
    print("\nFourteen-runner focus")
    print(f"  runner tournament size m={focus.m}; HYP-2091 lane={focus.lane}")
    print(f"  C=2n-1={focus.c}={focus.factors}; shell strata {fmt_strata(focus)}")
    print(
        "  round body: "
        f"{focus.round_classes} round classes, {focus.self_converse_round} self-converse, "
        f"{focus.converse_merged_round} converse-merged nodes"
    )
    print(
        "  D/U/N obligations: "
        f"D={focus.d_obligations}, U={focus.u_obligations}, "
        f"N={focus.n_obligations}, total={focus.total_obligations}"
    )
    print("  interpretation: clean geometry, composite clock burden.")
    print("  The hard n=14 work is not antipodal tie resolution; it is the")
    print("  gcd-3/gcd-9 nonunit shell descent plus owner-labelled D/U/N pivots.")

    print("\nClocks that matter")
    print("  pair-sum clocks: exact 1D maximin candidates and THM-397 endpoint owners")
    print("  n-clock: floor/tight equality and observer-source escape attempts")
    print("  2n-1 unit clock: antipodal summand witnesses visible to multiplication")
    print("  D clocks: small denominator divisibility witnesses")
    print("  labelled event clocks: runner owner, endpoint owner, pair-sum denominator")

    print("\nClocks or quotients to ignore until labels return")
    print("  primitive reset length: normalized integer rows all close at period 1")
    print("  binary cyclic lonely-word stabilizer: usually trivial in S571/S573")
    print("  unlabelled round class alone: loses D/U/N and endpoint ownership")
    print("  runner-vertex tournaments alone: lose proof-obligation incidence")

    print("\nLemma queue")
    print("  L1 clean seam: attach D/U/N labels to converse-merged round nodes.")
    print("  L2 clean exchange: full D/U/N cover has floor or pair-sum witness.")
    print("  L3 nonunit descent: gcd-stratum defects lift to a second clock or lower row.")
    print("  L4 wall discharge: antipodal tie resolutions either create a source gap")
    print("     or reduce to a neighboring clean-ladder labelled cover.")
    print("  L5 endpoint compatibility: THM-397 endpoint blockers are exactly the")
    print("     labelled owners that stop pair-pinch escape clocks.")

    print("\nTournament Analysis")
    print("  vertices: LRC n-ladder rows")
    print("  observable: (tie pairs, nonunit shells, D/U/N total, merged round nodes, D_ext)")
    print("  switch: orient toward larger remaining proof burden")
    print(f"  fingerprints: {tournament_fingerprint(rows)}")


if __name__ == "__main__":
    main()
