"""Digraph H=7 guardrail for the LRC14 forbidden-clique route.

The owner's prompt points at two facts:

* tournament arcs are binary orientation choices;
* H=7 is forbidden for tournaments.

This script audits the exact place where those facts are powerful and the
place where they are not enough.  Ordinary directed graphs also have binary
ordered-arc present/absent states, and they *do* realize seven Hamiltonian
paths.  Therefore the LRC14 disproof obstruction must prove that the LRC
apex over-cover lands in the tournament/OCF conflict-graph image, not merely
in a generic binary digraph shadow.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp[mask][last]
            if cur == 0:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += cur
    return sum(dp[(1 << n) - 1])


def tournament_from_bits(n: int, bits: int) -> list[list[bool]]:
    adj = [[False] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits & (1 << k):
                adj[i][j] = True
            else:
                adj[j][i] = True
            k += 1
    return adj


def digraph_from_bits(n: int, bits: int) -> list[list[bool]]:
    arcs = [(i, j) for i in range(n) for j in range(n) if i != j]
    adj = [[False] * n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if bits & (1 << k):
            adj[i][j] = True
    return adj


def oriented_from_states(n: int, states: tuple[int, ...]) -> list[list[bool]]:
    adj = [[False] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            state = states[k]
            if state == 1:
                adj[i][j] = True
            elif state == 2:
                adj[j][i] = True
            k += 1
    return adj


def pair_state_profile(adj: list[list[bool]]) -> Counter[str]:
    profile: Counter[str] = Counter()
    n = len(adj)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = adj[i][j], adj[j][i]
            if a and b:
                profile["both"] += 1
            elif a:
                profile["one_forward"] += 1
            elif b:
                profile["one_backward"] += 1
            else:
                profile["none"] += 1
    return profile


def arc_list(adj: list[list[bool]]) -> list[tuple[int, int]]:
    return [(i, j) for i, row in enumerate(adj) for j, val in enumerate(row) if val]


def tournament_spectrum(n: int) -> Counter[int]:
    m = n * (n - 1) // 2
    hist: Counter[int] = Counter()
    for bits in range(1 << m):
        hist[hamiltonian_paths(tournament_from_bits(n, bits))] += 1
    return hist


def first_general_digraph_with_h(n: int, target: int) -> tuple[list[list[bool]], int] | None:
    m = n * (n - 1)
    for bits in range(1 << m):
        adj = digraph_from_bits(n, bits)
        if hamiltonian_paths(adj) == target:
            return adj, bits
    return None


def oriented_h_count(n: int, target: int) -> tuple[int, list[list[bool]] | None]:
    count = 0
    first = None
    m = n * (n - 1) // 2
    for states in product((0, 1, 2), repeat=m):
        adj = oriented_from_states(n, states)
        if hamiltonian_paths(adj) == target:
            count += 1
            if first is None:
                first = adj
    return count, first


def winding_tournament(E: list[int], x: Fraction) -> list[list[bool]]:
    n = len(E)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            phase = ((E[i] - E[j]) * x) % 1
            adj[i][j] = Fraction(0) < phase < Fraction(1, 2)
    return adj


def winding_breakpoints(E: list[int]) -> list[Fraction]:
    bps = {Fraction(0), Fraction(1)}
    diffs = {abs(a - b) for a in E for b in E if a != b}
    for d in diffs:
        for k in range(d + 1):
            bps.add(Fraction(k, d))
        for k in range(d):
            bps.add(Fraction(2 * k + 1, 2 * d))
    return sorted(x for x in bps if 0 <= x <= 1)


def winding_cell_spectrum(E: list[int]) -> Counter[int]:
    bps = winding_breakpoints(E)
    hist: Counter[int] = Counter()
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        hist[hamiltonian_paths(winding_tournament(E, mid))] += 1
    return hist


def independence_clique_value(r: int, x: int = 2) -> int:
    return 1 + r * x


def main() -> None:
    print("DIGRAPH H=7 GUARDRAIL FOR LRC14")
    print("================================\n")

    print("1. Tournament orientation bits: H=7 is absent in exact small spectra.")
    for n in range(3, 7):
        hist = tournament_spectrum(n)
        values = sorted(hist)
        print(
            f"   tournaments n={n}: count={sum(hist.values())}, "
            f"H=7 count={hist.get(7, 0)}, all H odd={all(h % 2 for h in hist)}, "
            f"values={values[:12]}{'...' if len(values) > 12 else ''}"
        )

    print("\n2. Ordinary binary ordered-arc digraphs: H=7 is realizable.")
    for n in range(1, 4):
        max_h = 1
        for k in range(2, n + 1):
            max_h *= k
        print(f"   n={n}: max possible Hamiltonian paths is {max_h}; H=7 impossible by size.")
    found = first_general_digraph_with_h(4, 7)
    assert found is not None
    adj, bits = found
    print(f"   first n=4 simple digraph with H=7: bitmask={bits}")
    print(f"   arcs={arc_list(adj)}")
    print(f"   unordered pair-state profile={dict(pair_state_profile(adj))}")
    print("   This is not a tournament: the binary states are ordered-arc present/absent,")
    print("   not one-of-two orientations for every unordered pair.")

    print("\n3. Oriented but incomplete graphs also show the danger of forgetting completeness.")
    count4, first4 = oriented_h_count(4, 7)
    print(f"   oriented graphs on n=4 with H=7: {count4}")
    if first4 is not None:
        print(f"   first oriented arcs={arc_list(first4)}")
        print(f"   profile={dict(pair_state_profile(first4))}")
    count5, first5 = oriented_h_count(5, 7)
    print(f"   oriented graphs on n=5 with H=7: {count5}")
    if first5 is not None:
        print(f"   first n=5 oriented arcs={arc_list(first5)}")
        print(f"   profile={dict(pair_state_profile(first5))}")
        print("   These have at most one orientation per pair, but missing pairs break")
        print("   the tournament completeness constraint that drives THM-200.")

    print("\n4. Tie-free winding-tournament cells for AP7 avoid H=7.")
    E = list(range(7))
    whist = winding_cell_spectrum(E)
    print(f"   E={E}, cells={sum(whist.values())}, H=7 cells={whist.get(7, 0)}")
    print(f"   all cell H odd={all(h % 2 for h in whist)}")
    print(f"   cell H spectrum={sorted(whist.items())}")
    h_lonely = hamiltonian_paths(winding_tournament(E, Fraction(1, 7)))
    print(f"   H at exact lonely point x=1/7 is {h_lonely}")

    print("\n5. Clique law is an abstract preimage, not a generic digraph obstruction.")
    for r in range(1, 6):
        mark = " <-- K3 gives 7" if r == 3 else ""
        print(f"   I(K_{r},2) = {independence_clique_value(r)}{mark}")

    print("\nCONCLUSION")
    print("   The prompt's two facts are powerful only with the tournament/OCF")
    print("   realizability hypothesis included.  Generic binary digraphs realize H=7,")
    print("   while tournaments do not.  Therefore the LRC14 forbidden-clique route")
    print("   must prove:")
    print("       apex-7 over-cover => tournament conflict graph Omega = K3")
    print("   (or an equivalent labelled OCF packet), not merely a binary digraph with")
    print("   seven path/configuration states.  This works toward impossible-to-disprove:")
    print("   once that correspondence is exact, THM-200 blocks the counterexample.")


if __name__ == "__main__":
    main()
