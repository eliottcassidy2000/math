#!/usr/bin/env python3
"""LRC/tournament maximizer signal atlas for HYP-3100.

The goal is not another single cap computation.  It is to measure which
currency makes a row look maximal or minimal:

* pair-Pascal value shadow;
* live higher-order inclusion-exclusion anatomy;
* exchange-local traps in the one-swap improvement graph;
* tournament rigidity as Hamiltonian-path spectrum;
* controlled-forgetting transfer risk between these observers.

Tournament Analysis is included twice.  First we enumerate small tournament
Hamiltonian-path spectra.  Second we treat signal families themselves as
tournament vertices and orient A -> B when A retains more proof-relevant
coordinates than B under a majority gauge.
"""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb


KPLUS = 14


def fstr(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def norm(x: Fraction) -> Fraction:
    f = x - int(x)
    if f < 0:
        f += 1
    return min(f, 1 - f)


def joint_forbidden_measure(speeds: tuple[int, ...], kk1: int = KPLUS) -> Fraction:
    """Measure of x with ||p*x|| < 1/kk1 for every p in speeds."""
    speeds = tuple(sorted(p for p in speeds if p))
    if not speeds:
        return Fraction(1)
    breaks = {Fraction(0), Fraction(1)}
    for p in speeds:
        for a in range(p + 1):
            for sign in (-1, 1):
                x = Fraction(a, p) + sign * Fraction(1, kk1 * p)
                if 0 <= x <= 1:
                    breaks.add(x)
    cuts = sorted(breaks)
    total = Fraction(0)
    for left, right in zip(cuts, cuts[1:]):
        if right <= left:
            continue
        mid = (left + right) / 2
        if all(norm(p * mid) < Fraction(1, kk1) for p in speeds):
            total += right - left
    return total


def lonely_measure(speeds: tuple[int, ...], kk1: int = KPLUS) -> Fraction:
    """Measure of x with ||p*x|| >= 1/kk1 for every p in speeds."""
    speeds = tuple(sorted(set(p for p in speeds if p)))
    if not speeds:
        return Fraction(1)
    breaks = {Fraction(0), Fraction(1)}
    for p in speeds:
        for a in range(p + 1):
            for sign in (-1, 1):
                x = Fraction(a, p) + sign * Fraction(1, kk1 * p)
                if 0 <= x <= 1:
                    breaks.add(x)
    cuts = sorted(breaks)
    total = Fraction(0)
    for left, right in zip(cuts, cuts[1:]):
        if right <= left:
            continue
        mid = (left + right) / 2
        if all(norm(p * mid) >= Fraction(1, kk1) for p in speeds):
            total += right - left
    return total


def lonely_orders(speeds: tuple[int, ...]) -> list[Fraction]:
    speeds = tuple(speeds)
    out: list[Fraction] = []
    for r in range(len(speeds) + 1):
        raw = Fraction(0)
        for subset in itertools.combinations(speeds, r):
            raw += joint_forbidden_measure(subset)
        out.append(((-1) ** r) * raw)
    return out


def max_linear_run(speeds: tuple[int, ...]) -> int:
    if not speeds:
        return 0
    vals = sorted(set(speeds))
    best = cur = 1
    for a, b in zip(vals, vals[1:]):
        if b == a + 1:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def gap_vector(speeds: tuple[int, ...]) -> tuple[int, ...]:
    vals = sorted(speeds)
    return tuple(b - a for a, b in zip(vals, vals[1:]))


def print_cap_anatomy() -> None:
    minimizers = {
        1: (1,),
        2: (1, 13),
        3: (1, 12, 13),
        4: (1, 11, 12, 13),
        5: (1, 5, 7, 8, 9),
    }
    print("=" * 96)
    print("LRC CAP MAXIMIZER/MINIMIZER ANATOMY: PAIR SHADOW VS LIVE HIGHER ORDERS")
    print("=" * 96)
    for j, speeds in minimizers.items():
        k = 13 - j
        orders = lonely_orders(speeds)
        true = sum(orders)
        pascal = Fraction(comb(14 - j, 2), 91)
        pair_target = Fraction(comb(j, 2), 91)
        o2 = orders[2] if j >= 2 else Fraction(0)
        pair_excess = o2 - pair_target
        high_order = sum(orders[3:])
        first_live = next((r for r, val in enumerate(orders[3:], start=3) if val), None)
        false_positive = "yes" if pascal == true and (pair_excess or high_order) else "no"
        print(f"\nj={j} k={k} P={speeds}")
        print(f"  lonely={fstr(true)}  pair_pascal={fstr(pascal)}  final_dip={fstr(pascal - true)}")
        print("  orders=" + ", ".join(f"O{idx}={fstr(val)}" for idx, val in enumerate(orders)))
        print(
            "  pair_excess="
            f"{fstr(pair_excess)}  high_order_net={fstr(high_order)}"
            f"  first_live_order={first_live}"
        )
        print(
            "  shape="
            f"span={max(speeds)-min(speeds)} gaps={gap_vector(speeds)} "
            f"edge_charge={sum(1 for p in speeds if p <= 2 or p >= 12)} "
            f"max_linear_run={max_linear_run(speeds)}"
        )
        print(f"  pair_shadow_false_positive={false_positive}")
    print("\nReading: j=3 has zero final dip but nonzero O2 excess and O3 debt.")
    print("The value shadow is right there, but the anatomy already knows order-3 is live.")


def local_minima_one_swap(j: int, n: int) -> tuple[list[tuple[int, ...]], dict[tuple[int, ...], Fraction]]:
    pool = list(itertools.combinations(range(1, n + 1), j))
    values = {p: lonely_measure(p) for p in pool}
    universe = set(range(1, n + 1))
    locals_: list[tuple[int, ...]] = []
    for p in pool:
        pset = set(p)
        is_local = True
        for out in p:
            reduced = pset - {out}
            for inn in universe - pset:
                q = tuple(sorted(reduced | {inn}))
                if values[q] < values[p]:
                    is_local = False
                    break
            if not is_local:
                break
        if is_local:
            locals_.append(p)
    return locals_, values


def print_exchange_traps() -> None:
    print("\n" + "=" * 96)
    print("ONE-SWAP IMPROVEMENT GRAPH: LOCAL TRAPS ARE A SEPARATE MAXIMIZER SIGNAL")
    print("=" * 96)
    for n in (13, 16):
        print(f"\nsearch universe {{1..{n}}}")
        for j in range(2, 6):
            locals_, values = local_minima_one_swap(j, n)
            global_value = min(values.values())
            globals_ = [p for p, val in values.items() if val == global_value]
            traps = [p for p in locals_ if values[p] > global_value]
            trap_gap = min((values[p] - global_value for p in traps), default=Fraction(0))
            shown = sorted(locals_, key=lambda p: (values[p], p))[:8]
            print(
                f"  j={j}: configs={comb(n,j)} global_value={fstr(global_value)} "
                f"global_count={len(globals_)} local_minima={len(locals_)} "
                f"nonglobal_traps={len(traps)} min_trap_gap={fstr(trap_gap)}"
            )
            print("      leading_local_minima=" + "; ".join(f"{p}:{fstr(values[p])}" for p in shown))
    print("\nReading: exchange optimality is not equivalent to scalar optimality.")
    print("Even k=10/j=3 has nonglobal one-swap sinks, so tournament non-transitivity is live early.")


def hamiltonian_paths(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def tournament_from_mask(n: int, mask: int) -> list[list[bool]]:
    adj = [[False] * n for _ in range(n)]
    pairs = list(itertools.combinations(range(n), 2))
    for bit, (i, j) in enumerate(pairs):
        if mask & (1 << bit):
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def c3_count(adj: list[list[bool]]) -> int:
    total = 0
    for a, b, c in itertools.combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            total += 1
    return total


def score_sequence(adj: list[list[bool]]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def print_small_tournament_spectrum() -> None:
    print("\n" + "=" * 96)
    print("SMALL TOURNAMENT HAMILTONIAN-PATH SPECTRUM")
    print("=" * 96)
    for n in range(3, 7):
        pairs = n * (n - 1) // 2
        hist: Counter[int] = Counter()
        score_by_h: dict[int, Counter[tuple[int, ...]]] = defaultdict(Counter)
        c3_by_h: dict[int, Counter[int]] = defaultdict(Counter)
        for mask in range(1 << pairs):
            adj = tournament_from_mask(n, mask)
            h = hamiltonian_paths(adj)
            hist[h] += 1
            score_by_h[h][score_sequence(adj)] += 1
            c3_by_h[h][c3_count(adj)] += 1
        max_h = max(hist)
        missing_odd = [x for x in range(1, max_h + 1, 2) if x not in hist]
        print(f"\nn={n}: labeled_tournaments={1 << pairs} distinct_H={len(hist)} max_H={max_h}")
        print(f"  H_hist={dict(sorted(hist.items()))}")
        print(f"  missing_odd_H_below_max={missing_odd}")
        print(f"  max_H_score_sequences={dict(score_by_h[max_h])}")
        print(f"  max_H_c3_counts={dict(c3_by_h[max_h])}")
    print("\nReading: H maximization is a rigidity/balance problem, not an LRC value problem.")
    print("It should guard quotient transfer and forbidden coincidences rather than replace cap anatomy.")


def strongly_connected_components(adj: list[list[bool]]) -> list[list[int]]:
    n = len(adj)
    index = 0
    stack: list[int] = []
    on_stack = [False] * n
    indices = [-1] * n
    low = [0] * n
    comps: list[list[int]] = []

    def visit(v: int) -> None:
        nonlocal index
        indices[v] = low[v] = index
        index += 1
        stack.append(v)
        on_stack[v] = True
        for w in range(n):
            if not adj[v][w]:
                continue
            if indices[w] == -1:
                visit(w)
                low[v] = min(low[v], low[w])
            elif on_stack[w]:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack[w] = False
                comp.append(w)
                if w == v:
                    break
            comps.append(sorted(comp))

    for v in range(n):
        if indices[v] == -1:
            visit(v)
    return comps


def count_directed_3cycles(adj: list[list[bool]]) -> int:
    return c3_count(adj)


def print_signal_tournament() -> None:
    print("\n" + "=" * 96)
    print("TOURNAMENT ANALYSIS ON SIGNAL FAMILIES")
    print("=" * 96)
    print("vertices are signals/proof obligations, not runners or arcs")
    axes = [
        "predicate",
        "value",
        "high_order",
        "local_trap",
        "quotient_guard",
        "tournament_bridge",
        "computable",
        "hypothesis_yield",
    ]
    signals: dict[str, dict[str, int]] = {
        "pair_pascal_shadow": dict(predicate=1, value=3, high_order=0, local_trap=0, quotient_guard=1, tournament_bridge=0, computable=3, hypothesis_yield=1),
        "high_order_anatomy": dict(predicate=2, value=3, high_order=3, local_trap=1, quotient_guard=2, tournament_bridge=1, computable=3, hypothesis_yield=3),
        "local_exchange_traps": dict(predicate=2, value=2, high_order=1, local_trap=3, quotient_guard=2, tournament_bridge=3, computable=3, hypothesis_yield=3),
        "deep_tail_corner_signature": dict(predicate=3, value=2, high_order=2, local_trap=1, quotient_guard=2, tournament_bridge=0, computable=2, hypothesis_yield=2),
        "boundary_mass_charge": dict(predicate=2, value=2, high_order=2, local_trap=1, quotient_guard=1, tournament_bridge=0, computable=3, hypothesis_yield=2),
        "decorrelation_inward_flux": dict(predicate=2, value=1, high_order=2, local_trap=1, quotient_guard=2, tournament_bridge=1, computable=1, hypothesis_yield=3),
        "far_block_coherence": dict(predicate=3, value=2, high_order=2, local_trap=2, quotient_guard=2, tournament_bridge=0, computable=2, hypothesis_yield=3),
        "H_spectrum_rigidity": dict(predicate=1, value=1, high_order=0, local_trap=2, quotient_guard=3, tournament_bridge=3, computable=3, hypothesis_yield=2),
        "scissors_forgetting_cost": dict(predicate=3, value=1, high_order=2, local_trap=1, quotient_guard=3, tournament_bridge=2, computable=2, hypothesis_yield=3),
        "raw_scalar_value": dict(predicate=1, value=3, high_order=0, local_trap=0, quotient_guard=0, tournament_bridge=0, computable=3, hypothesis_yield=0),
    }
    names = list(signals)
    tie_order = {name: i for i, name in enumerate(names)}
    n = len(names)
    adj = [[False] * n for _ in range(n)]
    edge_votes: dict[tuple[str, str], tuple[int, int]] = {}
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i >= j:
                continue
            votes_a = 0
            votes_b = 0
            for axis in axes:
                av = signals[a][axis]
                bv = signals[b][axis]
                if av > bv:
                    votes_a += 1
                elif bv > av:
                    votes_b += 1
            if votes_a > votes_b or (votes_a == votes_b and tie_order[a] < tie_order[b]):
                adj[i][j] = True
                edge_votes[(a, b)] = (votes_a, votes_b)
            else:
                adj[j][i] = True
                edge_votes[(b, a)] = (votes_b, votes_a)
    scores = {names[i]: sum(adj[i]) for i in range(n)}
    hist = Counter(scores.values())
    comps = [[names[i] for i in comp] for comp in strongly_connected_components(adj)]
    path_count = hamiltonian_paths(adj)
    cycles = count_directed_3cycles(adj)
    print("pairwise observable: majority over retention axes")
    print("gauge axes=" + ", ".join(axes))
    print("score_hist=" + str(dict(sorted(hist.items()))))
    print("scores=" + str(dict(sorted(scores.items(), key=lambda kv: (-kv[1], kv[0])))))
    print(f"directed_3cycles={cycles}  SCCs={comps}  Hamiltonian_path_count={path_count}")
    sorted_vertices = sorted(names, key=lambda name: (-scores[name], name))
    print("one retention Hamiltonian path=" + " > ".join(sorted_vertices))
    print("top edge vote samples:")
    for a, b in list(edge_votes)[:12]:
        va, vb = edge_votes[(a, b)]
        print(f"  {a} -> {b} by {va}:{vb}")
    print("\nNew signals to measure next:")
    for signal in [
        "maximizer_currency_vector",
        "pascal_anatomy_residual",
        "first_live_interaction_order",
        "exchange_trap_index",
        "deep_tail_signature",
        "boundary_mass_charge",
        "decorrelation_inward_flux",
        "coarse_to_fine_transfer_risk",
        "scissors_split_count",
        "curvature_plateau_tail",
        "far_block_coherence",
        "pair_shadow_false_positive",
        "H_spectrum_transfer_risk",
    ]:
        print(f"  - {signal}")


def main() -> None:
    print_cap_anatomy()
    print_exchange_traps()
    print_small_tournament_spectrum()
    print_signal_tournament()
    print("\n" + "=" * 96)
    print("ASSUMPTION CHALLENGE")
    print("=" * 96)
    print("Do not force tournament vertices to be runners or arcs.")
    print("Useful vertices here include gaps, complements, interaction orders, local sinks,")
    print("observer charts, quotient debts, Hamiltonian spectra, and proof obligations.")
    print("The quotient preserves only the chosen predicate currency; it destroys the others.")


if __name__ == "__main__":
    main()
