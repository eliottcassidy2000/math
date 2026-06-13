#!/usr/bin/env python3
"""
lrc_obligation_hypergraph_s574.py

codex-2026-06-03 S574

Rearrange the S573 D/U/N clock-blocker ledger as a labelled obligation
hypergraph.

Speeds are hyperedges.  Clock obligations are vertices:
  D_q: some speed is divisible by q;
  U_a: unit antipodal shell {a,2n-1-a} is hit;
  N_j: the j/n clock is blocked at distance 0 or 1/n.

Question:
  What distinguishes floor rows, open-gap lifts, and nonunit-hole rows once the
  continuum has been replaced by this finite hitting-set object?

Tournament Analysis / assumption challenge:
- Vertices are proof obligations, not runners.
- Pair observable is (coverage count, layer, label), oriented toward the more
  fragile obligation (smaller coverage count).
- Tie Hamiltonian path orders D, then U, then N, by increasing label.
- This preserves possible strict-sub-edge status because failing any obligation
  gives an immediate clock witness at or above 2/(2n-1).  It destroys exact
  continuous interval geometry, so active witness data are printed separately.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


Obligation = tuple[str, int]


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm(x: Fraction) -> Fraction:
    r = frac_part(x)
    return min(r, 1 - r)


def score_time(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(norm(Fraction(v) * t) for v in speeds)


def exact_maximin(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    candidates: set[Fraction] = set()
    for i, a in enumerate(speeds):
        for m in range(a):
            t = Fraction(2 * m + 1, 2 * a)
            if 0 < t < 1:
                candidates.add(t)
        for b in speeds[i + 1 :]:
            for den in (a + b, abs(a - b)):
                if den <= 0:
                    continue
                for m in range(1, den):
                    candidates.add(Fraction(m, den))

    best = Fraction(0)
    best_t = Fraction(0)
    for t in candidates:
        score = score_time(speeds, t)
        if (score, -t) > (best, -best_t):
            best = score
            best_t = t
    return best, best_t


def fmt_frac(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)


def obligation_universe(k: int) -> list[Obligation]:
    n = k + 1
    c = 2 * k + 1
    out: list[Obligation] = []
    out.extend(("D", q) for q in range(2, n))
    out.extend(("U", a) for a in range(1, n) if gcd(a, c) == 1)
    out.extend(("N", j) for j in range(1, n))
    return out


def covers(v: int, obligation: Obligation, k: int) -> bool:
    layer, label = obligation
    n = k + 1
    c = 2 * k + 1
    if layer == "D":
        return v % label == 0
    if layer == "U":
        r = v % c
        return r == label or r == c - label
    if layer == "N":
        r = (v * label) % n
        return min(r, n - r) <= 1
    raise ValueError(obligation)


def incidence(speeds: tuple[int, ...]) -> dict[Obligation, tuple[int, ...]]:
    k = len(speeds)
    out: dict[Obligation, tuple[int, ...]] = {}
    for obligation in obligation_universe(k):
        out[obligation] = tuple(v for v in speeds if covers(v, obligation, k))
    return out


def layer_hist(items: list[Obligation] | tuple[Obligation, ...]) -> dict[str, int]:
    counts = Counter(layer for layer, _ in items)
    return dict(sorted(counts.items()))


def private_obligations(speeds: tuple[int, ...], inc: dict[Obligation, tuple[int, ...]]) -> dict[int, list[Obligation]]:
    out: dict[int, list[Obligation]] = {v: [] for v in speeds}
    for obligation, owners in inc.items():
        if len(owners) == 1:
            out[owners[0]].append(obligation)
    return out


def active_runners(speeds: tuple[int, ...], t: Fraction) -> tuple[int, ...]:
    score = score_time(speeds, t)
    return tuple(v for v in speeds if norm(Fraction(v) * t) == score)


def flipset(speeds: tuple[int, ...]) -> tuple[int, ...]:
    k = len(speeds)
    c = 2 * k + 1
    residues = {v % c for v in speeds}
    return tuple(a for a in range(1, k + 1) if c - a in residues)


def shell_defects(speeds: tuple[int, ...]) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    k = len(speeds)
    c = 2 * k + 1
    counts: Counter[int] = Counter()
    for v in speeds:
        r = v % c
        if r:
            counts[min(r, c - r)] += 1
    missed = tuple((a, c - a) for a in range(1, k + 1) if counts[a] == 0)
    doubled = tuple((a, c - a) for a in range(1, k + 1) if counts[a] > 1)
    return missed, doubled


def tournament_fingerprint(inc: dict[Obligation, tuple[int, ...]]) -> tuple[dict[int, int], int, list[int], int]:
    obligations = sorted(inc, key=lambda o: (len(inc[o]), {"D": 0, "U": 1, "N": 2}[o[0]], o[1]))
    idx = {o: i for i, o in enumerate(obligations)}
    n = len(obligations)
    adj = [[False] * n for _ in range(n)]
    scores = [0] * n
    for a, b in combinations(obligations, 2):
        ka = (len(inc[a]), {"D": 0, "U": 1, "N": 2}[a[0]], a[1])
        kb = (len(inc[b]), {"D": 0, "U": 1, "N": 2}[b[0]], b[1])
        winner, loser = (a, b) if ka < kb else (b, a)
        adj[idx[winner]][idx[loser]] = True
        scores[idx[winner]] += 1

    cycles = 0
    for i, j, k in combinations(range(n), 3):
        cyc = (
            (adj[i][j] and adj[j][k] and adj[k][i])
            or (adj[i][k] and adj[k][j] and adj[j][i])
        )
        cycles += int(cyc)

    scc_sizes = strongly_connected_component_sizes(adj)
    # The tie path is the sorted obligation order, and it is Hamiltonian by construction.
    ham_paths_at_least = 1
    return dict(sorted(Counter(scores).items())), cycles, scc_sizes, ham_paths_at_least


def strongly_connected_component_sizes(adj: list[list[bool]]) -> list[int]:
    n = len(adj)
    rev = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            rev[j][i] = adj[i][j]

    seen = [False] * n
    order: list[int] = []

    def dfs1(v: int) -> None:
        seen[v] = True
        for w, edge in enumerate(adj[v]):
            if edge and not seen[w]:
                dfs1(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs1(v)

    comps = []
    seen = [False] * n
    for start in reversed(order):
        if seen[start]:
            continue
        q = deque([start])
        seen[start] = True
        size = 0
        while q:
            v = q.popleft()
            size += 1
            for w, edge in enumerate(rev[v]):
                if edge and not seen[w]:
                    seen[w] = True
                    q.append(w)
        comps.append(size)
    return sorted(comps, reverse=True)


def print_row(row: Row) -> None:
    speeds = row.speeds
    k = len(speeds)
    n = k + 1
    c = 2 * k + 1
    floor = Fraction(1, n)
    edge = Fraction(2, c)
    maximin, witness = exact_maximin(speeds)
    inc = incidence(speeds)
    uncovered = [o for o, owners in inc.items() if not owners]
    critical = [o for o, owners in inc.items() if len(owners) == 1]
    counts_by_layer: dict[str, list[int]] = defaultdict(list)
    for obligation, owners in inc.items():
        counts_by_layer[obligation[0]].append(len(owners))
    private = private_obligations(speeds, inc)
    essential = {v: obs for v, obs in private.items() if obs}
    missed, doubled = shell_defects(speeds)
    score_hist, cycles, sccs, ham_paths = tournament_fingerprint(inc)

    print(f"[{row.label}]")
    print(f"  speeds={speeds}")
    print(f"  n={n} C={c} floor={fmt_frac(floor)} edge={fmt_frac(edge)}")
    print(
        f"  M={fmt_frac(maximin)} witness={fmt_frac(witness)} "
        f"active={active_runners(speeds, witness)}"
    )
    print(f"  flipset={flipset(speeds)} missed_shells={missed or '-'} doubled_shells={doubled or '-'}")
    print(f"  obligations={len(inc)} uncovered={uncovered or '-'} critical={len(critical)} {layer_hist(critical)}")
    layer_bits = []
    for layer in ("D", "U", "N"):
        vals = counts_by_layer[layer]
        layer_bits.append(f"{layer}:min={min(vals)},max={max(vals)},avg={sum(vals)/len(vals):.2f}")
    print(f"  coverage_by_layer {'; '.join(layer_bits)}")
    print("  private_obligations")
    for v, obs in sorted(essential.items()):
        preview = ", ".join(f"{a}{b}" for a, b in obs[:7])
        more = "" if len(obs) <= 7 else f", ... (+{len(obs)-7})"
        print(f"    v={v:<3} count={len(obs):<2} {layer_hist(obs)} {preview}{more}")
    if not essential:
        print("    -")
    print(
        "  obligation_tournament "
        f"score_hist={score_hist} directed_3_cycles={cycles} "
        f"sccs={sccs} hamiltonian_paths_at_least={ham_paths}"
    )
    print()


def main() -> None:
    print("LRC clock-obligation hypergraph (codex-2026-06-03 S574)")
    print("=" * 78)
    print("Rows are speed sets; columns are D/U/N clock obligations.")
    print("A strict sub-edge row must be a full hitting set for this hypergraph.")
    print("Private obligations are the candidate pivots for exchange/descent proofs.\n")

    rows = [
        Row("AP n=7 floor", (1, 2, 3, 4, 5, 6)),
        Row("lifted flip {2} n=7 open-gap", (1, 5, 6, 11, 16, 17)),
        Row("AP n=8 floor", (1, 2, 3, 4, 5, 6, 7)),
        Row("n=8 nonunit floor", (1, 2, 3, 4, 5, 7, 12)),
        Row("n=8 nonunit open A", (1, 2, 3, 4, 5, 7, 18)),
        Row("n=8 nonunit open B", (1, 3, 4, 5, 7, 13, 18)),
        Row("AP n=14 floor", tuple(range(1, 14))),
        Row("V* n=14 floor", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)),
    ]
    for row in rows:
        print_row(row)

    print("Synthesis")
    print("-" * 78)
    print("The D/U/N ledger is a labelled hitting-set hypergraph.")
    print("Open-gap lifts are not failures of the ledger; they are full covers whose")
    print("exact pair-sum witness rises above the floor but stays below the unit-shell edge.")
    print("Private obligations identify descent pivots: replacing or lowering a runner must")
    print("preserve its private D/U/N labels or immediately creates a clock witness.")
    print("This is the owner-labelled event word requested by the Burnside/Fourier framing,")
    print("but compressed to the obligations that are provably necessary for sub-edge rows.")


if __name__ == "__main__":
    main()
