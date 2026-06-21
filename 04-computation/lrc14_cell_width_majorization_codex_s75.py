#!/usr/bin/env python3
"""
codex-2026-06-21-S75

Cell-width majorization scout for the bounded LRC(14) AP/consec wall.

Incoming Huffer-Shepp work shows that per-cell inequalities are false:
some non-AP shapes beat consec on an individual W_a cell while losing in
the aggregate sum measS7=sum_a W_a.  This script asks whether the right
cell-level quotient is instead the sorted width vector.

Method:
  * enumerate anchored full-residue shapes in bounded spans;
  * compute the exact six cell widths W_a, a=1..6;
  * compare fixed cyclic windows against consec as a no-go baseline;
  * compare sorted top-L sums against consec's sorted W-vector.

The sorted comparison is compatible with the dilation action
W_a(cE)=W_{ca}(E): it forgets cell labels and keeps the aggregate shape of
the six survival widths.

Tournament Analysis:
  Vertices are proof lenses / aggregation choices, not runners or arcs.
  Pair observable is

      (affine_safe, pins_ap, finite_leaks, formal_ready, composes_with_wide).

  The switch is lexicographic comparison, with name as final tie-break.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import math


P = 7
CELLS = tuple(range(1, P))


def W_a(E, a):
    """Exact width in cell a where E covers all seven sectors."""
    E = tuple(sorted(set(E)))
    lo = Fraction(a, P) - Fraction(1, 2 * P)
    hi = Fraction(a, P) + Fraction(1, 2 * P)
    bps = {lo, hi}
    for e in E:
        if e == 0:
            continue
        d = P * abs(e)
        j0 = math.floor(lo * d)
        j1 = math.ceil(hi * d)
        for j in range(j0 - 1, j1 + 2):
            x = Fraction(j, d)
            if lo <= x <= hi:
                bps.add(x)
    bps = sorted(bps)
    total = Fraction(0)
    for left, right in zip(bps, bps[1:]):
        if right <= left:
            continue
        xm = (left + right) / 2
        hit = set()
        for e in E:
            v = e * xm
            v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * P) // v.denominator)
        if len(hit) == P:
            total += right - left
    return total


@lru_cache(None)
def Wvec(E):
    return tuple(W_a(E, a) for a in CELLS)


def full_residue_bank(k, span):
    for combo in combinations(range(1, span + 1), k - 1):
        E = (0,) + combo
        if set(e % P for e in E) == set(range(P)):
            yield E


def cyclic_window_sums(vec, length):
    n = len(vec)
    return tuple(sum(vec[(i + j) % n] for j in range(length)) for i in range(n))


def sorted_top_sums(vec):
    total = Fraction(0)
    out = []
    for x in sorted(vec, reverse=True):
        total += x
        out.append(total)
    return tuple(out)


def fmt(fr):
    return f"{fr} ({float(fr):.9f})"


def worst_violations(bank, base_vec, mode):
    base = (
        {L: cyclic_window_sums(base_vec, L) for L in range(1, 6)}
        if mode == "cyclic"
        else {L: sorted_top_sums(base_vec)[L - 1] for L in range(1, 6)}
    )
    out = {}
    for L in range(1, 6):
        count = 0
        worst = None
        for E, vec in bank:
            if mode == "cyclic":
                vals = cyclic_window_sums(vec, L)
                for start, value in enumerate(vals, start=1):
                    gain = value - base[L][start - 1]
                    if gain > 0:
                        count += 1
                        item = (gain, E, start, value, base[L][start - 1], sum(vec) - sum(base_vec))
                        if worst is None or item[0] > worst[0]:
                            worst = item
            else:
                value = sorted_top_sums(vec)[L - 1]
                gain = value - base[L]
                if gain > 0:
                    count += 1
                    item = (gain, E, None, value, base[L], sum(vec) - sum(base_vec))
                    if worst is None or item[0] > worst[0]:
                        worst = item
        out[L] = (count, worst)
    return out


def leak_detail(E, base_vec):
    vec = Wvec(E)
    base_sorted = sorted(base_vec, reverse=True)
    e_sorted = sorted(vec, reverse=True)
    rows = []
    running_base = Fraction(0)
    running_e = Fraction(0)
    for i, (b, e) in enumerate(zip(base_sorted, e_sorted), start=1):
        running_base += b
        running_e += e
        rows.append((i, e - b, running_e - running_base))
    return vec, rows


LENSES = {
    "sorted_cell_majorization": (3, 2, 2, 2, 3),
    "one_sink_compensation": (3, 2, 3, 1, 3),
    "relation_marginal_bound": (3, 3, 2, 2, 3),
    "conductance_bottleneck": (3, 2, 1, 2, 2),
    "fixed_cyclic_windows": (1, 1, 1, 2, 1),
    "win_disconnected_split": (2, 1, 1, 2, 1),
    "per_cell_reflection": (1, 0, 0, 3, 0),
    "residue_affine_moments": (3, 0, 0, 2, 1),
}


def tournament_from_scores(scores):
    names = list(scores)
    adj = {u: set() for u in names}
    for i, u in enumerate(names):
        for v in names[i + 1 :]:
            if (scores[u], u) >= (scores[v], v):
                adj[u].add(v)
            else:
                adj[v].add(u)
    return names, adj


def count_three_cycles(names, adj):
    count = 0
    for a, b, c in combinations(names, 3):
        edges = sum(1 for u, v in ((a, b), (b, c), (c, a)) if v in adj[u])
        if edges in (0, 3):
            count += 1
    return count


def hamiltonian_paths(names, adj):
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    outmask = [0] * n
    for u in names:
        for v in adj[u]:
            outmask[idx[u]] |= 1 << idx[v]

    @lru_cache(None)
    def dp(mask, last):
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if (prev_mask & (1 << prev)) and (outmask[prev] & (1 << last)):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def main():
    spans = {8: 15, 9: 14, 10: 14, 11: 14, 12: 14}
    print("LRC14 cell-width majorization scout (codex S75)")
    print()
    for k, span in spans.items():
        consec = tuple(range(k))
        base_vec = Wvec(consec)
        bank = [(E, Wvec(E)) for E in full_residue_bank(k, span)]
        print(f"== k={k}, span<={span}, full-residue shapes={len(bank)} ==")
        print(f"consec W sorted = {[str(x) for x in sorted(base_vec, reverse=True)]}")
        print(f"consec measS7 = {fmt(sum(base_vec))}")

        print("fixed cyclic cell-window violations (not dilation-invariant):")
        for L, (count, worst) in worst_violations(bank, base_vec, "cyclic").items():
            if worst is None:
                print(f"  L={L}: none")
            else:
                gain, E, start, value, ref, dtotal = worst
                print(
                    f"  L={L}: count={count}, worst_gain={gain}, start={start}, "
                    f"E={E}, dtotal={dtotal}"
                )

        print("sorted top-L width violations (dilation-safe multiset quotient):")
        leaks = []
        for L, (count, worst) in worst_violations(bank, base_vec, "sorted").items():
            if worst is None:
                print(f"  top{L}: none")
            else:
                gain, E, _start, value, ref, dtotal = worst
                leaks.append(E)
                print(
                    f"  top{L}: count={count}, worst_gain={gain}, "
                    f"E={E}, dtotal={dtotal}"
                )
        for E in sorted(set(leaks)):
            vec, rows = leak_detail(E, base_vec)
            print(f"  leak detail E={E}: W={[str(x) for x in vec]}")
            for i, point_delta, prefix_delta in rows:
                print(
                    f"    sorted index {i}: point_delta={point_delta}, "
                    f"prefix_delta={prefix_delta}"
                )
        print()

    print("Tournament Analysis over proof lenses")
    names, adj = tournament_from_scores(LENSES)
    scores = Counter(len(adj[name]) for name in names)
    ranking = sorted(names, key=lambda name: (LENSES[name], name), reverse=True)
    print("pair observable = (affine_safe, pins_ap, finite_leaks, formal_ready, composes_with_wide)")
    print("ranking:")
    for name in ranking:
        print(f"  {name}: score={len(adj[name])}, observable={LENSES[name]}")
    print(f"score histogram = {dict(sorted(scores.items()))}")
    print(f"directed 3-cycles = {count_three_cycles(names, adj)}")
    print(f"Hamiltonian path count = {hamiltonian_paths(names, adj)}")


if __name__ == "__main__":
    main()
