#!/usr/bin/env python3
"""
Signed LRC pairwise floor via the modulus-n cut obstruction.  (monad-explorer-S711)

Challenge the default vertex set: at the floor witness t = 1/n, the useful
quotient is not "runners" but the residue classes mod n.  This quotient
preserves exactly the divisibility obstruction

    n | (eps_i v_i - eps_j v_j),

and forgets the metric size data that controls off-grid improvements.

Tournament lens, stated honestly:
  - pairwise observable: for a cut eps, pair {i,j} is BAD on the n-grid iff
      eps_i v_i - eps_j v_j == 0 mod n
  - switch/gauge: the sign cut eps / bipartition (THM-426)
  - why no tournament fingerprints: badness is a symmetric incompatibility
    relation, not an orientation, so this script reports cut-feasibility counts
    and residue contradiction motifs instead of SCC/triangle data

Scan plan:
  1. enumerate gcd-1 speed sets V with |V| = r, max(V) <= B
  2. test whether some cut has NO n-multiple relative speed ("an n-cut")
     -> then t = 1/n immediately witnesses Gstar(V) >= 1/n
  3. run the exact maximin search only on the NO-N-CUT residue obstruction class

Outputs:
  - counts of n-cut-feasible vs no-n-cut sets
  - exact below/equal/above-floor split inside the no-n-cut class
  - new bounded minima and the residue motifs of all below-floor examples

This preserves the on-grid obstruction and exposes its limit: many no-n-cut
sets remain above floor, so the residue quotient is informative but incomplete.
"""

from collections import Counter
from fractions import Fraction as F
from functools import reduce
from itertools import combinations, product
from math import gcd


CASES = [
    (3, 11),
    (4, 9),
    (5, 8),
    (5, 9),
    (5, 10),
    (6, 7),
    (6, 8),
]

GAP_CACHE = {}


def pr(*args):
    print(*args, flush=True)


def norm(x):
    f = x - int(x)
    if f < 0:
        f += 1
    return min(f, 1 - f)


def candidate_times(W):
    ms = set(abs(w) for w in W if w)
    for a, b in combinations(W, 2):
        ms.add(abs(a + b))
        ms.add(abs(a - b))
    ms.discard(0)
    times = set()
    for m in ms:
        for a in range(1, 2 * m):
            times.add(F(a, 2 * m))
    return times


def maximin_gap(W):
    key = tuple(sorted(abs(w) for w in W if w))
    if key in GAP_CACHE:
        return GAP_CACHE[key]
    best = F(0)
    for t in candidate_times(key):
        m = min(norm(w * t) for w in key)
        if m > best:
            best = m
    GAP_CACHE[key] = best
    return best


def best_times(W, best):
    key = tuple(sorted(abs(w) for w in W if w))
    times = []
    for t in candidate_times(key):
        if min(norm(w * t) for w in key) == best:
            times.append(t)
    return sorted(times)


def signed_relative_speeds(V, eps):
    return [eps[i] * V[i] - eps[j] * V[j] for i, j in combinations(range(len(V)), 2)]


def exact_gstar(V, keep_data=False):
    best = F(0)
    data = []
    r = len(V)
    for tail in product([1, -1], repeat=r - 1):
        eps = (1,) + tail
        W = signed_relative_speeds(V, eps)
        g = maximin_gap(W)
        if g > best:
            best = g
            data = []
        if keep_data and g == best:
            data.append((eps, W))
    if keep_data:
        return best, data
    return best


def n_cut_data(V):
    n = len(V) + 1
    out = []
    r = len(V)
    for tail in product([1, -1], repeat=r - 1):
        eps = (1,) + tail
        W = signed_relative_speeds(V, eps)
        if all(w % n != 0 for w in W):
            out.append((eps, W))
    return out


def contradiction_motifs(V):
    n = len(V) + 1
    counts = Counter(v % n for v in V)
    motifs = []
    seen = set()
    for a in sorted(counts):
        b = (-a) % n
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        if counts[a] >= 2 and counts[b] >= 1:
            motifs.append((a, b, counts[a], counts[b]))
            seen.add(key)
        elif a == b and counts[a] >= 2:
            motifs.append((a, b, counts[a], counts[b]))
            seen.add(key)
    return motifs


def motif_label(motif):
    a, b, ca, cb = motif
    if a == b:
        return f"{a},{a} (count {ca})"
    x, y = sorted((a, b))
    return f"{x},{y} ({ca}/{cb})"


def enum_sets(r, B):
    return [c for c in combinations(range(1, B + 1), r) if reduce(gcd, c) == 1]


def cut_sides(V, eps):
    plus = [V[i] for i, e in enumerate(eps) if e == 1]
    minus = [V[i] for i, e in enumerate(eps) if e == -1]
    return plus, minus


def summarize_case(r, B):
    n = r + 1
    floor = F(1, n)
    sets = enum_sets(r, B)
    feasible = 0
    no_n_cut = []

    for V in sets:
        cuts = n_cut_data(V)
        if cuts:
            feasible += 1
        else:
            no_n_cut.append(V)

    below = []
    equal = []
    above_examples = []
    worst = None
    worst_sets = []

    for idx, V in enumerate(no_n_cut, start=1):
        g = exact_gstar(V)
        if worst is None or g < worst:
            worst = g
            worst_sets = [V]
        elif g == worst:
            worst_sets.append(V)

        if g < floor:
            below.append((V, g))
        elif g == floor:
            equal.append((V, g))
        elif len(above_examples) < 3:
            above_examples.append((V, g))

    detailed_below = []
    for V, g in below:
        g2, opt = exact_gstar(V, keep_data=True)
        assert g2 == g
        eps, W = opt[0]
        plus, minus = cut_sides(V, eps)
        detailed_below.append(
            {
                "V": V,
                "g": g,
                "residues": tuple(v % n for v in V),
                "motifs": contradiction_motifs(V),
                "cut": eps,
                "plus": plus,
                "minus": minus,
                "D": tuple(sorted(abs(w) for w in W)),
                "times": best_times(W, g),
            }
        )

    return {
        "r": r,
        "n": n,
        "B": B,
        "total": len(sets),
        "feasible": feasible,
        "no_n_cut": len(no_n_cut),
        "floor": floor,
        "below": below,
        "equal": equal,
        "above_examples": above_examples,
        "worst": worst,
        "worst_sets": worst_sets,
        "detailed_below": detailed_below,
    }


def main():
    pr("=" * 78)
    pr("SIGNED PAIRWISE FLOOR VIA THE MODULUS-n CUT OBSTRUCTION  monad-S711")
    pr("=" * 78)
    pr("A cut is n-grid FEASIBLE if none of its relative speeds is 0 mod n.")
    pr("If a feasible cut exists, t = 1/n certifies Gstar >= 1/n immediately.")
    pr("Exact search is run only on the NO-N-CUT obstruction class.\n")

    summaries = []
    global_motif_hist = Counter()

    for r, B in CASES:
        summary = summarize_case(r, B)
        summaries.append(summary)

        pr("-" * 78)
        pr(
            f"r={summary['r']} movers, n={summary['n']}, B<={summary['B']}, "
            f"{summary['total']} gcd-1 sets, floor 1/n={summary['floor']}"
        )
        pr(
            f"  n-cut feasible: {summary['feasible']}   "
            f"no-n-cut: {summary['no_n_cut']}"
        )
        pr(
            f"  exact split inside no-n-cut: below={len(summary['below'])}  "
            f"equal={len(summary['equal'])}  "
            f"strictly above={summary['no_n_cut'] - len(summary['below']) - len(summary['equal'])}"
        )
        if summary["no_n_cut"] == 0:
            pr("  exact phase skipped: every set already has an n-cut witness at t=1/n")
        else:
            pr(
                f"  worst exact Gstar on this bounded family: {summary['worst']}   "
                f"witness sets: {summary['worst_sets'][:8]}"
            )

        if summary["detailed_below"]:
            pr("  below-floor examples:")
            for row in summary["detailed_below"]:
                labels = [motif_label(m) for m in row["motifs"]]
                for label in labels:
                    global_motif_hist[(summary["n"], label)] += 1
                pr(
                    f"    V={row['V']} residues={row['residues']} Gstar={row['g']} "
                    f"motifs={labels}"
                )
                pr(
                    f"      best cut A={row['plus']}  B={row['minus']}  "
                    f"D={row['D']}  t*={row['times'][:6]}"
                )
        if summary["above_examples"]:
            pr("  no-n-cut but still ABOVE floor (shows residue quotient loses metric size):")
            for V, g in summary["above_examples"]:
                pr(f"    V={V}  Gstar={g}")
        if summary["equal"]:
            pr(f"  first no-n-cut equality example at the floor: {summary['equal'][0][0]}")

    pr("\n" + "=" * 78)
    pr("GLOBAL SUMMARY")
    pr("=" * 78)
    pr(f"{'r':>2} {'n':>2} {'B':>3} {'total':>6} {'n-cut':>6} {'no-cut':>6} {'below':>6} {'worst':>8}")
    for summary in summaries:
        worst = "-" if summary["worst"] is None else str(summary["worst"])
        pr(
            f"{summary['r']:>2} {summary['n']:>2} {summary['B']:>3} "
            f"{summary['total']:>6} {summary['feasible']:>6} {summary['no_n_cut']:>6} "
            f"{len(summary['below']):>6} {worst:>8}"
        )

    pr("\nResidue contradiction motifs among all below-floor examples:")
    for (n, label), count in sorted(global_motif_hist.items()):
        pr(f"  n={n}: {label}  ->  {count}")


if __name__ == "__main__":
    main()
