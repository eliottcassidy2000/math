"""Session opus-2026-10-06-S14: exact checks behind the typed bridges from five owner-supplied papers
(Mumford-Shah; Golod-Shafarevich towers; KLS O(1); the countable reals; subquadratic 3SUM/APSP) to the
Collatz window-alphabet thread.  Companion of syracuse_window_alphabet_general_m_20261005.py (imported).

  A. Collatz backward tree as a "tower of covers" (Golod-Shafarevich bridge): the compressed map T has
     preimages 2y (always) and (2y-1)/3 (iff y = 2 mod 3); the even child of a branching vertex is = 1 mod 3,
     so no vertex has a complete binary backward subtree of depth 2 (the section of an etale tower that
     "splits completely at every level" has no Collatz counterpart); layer sizes grow like (4/3)^n.
  B. Near-zero-sum reading of the gate crossings (3SUM / Exact-Triangle bridge): the threshold margin
     2^(t_k)/3^k - 1 and the gate radius R_k at the semiconvergent lags of log2 3; exceptions only at 17, 29, 41.
  C. The first pinched windows (countable-reals bridge, "uniform versus exact"): the 12 lag-17 starts on the
     plus sheet and the 3 lag-12 starts on the minus sheet; which chord differs from the word prediction
     (always the long chord), and whether the actual pattern is realisable by any word (Farkas / Bellman-Ford
     feasibility of the staircase difference constraints): exactly 1 of 12 (437) and 2 of 3 (165, 549).

Run:  python -X utf8 five_papers_bridges_20261006.py [--save]
"""
from __future__ import annotations

import os
import sys
import warnings
from fractions import Fraction

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from syracuse_window_alphabet_general_m_20261005 import threshold, pairs  # noqa: E402
from syracuse_window_alphabet_20261005 import v2, chord_bits  # noqa: E402

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


# ------------------------------------------------------------------ A
def preimages(y):
    out = [2 * y]
    if y % 3 == 2:
        out.append((2 * y - 1) // 3)
    return out


def layer_counts(root, depth):
    layer, seen, counts = [root], {root}, [1]
    for _ in range(depth):
        nxt = []
        for y in layer:
            for p in preimages(y):
                if p not in seen:
                    seen.add(p)
                    nxt.append(p)
        layer = nxt
        counts.append(len(layer))
    return counts


def complete_depth(y, maxd=4):
    layer = [y]
    for d in range(maxd):
        if not all(v % 3 == 2 for v in layer):
            return d
        layer = [p for v in layer for p in preimages(v)]
    return maxd


def section_A():
    say("## A. The Collatz backward tree never splits completely (Golod-Shafarevich bridge)")
    say()
    say("Compressed map T(x) = x/2 (x even), (3x+1)/2 (x odd); preimages of y: 2y always, (2y-1)/3 iff y = 2 mod 3.")
    assert all((2 * y) % 3 == 1 for y in range(2, 10 ** 5, 3))
    say("PROVED (congruence): if y = 2 mod 3 then 2y = 1 mod 3, so the even child of a branching vertex never branches;")
    say("hence no vertex has a complete binary backward subtree of depth 2, and no 'section splits completely at every level'.")
    md = max(complete_depth(y) for y in range(1, 10 ** 5))
    assert md == 1
    say("FINITE-EXACT: maximal complete-splitting depth over all starts below 10^5 = %d." % md)
    c = layer_counts(8, 40)
    say("Backward layer sizes from 8, depth 0..40: %s" % c)
    r = [c[i + 1] / c[i] for i in range(30, 40)]
    say("Successive ratios at depth 30..40: %s (mean branching 1 + P(y = 2 mod 3) = 4/3 = 1.3333)" % [round(x, 4) for x in r])
    assert all(abs(x - 4 / 3) < 0.01 for x in r)
    say()


# ------------------------------------------------------------------ B
def section_B():
    say("## B. Near-zero-sum reading of the gate crossings (3SUM / Exact-Triangle bridge)")
    say()
    say("A k-step word w with sum s closes a cycle iff (2^s - 3^k) y = c_w exactly (an exact zero sum of the two ledgers plus")
    say("the correction); a gate crossing is an approximate zero sum: y < c_w/(2^s - 3^k).  Margin 2^(t_k)/3^k - 1 and gate")
    say("radius R_k = c_max/(2^(t_k) - 3^k), c_max = 3^(k-1) + 2^(t_k-k+1)(3^(k-1) - 2^(k-1)), for the lags with margin < 0.06:")
    rows = []
    for k in range(2, 46):
        t = threshold(k)
        margin = 2 ** t / 3 ** k - 1
        c_max = 3 ** (k - 1) + 2 ** (t - k + 1) * (3 ** (k - 1) - 2 ** (k - 1))
        R = Fraction(c_max, 2 ** t - 3 ** k)
        if margin < 0.06:
            rows.append((k, t, margin, float(R)))
            say("  k=%2d t_k=%3d margin=%.4f R_k=%.1f%s" % (k, t, margin, float(R), {17: "  (12 exceptions)", 29: "  (56)", 41: "  (62)", 5: "  (none: no lattice point below the gate)"}.get(k, "")))
    assert [k for k, *_ in rows] == [5, 17, 29, 41]
    say("The small-margin lags below 45 are exactly 5, 17, 29, 41 (the upper semiconvergents 8/5, 27/17, 46/29, 65/41 of log2 3);")
    say("exceptions exist at 17, 29, 41 and not at 5, so the margin is necessary and the cylinder offset decides (S13).")
    say()


# ------------------------------------------------------------------ C
def syr(x, sheet):
    y = 3 * x + sheet
    a = v2(y)
    return y >> a, a


def window(x, m, sheet):
    vals, avs = [x], []
    for _ in range(m - 1):
        y, a = syr(vals[-1], sheet)
        vals.append(y)
        avs.append(a)
    return vals, avs


def generic(avs, m):
    return tuple(1 if sum(avs[i:j]) >= threshold(j - i) else 0 for i, j in pairs(m))


def farkas_feasible(m, pat):
    """Staircase s_0 < s_1 < ... < s_(m-1), steps >= 1; forward chord (i,j): s_j - s_i >= t_(j-i); backward: <= t_(j-i) - 1.
    Feasible iff the difference-constraint graph has no negative cycle (Bellman-Ford)."""
    edges = [(i + 1, i, -1) for i in range(m - 1)]
    for (i, j), b in zip(pairs(m), pat):
        t = threshold(j - i)
        edges.append((j, i, -t) if b else (i, j, t - 1))
    dist = [0] * m
    for _ in range(m):
        changed = False
        for u, v, w in edges:
            if dist[u] + w < dist[v]:
                dist[v] = dist[u] + w
                changed = True
        if not changed:
            return True
    return False


def section_C():
    say("## C. The first pinched windows: which chord separates, and whether the actual pattern is word-realisable")
    say()
    plus = [165, 171, 231, 257, 259, 387, 389, 391, 437, 581, 587, 589]
    real = []
    for x in plus:
        vals, avs = window(x, 18, +1)
        pat = chord_bits(vals, "DESC")
        diff = [(i, j) for (i, j), a, b in zip(pairs(18), pat, generic(avs, 18)) if a != b]
        feas = farkas_feasible(18, pat)
        real.append(feas)
        assert diff == [(0, 17)] and sum(avs) == 27
        say("  plus sheet x=%3d word %s: differing chord %s (sum 27 = t_17); actual pattern word-realisable: %s" % (x, avs, diff, feas))
    assert sum(real) == 1 and real[plus.index(437)]
    say("  -> at m = 18 the plus sheet gains exactly 11 non-generic chord patterns (only 437's actual pattern is word-realisable).")
    minus = [165, 309, 549]
    realm = []
    for x in minus:
        vals, avs = window(x, 13, -1)
        pat = chord_bits(vals, "DESC")
        diff = [(i, j) for (i, j), a, b in zip(pairs(13), pat, generic(avs, 13)) if a != b]
        feas = farkas_feasible(13, pat)
        realm.append(feas)
        assert diff == [(0, 12)] and sum(avs) == 19
        say("  minus sheet x=%3d word %s: differing chord %s (sum 19 = t_12); realisable: %s" % (x, avs, diff, feas))
    assert realm == [True, False, True]
    say("  -> at m = 13 the minus sheet gains exactly 1 non-generic pattern (309's).")
    say("In every crossing only the long chord (0, m-1) separates from the word prediction: the pinch is the single gate of the")
    say("whole window, and the sub-windows are generic (short-lag genericity).")
    say()


def main(save):
    say("# Five-paper bridge checks (S14, 2026-10-06)")
    say()
    section_A()
    section_B()
    section_C()
    say("ALL ASSERTIONS PASSED")
    if save:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.normpath(os.path.join(here, "..", "..", "05-knowledge", "results", "five_papers_bridges_20261006.out"))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(OUT) + "\n")
        print("saved", path)


if __name__ == "__main__":
    main("--save" in sys.argv)
