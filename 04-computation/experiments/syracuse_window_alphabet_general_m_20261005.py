"""The Syracuse m-window alphabet for m = 4..8: exact densities by a capped-sum DP, census support,
classes occurring, the all-ascent and all-descent words, score separation, and the exotic twins.

Session opus-2026-10-05-S13. Companion of syracuse_window_alphabet_20261005.py (imported) and
trees5_tournaments4_converse_20261005.py (class builder, n <= 8).

Chord (i, j), j >= i+2, of the m-window is a descent iff A(i,j) = a_(i+1)+...+a_j >= t_(j-i), where
t_k = least s with 2^s > 3^k (= floor(k log2 3) + 1). Gauge DESC throughout (forward chord = descent).

DP: a word is processed letter by letter; the state is (capped suffix sums, bits decided so far) with
cap = t_(m-1); letters a >= cap are lumped with mass 2^-(cap-1). The pattern of a word depends only on
the capped letters, so the DP gives exact natural densities of every chord pattern (the capped word is a
union of residue classes of odd x modulo a power of two).

Run:  python -X utf8 syracuse_window_alphabet_general_m_20261005.py [--N 20] [--mmax 8] [--minus] [--save]
"""
from __future__ import annotations

import itertools
import os
import sys
import warnings
from collections import Counter, defaultdict
from fractions import Fraction

warnings.filterwarnings("ignore")
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trees5_tournaments4_converse_20261005 import tournament_classes, scores, three_cycles, largest_strong  # noqa: E402
from syracuse_window_alphabet_20261005 import v2, chord_bits, tournament_from_bits  # noqa: E402

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


def threshold(k):
    s = 0
    while 2 ** s <= 3 ** k:
        s += 1
    return s


def pairs(m):
    return [(i, j) for i in range(m) for j in range(i + 2, m)]


def generic_pattern(avs, m):
    return tuple(1 if sum(avs[i:j]) >= threshold(j - i) else 0 for i, j in pairs(m))


def dp_densities(m):
    """Exact densities of chord patterns (gauge DESC) for the m-window; returns {pattern: Fraction}."""
    cap = threshold(m - 1)
    letters = [(a, Fraction(1, 2 ** a) if a < cap else Fraction(1, 2 ** (cap - 1))) for a in range(1, cap + 1)]
    # state: (sums, bits) ; sums[i] = min(A(i, j), cap) for the current j; bits = dict pairs -> bit (as tuple in (j,i) order)
    states = {((), ()): Fraction(1)}
    for j in range(1, m):
        new = defaultdict(Fraction)
        for (sums, bits), mass in states.items():
            for a, pa in letters:
                nsums = tuple(min(s + a, cap) for s in sums) + (min(a, cap),)
                nbits = bits
                for i in range(0, j - 1):  # chords (i, j) with j - i >= 2
                    nbits = nbits + ((1 if nsums[i] >= threshold(j - i) else 0),)
                new[(nsums, nbits)] += mass * pa
        states = new
    # convert bits (ordered by j then i) to the (i, j) order of pairs(m)
    order = [(i, j) for j in range(2, m) for i in range(0, j - 1)]
    pos = {p: k for k, p in enumerate(order)}
    dens = defaultdict(Fraction)
    for (sums, bits), mass in states.items():
        pat = tuple(bits[pos[p]] for p in pairs(m))
        dens[pat] += mass
    assert sum(dens.values()) == 1
    return dict(dens)


def syr_plus(x):
    y = 3 * x + 1
    a = v2(y)
    return y >> a, a


def syr_minus(x):
    y = 3 * x - 1
    a = v2(y)
    return y >> a, a


def census(m, N, step):
    cnt = Counter()
    degen = 0
    for x in range(1, 1 << N, 2):
        vals = [x]
        for _ in range(m - 1):
            y, _a = step(vals[-1])
            vals.append(y)
        if len(set(vals)) < m:
            degen += 1
            continue
        cnt[chord_bits(vals, "DESC")] += 1
    return cnt, degen


def main(N, mmax, do_minus, do_save):
    say("# The Syracuse m-window alphabet, m = 4..%d (gauge DESC: forward chord = descent)" % mmax)
    say()
    say("thresholds t_k (least s with 2^s > 3^k): " + ", ".join("t_%d = %d" % (k, threshold(k)) for k in range(1, 9)))
    say("Sturmian differences t_(k+1) - t_k: " + ", ".join(str(threshold(k + 1) - threshold(k)) for k in range(1, 8)))
    say()
    TC = tournament_classes(mmax)
    summary = []
    strong_occ = {3: {TC[3].find(tournament_from_bits(3, (0,)))}}
    for m in range(4, mmax + 1):
        C = TC[m]
        reps = C.reps
        conv = [C.find(G.reverse(copy=True)) for G in reps]
        dens = dp_densities(m)
        cnt, degen = census(m, N, syr_plus)
        tot = sum(cnt.values())
        assert set(cnt) == set(dens), (m, len(set(cnt) ^ set(dens)))
        # classes
        pat_class = {pat: C.find(tournament_from_bits(m, pat)) for pat in dens}
        cls_dens = defaultdict(Fraction)
        cls_cnt = Counter()
        for pat, d in dens.items():
            cls_dens[pat_class[pat]] += d
        for pat, c in cnt.items():
            cls_cnt[pat_class[pat]] += c
        occ = set(cls_dens)
        sc = sum(1 for i in occ if conv[i] == i)
        pairs_ = (len(occ) - sc) // 2
        assert all(conv[i] in occ for i in occ)
        missing = [i for i in range(len(reps)) if i not in occ]
        Locc = Counter(largest_strong(reps[i]) for i in occ)
        Lmiss = Counter(largest_strong(reps[i]) for i in missing)
        strong_missing = sum(1 for i in missing if largest_strong(reps[i]) == m)
        regular_occ = [i for i in occ if len(set(scores(reps[i]))) == 1]
        score_sep = len({scores(reps[i]) for i in occ})
        all_desc = tuple([1] * len(pairs(m)))
        all_asc = tuple([0] * len(pairs(m)))
        asc_cls = pat_class[all_asc]
        max_dev = max(abs(float(cls_dens[i]) - cls_cnt[i] / tot) for i in occ)
        say("## m = %d  (n = %d vertices, %d classes)" % (m, m, len(reps)))
        say()
        say("chord patterns realised: %d of %d (DP support = census support below 2^%d; %d degenerate windows)" % (
            len(dens), 2 ** len(pairs(m)), N, degen))
        say("classes occurring: %d of %d = %d self-converse + %d converse pairs -> %d modulo converse" % (
            len(occ), len(reps), sc, pairs_, sc + pairs_))
        say("occurring by largest strong component L: %s" % dict(sorted(Locc.items())))
        say("missing: %d classes, by L: %s; strong among the missing: %d; regular tournament occurs: %s" % (
            len(missing), dict(sorted(Lmiss.items())), strong_missing, bool(regular_occ)))
        say("score-separated (distinct score sequences among occurring classes / occurring classes): %d / %d" % (score_sep, len(occ)))
        say("TT_%d density (all chords descend): %s = %.5f" % (m, dens[all_desc], float(dens[all_desc])))
        say("all-ascent class (every chord backward): scores %s, c3 = %d, L = %d, density %s = %.5f" % (
            scores(reps[asc_cls]), three_cycles(reps[asc_cls]), largest_strong(reps[asc_cls]), cls_dens[asc_cls], float(cls_dens[asc_cls])))
        say("largest class density: %s (scores %s); max |density - census frequency| = %.1e" % (
            max(cls_dens.values()), scores(reps[max(cls_dens, key=cls_dens.get)]), max_dev))
        # expected scores of the all-ascent class: (1,1,2,...,m-2,m-2)
        exp = tuple(sorted([1] + [k for k in range(1, m - 1)] + [m - 2]))
        assert scores(reps[asc_cls]) == exp, (scores(reps[asc_cls]), exp)
        if m % 2 == 1:
            assert not regular_occ
        summary.append((m, len(dens), len(occ), len(reps), sc, pairs_, len(missing), strong_missing, score_sep))
        # ---- structure checks -------------------------------------------------
        # (i) the all-ascent class has exactly m-k+1 cycles of each length k (its cycles are the intervals): Moon-extremal
        cyc = Counter(len(c) for c in nx.simple_cycles(reps[asc_cls]))
        assert all(cyc[k] == m - k + 1 for k in range(3, m + 1)), cyc
        say("all-ascent class cycle counts by length: %s = Moon bound m-k+1 at every length (interval cycles)" % dict(sorted(cyc.items())))
        # nearly transitive tournament (TT_m with the arc 0 -> m-1 reversed): same scores, occurs?
        Nt = nx.DiGraph(); Nt.add_nodes_from(range(m))
        for i in range(m):
            for j in range(i + 1, m):
                Nt.add_edge(i, j)
        Nt.remove_edge(0, m - 1); Nt.add_edge(m - 1, 0)
        nt_idx = C.find(Nt)
        cyc_nt = Counter(len(c) for c in nx.simple_cycles(Nt))
        say("nearly transitive tournament (TT_%d with 0->%d reversed): class %d, scores %s, cycles %s, occurs: %s, same class as all-ascent: %s" % (
            m, m - 1, nt_idx, scores(Nt), dict(sorted(cyc_nt.items())), nt_idx in occ, nt_idx == asc_cls))
        # (ii) TT_m density = transfer matrix over capped letters with all consecutive sums >= 4
        Mt = [[Fraction(0)] * 4 for _ in range(4)]
        mass4 = [Fraction(1, 2), Fraction(1, 4), Fraction(1, 8), Fraction(1, 8)]
        for a in range(4):
            for b in range(4):
                if (a + 1) + (b + 1) >= 4:
                    Mt[a][b] = mass4[b]
        vec = mass4[:]
        for _ in range(m - 2):
            vec = [sum(vec[a] * Mt[a][b] for a in range(4)) for b in range(4)]
        assert sum(vec) == dens[all_desc], (sum(vec), dens[all_desc])
        say("TT_%d density equals the transfer-matrix value for 'all consecutive valuation pairs sum to >= 4' (char. poly. 32x^3 - 16x^2 - 4x + 1, top root 0.6203)" % m)
        # (iii) condensation multiplicativity: a class occurs iff every strong block occurs as a strong window class of its size
        strong_occ[m] = {i for i in occ if largest_strong(reps[i]) == m}
        if m >= 5:
            ok = True
            for i in range(len(reps)):
                G = reps[i]
                comps = list(nx.strongly_connected_components(G))
                blocks_ok = True
                for comp in comps:
                    if len(comp) == 1:
                        continue
                    H = G.subgraph(comp).copy()
                    k = len(comp)
                    kidx = TC[k].find(nx.relabel_nodes(H, {v: t for t, v in enumerate(sorted(comp))}))
                    if kidx not in strong_occ[k]:
                        blocks_ok = False
                        break
                if blocks_ok != (i in occ):
                    ok = False
                    say("  multiplicativity FAILS at class %d (scores %s)" % (i, scores(G)))
                    break
            assert ok
            say("condensation multiplicativity holds: a class occurs iff each of its strong blocks is an occurring strong window class of its size (all %d classes checked)" % len(reps))
        if m <= 6:
            say("densities by class (index, scores, L, self-converse, density, census frequency):")
            for i in sorted(occ, key=lambda i: (-cls_dens[i], i)):
                say("  %3d %s L=%d sc=%s  %s = %.5f  | %.5f" % (i, scores(reps[i]), largest_strong(reps[i]), conv[i] == i, cls_dens[i], float(cls_dens[i]), cls_cnt[i] / tot))
        else:
            top = sorted(occ, key=lambda i: (-cls_dens[i], i))[:12]
            say("twelve densest classes (index, scores, L, self-converse, density):")
            for i in top:
                say("  %3d %s L=%d sc=%s  %s = %.5f" % (i, scores(reps[i]), largest_strong(reps[i]), conv[i] == i, cls_dens[i], float(cls_dens[i])))
            say("density of the smallest occurring class: %s" % min(cls_dens.values()))
        # exotic twins at m = 5: same score sequence, one occurring, one missing; Ryser 3-cycle reversal
        if m == 5:
            say()
            say("Exotic twins at m = 5 (same score sequence, different class):")
            by_scores = defaultdict(list)
            for i in range(len(reps)):
                by_scores[scores(reps[i])].append(i)
            for s, lst in sorted(by_scores.items()):
                if len(lst) < 2:
                    continue
                occ_here = [i for i in lst if i in occ]
                mis_here = [i for i in lst if i not in occ]
                say("  scores %s: classes %s, occurring %s, missing %s" % (s, lst, occ_here, mis_here))
                for i in occ_here:
                    G = reps[i]
                    for (u, v, w) in itertools.permutations(G.nodes, 3):
                        if u < v and u < w and G.has_edge(u, v) and G.has_edge(v, w) and G.has_edge(w, u):
                            H = G.copy()
                            H.remove_edges_from([(u, v), (v, w), (w, u)])
                            H.add_edges_from([(v, u), (w, v), (u, w)])
                            k = C.find(H)
                            if k in mis_here:
                                say("    reversing the 3-cycle (%d,%d,%d) of the occurring class %d gives the missing class %d (Ryser interchange)" % (u, v, w, i, k))
                assert all(len(occ_here) == 1 for _ in [0])
        say()
    say("## Summary")
    say()
    say("| m | patterns realised | classes occurring / all | self-converse | pairs | missing | strong among missing | score-separated |")
    say("|---|---|---|---|---|---|---|---|")
    for (m, npat, nocc, nall, sc, pr, nmiss, smiss, ssep) in summary:
        say("| %d | %d | %d / %d | %d | %d | %d | %d | %d / %d |" % (m, npat, nocc, nall, sc, pr, nmiss, smiss, ssep, nocc))
    say()
    say("realised-pattern counts: " + ", ".join(str(s[1]) for s in summary))
    say("strong occurring classes by size: " + ", ".join("%d: %d of %d" % (k, len(strong_occ[k]), sum(1 for G in TC[k].reps if largest_strong(G) == k)) for k in sorted(strong_occ)))
    # minus sheet comparison
    if do_minus:
        say()
        say("## The 3x-1 sheet (odd x, S^-(x) = (3x-1)/2^v), same thresholds")
        say()
        for m in range(4, min(mmax, 6) + 1):
            dens = dp_densities(m)
            cnt, degen = census(m, N, syr_minus)
            tot = sum(cnt.values())
            same = set(cnt) == set(dens)
            maxdev = max(abs(float(dens[p]) - cnt.get(p, 0) / tot) for p in dens)
            extra = [p for p in cnt if p not in dens]
            say("m = %d: census patterns %d, plus-sheet patterns %d, supports equal: %s, extra minus-sheet patterns: %d, %d degenerate windows, max |density - frequency| = %.1e" % (
                m, len(cnt), len(dens), same, len(extra), degen, maxdev))
    say()
    say("ALL ASSERTIONS PASSED (N = %d, mmax = %d)" % (N, mmax))
    if do_save:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.normpath(os.path.join(here, "..", "..", "05-knowledge", "results", "syracuse_window_alphabet_general_m_20261005.out"))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(OUT) + "\n")
        print("saved", path)


if __name__ == "__main__":
    N = 20
    mmax = 8
    if "--N" in sys.argv:
        N = int(sys.argv[sys.argv.index("--N") + 1])
    if "--mmax" in sys.argv:
        mmax = int(sys.argv[sys.argv.index("--mmax") + 1])
    main(N, mmax, "--minus" in sys.argv, "--save" in sys.argv)
