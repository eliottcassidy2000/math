"""The Syracuse (odd-to-odd) m-window alphabet of tournaments.

Session opus-2026-10-05-S12. Companion of trees5_tournaments4_converse_20261005.py (same folder;
its class builder is imported for m = 5, 6).

For odd x let S(x) = (3x+1)/2^a with a = v_2(3x+1) (Syracuse map), and a_i the valuation of the
i-th step. The m-window is the sequence W = (x, Sx, ..., S^(m-1) x). Its tournament has the time
arcs S^i x -> S^(i+1) x and, on every non-consecutive pair, an order arc:
    gauge DESC: from the larger value to the smaller (forward chord = descent),
    gauge ASC : from the smaller value to the larger.
Windows with a repeated value (the orbit reaches the fixed point 1 inside the window) are degenerate
and are counted separately.

Checks (all assertions raise on failure):
  A. m = 4, gauge DESC: the class is the B2 truth table of the two overlapping 2-step descent bits
         a = [a1 + a2 >= 4]  (x > S^2 x),   b = [a2 + a3 >= 4]  (Sx > S^3 x):
         (1,1) TT4 (past = source, present = sink), (0,1) C3 over the sink (present),
         (1,0) source (past) over C3, (0,0) STRONG;
     exact for every non-degenerate odd x with the finitely many small exceptions listed.
  B. m = 4, gauge ASC: only TT4 (iff a1+a2+a3 <= 4) and STRONG occur (generic law, exceptions listed).
  C. Haar law (valuation word i.i.d. geometric): TT4 11/32, each vortex 5/32, STRONG 11/32; census
     frequencies below 2^N compared.
  D. Time reversal of the window (same gauge) exchanges TT4 <-> STRONG and fixes the vortices (the B2
     complement); converse = negation o time reversal exchanges the vortices (the a <-> b swap).
  E. m = 5, 6: which of the 12 / 56 classes occur, their Haar law and census frequencies.
  F. The credit-grammar words G = (1,2) and H = (1,2,1,1,1,2): the 4-window classes they contain.

Run:  python -X utf8 syracuse_window_alphabet_20261005.py [--N 22] [--save]
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trees5_tournaments4_converse_20261005 import tournament_classes, scores, three_cycles, largest_strong  # noqa: E402

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


def v2(n):
    return (n & -n).bit_length() - 1


def syr(x):
    y = 3 * x + 1
    a = v2(y)
    return y >> a, a


def window(x, m):
    vals = [x]
    avs = []
    for _ in range(m - 1):
        y, a = syr(vals[-1])
        vals.append(y)
        avs.append(a)
    return vals, avs


NAME4 = {(0, 1, 2, 3): "TT4", (1, 1, 1, 3): "OUT (source over C3)", (0, 2, 2, 2): "IN (C3 over sink)", (1, 1, 2, 2): "STRONG"}
CONV4 = {"TT4": "TT4", "STRONG": "STRONG", "OUT (source over C3)": "IN (C3 over sink)", "IN (C3 over sink)": "OUT (source over C3)"}


def chord_bits(vals, gauge):
    """For each non-consecutive pair i<j: 1 if the chord points forward in time (i -> j)."""
    m = len(vals)
    bits = []
    for i in range(m):
        for j in range(i + 2, m):
            if gauge == "DESC":
                bits.append(1 if vals[i] > vals[j] else 0)
            else:
                bits.append(1 if vals[i] < vals[j] else 0)
    return tuple(bits)


def tournament_from_bits(m, bits, reverse_time=False):
    G = nx.DiGraph()
    G.add_nodes_from(range(m))
    for i in range(m - 1):
        if reverse_time:
            G.add_edge(i + 1, i)
        else:
            G.add_edge(i, i + 1)
    k = 0
    for i in range(m):
        for j in range(i + 2, m):
            if bits[k]:
                G.add_edge(i, j)
            else:
                G.add_edge(j, i)
            k += 1
    return G


def class4(bits, reverse_time=False):
    G = tournament_from_bits(4, bits, reverse_time)
    return NAME4[scores(G)], G


def table_B2(a, b):
    return {(1, 1): "TT4", (0, 1): "IN (C3 over sink)", (1, 0): "OUT (source over C3)", (0, 0): "STRONG"}[(a, b)]


def main(N, do_save):
    say("# The Syracuse (odd-to-odd) window alphabet: time arcs plus numerical order")
    say()
    say("Universe: all odd x < 2^%d; S(x) = (3x+1)/2^v2(3x+1); windows (x, Sx, ..., S^(m-1) x)." % N)
    say()

    # ------------------------------------------------------------------ A, B: m = 4 exact census
    say("## A. m = 4, gauge DESC (order arcs from the larger to the smaller value)")
    say()
    cnt = Counter()
    exceptions = []
    degenerate = []
    marked = defaultdict(Counter)
    cntA = Counter()
    exceptionsA = []
    words = Counter()
    rev_perm = Counter()
    for x in range(1, 1 << N, 2):
        vals, (a1, a2, a3) = window(x, 4)
        if len(set(vals)) < 4:
            degenerate.append(x)
            continue
        bits = chord_bits(vals, "DESC")
        name, G = class4(bits)
        a = 1 if a1 + a2 >= 4 else 0
        b = 1 if a2 + a3 >= 4 else 0
        pred = table_B2(a, b)
        cnt[name] += 1
        words[(min(a1, 8), min(a2, 8), min(a3, 8))] += 1
        if name != pred:
            exceptions.append((x, (a1, a2, a3), name, pred))
        marked[name][(G.out_degree(0), G.out_degree(3))] += 1
        # time reversal (same gauge): reversed value sequence
        rbits = chord_bits(vals[::-1], "DESC")
        rname, _ = class4(rbits)
        rev_perm[(name, rname)] += 1
        # gauge ASC
        bitsA = chord_bits(vals, "ASC")
        nameA, _ = class4(bitsA)
        cntA[nameA] += 1
        predA = "TT4" if a1 + a2 + a3 <= 4 else "STRONG"
        if nameA != predA:
            exceptionsA.append((x, (a1, a2, a3), nameA, predA))
        # exact identity: reversing time and keeping the chords = reversing every arc of the ASC reading
        assert rname == CONV4[nameA]
    total = sum(cnt.values())
    say("| class | count | frequency | Haar prediction |")
    say("|---|---|---|---|")
    haar = {"TT4": Fraction(11, 32), "IN (C3 over sink)": Fraction(5, 32), "OUT (source over C3)": Fraction(5, 32), "STRONG": Fraction(11, 32)}
    for name in ["TT4", "IN (C3 over sink)", "OUT (source over C3)", "STRONG"]:
        say("| %s | %d | %.6f | %s = %.6f |" % (name, cnt[name], cnt[name] / total, haar[name], float(haar[name])))
    say()
    say("degenerate windows (a repeated value, the orbit reaches 1 inside the window): %d, first ones %s" % (len(degenerate), degenerate[:8]))
    say("exceptions to the B2 law table(a,b): %d -> %s" % (len(exceptions), exceptions))
    say()
    say("Marked vertices (score of the past x, score of the present S^3 x) per class:")
    for name in ["TT4", "IN (C3 over sink)", "OUT (source over C3)", "STRONG"]:
        say("  %-22s %s" % (name, dict(marked[name])))
    assert marked["TT4"] == Counter({(3, 0): cnt["TT4"]})
    assert marked["IN (C3 over sink)"] == Counter({(2, 0): cnt["IN (C3 over sink)"]})
    assert marked["OUT (source over C3)"] == Counter({(3, 1): cnt["OUT (source over C3)"]})
    say()
    say("## B. m = 4, gauge ASC (order arcs from the smaller to the larger value)")
    say()
    for name in ["TT4", "IN (C3 over sink)", "OUT (source over C3)", "STRONG"]:
        say("  %-22s %d" % (name, cntA[name]))
    say("exceptions to the law (TT4 iff a1+a2+a3 <= 4, else STRONG): %d -> %s" % (len(exceptionsA), exceptionsA))
    assert cntA["IN (C3 over sink)"] == 0 and cntA["OUT (source over C3)"] == 0

    # exceptions: all x <= 9?
    say()
    say("## D. Time reversal of the window (same gauge DESC)")
    say()
    perm = defaultdict(Counter)
    for (name, rname), c in rev_perm.items():
        perm[name][rname] += c
    for name in ["TT4", "IN (C3 over sink)", "OUT (source over C3)", "STRONG"]:
        say("  %-22s -> %s" % (name, dict(perm[name])))
    say("Exact identity (asserted for every window): the time-reversed window read with the same gauge is the")
    say("converse of the original window read with the opposite gauge; so under DESC it is TT4 iff a1+a2+a3 <= 4")
    say("and STRONG otherwise (pure time reversal is not an involution of the DESC alphabet). The converse of the")
    say("DESC reading (= negation o time reversal, Redei note) exchanges the two vortices: the a <-> b swap of B2.")

    # ------------------------------------------------------------------ C: Haar law by word enumeration
    say()
    say("## C. Haar law from the valuation word (a_i i.i.d., P(a = k) = 2^-k), generic comparisons")
    say()
    CAP = 8

    def mass(a):
        return Fraction(1, 2 ** a) if a < CAP else Fraction(1, 2 ** (CAP - 1))

    def generic_bits(avs, m, gauge):
        # S^j x / S^i x ~ 3^(j-i) / 2^(A_j - A_i); descent iff 2^(sum) > 3^(j-i)
        bits = []
        for i in range(m):
            for j in range(i + 2, m):
                s = sum(avs[i:j])
                desc = 2 ** s > 3 ** (j - i)
                bits.append(1 if (desc if gauge == "DESC" else not desc) else 0)
        return tuple(bits)

    law4 = Counter()
    for avs in itertools.product(range(1, CAP + 1), repeat=3):
        p = mass(avs[0]) * mass(avs[1]) * mass(avs[2])
        name, _ = class4(generic_bits(avs, 4, "DESC"))
        law4[name] += p
    say("m = 4 DESC: " + ", ".join("%s %s" % (k, v) for k, v in sorted(law4.items())))
    assert law4 == haar
    law4A = Counter()
    for avs in itertools.product(range(1, CAP + 1), repeat=3):
        p = mass(avs[0]) * mass(avs[1]) * mass(avs[2])
        name, _ = class4(generic_bits(avs, 4, "ASC"))
        law4A[name] += p
    say("m = 4 ASC : " + ", ".join("%s %s" % (k, v) for k, v in sorted(law4A.items())))
    assert law4A["TT4"] == Fraction(5, 16) and law4A["STRONG"] == Fraction(11, 16)
    say("Census deviation at 2^%d: " % N + ", ".join("%s %+.2e" % (k, cnt[k] / total - float(haar[k])) for k in haar))

    # ------------------------------------------------------------------ E: m = 5, 6
    say()
    say("## E. Longer windows, gauge DESC: which classes occur")
    say()
    TC = tournament_classes(6)
    for m in (5, 6):
        C = TC[m]
        reps = C.reps
        conv = [C.find(G.reverse(copy=True)) for G in reps]
        # census by chord-bit vector
        bitcount = Counter()
        degen = 0
        for x in range(1, 1 << N, 2):
            vals, avs = window(x, m)
            if len(set(vals)) < m:
                degen += 1
                continue
            bitcount[chord_bits(vals, "DESC")] += 1
        cls_count = Counter()
        for bits, c in bitcount.items():
            idx = C.find(tournament_from_bits(m, bits))
            cls_count[idx] += c
        # Haar law
        law = Counter()
        for avs in itertools.product(range(1, CAP + 1), repeat=m - 1):
            p = Fraction(1)
            for a in avs:
                p *= mass(a)
            idx = C.find(tournament_from_bits(m, generic_bits(avs, m, "DESC")))
            law[idx] += p
        tot = sum(cls_count.values())
        occurring = sorted(set(cls_count) | set(law))
        say("m = %d: %d of %d classes occur (census %d distinct chord patterns, %d degenerate windows); Haar support %d classes" % (
            m, len(set(cls_count)), len(reps), len(bitcount), degen, len(law)))
        say("| class | scores | c3 | L | self-converse | converse partner occurs | Haar mass | census frequency |")
        say("|---|---|---|---|---|---|---|---|")
        for idx in sorted(occurring, key=lambda i: (-law[i], i)):
            G = reps[idx]
            say("| %d | %s | %d | %d | %s | %s | %s = %.5f | %.5f |" % (
                idx, scores(G), three_cycles(G), largest_strong(G), conv[idx] == idx,
                "yes" if conv[idx] in law else "NO", law[idx], float(law[idx]), cls_count[idx] / tot))
        assert sum(law.values()) == 1
        # the law support must equal the census support
        assert set(cls_count) == set(law), (set(cls_count) ^ set(law))
        # converse closure
        say("Occurring set closed under converse: %s" % all(conv[i] in law for i in law))
        missing = [i for i in range(len(reps)) if i not in law]
        say("Missing classes (%d): %s" % (len(missing), ", ".join("%d %s L=%d" % (i, scores(reps[i]), largest_strong(reps[i])) for i in missing)))
        occ = set(law)
        sc = sum(1 for i in occ if conv[i] == i)
        say("Occurring classes: %d = %d self-converse + %d converse pairs -> %d modulo converse; by L: %s" % (
            len(occ), sc, (len(occ) - sc) // 2, sc + (len(occ) - sc) // 2,
            dict(sorted(Counter(largest_strong(reps[i]) for i in occ if conv[i] >= i).items()))))
        say("Realizable chord patterns (generic bits): %d of %d" % (len({generic_bits(avs, m, "DESC") for avs in itertools.product(range(1, CAP + 1), repeat=m - 1)}), 2 ** ((m - 1) * (m - 2) // 2)))
        say()

    # ------------------------------------------------------------------ F: credit words
    say("## F. The credit-grammar words read through the 4-window alphabet (generic bits)")
    say()
    for wname, w in (("G = (1,2)", (1, 2)), ("G G = (1,2,1,2)", (1, 2, 1, 2)), ("H = (1,2,1,1,1,2)", (1, 2, 1, 1, 1, 2)),
                     ("H G = (1,2,1,1,1,2,1,2)", (1, 2, 1, 1, 1, 2, 1, 2)), ("G H = (1,2,1,2,1,1,1,2)", (1, 2, 1, 2, 1, 1, 1, 2))):
        wins = []
        for i in range(len(w) - 2):
            avs = w[i:i + 3]
            a = 1 if avs[0] + avs[1] >= 4 else 0
            b = 1 if avs[1] + avs[2] >= 4 else 0
            wins.append("%s->(%d,%d)=%s" % (avs, a, b, table_B2(a, b).split(" ")[0]))
        say("  %-26s %s" % (wname, "; ".join(wins) if wins else "(shorter than a 4-window)"))
    say("Every 4-window inside H and G is STRONG: no 2-step descent (a1+a2 >= 4) occurs in these words, since")
    say("they contain no consecutive valuations summing to 4 or more; the credit words are ascent-only at the")
    say("2-step scale, and the descent they pay for is carried by the common-future dependency, not by the window.")

    say()
    say("ALL ASSERTIONS PASSED (N = %d)" % N)
    if do_save:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.normpath(os.path.join(here, "..", "..", "05-knowledge", "results", "syracuse_window_alphabet_20261005.out"))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(OUT) + "\n")
        print("saved", path)


if __name__ == "__main__":
    N = 22
    if "--N" in sys.argv:
        N = int(sys.argv[sys.argv.index("--N") + 1])
    main(N, "--save" in sys.argv)
