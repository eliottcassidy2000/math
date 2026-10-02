#!/usr/bin/env python3
"""Shaved tournaments II: the n = 9 case, an independent audit of n <= 8, the excess law, and the typed links to
the lonely-runner harmful relations and to Collatz rigidity (thread session, 2026-10-01).

A shaving of order n is a spanning oriented graph contained in every n-tournament; u(n) = max #arcs = C(n,2) - kappa(n).
Companion of 04-computation/experiments/shaved_tournaments_20261001.py (THM-4526).  Note:
05-knowledge/results/shaved_tournaments_proved_vs_conjectured_20261001.md.

Checks (nauty's gentourng/labelg/converseg/countg/pickg and gcc are required; without them the recorded results are
printed and the check is skipped):
  N1. class lists: T(n) = 456, 6880, 191536 for n = 7, 8, 9 (A000568); at n = 9 the largest automorphism group is 81,
      attained only by C3[C3] = Cay(Z9, {1,3,4,7})
  N2. Theorem A (Grunbaum 1971): the only classes without H_n (path + first->last arc) for 3 <= n <= 9 are C3 and
      C3[1,C3,1]  [--full: n = 10, all 9733056 classes]
  N3. u(7) = 9 and u(8) = 11 again, by an orderly search inside P7 and inside TT1 => P7 (a method independent of the
      first script): 54 and 2290 host orbits = 51 and 1617 classes; no 10-arc (n = 7) or 12-arc (n = 8) shaving
  N4. u(9) = 14, kappa(9) = 22: inside C3[C3] there are 64 host orbits of 14-arc shavings = 54 classes and no 15-arc
      shaving  [--full: also inside TT3[C3] (104 orbits) and C3[1,P7,1] (66 orbits): the same 54 classes, none at 15]
  N5. exhaustive embedding counts (plain DFS, no look-ahead pruning) re-verify every maximum shaving for n = 7, 8, 9
      against every class; odd-in-every-class census: 2 of 51 (a converse pair), 0 of 1617, 0 of 54
  N6. a second proof of u(9) <= 14 from the list of 54: every acyclic one-arc extension (1739 classes) is avoided by
      some class, found by brute force over all 9! bijections (shavings are hereditary); the blockers: C3[C3] avoids
      1018 of the 1739, the 15 most frequent blockers are all symmetric, and every extension is avoided by some
      symmetric class (|Aut| > 1) and by some rigid class; at n = 7, 8 every one-arc extension of a maximum shaving
      is avoided by a symmetric class, and some only by symmetric classes (at n = 8 one only by C3 => R5)
  N7. structure: every maximum shaving for n = 7, 8, 9 is rigid; the sets are converse-closed with 1, 15, 0
      self-converse classes; Hamiltonian paths 4, 5, 0; path + all span-3 arcs is avoided by 2, 17, 81 classes
  N8. the symmetry penalty |Aut S| * 2^u <= n! and the Burnside form; neither forces rigidity at n = 9
  N9. the excess law of HYP-3817, kappa(n) - ceil(log2 T(n)) = #{self-converse classes with |Aut| > n}: holds for
      n = 3..9 (n = 9: 22 - 18 = 4 = 4); self-converse counts by orbit counting agree with nauty for n <= 10; the
      substitution X -> X[C3] is injective, keeps self-converse, |Aut| >= 3^m; so the law FAILS for every n >= 21 with
      n = 0 or 2 (mod 3)
  N10. lazy caterer kappa = 1 + C(n-2,2): exact exactly for n = 3..6 and 9 among n <= 9; u(n) > 2n - 4 for all n >= 38
  N11. lonely-runner probe: primitive tight sets for k = 2..8 in boxes; the k = 4, 5 sporadic sets {1,3,4,7},
      {1,3,4,5,9} are the connection sets of C3[C3] and of the Paley tournament P11, the k = 7 ones are not circulant
Reproduce: python 04-computation/experiments/shaved_tournaments_n9_20261001.py [--full]
(default: 57 checks, about 8 minutes on 4 cores; --full: 62 checks, about 36 minutes.  The "full checks" counts
and candidate numbers in the output can depend on the parallel search order; the PASS lines do not.)
"""
import argparse
import math
import os
import shutil
import subprocess
import tempfile
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations

FAILS = []
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser()
ap.add_argument("--full", action="store_true", help="also n = 10 for Theorem A and two more hosts at n = 9")
ARGS = ap.parse_args()


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def tool(name):
    return shutil.which("nauty-" + name) or shutil.which(name)


GENT, LABELG, CONVG, COUNTG, PICKG = (tool(t) for t in ("gentourng", "labelg", "converseg", "countg", "pickg"))
GCC = shutil.which("gcc")
WORK = tempfile.mkdtemp(prefix="shv9_")
HAVE = all([GENT, LABELG, CONVG, COUNTG, PICKG, GCC])
EXE = os.path.join(WORK, "shv9")
if HAVE:
    subprocess.run([GCC, "-O2", "-fopenmp", "-o", EXE, os.path.join(HERE, "shaved_tournaments_n9_20261001.c")],
                   check=True)
else:
    print("SKIP: nauty (gentourng, labelg, converseg, countg, pickg) and gcc are needed; recorded results follow.")


def run(args, stdin=None):
    return subprocess.run(args, input=stdin, capture_output=True, text=True, check=True).stdout


# ------------------------------------------------------------------ tournaments as out-neighbour lists
def from_ascii(n, s):
    A = [[0] * n for _ in range(n)]
    b = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[b] == "1":
                A[i][j] = 1
            else:
                A[j][i] = 1
            b += 1
    return A


def to_ascii(A):
    n = len(A)
    return "".join("1" if A[i][j] else "0" for i in range(n) for j in range(i + 1, n))


def circulant(m, D):
    return [[1 if (j - i) % m in D else 0 for j in range(m)] for i in range(m)]


def lex(X, Y):
    """X[Y]: every vertex of X replaced by a copy of Y."""
    a, b = len(X), len(Y)
    M = [[0] * (a * b) for _ in range(a * b)]
    for x1 in range(a):
        for y1 in range(b):
            for x2 in range(a):
                for y2 in range(b):
                    if (x1, y1) != (x2, y2):
                        M[x1 * b + y1][x2 * b + y2] = X[x1][x2] if x1 != x2 else Y[y1][y2]
    return M


def blowup_middle(Y):
    """C3[1, Y, 1]: a => Y => b => a."""
    m = len(Y)
    n = m + 2
    M = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m):
            M[i + 1][j + 1] = Y[i][j]
        M[0][i + 1] = 1
        M[i + 1][n - 1] = 1
    M[n - 1][0] = 1
    return M


def source_sink(Y, source=True, sink=True):
    m = len(Y)
    n = m + source + sink
    M = [[0] * n for _ in range(n)]
    off = 1 if source else 0
    for i in range(m):
        for j in range(m):
            M[i + off][j + off] = Y[i][j]
    if source:
        for v in range(1, n):
            M[0][v] = 1
    if sink:
        for v in range(n - 1):
            M[v][n - 1] = 1
    return M


C3 = circulant(3, {1})
TT1 = [[0]]
TT3 = [[1 if j > i else 0 for j in range(3)] for i in range(3)]
P7 = circulant(7, {1, 2, 4})


def arcs_of(A):
    n = len(A)
    return [(u, v) for u in range(n) for v in range(n) if A[u][v]]


def digraph6(n, arcs):
    bits = [0] * (n * n)
    for u, v in arcs:
        bits[u * n + v] = 1
    while len(bits) % 6:
        bits.append(0)
    s = "&" + chr(63 + n)
    for k in range(0, len(bits), 6):
        val = 0
        for t in range(6):
            val = (val << 1) | bits[k + t]
        s += chr(63 + val)
    return s


def canon(n, arc_lists):
    if not arc_lists:
        return []
    return run([LABELG, "-q"], "\n".join(digraph6(n, a) for a in arc_lists) + "\n").split()


def group_sizes(d6_list):
    out = run([COUNTG, "-q", "--a"], "\n".join(d6_list) + "\n")
    hist = Counter()
    for line in out.splitlines():
        line = line.strip()
        if "groupsize=" in line:
            cnt = int(line.split()[0])
            hist[int(line.split("groupsize=")[1].split()[0])] += cnt
    return hist


def linext(n, arcs):
    pred = [0] * n
    for a, b in arcs:
        pred[b] |= 1 << a

    @lru_cache(None)
    def f(m):
        if m == (1 << n) - 1:
            return 1
        return sum(f(m | 1 << v) for v in range(n) if not m >> v & 1 and pred[v] & ~m == 0)
    return f(0)


def acyclic(n, arcs):
    return linext(n, arcs) > 0


def bipartite(n, arcs):
    adj = [set() for _ in range(n)]
    for a, b in arcs:
        adj[a].add(b)
        adj[b].add(a)
    col = [-1] * n
    for s in range(n):
        if col[s] >= 0:
            continue
        col[s] = 0
        st = [s]
        while st:
            v = st.pop()
            for w in adj[v]:
                if col[w] < 0:
                    col[w] = col[v] ^ 1
                    st.append(w)
                elif col[w] == col[v]:
                    return False
    return True


def weakly_connected(n, arcs):
    adj = [set() for _ in range(n)]
    for a, b in arcs:
        adj[a].add(b)
        adj[b].add(a)
    seen, st = {0}, [0]
    while st:
        v = st.pop()
        for w in adj[v] - seen:
            seen.add(w)
            st.append(w)
    return len(seen) == n


def write_lines(name, lines):
    p = os.path.join(WORK, name)
    with open(p, "w") as f:
        f.write("\n".join(lines) + "\n")
    return p


# ------------------------------------------------------------------ N1 class lists
print("\nN1. class lists")
CLASSES = {}
RECORDED_T = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536, 10: 9733056}
if HAVE:
    for n in range(3, 10):
        p = os.path.join(WORK, f"t{n}.txt")
        with open(p, "w") as f:
            subprocess.run([GENT, "-q", str(n)], stdout=f, check=True)
        CLASSES[n] = p
        cnt = sum(1 for _ in open(p))
        check(cnt == RECORDED_T[n], f"gentourng n = {n}: {cnt} classes (A000568: {RECORDED_T[n]})")
    d9 = run([GENT, "-q", "-z", "9"])
    hist = group_sizes(d9.split())
    print("   |Aut| distribution at n = 9:", dict(sorted(hist.items())))
    # |Aut| of every class, in class order (gentourng lists the same classes in the same order with or without -z)
    AUT = {}
    for n in (7, 8, 9):
        dz = d9 if n == 9 else run([GENT, "-q", "-z", str(n)])
        AUT[n] = [int(line.split("groupsize=")[1]) for line in run([COUNTG, "-q", "-V", "--a"], dz).splitlines()
                  if "groupsize=" in line]
    check(all(len(AUT[n]) == RECORDED_T[n] for n in (7, 8, 9)) and Counter(AUT[9]) == hist,
          "n = 7, 8, 9: |Aut| read for every class")
    big = run([PICKG, "-q", "-a81:"], d9).split()
    c3c3 = canon(9, [arcs_of(lex(C3, C3))])[0]
    cay = canon(9, [arcs_of(circulant(9, {1, 3, 4, 7}))])[0]
    check(max(hist) == 81 and hist[81] == 1 and len(big) == 1 and run([LABELG, "-q"], big[0] + "\n").strip() == c3c3
          and cay == c3c3, "n = 9: the largest group order is 81, attained only by C3[C3] = Cay(Z9, {1,3,4,7})")
else:
    print("   recorded: T(n) = 456, 6880, 191536 (n = 7, 8, 9); max |Aut| at n = 9 is 81, only C3[C3] = Cay(Z9,{1,3,4,7})")

# ------------------------------------------------------------------ N2 Theorem A
print("\nN2. Theorem A: H_n = Hamiltonian path + first->last arc")
if HAVE:
    for n in range(3, 10):
        out = run([EXE, "thma", str(n)], open(CLASSES[n]).read())
        last = out.strip().splitlines()[-1]
        bad = int(last.split("=")[-1])
        exc = [line.split(":")[1].strip() for line in out.splitlines() if line.startswith("  no copy")]
        expect = {3: 1, 5: 1}.get(n, 0)
        ok = bad == expect
        if n == 3:
            ok = ok and canon(3, [arcs_of(from_ascii(3, exc[0]))]) == canon(3, [arcs_of(C3)])
        if n == 5:
            ok = ok and canon(5, [arcs_of(from_ascii(5, exc[0]))]) == canon(5, [arcs_of(blowup_middle(C3))])
        check(ok, f"n = {n}: classes without H_n = {bad}" + (" (C3)" if n == 3 else " (C3[1,C3,1])" if n == 5 else ""))
    if ARGS.full:
        p = subprocess.Popen([GENT, "-q", "10"], stdout=subprocess.PIPE)
        out = subprocess.run([EXE, "thma", "10"], stdin=p.stdout, capture_output=True, text=True, check=True).stdout
        p.wait()
        last = out.strip().splitlines()[-1]
        check("tournaments=9733056" in last and last.endswith("= 0"), "n = 10: " + last)
    else:
        print("   n = 10 skipped (--full); recorded: all 9733056 classes contain H_10")
else:
    print("   recorded: only C3 (n = 3) and C3[1,C3,1] (n = 5) for n <= 10")


# ------------------------------------------------------------------ orderly searches
def orderly(n, host, e, killfrom=None, seedk=4):
    hp = write_lines(f"host_{n}_{e}.txt", [to_ascii(host)])
    op = os.path.join(WORK, f"orbits_{n}_{e}.txt")
    args = [EXE, "orderly", str(n), hp, CLASSES[n], str(e), op]
    if killfrom is not None:
        args += [str(killfrom), str(seedk)]
    line = run(args).strip().splitlines()[-1]
    HA = arcs_of(host)
    masks = [int(x) for x in open(op) if x.strip()]
    lists = [[HA[b] for b in range(len(HA)) if m >> b & 1] for m in masks]
    return line, lists


def classes_of(n, lists):
    rep = {}
    for c, l in zip(canon(n, lists), lists):
        rep.setdefault(c, l)
    return rep


MAX = {}
print("\nN3. u(7) = 9 and u(8) = 11 by orderly search inside a host")
if HAVE:
    for n, host, e, orbits, ncls, name in ((7, P7, 9, 54, 51, "P7"), (8, source_sink(P7, True, False), 11, 2290, 1617,
                                                                          "TT1 => P7")):
        line, lists = orderly(n, host, e)
        print("  ", line)
        rep = classes_of(n, lists)
        check(len(lists) == orbits and len(rep) == ncls, f"n = {n}, host {name}: {len(lists)} orbits of {e}-arc "
              f"shavings = {len(rep)} classes")
        MAX[n] = rep
        line, lists = orderly(n, host, e + 1)
        print("  ", line)
        check(len(lists) == 0, f"n = {n}: no {e + 1}-arc shaving, so u({n}) = {e} and kappa({n}) = {n * (n - 1) // 2 - e}")
else:
    print("   recorded: n = 7: 54 P7-orbits = 51 classes, none at 10; n = 8: 2290 orbits = 1617 classes, none at 12")

print("\nN4. u(9) = 14")
if HAVE:
    hosts = [("C3[C3]", lex(C3, C3), 64)]
    if ARGS.full:
        hosts += [("TT3[C3]", lex(TT3, C3), 104), ("C3[1,P7,1]", blowup_middle(P7), 66)]
    reps9 = []
    for name, host, orbits in hosts:
        line, lists = orderly(9, host, 14)
        print("  ", line)
        rep = classes_of(9, lists)
        check(len(lists) == orbits and len(rep) == 54, f"host {name}: {len(lists)} orbits of 14-arc shavings = "
              f"{len(rep)} classes")
        reps9.append(rep)
        line, lists = orderly(9, host, 15)
        print("  ", line)
        check(len(lists) == 0, f"host {name}: no 15-arc shaving")
    check(all(set(r) == set(reps9[0]) for r in reps9), ("all three hosts give the same 54 classes; " if len(reps9) > 1
          else "") + "u(9) = 14, kappa(9) = 22")
    MAX[9] = reps9[0]
    if not ARGS.full:
        print("   hosts TT3[C3] and C3[1,P7,1] skipped (--full); recorded: 104 and 66 orbits, the same 54 classes, "
              "none at 15")
else:
    print("   recorded: C3[C3] 64 orbits, TT3[C3] 104, C3[1,P7,1] 66: the same 54 classes; none at 15 arcs")

# ------------------------------------------------------------------ N5 exhaustive counts and parity
print("\nN5. exhaustive embedding counts (no pruning) and parity")
if HAVE:
    for n, odd_expect in ((7, 2), (8, 0), (9, 0)):
        reps = list(MAX[n].values())
        cp = write_lines(f"max{n}.txt", [" ".join(f"{a} {b}" for a, b in l) for l in reps])
        out = run([EXE, "count", str(n), CLASSES[n], cp])
        print("  ", out.strip().replace("\n", "\n   "))
        first = out.splitlines()[0]
        allin = int(first.split("embedded in every class = ")[1].split()[0])
        odd = [int(line.split()[-1]) for line in out.splitlines() if "odd in every class: candidate" in line]
        ok = allin == len(reps) and len(odd) == odd_expect
        if n == 7 and len(odd) == 2:
            a, b = reps[odd[0]], reps[odd[1]]
            ok = ok and canon(7, [[(y, x) for x, y in a]]) == canon(7, [b])
        check(ok, f"n = {n}: all {len(reps)} maximum shavings embed in every class; odd count in every class: "
              f"{len(odd)}" + (" (a converse pair)" if n == 7 else ""))
else:
    print("   recorded: every maximum shaving embeds in every class (n = 7, 8, 9); odd everywhere: 2 (a converse pair), 0, 0")

# ------------------------------------------------------------------ N6 one-arc extensions at n = 9
print("\nN6. one-arc extensions of the 54 maximum 9-shavings (brute force over 9! bijections)")
if HAVE:
    ext = []
    for l in MAX[9].values():
        s = set(l)
        for a in range(9):
            for b in range(9):
                if a != b and (a, b) not in s and (b, a) not in s and acyclic(9, l + [(a, b)]):
                    ext.append(l + [(a, b)])
    erep = classes_of(9, ext)
    cp = write_lines("ext15.txt", [" ".join(f"{a} {b}" for a, b in l) for l in erep.values()])
    out = run([EXE, "bruteext", "9", CLASSES[9], cp]).strip().splitlines()[-1]
    print("  ", out)
    check(len(erep) == 1739 and out.endswith("= 0"), f"{len(erep)} acyclic 15-arc extension classes, none in every "
          "class: u(9) <= 14 again (shavings are hereditary)")
    ap9 = write_lines("aut9.txt", [str(a) for a in AUT[9]])
    out = run([EXE, "blockers", "9", CLASSES[9], cp, ap9])
    print("  ", out.strip().replace("\n", "\n   "))
    lines9 = [l.strip() for l in open(CLASSES[9]) if l.strip()]
    top = [int(l.split("class ")[1].split()[0]) for l in out.splitlines() if "top blocker" in l]
    greedy = [int(l.split("class ")[1].split()[0]) for l in out.splitlines() if "greedy" in l and "class" in l]
    d6_of = lambda i: digraph6(9, arcs_of(from_ascii(9, lines9[i])))
    auts = lambda idx: [max(group_sizes([d6_of(i)])) for i in idx]   # recomputed per class, not read from AUT[9]
    top_aut, greedy_aut = auts(top), auts(greedy)
    print("   |Aut| of the 15 most frequent blockers:", top_aut)
    print("   |Aut| along the greedy blocking set:", greedy_aut)
    first = [int(l.split("avoids")[1]) for l in out.splitlines() if "top blocker 1:" in l][0]
    num = lambda text, key: int(text.split(key)[1].split()[0])
    check("candidates avoided by no class = 0" in out and canon(9, [arcs_of(from_ascii(9, lines9[top[0]]))])[0] == c3c3
          and first == 1018 and min(top_aut) > 1 and [AUT[9][i] for i in top + greedy] == top_aut + greedy_aut,
          "every extension has a blocker; the most frequent one is C3[C3] (it avoids 1018 of the 1739); the 15 most "
          "frequent blockers all have |Aut| > 1 (188337 of the 191536 classes are rigid)")
    check(num(out, "only by rigid classes = ") == 0 and num(out, "only by symmetric classes = ") == 0
          and num(out, "fewest classes avoiding one candidate = ") == 27,
          "n = 9: every extension is avoided by at least 27 classes, among them a symmetric one and a rigid one: "
          "symmetric classes alone block all 1739, and so do rigid classes alone")
    # the same question at n = 7 and n = 8: the one-arc extensions of the 51 and of the 1617 maximum shavings
    C3_R5 = [[0] * 8 for _ in range(8)]                     # C3 => R5: a 3-cycle beating a rotational 5-tournament
    for a in range(8):
        for b in range(8):
            if a < 3 and b < 3:
                C3_R5[a][b] = int((b - a) % 3 == 1)
            elif a >= 3 and b >= 3:
                C3_R5[a][b] = int((b - a) % 5 in (1, 2))
            else:
                C3_R5[a][b] = int(a < 3)
    for n in (7, 8):
        extn = []
        for l in MAX[n].values():
            sn = set(l)
            for a in range(n):
                for b in range(n):
                    if a != b and (a, b) not in sn and (b, a) not in sn and acyclic(n, l + [(a, b)]):
                        extn.append(l + [(a, b)])
        erepn = classes_of(n, extn)
        cpn = write_lines(f"ext_{n}.txt", [" ".join(f"{a} {b}" for a, b in l) for l in erepn.values()])
        apn = write_lines(f"aut{n}.txt", [str(a) for a in AUT[n]])
        outn = run([EXE, "blockers", str(n), CLASSES[n], cpn, apn])
        print("  ", "\n   ".join(outn.strip().splitlines()[:4]))
        ok = ("candidates avoided by no class = 0" in outn and num(outn, "only by rigid classes = ") == 0
              and num(outn, "only by symmetric classes = ") > 0)
        if n == 8:   # one extension is avoided by a single class, C3 => R5 (|Aut| = 15)
            linesn = [l.strip() for l in open(CLASSES[n]) if l.strip()]
            sole = [int(x) for x in outn.split("avoided by class")[1].split("\n")[0].split()]
            ok = ok and num(outn, "fewest classes avoiding one candidate = ") == 1 and len(sole) == 1 \
                and canon(8, [arcs_of(from_ascii(8, linesn[sole[0]]))]) == canon(8, [arcs_of(C3_R5)]) \
                and AUT[8][sole[0]] == 15
        check(ok, f"n = {n}: every one of the {len(erepn)} one-arc extension classes is avoided by a symmetric class, "
              "and some are avoided only by symmetric classes" + (" (one only by C3 => R5)" if n == 8 else ""))
else:
    print("   recorded: 1739 extension classes, none contained in every class; the most frequent blocker is C3[C3] "
          "(1018 of 1739), the top 15 blockers are all symmetric; every extension is avoided by a symmetric class and "
          "by a rigid class; at n = 7, 8 every extension is avoided by a symmetric class and some only by symmetric "
          "classes (at n = 8 one only by C3 => R5); a greedy blocking set has 10 classes with "
          "|Aut| = 81, 21, 21, 21, 9, 15, 9, 9, 7, 1")

# ------------------------------------------------------------------ N7 structure
print("\nN7. structure of the maximum shavings")
if HAVE:
    expect = {7: (1, 4), 8: (15, 5), 9: (0, 0)}
    for n in (7, 8, 9):
        rep = MAX[n]
        keys = list(rep)
        hist = group_sizes(keys)
        conv = canon(n, [[(b, a) for a, b in l] for l in rep.values()])
        closed = all(c in rep for c in conv)
        sc = sum(1 for k, c in zip(keys, conv) if k == c)
        hp = sum(1 for l in rep.values() if linext(n, l) == 1)
        bip = sum(1 for l in rep.values() if bipartite(n, l))
        con = sum(1 for l in rep.values() if weakly_connected(n, l))
        check(hist == Counter({1: len(rep)}) and closed and (sc, hp) == expect[n],
              f"n = {n}: {len(rep)} classes, all rigid (|Aut| = 1), converse-closed, self-converse {sc}, "
              f"with a Hamiltonian path {hp}, bipartite {bip}, weakly connected {con}")
    for n, avoid in ((7, 2), (8, 17), (9, 81)):
        ref = [(i, i + 1) for i in range(n - 1)] + [(i, i + 3) for i in range(n - 3)]
        cp = write_lines(f"span3_{n}.txt", [" ".join(f"{a} {b}" for a, b in ref)])
        out = run([EXE, "count", str(n), CLASSES[n], cp])
        z = int(out.strip().splitlines()[-1].split("classes avoiding it")[1])
        check(z == avoid, f"n = {n}: path + all span-3 arcs ({len(ref)} arcs) is avoided "
              f"by {z} classes, so it is not a shaving")
else:
    print("   recorded: all rigid; self-converse 1, 15, 0; Hamiltonian path 4, 5, 0; bipartite 8, 27, 0; path + "
          "span-3 arcs avoided by 2, 17, 81 classes")

# ------------------------------------------------------------------ N8 symmetry penalty
print("\nN8. the symmetry penalty")
# n!/|Aut S| labelled copies, each in 2^(C(n,2)-u) labelled tournaments, must cover all 2^C(n,2) of them
for n, u in ((7, 9), (8, 11), (9, 14)):
    print(f"   n = {n}, u = {u}: a {u}-arc shaving has |Aut S| <= n!/2^u = {math.factorial(n) / 2 ** u:.2f}")
k9 = 36 - 14
check(2 ** k9 // 2 > 191536, "Burnside: a 9-shaving with |Aut S| = 2 still has >= 2^21 > T(9) completion orbits, so "
      "counting alone does not force rigidity")


# ------------------------------------------------------------------ N9 the excess law
def odd_partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n if n % 2 else n - 1
    if n == 0:
        yield []
        return
    for k in range(maxpart, 0, -2):
        for m in range(1, n // k + 1):
            for rest in odd_partitions(n - m * k, k - 2):
                yield [(k, m)] + rest


def T_count(n):
    """Davis-Polya: the number of n-tournaments (A000568)."""
    if n == 0:
        return 1
    tot = Fraction(0)
    for lam in odd_partitions(n):
        e = sum(k * m * (m - 1) // 2 + m * (k - 1) // 2 for k, m in lam)
        e += sum(m1 * m2 * math.gcd(k1, k2) for (k1, m1), (k2, m2) in combinations(lam, 2))
        z = 1
        for k, m in lam:
            z *= k ** m * math.factorial(m)
        tot += Fraction(2 ** e, z)
    assert tot.denominator == 1
    return int(tot)


def partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield []
        return
    for k in range(min(n, maxp), 0, -1):
        for rest in partitions(n - k, k):
            yield [k] + rest


def anti_fix(lam):
    """#labelled tournaments with a fixed permutation of cycle type lam as anti-automorphism (u->v iff s(v)->s(u))."""
    n = sum(lam)
    sigma, start = [], 0
    for c in lam:
        sigma += [start + (i + 1) % c for i in range(c)]
        start += c
    seen, orbits = set(), 0
    for u in range(n):
        for v in range(u + 1, n):
            if (u, v) in seen:
                continue
            a, b = u, v
            while True:
                a, b = sigma[b], sigma[a]
                seen.add((min(a, b), max(a, b)))
                if (min(a, b), max(a, b)) == (u, v):
                    if (a, b) != (u, v):
                        return 0
                    break
            orbits += 1
    return 2 ** orbits


def SC(n):
    """number of self-converse n-tournaments (A002785) by orbit counting."""
    tot = Fraction(0)
    for lam in partitions(n):
        f = anti_fix(lam)
        if f:
            z = 1
            for k in set(lam):
                z *= k ** lam.count(k) * math.factorial(lam.count(k))
            tot += Fraction(f, z)
    assert tot.denominator == 1
    return int(tot)


print("\nN9. the excess law (HYP-3817): kappa(n) - ceil(log2 T(n)) = #{SC classes with |Aut| > n}")
NMAX = 60
T = {n: T_count(n) for n in range(1, NMAX + 1)}
check([T[n] for n in range(3, 11)] == [RECORDED_T[n] for n in range(3, 11)], "Davis-Polya T(n) = A000568 for n <= 10")
clog2 = lambda x: (x - 1).bit_length()
U_ = {n: n * (n - 1) // 2 - clog2(T[n]) for n in T}
u_exact = {1: 0, 2: 1, 3: 2, 4: 4, 5: 6, 6: 8, 7: 9, 8: 11, 9: 14}
excess = {n: n * (n - 1) // 2 - u_exact[n] - clog2(T[n]) for n in range(3, 10)}
symsc_rec = {3: 0, 4: 0, 5: 0, 6: 1, 7: 3, 8: 4, 9: 4, 10: 1}
if HAVE:
    symsc = {}
    for n in range(3, 10):
        d6 = run([GENT, "-q", "-z", str(n)])
        big = run([PICKG, "-q", f"-a{n + 1}:"], d6).split()
        if big:
            c1 = run([LABELG, "-q"], "\n".join(big) + "\n").split()
            c2 = run([LABELG, "-q"], run([CONVG, "-q"], "\n".join(big) + "\n")).split()
            symsc[n] = sum(1 for x, y in zip(c1, c2) if x == y)
        else:
            symsc[n] = 0
else:
    symsc = {n: symsc_rec[n] for n in range(3, 10)}
print("   excess(n), n = 3..9:", [excess[n] for n in range(3, 10)])
print("   #{SC, |Aut| > n}:  ", [symsc[n] for n in range(3, 10)])
check(all(excess[n] == symsc[n] for n in range(3, 10)), "the law holds for n = 3..9 (n = 9: 22 - 18 = 4)")
print(f"   prediction at n = 10 (recorded #{{SC, |Aut| > 10}} = {symsc_rec[10]}): kappa(10) = "
      f"{clog2(T[10]) + symsc_rec[10]}, u(10) = {45 - clog2(T[10]) - symsc_rec[10]} (counting bound U(10) = {U_[10]})")
sc_known = [1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784]
check([SC(n) for n in range(1, 11)] == sc_known, "self-converse counts by orbit counting = A002785 for n <= 10")
if HAVE:
    ok = True
    for m in range(2, 7):
        Xs = [from_ascii(m, s) for s in run([GENT, "-q", str(m)]).split()]
        cX = canon(m, [arcs_of(X) for X in Xs])
        cXc = canon(m, [[(b, a) for a, b in arcs_of(X)] for X in Xs])
        sub = [lex(X, C3) for X in Xs]
        cs = canon(3 * m, [arcs_of(Y) for Y in sub])
        csc = canon(3 * m, [[(b, a) for a, b in arcs_of(Y)] for Y in sub])
        ok = ok and len(set(cs)) == len(Xs)                                       # injective on classes
        ok = ok and all((x == xc) <= (y == yc) for x, xc, y, yc in zip(cX, cXc, cs, csc))   # SC is kept
        ok = ok and min(group_sizes(cs)) >= 3 ** m                                # |Aut| >= 3^m
        if m <= 5:
            ss = [source_sink(Y) for Y in sub]
            css = canon(3 * m + 2, [arcs_of(Y) for Y in ss])
            cssc = canon(3 * m + 2, [[(b, a) for a, b in arcs_of(Y)] for Y in ss])
            ok = ok and len(set(css)) == len(Xs) and min(group_sizes(css)) >= 3 ** m
            ok = ok and all((x == xc) <= (y == yc) for x, xc, y, yc in zip(cX, cXc, css, cssc))
    check(ok, "X -> X[C3] (m <= 6) and X -> TT1 => X[C3] => TT1 (m <= 5): injective on classes, keep self-converse, "
          "|Aut| >= 3^m")
L_ = dict(u_exact)
for n in range(10, NMAX + 1):
    best = max(L_[a] + L_[n - a] for a in range(1, n))                 # superadditivity
    kk = 1 + int(math.floor(math.log2(n)))                             # Erdos-Moser: TT_kk in every n-tournament
    for k in range(2, kk + 1):
        best = max(best, k * (k - 1) // 2 + L_[n - k])
    L_[n] = best
fails = []
print("    n  U(n)  L(n)  U-L  lower bound on #{SC, |Aut| > n}")
for n in range(10, NMAX + 1):
    lb = SC(n // 3) if n % 3 == 0 else SC((n - 2) // 3) if n % 3 == 2 else None
    bad = lb is not None and lb > U_[n] - L_[n]
    if bad:
        fails.append(n)
    if n <= 33 or bad and n % 9 == 0:
        print(f"   {n:3d} {U_[n]:5d} {L_[n]:5d} {U_[n] - L_[n]:4d}  {lb if lb is not None else '-'}"
              + ("   law FAILS" if bad else ""))
check(fails == [n for n in range(21, NMAX + 1) if n % 3 != 1],
      f"the law fails at every n in [21, {NMAX}] with n = 0, 2 (mod 3), and at no smaller n by this construction")
# general m: an involutive anti-automorphism gives #SC(m) >= 2^((C(m,2) + floor(m/2))/2) / m!
check(all(anti_fix([2] * (m // 2) + [1] * (m % 2)) == 2 ** ((m * (m - 1) // 2 + m // 2) // 2) for m in range(2, 13)),
      "an involution with <= 1 fixed point is an anti-automorphism of 2^((C(m,2) + floor(m/2))/2) labelled tournaments")
check(all(Fraction(2 ** ((m * (m - 1) // 2 + m // 2) // 2), math.factorial(m)) > math.log2(math.factorial(3 * m + 2))
          for m in range(15, 400)),
      "for 15 <= m < 400 that bound exceeds log2((3m+2)!) >= U(3m+2) >= U(3m), so the law fails at every n >= 21 with "
      "n = 0, 2 (mod 3) (beyond 400 the gap only grows)")

# ------------------------------------------------------------------ N10 lazy caterer
print("\nN10. the lazy-caterer formula kappa(n) = 1 + C(n-2, 2)")
lc = [n for n in range(3, 10) if n * (n - 1) // 2 - u_exact[n] == 1 + (n - 2) * (n - 3) // 2]
check(lc == [3, 4, 5, 6, 9], f"exact for n in {lc} among 3 <= n <= 9 (above it at n = 7, 8)")
check(all(L_[n] > 2 * n - 4 for n in range(38, NMAX + 1)), "u(n) >= L(n) > 2n - 4 for 38 <= n <= 60; with "
      "u(n) >= C(k,2) + u(n-k), k = 1 + floor(log2 n) >= 6, C(k,2) >= 2k, this holds for all n >= 38")

# ------------------------------------------------------------------ N11 lonely-runner probe
print("\nN11. lonely-runner tight sets and circulant tournaments on Z_(2k+1)")
if HAVE:
    boxes = ((2, 30), (3, 30), (4, 40), (5, 30), (6, 26), (7, 26), (8, 22))
    rows = {}
    for k, B in boxes:
        out = run([EXE, "lrctight", str(k), str(B)])
        print("  ", out.strip().replace("\n", "\n   "))
        rows[k] = [line for line in out.splitlines() if line.startswith("  tight")]
        last = out.strip().splitlines()[-1]
        check(last.endswith("below 1/(k+1)=0"), f"k = {k}, speeds <= {B}: no set below 1/(k+1)")
    sporadic = {k: [r.split(":")[1].split("circulant")[0].strip() for r in rows[k]
                    if r.split(":")[1].split("circulant")[0].strip() != " ".join(map(str, range(1, k + 1)))]
                for k in rows}
    check(sporadic == {2: [], 3: [], 4: ["1 3 4 7"], 5: ["1 3 4 5 9"], 6: [], 7: ["1 2 3 4 5 7 12", "1 4 5 6 7 11 13"],
                       8: []}, "primitive tight sets besides {1..k}: {1,3,4,7}, {1,3,4,5,9}, and two at k = 7")
    check(all(r.endswith("yes") for k in (4, 5) for r in rows[k]) and all(r.endswith("no") for r in rows[7][1:]),
          "{1,3,4,7} and {1,3,4,5,9} are circulant connection sets mod 9 and 11 (C3[C3], Paley P11); the k = 7 "
          "sporadic sets are not (3 + 12 = 4 + 11 = 15)")
    check(canon(11, [arcs_of(circulant(11, {1, 3, 4, 5, 9}))]) == canon(11, [arcs_of(circulant(11, {
        (x * x) % 11 for x in range(1, 11)}))]), "{1,3,4,5,9} = QR_11")
else:
    print("   recorded: sporadic tight sets {1,3,4,7} (k = 4, = C3[C3]), {1,3,4,5,9} (k = 5, = P11), and at k = 7 "
          "{1,2,3,4,5,7,12}, {1,4,5,6,7,11,13} (not circulant)")

shutil.rmtree(WORK, ignore_errors=True)
print("\nALL CHECKS PASSED" if not FAILS else f"\n{len(FAILS)} FAILURES: {FAILS}")
