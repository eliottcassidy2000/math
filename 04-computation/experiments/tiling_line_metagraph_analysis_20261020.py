#!/usr/bin/env python3
"""tiling_line_metagraph_analysis_20261020.py -- opus-2026-10-20-S46
Reads the per-class and line TSV files written by tiling_line_metagraph_census_20261020.cpp and
reports node-level laws: degrees by colour, the fundamental black graph B (sigma-twins identified)
and its degree parities, self-line classes, the transitive class's partner, depth profiles by node
type, multiplicity structure, and the handshake constraint.  Usage: python3 ... n [n ...]
"""
import sys, os
from collections import Counter, defaultdict

def load(n):
    base = os.path.dirname(os.path.abspath(__file__))
    cls = {}
    with open(os.path.join(base, f"tiling_line_metagraph_census_20261020_classes_n{n}.tsv"), encoding="utf-8") as f:
        hdr = f.readline().strip().split("\t")
        for line in f:
            x = line.rstrip("\n").split("\t"); d = dict(zip(hdr, x))
            d["class"] = int(d["class"]); d["H"] = int(d["H"]); d["Aut"] = int(d["Aut"]); d["tc"] = int(d["tc"]); d["bm"] = int(d["bm"])
            d["SC"] = int(d["SC"]); d["sigmaClass"] = int(d["sigmaClass"]); d["scores"] = tuple(int(v) for v in d["scores"].split(","))
            d["depthProfile"] = tuple(int(v) for v in d["depthProfile"].split(","))
            cls[d["class"]] = d
    lines = []
    with open(os.path.join(base, f"tiling_line_metagraph_census_20261020_lines_n{n}.tsv"), encoding="utf-8") as f:
        f.readline()
        for line in f:
            a, b, bl, bk = (int(v) for v in line.split("\t")); lines.append((a, b, bl, bk))
    return cls, lines

def analyse(n):
    cls, lines = load(n)
    C = len(cls); print(f"=== n = {n}: {C} classes, {len(lines)} joined class pairs")
    blueDeg = Counter(); blackDeg = Counter(); selfBlue = Counter(); selfBlack = Counter()
    for a, b, bl, bk in lines:
        if a == b: selfBlue[a] += bl; selfBlack[a] += bk; blueDeg[a] += 2 * bl; blackDeg[a] += 2 * bk
        else: blueDeg[a] += bl; blueDeg[b] += bl; blackDeg[a] += bk; blackDeg[b] += bk
    # degree = tiling counts
    bad = [c for c in cls if blueDeg[c] != cls[c]["bm"] or blackDeg[c] != cls[c]["tc"] - cls[c]["bm"]]
    print(f"degree bookkeeping (blue degree = bm, black degree = tc - bm): violations {len(bad)}")
    # node types
    types = {c: ("PK" if not cls[c]["SC"] else ("PB" if cls[c]["bm"] == cls[c]["tc"] else "MX")) for c in cls}
    # transitive class = the one with H = 1
    trans = [c for c in cls if cls[c]["H"] == 1]
    for c in trans:
        partners = [(b if a == c else a, bl, bk) for a, b, bl, bk in lines if a == c or b == c]
        for p, bl, bk in partners:
            print(f"transitive class {c} (scores {cls[c]['scores']}) -> partner class {p}: H = {cls[p]['H']}, |Aut| = {cls[p]['Aut']}, tc = {cls[p]['tc']}, bm = {cls[p]['bm']}, type {types[p]}, scores {cls[p]['scores']}, line colour {'blue' if bl else 'black'}")
    # pure-blue classes and their partners
    for c in cls:
        if types[c] == "PB":
            partners = sorted(set((b if a == c else a) for a, b, bl, bk in lines if a == c or b == c))
            print(f"pure-blue class {c}: H = {cls[c]['H']}, |Aut| = {cls[c]['Aut']}, tc = {cls[c]['tc']}, scores {cls[c]['scores']}, partners {[(p, types[p], cls[p]['H']) for p in partners]}, self-lines {selfBlue[c]}")
    # fundamental black graph B: identify sigma-twins at the converse-merged level
    mid = lambda c: min(c, cls[c]["sigmaClass"])
    mergedBlack = Counter(); mergedBlue = Counter()
    for a, b, bl, bk in lines:
        u, v = sorted((mid(a), mid(b))); mergedBlack[(u, v)] += bk; mergedBlue[(u, v)] += bl
    oddm = sum(1 for v in mergedBlack.values() if v % 2); print(f"merged black multiplicities odd: {oddm} of {len(mergedBlack)}")
    Bdeg = Counter()
    for (u, v), k in mergedBlack.items():
        h = k // 2
        if u == v: Bdeg[u] += 2 * h
        else: Bdeg[u] += h; Bdeg[v] += h
    # B-degree law: NS merged node: tc (odd); SC node: (tc - bm)/2
    badB = 0
    for c in cls:
        if mid(c) != c: continue
        if cls[c]["SC"]: expect = (cls[c]["tc"] - cls[c]["bm"]) // 2
        else: expect = cls[c]["tc"]
        if Bdeg[c] != expect: badB += 1
    print(f"fundamental black graph B: {sum(mergedBlack.values())//2} edges; B-degree law violations {badB}; odd-degree B nodes {sum(1 for c in cls if mid(c)==c and Bdeg[c]%2)} (NS merged nodes {sum(1 for c in cls if mid(c)==c and not cls[c]['SC'])})")
    # self-line classes
    sb = {c: v for c, v in selfBlue.items() if v}; sk = {c: v for c, v in selfBlack.items() if v}
    print(f"classes with blue self-lines: {len(sb)} (lines {sum(sb.values())}); with black self-lines: {len(sk)} (lines {sum(sk.values())}); black self-line class types {Counter(types[c] for c in sk)}; max black self-lines at a class {max(sk.values()) if sk else 0}")
    hist = Counter();
    for c, v in sk.items(): hist[(cls[c]["H"], cls[c]["Aut"], v)] += 1
    print(f"black self-line classes by (H, |Aut|, self-lines): {dict(sorted(hist.items()))}" if n <= 7 else f"black self-line classes: |Aut| histogram {Counter(cls[c]['Aut'] for c in sk)}, SC among them {sum(cls[c]['SC'] for c in sk)}, s_0-vs-s_(n-1) not applicable at class level")
    # depth profiles by node type (sum over classes)
    agg = defaultdict(lambda: [0] * 8)
    for c in cls:
        for d, v in enumerate(cls[c]["depthProfile"]): agg[types[c]][d] += v
    print("depth profile by node type:", {k: tuple(v[:len(cls[next(iter(cls))]['depthProfile'])]) for k, v in agg.items()})
    # multiplicity vs fibre product: top deviations
    tot = 2 ** ((n - 1) * (n - 2) // 2)
    ratios = []
    for a, b, bl, bk in lines:
        if a == b: continue
        expect = cls[a]["tc"] * cls[b]["tc"] / tot
        ratios.append(((bl + bk) / expect, a, b, bl + bk))
    ratios.sort()
    print(f"actual/fibre-product ratio over joined pairs: min {ratios[0][0]:.3f}, median {ratios[len(ratios)//2][0]:.3f}, max {ratios[-1][0]:.1f} (classes {ratios[-1][1]},{ratios[-1][2]} with {ratios[-1][3]} lines)")
    # Klein / handshake
    oddSC = sum(1 for c in cls if cls[c]["SC"] and ((cls[c]["tc"] - cls[c]["bm"]) // 2) % 2 == 1)
    ns = sum(1 for c in cls if not cls[c]["SC"]) // 2
    print(f"handshake: NS merged {ns} + SC-with-odd-half {oddSC} = {ns + oddSC} (even: {(ns + oddSC) % 2 == 0})")
    sys.stdout.flush()

if __name__ == "__main__":
    for n in [int(a) for a in sys.argv[1:]] or [5, 6, 7, 8]: analyse(n)
