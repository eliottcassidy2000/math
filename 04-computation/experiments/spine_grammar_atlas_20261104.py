#!/usr/bin/env python3
"""spine_grammar_atlas_20261104.py -- opus-2026-11-04-S79

Minimal-grade atlas from census outputs of spine_grammar_census_20261104.cpp: for each odd N <= X, which grades
G(k, w) attain N (directly, or as an ordered-join product of directly attained values), the miss counts per grammar in
[1, 300], [1, 511], [1, X] and per dyadic scale, and the pre-registered prediction checks.
Usage: python3 spine_grammar_atlas_20261104.py X file1 file2 ...   (file names contain k<k>_w<w>)
"""
import sys, re

def load(path):
    att = set()
    for line in open(path):
        if line.startswith('ATTAINED_DIRECT:'):
            att = set(int(x) for x in line.split()[1:])
    return att

def closure(att, X):
    cl = set(att)
    vals = sorted(a for a in att if a >= 3)
    changed = True
    while changed:
        changed = False
        for a in sorted(cl):
            if a < 3: continue
            for b in vals:
                if a * b > X: break
                if a * b not in cl: cl.add(a * b); changed = True
    return cl

def main():
    X = int(sys.argv[1]); files = sys.argv[2:]
    grades = {}
    for f in files:
        m = re.search(r'k(\d)_w(\d)', f); k, w = int(m.group(1)), int(m.group(2))
        att = load(f); grades[(k, w)] = (att, closure(att, X))
    odds = list(range(1, X + 1, 2))
    print("grammar | misses<=300 (besides 7,21) | misses<=511 | misses<=X | per dyadic scale [2^j,2^(j+1)) misses/total")
    for (k, w) in sorted(grades):
        att, cl = grades[(k, w)]
        m300 = [n for n in odds if n <= 300 and n not in cl and n not in (7, 21)]
        m511 = [n for n in odds if n <= 511 and n not in cl and n not in (7, 21)]
        mX = [n for n in odds if n not in cl and n not in (7, 21)]
        scales = []
        j = 2
        while (1 << j) <= X:
            lo, hi = 1 << j, min(X, (1 << (j + 1)) - 1)
            tot = [n for n in odds if lo <= n <= hi]; mis = [n for n in tot if n not in cl]
            scales.append(f"2^{j}:{len(mis)}/{len(tot)}"); j += 1
        print(f"G({k},{w}) | {len(m300)}: {m300[:30]}{'...' if len(m300) > 30 else ''} | {len(m511)} | {len(mX)} | {' '.join(scales)}")
    # minimal grades per odd number
    print("\nodd N <= 300 with the set of grades attaining it (closure under products):")
    hard = []
    for n in odds:
        if n > 300: break
        gs = [f"({k},{w})" for (k, w) in sorted(grades) if n in grades[(k, w)][1]]
        if len(gs) <= 2 or n in (13, 67, 173): print(f"  {n}: {' '.join(gs) if gs else 'NONE'}")
    print("\nodd N <= X attained by no grade in the census (besides 7, 21):")
    none = [n for n in odds if n not in (7, 21) and not any(n in grades[g][1] for g in grades)]
    print(f"  {len(none)}: {none[:60]}")

main()
