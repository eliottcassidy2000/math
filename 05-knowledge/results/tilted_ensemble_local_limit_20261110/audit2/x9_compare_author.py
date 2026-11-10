#!/usr/bin/env python3
"""AUDIT 2.  Field-by-field comparison of the packet's raw totals (out/biasweights_2_9.out, out/biasweights_10_parts.out,
out/oddedge_2_9.out, out/oddedge_10_parts.out) with the auditor's raw totals (x1_raw_*.out).  Nothing is imported from
the packet; its output files are only read.
Usage: python3 x9_compare_author.py <packet out dir> x1_raw_direct_2_8.out x1_raw_direct_9.out x1_raw_cut_10.out"""
import sys, re
from collections import defaultdict
out = sys.argv[1]
def add(d, k, v): d[k] = d.get(k, 0) + v
auth = defaultdict(dict); authE = defaultdict(dict)
for fn in ('biasweights_2_9.out', 'biasweights_10_parts.out'):
    for line in open(out + '/' + fn):
        g = re.match(r'order (\d+) classes (\d+) weight (\d+) failures\(mu odd, a odd, C3, aut\) (\d+) (\d+) (\d+) (\d+) S (-?\d+) U (-?\d+) A (\d+)(.*)', line)
        if not g: continue
        m = int(g.group(1)); d = auth[m]
        for k, i in (('classes', 2), ('W', 3), ('S', 8), ('U', 9), ('A', 10)): add(d, k, int(g.group(i)))
        add(d, 'fails', sum(int(g.group(i)) for i in (4, 5, 6, 7)))
        for h in re.finditer(r' r(\d+) (\d+) (\d+) (\d+) (\d+)', g.group(11)):
            for j, nm in enumerate(('c0', 'c1', 'l0', 'l1')): add(d, (int(h.group(1)), nm), int(h.group(2 + j)))
for fn in ('oddedge_2_9.out', 'oddedge_10_parts.out'):
    for line in open(out + '/' + fn):
        g = re.match(r'order (\d+) classes (\d+) failures: sum mu != sum a (\d+) ; sum a != sum_S H\(S\)H\(S\^c\) (\d+) ; sum mu != 2 mod 4 (\d+) ; even number of beta-edges (\d+) ; even weight of ell (\d+) ; edges:(.*)', line)
        if not g: continue
        m = int(g.group(1)); d = authE[m]; add(d, 'classes', int(g.group(2))); add(d, 'fails', sum(int(g.group(i)) for i in (3, 4, 5, 6, 7)))
        for x in g.group(8).split(): e, c = x.split(':'); add(d, ('e', int(e)), int(c))
mine = defaultdict(dict); mineE = defaultdict(dict)
for fn in sys.argv[2:]:
    for line in open(fn):
        g = re.match(r'order (\d+) mode (\S+) classes (\d+) W (-?\d+) S1 (-?\d+) U (-?\d+) A (-?\d+) ', line)
        if not g or g.group(2) == 'labelled': continue
        m = int(g.group(1)); d = mine[m]
        for k, i in (('classes', 3), ('W', 4), ('S', 5), ('U', 6), ('A', 7)): add(d, k, int(g.group(i)))
        add(d, 'fails', 0)
        for h in re.finditer(r' rank (\d+) (\d+) (\d+) (\d+) (\d+)', line):
            for j, nm in enumerate(('c0', 'c1', 'l0', 'l1')): add(d, (int(h.group(1)), nm), int(h.group(2 + j)))
        e = mineE[m]; add(e, 'classes', int(g.group(3))); add(e, 'fails', 0)
        h = re.search(r' edges((?: \d+:\d+)*)', line)
        for x in h.group(1).split(): ee, c = x.split(':'); add(e, ('e', int(ee)), int(c))
def nz(d): return {k: v for k, v in d.items() if v != 0 or k in ('classes', 'W', 'S', 'U', 'A', 'fails')}
for m in sorted(set(auth) | set(mine)):
    a, b = nz(auth.get(m, {})), nz(mine.get(m, {})); ea, eb = nz(authE.get(m, {})), nz(mineE.get(m, {}))
    print("order %2d: bias totals (classes, weight, S, U, A, classes and labelled weight per rank and zero/non-zero bias): %s (%d fields) ; edge histogram of beta_R: %s (%d fields) ; packet failure counters: %d + %d"
          % (m, "IDENTICAL" if a == b else "DIFFERENT", len(a), "IDENTICAL" if ea == eb else "DIFFERENT", len(ea), a.get('fails', -1), ea.get('fails', -1)))
    if a != b: print("   ", {k: (a.get(k), b.get(k)) for k in set(a) | set(b) if a.get(k) != b.get(k)})
    if ea != eb: print("   ", {k: (ea.get(k), eb.get(k)) for k in set(ea) | set(eb) if ea.get(k) != eb.get(k)})
