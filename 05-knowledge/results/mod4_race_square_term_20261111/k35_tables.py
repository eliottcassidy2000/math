#!/usr/bin/env python3
"""Summary of the outputs of k35cube.c: (1) the censuses of all classes of orders m = 1..10 (exact vertex-cube sums
for K = {3}, {5}, {3,5} with labelled weights, hence the exact E_{m+1}[K]); (2) the Monte Carlo runs at orders
m = 11..40 (4,000,000 uniformly random labelled bases each).  Compares every number with the closed forms of
Theorem P (pentagon_closed_forms.py) and reports the sign statistics of the cube sum G_{3,5}(R).
Usage: python3 k35_tables.py   (reads out/k35_census_*.txt, out/k35mc_*.out, out/*.err)"""
import re, glob, math, os
from fractions import Fraction as Fr
from math import comb
from pentagon_closed_forms import E3, E5, E35
here = os.path.dirname(os.path.abspath(__file__)); out = os.path.join(here, 'out')
exact = {'{3}': E3, '{5}': E5, '{3,5}': E35}

def parse(line):
    g = re.match(r'(census|mc) m (\d+) K (\S+) (.*)', line)
    if not g: return None
    d = dict(mode=g.group(1), m=int(g.group(2)), K=g.group(3)); rest = g.group(4)
    for key in ('classes', 'weight', 'signedsum', 'N', 'neg', 'zero', 'pos'):
        h = re.search(r'\b%s (-?\d+)' % key, rest)
        if h: d[key] = int(h.group(1))
    for key in ('S1', 'S2'):
        h = re.search(r'\b%s (\S+)' % key, rest)
        if h: d[key] = float(h.group(1))
    d['rank'] = {int(a): int(b) for a, b in re.findall(r' r(\d+) (\d+)', rest)}
    return d

print("(1) Censuses: every class of order m, labelled weights m!/|Aut|.  E_{m+1}[K] = (signed sum)/(2^C(m,2) 2^m).")
print(" m   classes     K       E_(m+1)[K] exact from the census      = closed form?   G<0  G=0      G>0")
tot_classes = 0; mism = 0; negs = 0
for m in range(1, 11):
    files = sorted(glob.glob(os.path.join(out, 'k35_census_%d.txt' % m)) + glob.glob(os.path.join(out, 'k35_census_%d_part*.txt' % m)))
    acc = {}
    for fn in files:
        for line in open(fn):
            d = parse(line)
            if not d: continue
            a = acc.setdefault(d['K'], dict(classes=0, weight=0, signedsum=0, neg=0, zero=0, pos=0, rank={}))
            for key in ('classes', 'weight', 'signedsum', 'neg', 'zero', 'pos'): a[key] += d[key]
            for r, c in d['rank'].items(): a['rank'][r] = a['rank'].get(r, 0) + c
    for K in ('{3}', '{5}', '{3,5}'):
        if K not in acc: continue
        a = acc[K]
        if m >= 2: assert a['weight'] == 2 ** comb(m, 2), (m, a['weight'])
        if a['classes'] == 0: continue
        E = Fr(a['signedsum'], a['weight'] * 2 ** m); ok = (E == exact[K](m + 1)); mism += (not ok)
        print("%2d  %9d   %-6s  %-36s  %-14s  %4d %6d %9d" % (m, a['classes'], K, E, "yes" if ok else "NO", a['neg'], a['zero'], a['pos']))
        if K == '{3,5}': tot_classes += a['classes']; negs += a['neg']
    if '{3,5}' in acc and acc['{3,5}']['classes']:
        print("      rank of the polar form of {3,5} (classes): " + ", ".join("%d: %d" % kv for kv in sorted(acc['{3,5}']['rank'].items())))
print("classes of orders 1-10 examined: %d; with G_{3,5} < 0: %d; mismatches with the closed forms: %d" % (tot_classes, negs, mism))
print()
print("(2) Monte Carlo, N = 4,000,000 random labelled bases per order m; estimator = exact cube sum G_K(R)/2^m.")
print(" m   n    E_n[{3,5}] closed form   estimate        s.e.       z      E_n[{5}] closed form   estimate        z      G_{3,5}<0   G_{3,5}=0")
totN = 0; negs = 0; zmax = 0.0; nz = 0
for m in range(11, 41):
    fn = os.path.join(out, 'k35mc_%d.out' % m)
    if not os.path.exists(fn): continue
    rows = {}
    for line in open(fn):
        d = parse(line)
        if d: rows[d['K']] = d
    n = m + 1; a = rows['{3,5}']; b = rows['{5}']; N = a['N']; totN += N; negs += a['neg']
    e = a['S1'] / N; se = math.sqrt(max(a['S2'] / N - e * e, 0) / N); x = float(E35(n)); z = (e - x) / se
    e5 = b['S1'] / N; se5 = math.sqrt(max(b['S2'] / N - e5 * e5, 0) / N); x5 = float(E5(n)); z5 = (e5 - x5) / se5
    zmax = max(zmax, abs(z), abs(z5)); nz += 2
    assert rows['{3}']['S1'] == N * 2.0 ** ((m + 1) // 2 - m)   # Theorem T holds sample by sample
    print("%2d  %2d   %.10e         %.10e  %.2e  %+5.2f    %.6e           %+.6e  %+5.2f   %9d   %9d" % (m, n, x, e, se, z, x5, e5, z5, a['neg'], a['zero']))
print("random bases examined: %d; with G_{3,5} < 0: %d; largest |z| among %d comparisons: %.2f" % (totN, negs, nz, zmax))
print()
print("(3) Check lines of the program (must all be zero):")
bad = 0; lines = 0
for fn in sorted(glob.glob(os.path.join(out, 'k35_census_*.err')) + glob.glob(os.path.join(out, 'k35mc_checks.err'))):
    for line in open(fn):
        lines += 1
        nums = re.findall(r'(?:trace formulas|direct summation|Theorem T|score formulas|closed formula|kernel weights) (\d+)', line)
        if len(nums) != 6 or any(int(v) for v in nums): bad += 1; print("  FAILURE: " + fn + ": " + line.strip())
print("  %d check lines read, %d with a failure (checks: trace formulas, direct summation for m <= 9, Theorem T, (P1) score formulas, (P2) closed formula for the cube sum, kernel weights mod 4)" % (lines, bad))
