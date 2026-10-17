#!/usr/bin/env python3
"""Audit of THM-4665: compare the S39 outputs with the independent recounts (a2_*.out), and the base rates.

Reads 05-knowledge/results/ei17_finite_field_realizations_20261017.out and ei17_twisted_unit_graph_embeddings_20261017.out
(faithful column) and this folder's a2_*.out files; prints every prime where both exist and any disagreement; checks
the lists of statements 1, 2, 4; computes the control rates and the binomial probability of the zero counts.
"""
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RES = os.path.join(ROOT, '05-knowledge', 'results')


def parse_s39(name):
    d = {}
    for line in open(os.path.join(RES, name), encoding='utf-8'):
        m = re.match(r'^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$', line)
        if m:
            p, f = int(m.group(1)), int(m.group(3))
            if p in d:
                assert d[p] == f, ('inconsistent duplicate in S39 output', name, p)
            d[p] = f
    return d


def parse_a2(pattern):
    d = {}
    for fn in sorted(glob.glob(os.path.join(HERE, pattern))):
        for line in open(fn, encoding='utf-8'):
            if line.startswith('#'):
                continue
            m = re.match(r'^\s*(\d+)\s+(\d+)\s+(\d+)\s', line)
            if m:
                p, f = int(m.group(1)), int(m.group(3))
                if p in d:
                    assert d[p] == f, ('methods disagree', fn, p, d[p], f)
                d[p] = f
    return d


s39s = parse_s39('ei17_finite_field_realizations_20261017.out')
s39t = parse_s39('ei17_twisted_unit_graph_embeddings_20261017.out')
mys = parse_a2('a2_split_*.out')
myt = parse_a2('a2_twisted_*.out')
for lab, a, b in (('split', s39s, mys), ('twisted', s39t, myt)):
    common = sorted(set(a) & set(b))
    bad = [(p, a[p], b[p]) for p in common if a[p] != b[p]]
    only_mine = sorted(set(b) - set(a))
    print('%s: %d primes in both, disagreements %s; primes only in the audit runs: %s'
          % (lab, len(common), bad, [(p, b[p]) for p in only_mine]))

S1_ZERO = [331, 661, 881, 991, 1321, 1871, 2311, 2531, 2971, 3301, 3631, 3851, 5171, 5281, 5501, 6271, 6491, 7151, 7591, 8581,
           9241, 9461, 9791, 9901, 10891, 11551]
S1_POS = [2861, 3191, 4621, 4951, 7481, 10781]
T54 = [109, 439, 659, 769, 1319, 1429, 1759, 1979, 2089, 2309, 2749, 2969, 3079, 3299, 3739, 4289, 4729, 5059, 5279, 5939]
CS = [103, 139, 281, 313, 401, 499, 521, 541, 631, 641, 677]
CT = [127, 149, 227, 277, 283, 311, 317, 349, 383, 397, 443, 463, 547, 587]
print('statement 1: zero at the 26 listed primes: %s; 8 or 16 at the six: %s'
      % (all(mys.get(p) == 0 for p in S1_ZERO), [(p, mys.get(p)) for p in S1_POS]))
print('statement 2: zero at the 20 conjugation-fixed primes: %s' % all(myt.get(p) == 0 for p in T54))
ctrl_s = sorted(p for p in mys if p < 700)
ctrl_t = sorted(p for p in myt if p < 600)
print('statement 4: split primes < 700 tested %d (from %d), positive %s = list: %s'
      % (len(ctrl_s), ctrl_s[0], [p for p in ctrl_s if mys[p] > 0], [p for p in ctrl_s if mys[p] > 0] == CS))
print('             twisted primes < 600 tested %d (from %d), positive %s = list: %s; 2441 -> %s, 2551 -> %s'
      % (len(ctrl_t), ctrl_t[0], [p for p in ctrl_t if myt[p] > 0], [p for p in ctrl_t if myt[p] > 0] == CT,
         myt.get(2441), myt.get(2551)))
allc = [v for v in list(mys.values()) + list(myt.values()) if v > 0]
print('             all positive counts: %s (multiples of 4: %s; of 8: %d of %d)'
      % (sorted(set(allc)), all(v % 4 == 0 for v in allc), sum(v % 8 == 0 for v in allc), len(allc)))
# base rates (controls exclude the special classes)
cs_gen = [p for p in ctrl_s if p % 55 not in (1, 54)]
ct_gen = [p for p in ctrl_t if p % 55 not in (1, 54)]
rs = sum(mys[p] > 0 for p in cs_gen) / len(cs_gen)
rt = sum(myt[p] > 0 for p in ct_gen) / len(ct_gen)
s1 = [p for p in mys if p % 55 == 1]
t54 = [p for p in myt if p % 55 == 54]
ks = sum(mys[p] > 0 for p in s1)
kt = sum(myt[p] > 0 for p in t54)
from math import comb
binom_le = lambda k, n, r: sum(comb(n, j) * r ** j * (1 - r) ** (n - j) for j in range(k + 1))
binom_ge = lambda k, n, r: 1 - binom_le(k - 1, n, r) if k > 0 else 1.0
print('rates: split controls (classes != +-1 mod 55) %d/%d = %.3f; class 1 mod 55: %d/%d = %.3f; P(X >= %d) = %.3f'
      % (sum(mys[p] > 0 for p in cs_gen), len(cs_gen), rs, ks, len(s1), ks / len(s1), ks, binom_ge(ks, len(s1), rs)))
print('       twisted controls %d/%d = %.3f; class 54 mod 55: %d/%d; P(X = 0) = %.3f'
      % (sum(myt[p] > 0 for p in ct_gen), len(ct_gen), rt, kt, len(t54), binom_le(0, len(t54), rt)))
m10 = [p for p in mys if p % 11 == 10 and p < 2600]
print('note observation: split, primes = 10 mod 11 below 2600: %d tested, positive %s; P(0) at the split rate = %.3f'
      % (len(m10), [p for p in m10 if mys[p] > 0], binom_le(0, len(m10), rs)))
m10t = [p for p in myt if p % 11 == 10 and p < 2600]
print('                  twisted, primes = 10 mod 11 below 2600: %d tested, positive %s'
      % (len(m10t), [(p, myt[p]) for p in m10t if myt[p] > 0]))
