"""pycheck.py -- independent pure-Python edge-multiset resolvability checker (no numpy, no shared code).

Straight from the definition: Q_d has vertices 0..2^d-1, edges {u, u + 2^i} (bit i of u = 0);
d(e, s) = min(d_H(u, s), d_H(v, s)) for e = {u, v}; the histogram of e w.r.t. S is
H_e(r) = #{s in S : d(e, s) = r}, r = 0..d-1;  S is edge-multiset resolving iff all H_e are distinct.
Uses NO projection lemma, NO keys/hashing: histograms are compared as tuples.

usage: python3 pycheck.py d < file_with_sets      (sets: one per line, vertices as integers; any separators)
       python3 pycheck.py selftest
"""
import sys, re, random
from collections import Counter

def edges(d):
    return [(u, u | (1 << i)) for u in range(1 << d) for i in range(d) if not (u >> i) & 1]

def hamming(x, y):
    return bin(x ^ y).count('1')

def histograms(d, S):
    out = []
    for (u, v) in edges(d):
        H = [0] * d
        for s in S:
            H[min(hamming(u, s), hamming(v, s))] += 1
        out.append(tuple(H))
    return out

def defect(d, S):
    H = histograms(d, S)
    return len(H) - len(set(H))

def is_resolving(d, S):
    S = list(S)
    assert len(set(S)) == len(S) and all(0 <= s < (1 << d) for s in S)
    return len(S) > 0 and defect(d, S) == 0

def parse_sets(text):
    for line in text.splitlines():
        if 'set=' in line:
            line = line.split('set=')[1]
        nums = [int(x) for x in re.findall(r'\d+', line)]
        if nums:
            yield nums

if __name__ == '__main__':
    if sys.argv[1] == 'selftest':
        # the explicit resolving 19-set of Q7 from procgen_edim_20261001_run.py (S10) and the paper's Q6 15-set
        S19 = [5, 57, 54, 104, 35, 109, 115, 49, 6, 39, 55, 102, 15, 32, 85, 97, 21, 113, 41]
        print('Q7 19-set resolving:', is_resolving(7, S19), ' defect', defect(7, S19))
        P6 = [v for v in range(64) if (0x02283022a042a00a >> v) & 1]
        print('Q6 paper 15-set resolving:', is_resolving(6, P6), ' |S| =', len(P6))
        # each 18-subset of the 19-set: defects
        print('Q7 19-set minus one landmark, defects:', [defect(7, S19[:j] + S19[j + 1:]) for j in range(19)])
        sys.exit(0)
    d = int(sys.argv[1])
    for S in parse_sets(sys.stdin.read()):
        dd = defect(d, S)
        print('%s k=%d defect=%d' % ('RESOLVING' if dd == 0 and len(set(S)) == len(S) else 'NOT', len(S), dd))
