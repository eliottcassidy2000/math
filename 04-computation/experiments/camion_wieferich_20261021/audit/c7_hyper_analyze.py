# c7_hyper_analyze.py -- second-pass audit of THM-4750 section 6b: re-derive the r = 3 and r = 4 root tables.
#  (1) my own parser of the session outputs (hyper4_second_order_5_2200.out; hyper3_second_order_*.out);
#  (2) independent recomputation from c7_hyper_r.c (different exponent method and different recurrence), with the root count
#      of g(l) = sum_i c_i C(l,i) over F_p obtained as deg gcd(g, x^p - x) (FLINT), not by evaluation at every l;
#  (3) per-prime, per-class comparison with the session's root counts;
#  (4) the forced roots: classes 2 (r = 3, 4) and 3 (r = 4) vanish at l = -1 (n = 2, 3), and the low bases m = 2, 3 obey
#      the degree-r law below the proved threshold r + 2;
#  (5) exact S_r fixed-point laws and "forced root + S_(r-1)";
#  (6) the anomalous class-2 count 3 at r = 4 (p = 71) and the non-genuine classes among zero classes.
# Usage: python3 c7_hyper_analyze.py RAW_R4 RAW_R3      (RAW = stdout of c7_hyper_r.c for r = 4 / r = 3 on [7, 2200])
# Second-pass auditor, 2026-10-21.
import sys, os, re, flint
from fractions import Fraction
from itertools import permutations
from collections import Counter, defaultdict
from math import factorial
here = os.path.dirname(os.path.abspath(__file__)); exp = os.path.join(here, '..')
def parse_session(files):
    d = {}
    for f in files:
        for line in open(os.path.join(exp, f)):
            m = re.match(r'p=(\d+) ord2=(\d+) zero-classes:(.*)', line)
            if not m: continue
            d[int(m.group(1))] = [(int(a), int(b)) for a, b in re.findall(r'(\d+):(\d+)', m.group(3))]
    return d
def parse_raw(path):
    P = {}; cur = None
    for line in open(path):
        t = line.split()
        if t[0] == 'P': cur = int(t[1]); P[cur] = {'ord': int(t[2]), 'lawfail': int(t[3]), 'notfull': int(t[4]), 'Z': [], 'B': {}}
        elif t[0] == 'Z': P[cur]['Z'].append((int(t[1]), int(t[2]), [int(x) for x in t[3:]]))
        elif t[0] == 'B':
            kv = dict(x.split('=') for x in t[2:]); P[cur]['B'][int(t[1])] = {k: int(v) for k, v in kv.items()}
    return P
def gpoly(c, p):
    x = flint.nmod_poly([0, 1], p); g = flint.nmod_poly([0], p); basis = flint.nmod_poly([1], p)
    for i, ci in enumerate(c):
        g += basis * (ci * pow(factorial(i), -1, p) % p)
        basis *= (x - i)
    return g
def nroots(g, p):
    if g.is_zero(): return p
    if g.degree() == 0: return 0
    x = flint.nmod_poly([0, 1], p)
    h = g.gcd(x.pow_mod(p, g) - x)
    return h.degree()
def fixdist(k):
    c = Counter(sum(1 for i, x in enumerate(s) if i == x) for s in permutations(range(k)))
    return {j: Fraction(c[j], factorial(k)) for j in range(k + 1)}
def typ(label): return label if label in (1, 2, 3) else 'sporadic'
for r, raw, sfiles in ((4, sys.argv[1], ['hyper4_second_order_5_2200.out']),
                       (3, sys.argv[2], ['hyper3_second_order_5_1400.out', 'hyper3_second_order_1401_2200.out'])):
    S = parse_session(sfiles); A = parse_raw(raw)
    print(f"===== r = {r} =====")
    print(f"primes: session {len(S)}, audit {len(A)}, same set: {sorted(S) == sorted(A)}")
    tot_law = sum(v['lawfail'] for v in A.values()); tot_nf = sum(v['notfull'] for v in A.values())
    print(f"audit: degree-{r} law failures over all classes m in [r+3, r+p+1]: {tot_law}; classes with p^2 | a_r (not genuine degree {r}): {tot_nf}")
    sess_tab = Counter((typ(l), k) for p in S for (l, k) in S[p])
    aud_tab = Counter(); mism = []; nongen = []; forced = Counter(); details = {}
    for p, v in sorted(A.items()):
        mine = {}
        for (m, label, c) in v['Z']:
            g = gpoly(c, p); k = nroots(g, p); mine[label] = k; aud_tab[(typ(label), k)] += 1
            if c[r] % p == 0: nongen.append((p, label, k, g.degree() if not g.is_zero() else -1))
            if label in ((2,) if r == 3 else (2, 3)):
                forced[(label, int(g(p - 1)) == 0)] += 1                    # root at l = -1 ?
            details[(p, label)] = (c, g)
        sm = dict(S.get(p, []))
        if sm != mine: mism.append((p, sm, mine))
    print(f"per-prime per-class root counts, session vs audit: {len(mism)} primes differ" + (f": {mism[:3]}" if mism else ""))
    print("table (class: counts of 0..r roots)   session  |  audit")
    for t in (1, 2, 3, 'sporadic'):
        sv = [sess_tab[(t, k)] for k in range(r + 1)]; av = [aud_tab[(t, k)] for k in range(r + 1)]
        if sum(sv) or sum(av): print(f"   class {str(t):9s}: {' / '.join(map(str, sv))}   |  {' / '.join(map(str, av))}   (n = {sum(av)})")
    print(f"forced root at l = -1 (n = 2{', 3' if r == 4 else ''}): {dict(forced)}")
    lowlaw = Counter((m, v['B'][m]['law'], v['B'][m]['zeroclass'], v['B'][m]['value_mod_p2'] == 0) for v in A.values() for m in v['B'])
    print(f"low bases (m, law holds at base m, class is a zero class, class polynomial of base m+(p-1) vanishes at l=-1 mod p^2): {dict(lowlaw)}")
    print(f"zero classes not genuinely of degree {r} (p, class, roots, degree of g): {nongen}")
    fd = fixdist(r); ff = {j + 1: v for j, v in fixdist(r - 1).items()}
    print(f"S_{r} fixed points: {[str(fd[j]) for j in range(r + 1)]};  forced root + S_{r-1}: {[str(ff.get(j, 0)) for j in range(r + 1)]}")
    if r == 4:
        for (p, label), (c, g) in details.items():
            if label == 2 and nroots(g, p) == 3:
                print(f"   anomalous class 2 at p={p}: c = {c}, g = {g}, factorisation {g.factor()}")
