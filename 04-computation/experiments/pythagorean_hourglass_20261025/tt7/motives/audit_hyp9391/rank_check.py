# Audit of HYP-9391: exact ranks.  (1) rank of [base | Y02 Y05 Y09 Y18] over Q at the primes q = 3 mod 4 in [7, 4999]
# (auditor's Y sums and PARI coefficients); (2) rank of the 16 'new' genuine sums modulo the base with the auditor's
# brute-force S_H (all primes available); (3) coincidences among the genuine sums (S_H36 = S_H05 = S_H06 = -S_H21 ...).
# Usage: rank_check.py ysums forms brute_files... ; writes rank_check.out
import sys, os, json
from sympy import Matrix
HERE = os.path.dirname(os.path.abspath(__file__))
Yf, Ff, bfiles = sys.argv[1], sys.argv[2], sys.argv[3:]
forms = {}
for l in open(Ff):
    p, a, c, d, b = map(int, l.split()); forms[p] = (a, c, d, b)
Y = {}
for l in open(Yf):
    r = list(map(int, l.split())); Y[r[0]] = (r[1], r[3], r[5], r[6])
def base(q):
    a, c, d, b = forms[q]
    return [1, q, q * q, q ** 3, a, q * a, c, d, b]
out = []
qs = sorted(Y)[:60]
Mb = Matrix([base(q) for q in qs]); MbY = Matrix([base(q) + list(Y[q]) for q in qs])
out.append(f"(1) {len(qs)} primes {qs[0]}..{qs[-1]}: rank(base) = {Mb.rank()}, rank(base + Y02,Y05,Y09,Y18) = {MbY.rank()}")
S = {}
for fn in bfiles:
    for l in open(fn):
        r = list(map(int, l.split()))
        if r: S[r[0]] = r[1:]
G = json.load(open(os.path.join(HERE, '..', 'genuine37.json')))
NEW = [2, 4, 5, 6, 9, 10, 16, 17, 18, 20, 21, 24, 28, 30, 31, 36]
qb = sorted(S)
Mb2 = Matrix([base(q) for q in qb]); MbN = Matrix([base(q) + [S[q][i] for i in NEW] for q in qb])
MbNY = Matrix([base(q) + [S[q][i] for i in NEW] + list(Y[q]) for q in qb])
MbYonly = Matrix([base(q) + list(Y[q]) for q in qb])
out.append(f"(2) {len(qb)} primes {qb[0]}..{qb[-1]}: rank(base) = {Mb2.rank()}, rank(base + 16 new S_H) = {MbN.rank()}, "
           f"rank(base + Y's) = {MbYonly.rank()}, rank(base + new S_H + Y's) = {MbNY.rank()}")
gen = [i for i, H in enumerate(G) if not H['all_even']]
groups = {}
for i in gen:
    key = tuple(S[q][i] for q in qb)
    neg = tuple(-x for x in key)
    if neg in groups: groups[neg].append(-i if i else 'neg0')
    else: groups.setdefault(key, []).append(i)
coinc = [v for v in groups.values() if len(v) > 1]
out.append(f"(3) genuine sums equal up to sign at all {len(qb)} primes (negative index = opposite sign): {coinc}")
print('\n'.join(out))
open(os.path.join(HERE, 'rank_check.out'), 'w').write('\n'.join(out) + '\n')
