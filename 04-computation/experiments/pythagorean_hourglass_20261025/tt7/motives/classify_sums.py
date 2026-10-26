# Classify the 37 genuine seven-vertex Paley sums: for each S_H, exact membership in the span of
#   1, q, q^2, q^3, a(8.4.a.a), q a(8.4.a.a), a(8.6.a.a), a(4.6.a.a), |b(q)|^2 (32.3.c.a)
# over all computed primes; then the rank of the 37 sequences modulo that span (the number of independent new
# motive combinations).  Usage: classify_sums.py sums.jsonl [...].  mac-mini-2026-10-25-pyth.
import json, sys, os, subprocess, tempfile
from fractions import Fraction
G = json.load(open('genuine37.json'))
data = {}
for fn in sys.argv[1:]:
    for l in open(fn):
        d = json.loads(l); data[d['q']] = d['S']
qs = sorted(data)
P = max(qs) + 2
gpf = os.path.join(tempfile.gettempdir(), 'forms_cls.gp')
open(gpf, 'w').write(f'default(parisize,400000000);\nf=mfeigenbasis(mfinit([8,4,1],0))[1];g=mfeigenbasis(mfinit([8,6,1],0))[1];h=mfeigenbasis(mfinit([4,6,1],0))[1];b=mfeigenbasis(mfinit([32,3,-4],0))[1];\n'
                     f'print(mfcoefs(f,{P}));print(mfcoefs(g,{P}));print(mfcoefs(h,{P}));print(apply(x->lift(x),mfcoefs(b,{P})));\nquit;\n')
lines = subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()
f84, f86, f46 = (eval(l) for l in lines[-4:-1])
braw = lines[-1].replace('y', '*1j').replace('[', '[').replace('*1j^', '**')
# b coefficients are a + c*y with y = i; parse via sympy for safety
import sympy as sp
y = sp.symbols('y')
bl = [sp.sympify(t) for t in lines[-1].strip('[]').split(',')]
def absb2(q):
    v = sp.expand(bl[q]); re_ = v.subs(y, 0); im_ = sp.expand((v - re_) / y) if v != re_ else 0
    return int(re_ ** 2 + im_ ** 2)
B2 = {q: absb2(q) for q in qs}
assert absb2(3) == 16 and absb2(7) == 64 and absb2(11) == 16
base = {'1': lambda q: 1, 'q': lambda q: q, 'q^2': lambda q: q * q, 'q^3': lambda q: q ** 3,
        'a84': lambda q: f84[q], 'q a84': lambda q: q * f84[q], 'a86': lambda q: f86[q], 'a46': lambda q: f46[q],
        '|b|^2': lambda q: B2[q]}
names = list(base)
def solve_span(target):
    A = [[Fraction(base[n](q)) for n in names] + [Fraction(target[i])] for i, q in enumerate(qs)]
    m, n = len(A), len(names); r = 0; piv = []
    for c in range(n + 1):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                fct = A[i][c]; A[i] = [x - fct * y_ for x, y_ in zip(A[i], A[r])]
        piv.append(c); r += 1
    if n in piv: return None
    return {names[c]: A[i][n] for i, c in enumerate(piv)}
print(f"{len(qs)} primes, {qs[0]}..{qs[-1]}")
explained, new = [], []
for idx, H in enumerate(G):
    if H['all_even']: continue
    target = [data[q][idx] for q in qs]
    sol = solve_span(target)
    comp = sorted({(i, j) for i in range(7) for j in range(i + 1, 7)} - {tuple(e) for e in H['edges']})
    if sol:
        explained.append(idx)
        print(f"H{idx:02d} m={H['m']:+4d} deg={sorted(H['degrees'])} EXPLAINED: " + ' + '.join(f'({v})*{k}' for k, v in sol.items() if v))
    else:
        new.append(idx)
        print(f"H{idx:02d} m={H['m']:+4d} deg={sorted(H['degrees'])} NEW  complement={comp}")
print(f"explained {len(explained)}, new {len(new)}")
# rank of the new sequences modulo the base span (project out the span by solving least squares in exact arithmetic)
import numpy as np
from sympy import Matrix
Bm = Matrix([[base[n](q) for n in names] for q in qs])
cols = [Matrix([data[q][i] for q in qs]) for i in new]
Mfull = Bm.row_join(Matrix.hstack(*cols)) if cols else Bm
print("rank(base) =", Bm.rank(), "; rank(base + new sums) =", Mfull.rank(), "-> independent new directions:", Mfull.rank() - Bm.rank())
