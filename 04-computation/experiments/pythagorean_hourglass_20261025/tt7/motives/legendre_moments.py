# chi-twisted moments of the Legendre family, M_k(q) = sum_{l != 0,1} chi(l) a_l^k, a_l = -sum_t chi(t(t-1)(t-l)),
# at odd primes q < P.  The local system Sym^(2j) H^1 (x) L_chi on P^1 - {0,1,inf} has parabolic H^1 of dimension
# 2j (Euler-Poincare: rank 2j+1 times -1, minus one invariant at lambda = 1), weight 2j+1: j = 1 gives the Apery form
# 8.4.a.a (the wheel); j = 2 a 4-dimensional weight-5 motive.  Test M_4 against span{q^i, q^i a(8.4.a.a), 8.6.a.a,
# 4.6.a.a}.  mac-mini-2026-10-25-pyth.
import sys, os, subprocess, tempfile
import numpy as np
from fractions import Fraction
from sympy import primerange
P = int(sys.argv[1]) if len(sys.argv) > 1 else 600
def moments(q, ks=(2, 4, 6)):
    chi = np.full(q, -1, dtype=np.int64); chi[0] = 0; chi[[(x * x) % q for x in range(1, q)]] = 1
    t = np.arange(q)
    base = chi[(t * (t - 1)) % q]                         # chi(t(t-1))
    a = np.zeros(q, dtype=np.int64)
    for l in range(2, q):
        a[l] = -int((base * chi[(t - l) % q]).sum())
    lam = np.arange(2, q)
    return {k: int((chi[lam] * a[lam] ** k).sum()) for k in ks}, a
gpf = os.path.join(tempfile.gettempdir(), 'forms_8_4_6.gp')
open(gpf, 'w').write(f'default(parisize,400000000);\nf=mfeigenbasis(mfinit([8,4,1],0))[1];g=mfeigenbasis(mfinit([8,6,1],0))[1];h=mfeigenbasis(mfinit([4,6,1],0))[1];\n'
                     f'print(mfcoefs(f,{P}));print(mfcoefs(g,{P}));print(mfcoefs(h,{P}));\nquit;\n')
lines = subprocess.run(['gp', '-q', gpf], capture_output=True, text=True).stdout.strip().splitlines()
f84, f86, f46 = (eval(l) for l in lines[-3:])
assert f84[3] == -4 and f86[3] == 20 and f46[3] == -12
rows = []
for q in primerange(5, P):
    M, _ = moments(q)
    rows.append((q, M))
print("M_2 = -q - 1 - a_q(8.4.a.a) (the wheel / Ahlgren-Ono):", all(M[2] == -q - 1 - f84[q] for q, M in rows))
def in_span(colfun, names, target):
    A = [[Fraction(colfun[n](q)) for n in names] + [Fraction(target(q, M))] for q, M in rows]
    m, n = len(A), len(names); r = 0; piv = []
    for c in range(n + 1):
        pr = next((i for i in range(r, m) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]; inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] != 0:
                fct = A[i][c]; A[i] = [x - fct * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    if n in piv: return None
    return {names[c]: A[i][n] for i, c in enumerate(piv)}
cols = {'1': lambda q: 1, 'q': lambda q: q, 'q^2': lambda q: q * q, 'q^3': lambda q: q ** 3,
        'a84': lambda q: f84[q], 'q a84': lambda q: q * f84[q], 'a86': lambda q: f86[q], 'a46': lambda q: f46[q],
        'chi-4': lambda q: 1 if q % 4 == 1 else -1, 'q chi-4': lambda q: q * (1 if q % 4 == 1 else -1),
        'q^2 chi-4': lambda q: q * q * (1 if q % 4 == 1 else -1)}
for label, names in [('poly + a84 terms', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84']),
                     ('+ 8.6.a.a', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86']),
                     ('+ 4.6.a.a', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a46']),
                     ('+ both weight 6', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86', 'a46']),
                     ('+ both weight 6 + chi_-4 polys', ['1', 'q', 'q^2', 'q^3', 'a84', 'q a84', 'a86', 'a46', 'chi-4', 'q chi-4', 'q^2 chi-4'])]:
    sol = in_span(cols, names, lambda q, M: M[4])
    print(f"M_4 in span({label}) over {len(rows)} primes:", ('EXACT ' + str({k: str(v) for k, v in sol.items() if v})) if sol else 'NOT')
