# rho over all 1223 regular 11-tournaments (gentourng -d5 -D5 11, file argument); numeric, then exact for the top 3.
import sys
import numpy as np
import sympy as sp
d = 11
lines = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith('>')]
def adj(line):
    A = np.zeros((d, d), dtype=int); b = 0
    for i in range(d):
        for j in range(i + 1, d):
            if line[b] == '1': A[i, j] = 1
            else: A[j, i] = 1
            b += 1
    return A
vals = []
for l in lines:
    A = adj(l); S = A - A.T
    ev = np.linalg.eigvals(S.astype(float))
    s2 = np.sort(np.abs(ev.imag))[1::2] ** 2 if False else None
    sig = np.abs(ev.imag)
    r = np.prod([(d * d + s * s) / (d * d - s * s) for s in sig if s > 1e-9]) ** 0.5
    vals.append((r, l))
vals.sort(key=lambda t: -t[0])
print("d=11 regular classes:", len(vals))
for r, l in vals[:3]:
    A = adj(l); S = sp.Matrix(A - A.T); x = sp.symbols('x'); cp = S.charpoly(x).as_expr()
    ex = sp.Rational(cp.subs(x, d)) / sp.Abs(sp.expand(cp.subs(x, sp.I * d)))
    print(f"   rho = {ex} = {float(ex):.10f}")
print(f"   min rho = {vals[-1][0]:.10f}; (12/10)^5 = {(12/10)**5:.10f}")
