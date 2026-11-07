# mac-mini-2026-11-06-microcosm. Data (note): minimum of phi over strong tournaments of orders 5..9 via nauty gentourng (no conjecture).
import sys, subprocess, numpy as np
from fractions import Fraction
def parse(line, n):
    # gentourng ascii output: upper triangle row by row, '1' means i -> j
    bits = line.strip(); A = np.zeros((n, n), dtype=int); k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[k] == '1': A[i, j] = 1
            else: A[j, i] = 1
            k += 1
    return A
n = int(sys.argv[1])
out = subprocess.run(['/opt/homebrew/bin/gentourng', '-c', '-q', str(n)], capture_output=True, text=True).stdout.split('\n')
lines = [l for l in out if l and not l.startswith('>')]
best = None; vals = []
u = np.ones(n)
for l in lines:
    A = parse(l, n); M = (A - A.T).astype(float)
    v = float(u @ np.linalg.solve(np.eye(n) - M, u)); vals.append(v)
    if best is None or v < best[0] - 1e-12: best = (v, A)
v, A = best
M = (A - A.T)
# exact value via fractions (sympy-free Gaussian elimination)
import sympy as sp
ex = (sp.ones(1, n) * (sp.eye(n) - sp.Matrix(M)).inv() * sp.ones(n, 1))[0, 0]
print('n = %d: %d strong classes; min phi = %.6f (exact %s); max %.4f; minimiser scores %s'
      % (n, len(lines), v, sp.nsimplify(ex), max(vals), sorted(A.sum(axis=1).tolist())))
cnt = sum(1 for x in vals if abs(x - v) < 1e-9)
print('   number of strong classes attaining the minimum:', cnt)
