# Audit of the (uncommitted) note section 10 table: rho over all regular d-tournaments (gentourng -d -D), d = 7, 9.
import sys, subprocess
from fractions import Fraction
import sympy as sp
def rho_exact(A):
    d = len(A)
    S = sp.Matrix(d, d, lambda i, j: A[i][j] - A[j][i])
    x = sp.symbols('x')
    cp = S.charpoly(x).as_expr()
    num = cp.subs(x, d)
    den = sp.Abs(sp.expand(cp.subs(x, sp.I * d)))
    return sp.Rational(num) / sp.Rational(den)
def read(d):
    out = subprocess.run(['gentourng', '-q', f'-d{(d-1)//2}', f'-D{(d-1)//2}', str(d)], capture_output=True, text=True).stdout.split()
    for line in out:
        A = [[0] * d for _ in range(d)]
        b = 0
        for i in range(d):
            for j in range(i + 1, d):
                if line[b] == '1': A[i][j] = 1
                else: A[j][i] = 1
                b += 1
        yield A
def is_c3c3(A):
    # C3[C3] test: vertex set splits into 3 cyclic triangles that are modules
    d = len(A)
    import itertools
    for blocks in [None]:
        pass
    # module test: find partition into 3 modules of size 3 each inducing C3
    verts = range(d)
    def is_module(M):
        out = [v for v in verts if v not in M]
        return all(len({A[m][v] for m in M}) == 1 for v in out)
    mods = [M for M in itertools.combinations(verts, 3) if is_module(M) and sum(A[a][b] for a in M for b in M) == 3 and all(sum(A[a][b] for b in M) == 1 for a in M)]
    for M1, M2, M3 in itertools.combinations(mods, 3):
        if len(set(M1) | set(M2) | set(M3)) == 9: return True
    return False
for d in (7, 9):
    vals = []
    for A in read(d):
        vals.append((rho_exact(A), is_c3c3(A) if d == 9 else False))
    vals.sort(key=lambda t: -t[0])
    print(f"d={d}: regular classes = {len(vals)}")
    for r, c in vals[:4]:
        print(f"   rho = {r} = {float(r):.10f}{'   <- C3[C3]' if c else ''}")
    print(f"   min rho = {vals[-1][0]} = {float(vals[-1][0]):.10f}; doubly-regular bound ((d+1)/(d-1))^((d-1)/2) = {float(Fraction(d+1, d-1)**((d-1)//2)):.10f}")
