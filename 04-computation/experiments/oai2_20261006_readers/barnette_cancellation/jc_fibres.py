#!/usr/bin/env python3
"""Cancellation-type invariants of the THM-1300 Keller map F: C^3 -> C^3 (det JF = -2).

F1 = u^3 z + y^2 u (4+3xy), F2 = y + 3x u^2 z + 3x y^2 (4+3xy), F3 = 2x - 3x^2 y - x^3 z, u = 1+xy.
Checks: (1) det JF = -2; factorizations of F1, F2, F3; (2) brute-force point counts of the fibres
F_i = c over F_p for c != 0 against the closed forms q^2-q+1, q^2+1, q^2-q; (3) the zero fibres;
(4) the dual commuting fields V_j (V_j F_i = delta_ij) are not locally nilpotent: degree growth of
V_j^k(x), and the slice criterion (V_j LND => C^3 = ker(V_j)[F_j] => #fibre = q^2)."""
from sympy import symbols, expand, factor_list, Matrix, Poly, simplify, diff, factor
x, y, z = symbols('x y z')
u = 1 + x*y
F1 = expand(u**3*z + y**2*u*(4 + 3*x*y))
F2 = expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y))
F3 = expand(2*x - 3*x**2*y - x**3*z)
Fs = [F1, F2, F3]
JF = Matrix([[diff(f, v) for v in (x, y, z)] for f in Fs])
print("det JF =", simplify(JF.det()))
for i, f in enumerate(Fs, 1):
    print(f"F{i} factor_list:", factor_list(f))
# C*-weights: (x,y,z) -> (lx, y/l, z/l^2)
import sympy
l = symbols('l')
for i, f in enumerate(Fs, 1):
    g = expand(f.subs({x: l*x, y: y/l, z: z/l**2}, simultaneous=True))
    for w in range(-3, 4):
        if expand(g - l**w*f) == 0:
            print(f"F{i} has C*-weight {w}")

def count(p, f, c):
    fl = sympy.lambdify((x, y, z), f, 'math')
    cnt = 0
    for a in range(p):
        for b in range(p):
            for cc in range(p):
                if int(fl(a, b, cc)) % p == c % p:
                    cnt += 1
    return cnt

# faster: use the z-affine structure F_i = A_i(x,y) z + B_i(x,y), but brute force for small p as a check
for p in (5, 7, 11, 13):
    q = p
    res = []
    for i, f in enumerate(Fs, 1):
        cs = [count(p, f, c) for c in (1, 2, 3)]
        res.append((i, cs, count(p, f, 0)))
    print(f"p={p}: fibre counts c=1,2,3 / c=0:", [(f'F{i}', cs, c0) for i, cs, c0 in res],
          " closed forms q^2-q+1, q^2+1, q^2-q =", (q*q - q + 1, q*q + 1, q*q - q))

# dual fields: V = (JF^T)^{-1} applied to gradients; V_j = sum_k B_jk d/dx_k with B = (JF^T)^{-1}
B = (JF.T).inv()
B = B.applyfunc(lambda e: expand(simplify(e)))
def V(j, e):
    return expand(sum(B[j, k]*diff(e, v) for k, v in enumerate((x, y, z))))
for j in range(3):
    chk = [expand(V(j, f)) for f in Fs]
    degs = []
    e = x
    for k in range(5):
        e = V(j, e)
        degs.append(Poly(e, x, y, z).total_degree() if e != 0 else -1)
        if e == 0:
            break
    ey = y; degy = []
    for k in range(4):
        ey = V(j, ey)
        degy.append(Poly(ey, x, y, z).total_degree() if ey != 0 else -1)
        if ey == 0:
            break
    print(f"V{j+1}: V(F1,F2,F3) = {chk}; deg V^k(x), k=1..: {degs}; deg V^k(y): {degy}")
