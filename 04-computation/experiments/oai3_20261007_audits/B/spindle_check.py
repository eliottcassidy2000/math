# Generalized Moser spindle (THM-4558 construction) for N = 3n, n Loeschian:
#   v = sqrt(-3)*w in Z[omega], |v|^2 = N; omega_N = ((2N-1) + sqrt(-(4N-1)))/(2N) (unit, Re = 1 - 1/(2N)).
# Proof graph: triangulated Eisenstein parallelogram P containing 0 and v, its rotation omega_N*P (sharing 0),
# and the cross edge v ~ omega_N v.  Check: |v - omega_N v| = 1 exactly; 3-colourability by SAT (expect UNSAT).
import sympy as sp
from pysat.solvers import Cadical153
from math import isqrt
def loeschian_vec(n):          # a,b with a^2 - a b + b^2 = n  (omega = e^{2 pi i/3})
    for a in range(0, 2*isqrt(n)+3):
        for b in range(0, 2*isqrt(n)+3):
            if a*a - a*b + b*b == n: return a, b
    return None
def check(N):
    n = N // 3
    w = loeschian_vec(n)
    if w is None: return None
    a, b = w
    # sqrt(-3) = 1 + 2 omega ; v = (1+2w)(a + b w) with w^2 = -1 - w
    # (1 + 2w)(a + b w) = a + b w + 2a w + 2b w^2 = a + (b+2a) w + 2b(-1-w) = (a - 2b) + (2a - b) w
    va, vb = a - 2*b, 2*a - b
    assert va*va - va*vb + vb*vb == N and (va + vb) % 3 == 0
    om = sp.Rational(-1, 2) + sp.sqrt(-3)/2
    v = va + vb*om
    wN = (sp.Integer(2*N - 1) + sp.sqrt(-(4*N - 1))) / (2*N)
    unit = sp.simplify(sp.expand(wN*sp.conjugate(wN)))
    d2 = sp.simplify(sp.expand((v - wN*v)*sp.conjugate(v - wN*v)))
    # proof graph: parallelogram in (a,b) lattice coords covering 0 and v
    a0, a1 = min(0, va) - 1, max(0, va) + 1; b0, b1 = min(0, vb) - 1, max(0, vb) + 1
    pts = [(x, y) for x in range(a0, a1+1) for y in range(b0, b1+1)]
    idx = {}
    for p in pts: idx[('P',) + p] = len(idx)
    for p in pts:
        if p == (0, 0): idx[('R', 0, 0)] = idx[('P', 0, 0)]
        else: idx[('R',) + p] = len(idx)
    edges = set()
    for tag in 'PR':
        for (x, y) in pts:
            for dx, dy in ((1,0),(0,1),(1,1)):
                q = (x+dx, y+dy)
                if q in set(pts) if False else (a0 <= q[0] <= a1 and b0 <= q[1] <= b1):
                    edges.add((idx[(tag, x, y)], idx[(tag,)+q]))
    edges.add((idx[('P', va, vb)], idx[('R', va, vb)]))
    nv = len(set(idx.values()))
    s = Cadical153()
    var = lambda i, c: 3*i + c + 1
    for i in range(nv): s.add_clause([var(i, c) for c in range(3)])
    for (i, j) in edges:
        for c in range(3): s.add_clause([-var(i, c), -var(j, c)])
    s.add_clause([var(0, 0)])
    sat = s.solve()
    return dict(N=N, n=n, v=(va, vb), sqrt_part=4*N-1, omega_N_unit=unit, dist2=d2, vertices=nv, edges=len(edges), three_colourable=sat)
for N in (3, 9, 144, 723, 3600):
    print(check(N))
