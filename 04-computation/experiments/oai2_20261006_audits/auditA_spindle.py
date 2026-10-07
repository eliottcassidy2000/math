"""Audit A: generalized spindle in L_N = Q(sqrt-3, sqrt m), m = sqfree(-(4N-1)), N = 3n (n Loeschian).
Exact arithmetic: elements a + b s3 + c sm + d s3 sm (s3^2 = -3, sm^2 = m), rationals. Complex conj negates s3, sm.
Patch = Eisenstein points (w = (1+s3)/2) within radius R of 0; tip v with |v|^2 = N, v in sqrt(-3) Z[w].
Graph = patch U rho*patch (rho = omega_N = (2N-1 + sqrt(-(4N-1)))/(2N)) with ALL exact unit-distance edges; test 3-colourability."""
import sys, itertools
from fractions import Fraction as Fr
from pysat.solvers import Solver
from sympy import factorint

def sqf(n):
    s = -1 if n < 0 else 1; r = 1
    for p, e in factorint(abs(n)).items():
        if e % 2: r *= p
    return s * r

class E:
    __slots__ = ('v',)
    m = None
    def __init__(self, a, b=0, c=0, d=0): self.v = (Fr(a), Fr(b), Fr(c), Fr(d))
    def __add__(s, o): return E(*[x + y for x, y in zip(s.v, o.v)])
    def __sub__(s, o): return E(*[x - y for x, y in zip(s.v, o.v)])
    def __mul__(s, o):
        a, b, c, d = s.v; e, f, g, h = o.v; m = E.m
        # s3^2=-3, sm^2=m, (s3 sm)^2 = -3m ; s3*(s3 sm) = -3 sm ; sm*(s3 sm) = m s3
        return E(a*e - 3*b*f + m*c*g - 3*m*d*h,
                 a*f + b*e + m*(c*h + d*g),
                 a*g + c*e - 3*(b*h + d*f),
                 a*h + d*e + b*g + c*f)
    def conj(s): a, b, c, d = s.v; return E(a, -b, -c, d)
    def key(s): return s.v
    def __eq__(s, o): return s.v == o.v
    def __hash__(s): return hash(s.v)

def run(N, R):
    m = sqf(-(4 * N - 1)); E.m = m
    k = Fr(4 * N - 1, 1)
    # sqrt(-(4N-1)) = t * sm with t^2 * m = -(4N-1)
    t2 = Fr(-(4 * N - 1), m); t = Fr(int(round(float(t2) ** 0.5)))
    assert t * t == t2
    rho = E(Fr(2 * N - 1, 2 * N), 0, t / (2 * N), 0)
    one = E(1)
    assert (rho * rho.conj()) == one
    w = E(Fr(1, 2), Fr(1, 2))                         # e^{i pi/3}
    pts = []
    for a in range(-R - 2, R + 3):
        for b in range(-R - 2, R + 3):
            if a * a + a * b + b * b <= R * R:
                pts.append(E(a) + E(b) * w)
    tips = [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1) if a * a + a * b + b * b == N and (a - b) % 3 == 0]
    assert tips, "no tip"
    V = list(dict.fromkeys(pts + [rho * p for p in pts]))
    idx = {p: i for i, p in enumerate(V)}
    # exact edges: brute force O(n^2) (n small)
    Ed = []
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            dz = V[i] - V[j]
            if dz * dz.conj() == one:
                Ed.append((i, j))
    s = Solver(name='glucose4')
    var = lambda i, c: 3 * i + c + 1
    for i in range(len(V)): s.add_clause([var(i, c) for c in range(3)])
    for i, j in Ed:
        for c in range(3): s.add_clause([-var(i, c), -var(j, c)])
    s.add_clause([var(idx[E(0)], 0)])
    r = s.solve(); s.delete()
    print(f"N={N}: field Q(sqrt-3, sqrt{m}); tips {tips[:3]}...; |V|={len(V)} |E|={len(Ed)}; 3-colourable: {r}", flush=True)
    return r

for N, R in [(3, 2), (9, 3), (21, 5), (27, 6)]:
    run(N, R)
