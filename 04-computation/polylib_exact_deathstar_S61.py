#!/usr/bin/env python3
"""
death-star-2026-07-20-S61 (HYP-8265) -- exact multivariate polynomial arithmetic over
Gaussian rationals Q(i), pure Python (no sympy in sandbox). Supports the JC-counterexample
transport: Jacobians, determinants, matrix powers (nilpotency), evaluation at rational
points, Laplacian. Coefficients are a+bi with a,b in Fraction so det=const / nilpotency /
collision are all decided EXACTLY.
"""
from fractions import Fraction as Fr

class QI:
    __slots__ = ('a', 'b')
    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, Fr) else Fr(a)
        self.b = b if isinstance(b, Fr) else Fr(b)
    def __add__(s, o):
        o = asQI(o); return QI(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __sub__(s, o):
        o = asQI(o); return QI(s.a - o.a, s.b - o.b)
    def __neg__(s): return QI(-s.a, -s.b)
    def __mul__(s, o):
        o = asQI(o); return QI(s.a*o.a - s.b*o.b, s.a*o.b + s.b*o.a)
    __rmul__ = __mul__
    def conj(s): return QI(s.a, -s.b)
    def norm2(s): return s.a*s.a + s.b*s.b     # a^2+b^2 (rational)
    def inv(s):
        n = s.norm2()
        if n == 0: raise ZeroDivisionError
        return QI(s.a/n, -s.b/n)
    def __truediv__(s, o):
        o = asQI(o); return s * o.inv()
    def is_zero(s): return s.a == 0 and s.b == 0
    def __eq__(s, o):
        o = asQI(o); return s.a == o.a and s.b == o.b
    def __repr__(s):
        if s.b == 0: return f"{s.a}"
        if s.a == 0: return f"{s.b}i"
        return f"({s.a}{'+' if s.b>0 else ''}{s.b}i)"

def asQI(o):
    if isinstance(o, QI): return o
    return QI(o, 0)

I = QI(0, 1)

# ---- polynomials: dict {exponent-tuple (len NV): QI}, zero coeffs dropped ----
def pzero(): return {}
def pconst(c, NV):
    c = asQI(c)
    return {} if c.is_zero() else {tuple([0]*NV): c}
def pvar(i, NV):
    e = [0]*NV; e[i] = 1; return {tuple(e): QI(1, 0)}
def pclean(p):
    return {e: c for e, c in p.items() if not c.is_zero()}
def padd(p, q):
    r = dict(p)
    for e, c in q.items():
        r[e] = (r.get(e, QI(0,0)) + c)
    return pclean(r)
def pneg(p): return {e: -c for e, c in p.items()}
def psub(p, q): return padd(p, pneg(q))
def pscale(p, c):
    c = asQI(c)
    if c.is_zero(): return {}
    return {e: c*ce for e, ce in p.items()}
def pmul(p, q):
    r = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple(a+b for a, b in zip(e1, e2))
            r[e] = r.get(e, QI(0,0)) + c1*c2
    return pclean(r)
def ppow(p, n, NV):
    r = pconst(1, NV)
    for _ in range(n): r = pmul(r, p)
    return r
def pderiv(p, i):
    r = {}
    for e, c in p.items():
        if e[i] > 0:
            e2 = list(e); k = e2[i]; e2[i] -= 1
            r[tuple(e2)] = r.get(tuple(e2), QI(0,0)) + c*QI(k,0)
    return pclean(r)
def peval(p, pt):     # pt: list of QI (len NV)
    s = QI(0, 0)
    for e, c in p.items():
        t = c
        for i, k in enumerate(e):
            if k: t = t * (pt[i] if isinstance(pt[i], QI) else QI(pt[i],0))
            for _ in range(k-1) if False else []: pass
        # proper power:
        s = s  # placeholder
    # redo cleanly:
    s = QI(0,0)
    for e, c in p.items():
        t = c
        for i, k in enumerate(e):
            if k:
                base = pt[i] if isinstance(pt[i], QI) else QI(pt[i],0)
                for _ in range(k): t = t*base
        s = s + t
    return s
def pdeg(p):
    return max((sum(e) for e in p), default=-1)
def pis_zero(p): return len(pclean(p)) == 0
def pequal(p, q): return pis_zero(psub(p, q))

# ---- matrices of polynomials (list of lists) ----
def mzero(n, m, NV): return [[pzero() for _ in range(m)] for _ in range(n)]
def mmul(A, B):
    n = len(A); k = len(B); m = len(B[0])
    C = [[pzero() for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            s = pzero()
            for t in range(k):
                if A[i][t] and B[t][j]:
                    s = padd(s, pmul(A[i][t], B[t][j]))
            C[i][j] = s
    return C
def mis_zero(A): return all(pis_zero(A[i][j]) for i in range(len(A)) for j in range(len(A[0])))
def jacobian(F, NV):     # F: list of polys (len n); returns n x NV matrix
    return [[pderiv(F[i], j) for j in range(NV)] for i in range(len(F))]
def det3(M):
    a,b,c = M[0]; d,e,f = M[1]; g,h,i_ = M[2]
    return psub(padd(pmul(a, psub(pmul(e,i_), pmul(f,h))),
                     pmul(c, psub(pmul(d,h), pmul(e,g)))),
                pmul(b, psub(pmul(d,i_), pmul(f,g))))

if __name__ == "__main__":
    NV = 3
    x, y, z = pvar(0,NV), pvar(1,NV), pvar(2,NV)
    one = pconst(1,NV)
    u = padd(one, pmul(x,y))                       # 1+xy
    # F1 = u^3 z + y^2 u (4+3xy)
    F1 = padd(pmul(ppow(u,3,NV), z),
              pmul(pmul(ppow(y,2,NV), u), padd(pconst(4,NV), pscale(pmul(x,y),3))))
    # F2 = y + 3 x u^2 z + 3 x y^2 (4+3xy)
    F2 = padd(padd(y, pscale(pmul(pmul(x, ppow(u,2,NV)), z), 3)),
              pscale(pmul(pmul(x, ppow(y,2,NV)), padd(pconst(4,NV), pscale(pmul(x,y),3))), 3))
    # F3 = 2x - 3 x^2 y - x^3 z
    F3 = psub(psub(pscale(x,2), pscale(pmul(ppow(x,2,NV), y),3)), pmul(ppow(x,3,NV), z))
    F = [F1, F2, F3]
    print("degrees of F:", [pdeg(f) for f in F])

    JF = jacobian(F, NV)
    dJF = det3(JF)
    print("det JF =", dJF, " (expect -2)")

    pts = [[QI(0),QI(0),QI(Fr(-1,4))],
           [QI(1),QI(Fr(-3,2)),QI(Fr(13,2))],
           [QI(-1),QI(Fr(3,2)),QI(Fr(13,2))]]
    for p in pts:
        print("F(",[str(t) for t in p],") =", [str(peval(f,p)) for f in F])

    # ---- Yagzhev normalize: L = JF(0), G = L^{-1} F, H = G - x ----
    zero_pt = [QI(0),QI(0),QI(0)]
    L = [[peval(JF[i][j], zero_pt) for j in range(3)] for i in range(3)]
    print("L = JF(0) =", [[str(L[i][j]) for j in range(3)] for i in range(3)])
    # invert L (3x3 numeric QI)
    Linv = [[QI(0),QI(0),QI(Fr(1,2))],[QI(0),QI(1),QI(0)],[QI(1),QI(0),QI(0)]]
    # verify Linv*L = I
    prod = [[sum((Linv[i][k]*L[k][j] for k in range(3)), QI(0)) for j in range(3)] for i in range(3)]
    print("Linv*L =", [[str(prod[i][j]) for j in range(3)] for i in range(3)])
    G = [padd(padd(pscale(F[0], Linv[i][0].a if False else 0), pzero()), pzero()) for i in range(3)]
    G = []
    for i in range(3):
        gi = pzero()
        for k in range(3):
            gi = padd(gi, pscale(F[k], Linv[i][k]))
        G.append(gi)
    Xvars = [x,y,z]
    H = [psub(G[i], Xvars[i]) for i in range(3)]
    print("degrees of H (Yagzhev nonlinear part):", [pdeg(h) for h in H])
    # collision transports through G?
    print("G at collision pts (expect all equal):")
    for p in pts:
        print("  ", [str(peval(g,p)) for g in G])

    # nilpotency of JH
    JH = jacobian(H, NV)
    JH2 = mmul(JH, JH); JH3 = mmul(JH2, JH); JH4 = mmul(JH3, JH)
    print("JH nilpotent? JH^2==0:", mis_zero(JH2), " JH^3==0:", mis_zero(JH3), " JH^4==0:", mis_zero(JH4))
    # report char-poly-ish invariants of JH: trace, and whether entries nilpotent
    trJH = pzero()
    for i in range(3): trJH = padd(trJH, JH[i][i])
    print("trace(JH) =", trJH, " (nilpotent needs trace 0)")
