#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- orbit-tangent computation.
Question: of the 9-dim kernel of the det-constancy linearization at F (within
the 22-param six-function ansatz), how much is the tangent of the
(equivariant Aut_target) o F o (equivariant Aut_source) orbit?
Moduli tangent dim = 9 - dim(orbit tangent within the ansatz space).

Source equivariant fields: W = x p(t,s) dx + (y q + xz m) dy + (z r + y^2 n) dz
Target equivariant fields: V = (a ph + b^2 uh) da + (b qh + ac mh) db + (c rh) dc
Orbit tangents: JF.W  and  V(F).
Each tangent is itself of the six-function shape; we extract its
(A..E1)-coefficient vector, keep those inside the ansatz degree caps,
and compute the rank.
"""
from fractions import Fraction as Fr

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = (ka[0]+kb[0], ka[1]+kb[1], ka[2]+kb[2])
            v = r.get(k, 0) + ca*cb
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def padd(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            v = r.get(k, 0) + c
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def psc(p, s): return {k: c*s for k, c in p.items() if c*s != 0}
def pdiff(p, i):
    r = {}
    for k, c in p.items():
        if k[i] > 0:
            k2 = list(k); k2[i] -= 1
            r[tuple(k2)] = c*k[i]
    return r

X = {(1,0,0):1}; Y = {(0,1,0):1}; Z = {(0,0,1):1}; ONE = {(0,0,0):1}
t = pmul(X, Y); s = pmul(pmul(X,X), Z)
u = padd(ONE, t)
A  = pmul(pmul(u,u),u)
B  = pmul(u, padd(psc(ONE,4), psc(t,3)))
C  = padd(ONE, psc(t,12), psc(pmul(t,t),9))
D  = psc(pmul(u,u),3)
E0 = padd(psc(ONE,2), psc(t,-3))
F1 = padd(pmul(Z, A), pmul(pmul(Y,Y), B))
F2 = padd(pmul(Y, C), pmul(pmul(X,Z), D))
F3 = pmul(X, padd(E0, psc(s,-1)))
F = [F1, F2, F3]
J = [[pdiff(Fi, i) for i in range(3)] for Fi in F]

# ---- basis of (t,s)-multiplier polys up to the caps
TS = []
for i in range(0, 4):
    for j in range(0, 2):
        if i + 2*j <= 3:
            m = ONE
            for _ in range(i): m = pmul(m, t)
            for _ in range(j): m = pmul(m, s)
            TS.append((f"t^{i}s^{j}", m))

tangents = []
labels = []
# source fields
for nm, mult in TS:
    fields = [
        ("xp.dx",  [pmul(X, mult), {}, {}]),
        ("yq.dy",  [{}, pmul(Y, mult), {}]),
        ("zr.dz",  [{}, {}, pmul(Z, mult)]),
        ("xzm.dy", [{}, pmul(pmul(X,Z), mult), {}]),
        ("y2n.dz", [{}, {}, pmul(pmul(Y,Y), mult)]),
    ]
    for fn, W in fields:
        dF = [padd(*[pmul(J[i][j], W[j]) for j in range(3)]) for i in range(3)]
        tangents.append(dF); labels.append(f"src:{fn}*{nm}")
# target fields (compose with F): need target invariants at F
th = pmul(F2, F3); sh = pmul(F1, pmul(F3, F3))
TSh = []
for i in range(0, 4):
    for j in range(0, 2):
        if i + 2*j <= 3:
            m = ONE
            for _ in range(i): m = pmul(m, th)
            for _ in range(j): m = pmul(m, sh)
            TSh.append((f"th^{i}sh^{j}", m))
for nm, mult in TSh:
    fields = [
        ("a.da",   [pmul(F1, mult), {}, {}]),
        ("b2.da",  [pmul(pmul(F2,F2), mult), {}, {}]),
        ("b.db",   [{}, pmul(F2, mult), {}]),
        ("ac.db",  [{}, pmul(pmul(F1,F3), mult), {}]),
        ("c.dc",   [{}, {}, pmul(F3, mult)]),
    ]
    for fn, V in fields:
        tangents.append(V); labels.append(f"tgt:{fn}*{nm}")

# ---- extract six-function coefficients (or reject if outside ansatz caps)
# caps: A deg<=4, B<=3, C<=3, D<=3, E0<=2, E1<=1 (t-degree), E has no s^2 term
def extract(dF):
    vec = {}
    ok = True
    used = [dict(d) for d in dF]
    def take(comp, xe_off, ye_off, ze_off, name, cap):
        nonlocal ok
        # monomials x^(k+xe_off) y^(k+ye_off) z^(ze_off) -> coeff of t^k
        for k in list(used[comp].keys()):
            pass
        for kk in range(0, cap+1):
            key = (kk+xe_off, kk+ye_off, ze_off)
            c = used[comp].pop(key, 0)
            if c: vec[f"{name}{kk}"] = c
    take(0, 0, 0, 1, "A", 4)     # z t^k
    take(0, 0, 2, 0, "B", 3)     # y^2 t^k
    take(1, 0, 1, 0, "C", 3)     # y t^k
    take(1, 1, 0, 1, "D", 3)     # xz t^k
    take(2, 1, 0, 0, "E0", 2)    # x t^k
    # E1: -x s t^k = -x^(3+k) y^k z
    for kk in range(0, 2):
        key = (3+kk, kk, 1)
        c = used[2].pop(key, 0)
        if c: vec[f"E1{kk}"] = -c
    if any(used[i] for i in range(3)):
        ok = False
    return ok, vec

names = [f"A{k}" for k in range(5)] + [f"B{k}" for k in range(4)] + \
        [f"C{k}" for k in range(4)] + [f"D{k}" for k in range(4)] + \
        [f"E0{k}" for k in range(3)] + [f"E1{k}" for k in range(2)]
rowsM = []
kept = []
for dF, lb in zip(tangents, labels):
    ok, vec = extract(dF)
    if ok and vec:
        rowsM.append([Fr(vec.get(nm, 0)) for nm in names])
        kept.append(lb)

def rank_frac(Mat):
    Mat = [row[:] for row in Mat]
    r = 0; cols = len(Mat[0]) if Mat else 0
    for c in range(cols):
        piv = next((i for i in range(r, len(Mat)) if Mat[i][c] != 0), None)
        if piv is None: continue
        Mat[r], Mat[piv] = Mat[piv], Mat[r]
        pv = Mat[r][c]
        Mat[r] = [x/pv for x in Mat[r]]
        for i in range(len(Mat)):
            if i != r and Mat[i][c] != 0:
                f = Mat[i][c]
                Mat[i] = [a - f*b for a, b in zip(Mat[i], Mat[r])]
        r += 1
    return r

rk = rank_frac(rowsM)
print(f"orbit tangents generated: {len(tangents)}, inside ansatz caps: {len(rowsM)}")
print(f"ORBIT TANGENT DIM (within ansatz) = {rk}")
print(f"=> MODULI TANGENT >= 9 - {rk} = {9 - rk}")
print("kept generators:", len(kept))
