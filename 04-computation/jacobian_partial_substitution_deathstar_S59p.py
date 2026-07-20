#!/usr/bin/env python3
"""
death-star-S59p (HYP-8120) -- D3 retry: PARTIAL substitutions u -> (1+w) in
selected occurrences; for each pattern solve the 4th row P exactly (linear via
row-multilinearity) and report whether det == nonzero const is achievable,
plus collision transport.  Also allows rows 1-3 corrections (w - xy)*R_i with
R_i in a small basis, solved jointly (still linear: R_i sit in rows 1-3 ONE
row at a time? no -- corrections make det bilinear; we keep R = 0 and vary
patterns instead)."""
from fractions import Fraction as Fr
from itertools import product

def pmul(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(p+q for p, q in zip(ka, kb))
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
def peval(p, pt):
    s = Fr(0)
    for k, c in p.items():
        m = Fr(c)
        for xi, e in zip(pt, k):
            m *= Fr(xi)**e
        s += m
    return s

def V4(i):
    k = [0,0,0,0]; k[i] = 1
    return {tuple(k): 1}
x, y, z, w = (V4(i) for i in range(4))
ONE4 = {(0,0,0,0): 1}
t4 = pmul(x, y)
u_t = padd(ONE4, t4)      # 1 + xy
u_w = padd(ONE4, w)       # 1 + w
W_t = padd(psc(ONE4,4), psc(t4,3))
W_w = padd(psc(ONE4,4), psc(w,3))

def build(pattern):
    """pattern = (n1, useWB, n2, useW2) : n1 of the three u's in the F1 cube
    are 1+w; useWB: the u and W in F1's B-part use w; n2 of F2's two u's use w;
    useW2: F2's W uses w."""
    n1, useWB, n2, useW2 = pattern
    cube = ONE4
    for i in range(3):
        cube = pmul(cube, u_w if i < n1 else u_t)
    F1 = padd(pmul(cube, z),
              pmul(pmul(pmul(y,y), u_w if useWB else u_t), W_w if useWB else W_t))
    u2 = pmul(u_w if n2 >= 1 else u_t, u_w if n2 >= 2 else u_t)
    F2 = padd(y, psc(pmul(pmul(x, u2), z), 3),
              psc(pmul(pmul(x, pmul(y,y)), W_w if useW2 else W_t), 3))
    F3 = padd(psc(x,2), psc(pmul(pmul(x,x),y),-3), psc(pmul(pmul(pmul(x,x),x),z),-1))
    return F1, F2, F3

monos = [k for k in product(range(6), repeat=4) if sum(k) <= 5]
NU = len(monos)

def try_pattern(pattern, verbose=False):
    F1, F2, F3 = build(pattern)
    Jt = [[pdiff(Fi, j) for j in range(4)] for Fi in (F1, F2, F3)]
    def minor3(cols):
        (c0, c1, c2) = cols
        return padd(
            pmul(Jt[0][c0], padd(pmul(Jt[1][c1], Jt[2][c2]), psc(pmul(Jt[1][c2], Jt[2][c1]), -1))),
            psc(pmul(Jt[0][c1], padd(pmul(Jt[1][c0], Jt[2][c2]), psc(pmul(Jt[1][c2], Jt[2][c0]), -1))), -1),
            pmul(Jt[0][c2], padd(pmul(Jt[1][c0], Jt[2][c1]), psc(pmul(Jt[1][c1], Jt[2][c0]), -1))))
    cof = []
    for j in range(4):
        cols = [cc for cc in range(4) if cc != j]
        cof.append(psc(minor3(cols), (-1)**(4 + (j+1))))
    rows = {}
    for ui, m in enumerate(monos):
        pm_ = {m: 1}
        contrib = padd(*[pmul(pdiff(pm_, j), cof[j]) for j in range(4)])
        for k, cc in contrib.items():
            rows.setdefault(k, [0]*NU)
            rows[k][ui] += cc
    M = []
    for k, rr in rows.items():
        if k == (0,0,0,0): continue
        M.append([Fr(v) for v in rr])
    aug = [row[:] for row in M]
    r = 0; piv = []
    for cc in range(NU):
        pv = next((i for i in range(r, len(aug)) if aug[i][cc] != 0), None)
        if pv is None: continue
        aug[r], aug[pv] = aug[pv], aug[r]
        inv = aug[r][cc]
        aug[r] = [xx/inv for xx in aug[r]]
        for i in range(len(aug)):
            if i != r and aug[i][cc] != 0:
                f = aug[i][cc]
                aug[i] = [aa - f*bb for aa, bb in zip(aug[i], aug[r])]
        r += 1; piv.append(cc)
    const_row = rows.get((0,0,0,0), [0]*NU)
    best = None
    for fc in [cc for cc in range(NU) if cc not in piv]:
        v = [Fr(0)]*NU; v[fc] = Fr(1)
        for i, cc in enumerate(piv): v[cc] = -aug[i][fc]
        cv = sum(Fr(const_row[i])*v[i] for i in range(NU))
        if cv != 0:
            best = (v, cv); break
    if best is None:
        return None
    v, cv = best
    Psol = {monos[i]: v[i] for i in range(NU) if v[i] != 0}
    # verify + collisions
    Ft4 = Psol
    Jt4 = [pdiff(Ft4, j) for j in range(4)]
    det4 = padd(*[pmul(Jt4[j], cof[j]) for j in range(4)])
    nonconst = {k: v_ for k, v_ in det4.items() if k != (0,0,0,0)}
    Pt = (1, Fr(-3,2), Fr(13,2), Fr(-3,2))
    Qt = (-1, Fr(3,2), Fr(13,2), Fr(-3,2))
    imP = tuple(peval(Fi, Pt) for Fi in (F1, F2, F3, Ft4))
    imQ = tuple(peval(Fi, Qt) for Fi in (F1, F2, F3, Ft4))
    deg = max(sum(k) for Fi in (F1,F2,F3,Ft4) for k in Fi)
    return dict(pattern=pattern, det_const=cv, det_clean=not nonconst,
                collide=(imP == imQ), image=imP, deg=deg, Pterms=len(Psol))

patterns = [(3,1,2,1), (3,0,0,0), (1,0,0,0), (2,0,1,0), (3,1,0,0), (0,1,0,1),
            (2,1,1,1), (1,1,1,1), (3,0,2,0), (0,0,1,0)]
for pat in patterns:
    res = try_pattern(pat)
    if res:
        print(f"  pattern {pat}: det = {res['det_const']} clean={res['det_clean']} "
              f"collide={res['collide']} deg={res['deg']} |P|={res['Pterms']}")
    else:
        print(f"  pattern {pat}: NO Keller completion (det forced 0)")
