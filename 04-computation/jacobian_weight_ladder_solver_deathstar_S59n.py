#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- weight-ladder solver, v2 (row-block).

General equivariant map for weights (1,-1,-k):
  Row1: F1 = z A(t) + y^k B(t)      (weight of F1 = weight of z = -k)
  Row2: F2 = y C(t) + x^{k-1} z D(t)
  Row3: F3 = x (E0(t) - E1(t) s_k),   t = xy, s_k = x^k z.
det J is multilinear in rows => fixing TWO rows makes det == const LINEAR in
the third row's coefficients.  Exact row solves + block iteration.
Observed anchor law from the k=2 witness: D = A' (3u^2 = (u^3)').
"""
from fractions import Fraction as Fr

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

def make_F(k, NV, rowpolys):
    """rowpolys: dict of engine polys A,B,C,D,E0,E1 (already lifted, NV vars)."""
    def var(i):
        kk = [0]*NV; kk[i] = 1
        return {tuple(kk): 1}
    X, Y, Z = var(0), var(1), var(2)
    A,B,C,D,E0,E1 = (rowpolys[n] for n in ("A","B","C","D","E0","E1"))
    Yk1 = {tuple([0]*NV):1}
    for _ in range(k): Yk1 = pmul(Yk1, Y)
    Xk_1 = {tuple([0]*NV):1}
    for _ in range(k-1): Xk_1 = pmul(Xk_1, X)
    sk = pmul(pmul(Xk_1, X), Z)
    F1 = padd(pmul(Z, A), pmul(Yk1, B))
    F2 = padd(pmul(Y, C), pmul(pmul(Xk_1, Z), D))
    F3 = pmul(X, padd(E0, psc(pmul(sk, E1), -1)))
    return F1, F2, F3

def detJ(F1, F2, F3):
    J = [[pdiff(F, i) for i in range(3)] for F in (F1, F2, F3)]
    return padd(
        pmul(J[0][0], padd(pmul(J[1][1], J[2][2]), psc(pmul(J[1][2], J[2][1]), -1))),
        psc(pmul(J[0][1], padd(pmul(J[1][0], J[2][2]), psc(pmul(J[1][2], J[2][0]), -1))), -1),
        pmul(J[0][2], padd(pmul(J[1][0], J[2][1]), psc(pmul(J[1][1], J[2][0]), -1))))

def lift_known(cl, NV):
    out = {}
    for j, c in enumerate(cl):
        if c:
            key = [0]*NV; key[0] = j; key[1] = j
            out[tuple(key)] = out.get(tuple(key), 0) + c
    return out  # t^j = x^j y^j

def solve_row(k, row, known, degs):
    """row in {1,2,3}; known: dict name->coeff list for the OTHER rows' four
    functions; degs: (d1,d2) degrees for this row's two unknowns.
    Returns list of solutions: (particular, basis, names)."""
    pair = {1: ("A","B"), 2: ("C","D"), 3: ("E0","E1")}[row]
    n1, n2 = degs[0]+1, degs[1]+1
    NU = n1 + n2
    NV = 3 + NU
    def var(i):
        kk = [0]*NV; kk[i] = 1
        return {tuple(kk): 1}
    polys = {}
    for nm in ("A","B","C","D","E0","E1"):
        if nm in known:
            polys[nm] = lift_known(known[nm], NV)
    t = pmul(var(0), var(1))
    for idx, nm in enumerate(pair):
        out = {}
        tp = {tuple([0]*NV): 1}
        base = 3 if idx == 0 else 3 + n1
        for j in range((n1 if idx == 0 else n2)):
            out = padd(out, pmul(var(base+j), tp))
            tp = pmul(tp, t)
        polys[nm] = out
    F1, F2, F3 = make_F(k, NV, polys)
    det = detJ(F1, F2, F3)
    rows = {}
    for key, c in det.items():
        xyz = key[:3]; unk = key[3:]
        tot = sum(unk)
        assert tot <= 1, f"nonlinear: {key}"
        rows.setdefault(xyz, [0]*(NU+1))
        if tot == 0: rows[xyz][NU] += c
        else: rows[xyz][unk.index(1)] += c
    M = []; rhs = []
    for xyz, rr in rows.items():
        if xyz == (0,0,0): continue
        M.append([Fr(c) for c in rr[:NU]]); rhs.append(Fr(-rr[NU]))
    aug = [M[i] + [rhs[i]] for i in range(len(M))]
    r = 0; piv_cols = []
    for c in range(NU):
        piv = next((i for i in range(r, len(aug)) if aug[i][c] != 0), None)
        if piv is None: continue
        aug[r], aug[piv] = aug[piv], aug[r]
        pv = aug[r][c]
        aug[r] = [x/pv for x in aug[r]]
        for i in range(len(aug)):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [a - f*b for a, b in zip(aug[i], aug[r])]
        r += 1; piv_cols.append(c)
    for i in range(r, len(aug)):
        if aug[i][NU] != 0:
            return None
    part = [Fr(0)]*NU
    for i, c in enumerate(piv_cols): part[c] = aug[i][NU]
    free = [c for c in range(NU) if c not in piv_cols]
    basis = []
    for fc in free:
        v = [Fr(0)]*NU; v[fc] = Fr(1)
        for i, c in enumerate(piv_cols): v[c] = -aug[i][fc]
        basis.append(v)
    return part, basis, pair, (n1, n2)

def show(sol, tag):
    if sol is None:
        print(f"  {tag}: INCONSISTENT"); return
    part, basis, pair, (n1, n2) = sol
    p1, p2 = part[:n1], part[n1:]
    print(f"  {tag}: dim {len(basis)}; particular {pair[0]} = {[str(x) for x in p1]}, "
          f"{pair[1]} = {[str(x) for x in p2]}")
    for i, b in enumerate(basis):
        print(f"     dir{i}: {pair[0]} += {[str(x) for x in b[:n1]]}, "
              f"{pair[1]} += {[str(x) for x in b[n1:]]}")

print("=== k=2 per-row solution spaces at the known witness ===")
known_all = {"A": [1,3,3,1], "B": [4,7,3], "C": [1,12,9], "D": [3,6,3],
             "E0": [2,-3], "E1": [1]}
show(solve_row(2, 3, {n: known_all[n] for n in ("A","B","C","D")}, (2,1)), "row3 (E0,E1) given rows 1,2")
show(solve_row(2, 2, {n: known_all[n] for n in ("A","B","E0","E1")}, (3,3)), "row2 (C,D) given rows 1,3")
show(solve_row(2, 1, {n: known_all[n] for n in ("C","D","E0","E1")}, (4,3)), "row1 (A,B) given rows 2,3")

print("\n=== k=3 block iteration from the k-pattern seed ===")
# pattern guesses (u = 1+t):  A = u^4, B = 4u^3 + u, C = 16u^2-8u-3, D = 4u^3,
#                             E0 = 3-4t, E1 = 1
seed = {"A": [1,4,6,4,1], "B": [5,13,12,4], "C": [5,24,16], "D": [4,12,12,4],
        "E0": [3,-4], "E1": [1]}
cur = {k: list(v) for k, v in seed.items()}
def to_list(part, basis, n1, n2, prefer):
    """pick the affine-space element closest to prefer (project onto particular
    + span); here spaces are small: if dim==1 scale the direction to match
    prefer's leading nonzero; if dim==0 take particular."""
    from fractions import Fraction as F2
    if not basis:
        v = part
    else:
        b = basis[0]
        idx = next((i for i, x in enumerate(b) if x != 0), None)
        scale = F2(prefer[idx]) - part[idx] if idx is not None and idx < len(prefer) else F2(1)
        if idx is not None and b[idx] != 0:
            scale = scale / b[idx]
        v = [p + scale*bb for p, bb in zip(part, b)]
    return v[:n1], v[n1:]

for it in range(6):
    changed = False
    # row 2 given rows 1,3
    r = solve_row(3, 2, {n: cur[n] for n in ("A","B","E0","E1")}, (3,4))
    if r is None:
        print(f"  iter {it}: row2 INCONSISTENT"); break
    part, basis, pair, (n1, n2) = r
    pref = cur["C"] + [0]*(n1-len(cur["C"])) + cur["D"] + [0]*(n2-len(cur["D"]))
    c_new, d_new = to_list(part, basis, n1, n2, pref)
    print(f"  iter {it}: row2 dim {len(basis)}; C = {[str(x) for x in c_new]}, D = {[str(x) for x in d_new]}")
    # row 1 given rows 2,3
    cur["C"], cur["D"] = c_new, d_new
    r = solve_row(3, 1, {"C": [str(x) for x in []] or cur["C"], "D": cur["D"], "E0": cur["E0"], "E1": cur["E1"]}, (5,4))
    if r is None:
        print(f"  iter {it}: row1 INCONSISTENT"); break
    part, basis, pair, (n1, n2) = r
    pref = cur["A"] + [0]*(n1-len(cur["A"])) + cur["B"] + [0]*(n2-len(cur["B"]))
    a_new, b_new = to_list(part, basis, n1, n2, pref)
    print(f"           row1 dim {len(basis)}; A = {[str(x) for x in a_new]}, B = {[str(x) for x in b_new]}")
    cur["A"], cur["B"] = a_new, b_new
    # row 3 given rows 1,2
    r = solve_row(3, 3, {n: cur[n] for n in ("A","B","C","D")}, (2,1))
    if r is None:
        print(f"           row3 INCONSISTENT"); break
    part, basis, pair, (n1, n2) = r
    pref = cur["E0"] + [0]*(n1-len(cur["E0"])) + cur["E1"] + [0]*(n2-len(cur["E1"]))
    e0_new, e1_new = to_list(part, basis, n1, n2, pref)
    print(f"           row3 dim {len(basis)}; E0 = {[str(x) for x in e0_new]}, E1 = {[str(x) for x in e1_new]}")
    if e0_new == cur["E0"] + [0]*(n1-len(cur["E0"])) and cur["A"] == a_new:
        cur["E0"], cur["E1"] = e0_new, e1_new
        print("  FIXED POINT reached")
        break
    cur["E0"], cur["E1"] = e0_new, e1_new
