#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- k=3 discovery via subset-linear solves.

General-k s-grading of det (derived by hand, verified by engine at k=2):
  c2 = -2A'D + (k+1)AD'  == 0   =>  A = v^{k+1}, D = delta v^2.
Iteration: A = (1+t)^{k+1} FIXED, E1 = 1 FIXED; alternate exact linear solves
of {B}, {C,D}, {E0} until fixed point; then verify det == const and build the
collision certificate (fixed branch + orbit branch at the root of
Phi_k = tC + E0 D).
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

def build_and_det(k, funcs, NV):
    def var(i):
        kk = [0]*NV; kk[i] = 1
        return {tuple(kk): 1}
    X, Y, Z = var(0), var(1), var(2)
    A,B,C,D,E0,E1 = (funcs[n] for n in ("A","B","C","D","E0","E1"))
    Yk = {tuple([0]*NV):1}
    for _ in range(k): Yk = pmul(Yk, Y)
    Xk_1 = {tuple([0]*NV):1}
    for _ in range(k-1): Xk_1 = pmul(Xk_1, X)
    sk = pmul(pmul(Xk_1, X), Z)
    F1 = padd(pmul(Z, A), pmul(Yk, B))
    F2 = padd(pmul(Y, C), pmul(pmul(Xk_1, Z), D))
    F3 = pmul(X, padd(E0, psc(pmul(sk, E1), -1)))
    J = [[pdiff(F, i) for i in range(3)] for F in (F1, F2, F3)]
    det = padd(
        pmul(J[0][0], padd(pmul(J[1][1], J[2][2]), psc(pmul(J[1][2], J[2][1]), -1))),
        psc(pmul(J[0][1], padd(pmul(J[1][0], J[2][2]), psc(pmul(J[1][2], J[2][0]), -1))), -1),
        pmul(J[0][2], padd(pmul(J[1][0], J[2][1]), psc(pmul(J[1][1], J[2][0]), -1))))
    return (F1, F2, F3), det

def lift_known(cl, NV):
    out = {}
    for j, c in enumerate(cl):
        if c:
            key = [0]*NV; key[0] = j; key[1] = j
            out[tuple(key)] = out.get(tuple(key), 0) + c
    return out

def solve_subset(k, free, fixed, degs):
    """free: list of function names (all in one det-row); fixed: name->list.
       degs: name->max t-degree. Exact affine solve of det==const."""
    NU = sum(degs[n]+1 for n in free)
    NV = 3 + NU
    def var(i):
        kk = [0]*NV; kk[i] = 1
        return {tuple(kk): 1}
    t = pmul(var(0), var(1))
    funcs = {}
    base = 3
    for n in ("A","B","C","D","E0","E1"):
        if n in free:
            out = {}; tp = {tuple([0]*NV):1}
            for j in range(degs[n]+1):
                out = padd(out, pmul(var(base+j), tp))
                tp = pmul(tp, t)
            funcs[n] = out; base += degs[n]+1
        else:
            funcs[n] = lift_known(fixed[n], NV)
    _, det = build_and_det(k, funcs, NV)
    rows = {}
    for key, c in det.items():
        xyz = key[:3]; unk = key[3:]
        tot = sum(unk)
        assert tot <= 1, f"nonlinear in frees {free}: {key}"
        rows.setdefault(xyz, [0]*(NU+1))
        if tot == 0: rows[xyz][NU] += c
        else: rows[xyz][unk.index(1)] += c
    aug = []
    for xyz, rr in rows.items():
        if xyz == (0,0,0): continue
        aug.append([Fr(c) for c in rr[:NU]] + [Fr(-rr[NU])])
    r = 0; piv = []
    for c in range(NU):
        pv = next((i for i in range(r, len(aug)) if aug[i][c] != 0), None)
        if pv is None: continue
        aug[r], aug[pv] = aug[pv], aug[r]
        p0 = aug[r][c]
        aug[r] = [x/p0 for x in aug[r]]
        for i in range(len(aug)):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [a - f*b for a, b in zip(aug[i], aug[r])]
        r += 1; piv.append(c)
    for i in range(r, len(aug)):
        if aug[i][NU] != 0: return None
    part = [Fr(0)]*NU
    for i, c in enumerate(piv): part[c] = aug[i][NU]
    basis = []
    for fc in [c for c in range(NU) if c not in piv]:
        v = [Fr(0)]*NU; v[fc] = Fr(1)
        for i, c in enumerate(piv): v[c] = -aug[i][fc]
        basis.append(v)
    return part, basis

def run_k(k, delta, iters=12):
    # A = (1+t)^{k+1}
    from math import comb
    A = [comb(k+1, j) for j in range(k+2)]
    state = {"A": A, "E1": [1],
             "B": [Fr(delta) + 1] + [0]*3,       # loose seed
             "C": [1] + [0]*3,
             "D": [Fr(delta), 2*Fr(delta), Fr(delta)],   # delta*(1+t)^2
             "E0": [k+3, -(k+1)] + [0]*1}
    degs = {"B": 3, "C": 3, "D": 2, "E0": 2}
    hist = []
    for it in range(iters):
        snap = {n: list(state[n]) for n in ("B","C","D","E0")}
        # {B}
        r = solve_subset(k, ["B"], {n: state[n] for n in ("A","C","D","E0","E1")}, degs)
        if r is None: return None, f"B inconsistent at iter {it}"
        part, basis = r
        # choose element closest to current B (project): particular + span; take
        # particular unless zero and basis exists w/ seed match
        Bv = part
        if basis:
            # align to previous B leading coeff if possible
            b0 = basis[0]
            idx = next((i for i, x in enumerate(b0) if x != 0), None)
            if idx is not None:
                sc = (Fr(state["B"][idx]) - part[idx]) / b0[idx]
                Bv = [p + sc*x for p, x in zip(part, b0)]
        state["B"] = Bv
        # {C,D}
        r = solve_subset(k, ["C","D"], {n: state[n] for n in ("A","B","E0","E1")}, degs)
        if r is None: return None, f"CD inconsistent at iter {it}"
        part, basis = r
        v = part
        if basis:
            b0 = basis[0]
            idx = next((i for i, x in enumerate(b0) if x != 0), None)
            if idx is not None:
                prev = state["C"] + state["D"]
                sc = (Fr(prev[idx]) - part[idx]) / b0[idx]
                v = [p + sc*x for p, x in zip(part, b0)]
        nC = degs["C"]+1
        state["C"], state["D"] = v[:nC], v[nC:]
        # {E0}
        r = solve_subset(k, ["E0"], {n: state[n] for n in ("A","B","C","D","E1")}, degs)
        if r is None: return None, f"E0 inconsistent at iter {it}"
        part, basis = r
        v = part
        if basis:
            b0 = basis[0]
            idx = next((i for i, x in enumerate(b0) if x != 0), None)
            if idx is not None:
                sc = (Fr(state["E0"][idx]) - part[idx]) / b0[idx]
                v = [p + sc*x for p, x in zip(part, b0)]
        state["E0"] = v
        if all(state[n] == snap[n] for n in snap):
            return state, f"fixed point at iter {it}"
        hist.append({n: list(state[n]) for n in snap})
    return state, "no fixed point (last state returned)"

for delta in [1, 2, 3, 4, 6]:
    st, msg = run_k(3, delta)
    print(f"delta={delta}: {msg}")
    if st:
        # verify det
        NV = 3
        funcs = {n: lift_known([Fr(x) for x in st[n]], NV) for n in st}
        _, det = build_and_det(3, funcs, NV)
        nonzero = {k: v for k, v in det.items() if k != (0,0,0)}
        c = det.get((0,0,0), 0)
        print(f"   det == const: {not nonzero} (const = {c})")
        if not nonzero and c != 0:
            print("   SOLUTION:", {n: [str(x) for x in st[n]] for n in st})
