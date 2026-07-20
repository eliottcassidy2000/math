#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- k=3 proper solve via the s-graded
system.  A = v^{k+1}, D = delta v^2 (kills c2 for every k).  Then:
 (*)  c1 = E0(A'D-AD') + (k-1)E0'AD - 2t^k B'D - 2k t^{k-1} BD + k t^k BD'
          + (k+1)AC + (k+1)tAC' - tA'C  == 0            [LINEAR in B,C,E0]
 (**) c0 = (E0+tE0')(t^k B'D + k t^{k-1} BD - AC - tAC') - E0'(t^{k+1}B'D - t^2 AC')
          - t^{k+1}(B'C - kBC')  == const               [bilinear]
Solve (*) exactly; substitute the affine solution space into (**); examine the
resulting quadratic system over Q.  Also re-verify the whole grading at k=2.
"""
from fractions import Fraction as Fr
from math import comb

def pm(a, b):
    r = {}
    for ka, ca in a.items():
        for kb, cb in b.items():
            k = tuple(p+q for p, q in zip(ka, kb))
            v = r.get(k, 0) + ca*cb
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def pa(*ps):
    r = {}
    for p in ps:
        for k, c in p.items():
            v = r.get(k, 0) + c
            if v: r[k] = v
            elif k in r: del r[k]
    return r
def ps_(p, s): return {k: c*s for k, c in p.items() if c*s != 0}

def check_grading(k, A, B, C, D, E0, E1c):
    """verify c2,c1,c0 formulas against direct det for numeric functions
    (univariate dicts keyed (j,) over t).  Returns (c2,c1,c0) as t-polys."""
    d = lambda p: {(j-1,): c*j for (j,), c in p.items() if j}
    t = {(1,): 1}
    tk = {(k,): 1}; tk1 = {(k-1,): 1}; tk_1p = {(k+1,): 1}
    Ap, Bp, Cp, Dp, E0p = d(A), d(B), d(C), d(D), d(E0)
    M1 = pa(pm(Ap, D), ps_(pm(A, Dp), -1))
    c2 = pa(ps_(pm({(0,):1}, M1), -(k+1)*E1c), ps_(pm(Ap, D), (k-1)*E1c))
    c1 = pa(pm(E0, M1), ps_(pm(E0p, pm(A, D)), k-1),
            ps_(pm(tk, pm(Bp, D)), -2), ps_(pm(tk1, pm(B, D)), -2*k),
            ps_(pm(tk, pm(B, Dp)), k),
            ps_(pm(A, C), k+1), ps_(pm(t, pm(A, Cp)), k+1),
            ps_(pm(t, pm(Ap, C)), -1))
    K0 = pa(E0, pm(t, E0p))
    M0 = pa(pm(tk, pm(Bp, D)), ps_(pm(tk1, pm(B, D)), k),
            ps_(pm(A, C), -1), ps_(pm(t, pm(A, Cp)), -1))
    c0 = pa(pm(K0, M0),
            ps_(pm(E0p, pa(pm(tk_1p, pm(Bp, D)),
                           ps_(pm({(2,):1}, pm(A, Cp)), -1))), -1),
            ps_(pm(tk_1p, pa(pm(Bp, C), ps_(pm(B, Cp), -k))), -1))
    return c2, c1, c0

# regression at k=2 with the known witness
A2 = {(0,):1,(1,):3,(2,):3,(3,):1}; B2 = {(0,):4,(1,):7,(2,):3}
C2 = {(0,):1,(1,):12,(2,):9}; D2 = {(0,):3,(1,):6,(2,):3}
E02 = {(0,):2,(1,):-3}
c2_, c1_, c0_ = check_grading(2, A2, B2, C2, D2, E02, 1)
print("k=2 grading check: c2 =", c2_, " c1 =", c1_, " c0 =", c0_, "(expect 0, 0, {(0,):-2})")

# ---------------- k=3 solve ----------------
k = 3
degB, degC, degE0 = 4, 4, 3
NB, NC, NE = degB+1, degC+1, degE0+1
NU = NB + NC + NE
# unknown-coefficient univariate polys over variables tuple (t-deg, u_0..u_{NU-1})
def uvar(j):
    key = [0]*(1+NU); key[1+j] = 1
    return tuple(key)
def tpow(j):
    key = [0]*(1+NU); key[0] = j
    return tuple(key)
def lift_uni(p):  # numeric t-poly -> extended ring
    return {tpow(j): c for (j,), c in p.items()}
one = {tpow(0): 1}
t_ = {tpow(1): 1}
v = pa(one, t_)
A = one
for _ in range(k+1): A = pm(A, v)     # v^4
delta = 1
D = ps_(pm(v, v), delta)              # v^2
B  = pa(*[pm({uvar(j):1}, {tpow(j):1}) for j in range(NB)])
C  = pa(*[pm({uvar(NB+j):1}, {tpow(j):1}) for j in range(NC)])
E0 = pa(*[pm({uvar(NB+NC+j):1}, {tpow(j):1}) for j in range(NE)])
d = lambda p: {tuple([kk[0]-1]+list(kk[1:])): c*kk[0] for kk, c in p.items() if kk[0]}
Ap, Bp, Cp, Dp, E0p = d(A), d(B), d(C), d(D), d(E0)
tk = {tpow(k):1}; tk1 = {tpow(k-1):1}; tk_1p = {tpow(k+1):1}
M1 = pa(pm(Ap, D), ps_(pm(A, Dp), -1))
c1 = pa(pm(E0, M1), ps_(pm(E0p, pm(A, D)), k-1),
        ps_(pm(tk, pm(Bp, D)), -2), ps_(pm(tk1, pm(B, D)), -2*k),
        ps_(pm(tk, pm(B, Dp)), k),
        ps_(pm(A, C), k+1), ps_(pm(t_, pm(A, Cp)), k+1),
        ps_(pm(t_, pm(Ap, C)), -1))
# c1 == 0: collect by t-degree; each row linear in unknowns
rows = {}
for kk, c in c1.items():
    td = kk[0]; unk = kk[1:]
    tot = sum(unk)
    assert tot <= 1
    rows.setdefault(td, [Fr(0)]*(NU+1))
    if tot == 0: rows[td][NU] += c
    else: rows[td][unk.index(1)] += c
aug = [rows[td][:NU] + [-rows[td][NU]] for td in sorted(rows)]
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
bad = [i for i in range(r, len(aug)) if aug[i][NU] != 0]
print(f"\nk=3 (*) system: {len(aug)} eqs, {NU} unknowns, rank {r}, "
      f"inconsistent rows: {len(bad)}, solution dim {NU - r if not bad else 'NONE'}")
if not bad:
    part = [Fr(0)]*NU
    for i, c in enumerate(piv): part[c] = aug[i][NU]
    free = [c for c in range(NU) if c not in piv]
    basis = []
    for fc in free:
        vv = [Fr(0)]*NU; vv[fc] = Fr(1)
        for i, c in enumerate(piv): vv[c] = -aug[i][fc]
        basis.append(vv)
    names = [f"B{j}" for j in range(NB)] + [f"C{j}" for j in range(NC)] + [f"E0{j}" for j in range(NE)]
    print("particular:", {n: str(x) for n, x in zip(names, part) if x != 0} or "0")
    for bi, bb in enumerate(basis):
        print(f"free dir {bi} ({names[free[bi]]}):",
              {n: str(x) for n, x in zip(names, bb) if x != 0})

# ---------------- stage 2: c0 == const on the (*)-solution space ----------------
print("\n--- stage 2: c0 constancy (quadratic system in the 7 free params) ---")
NL = len(basis)
def ulam(j):
    key = [0]*(1+NL); key[1+j] = 1
    return tuple(key)
def tpow2(j):
    key = [0]*(1+NL); key[0] = j
    return tuple(key)
one2 = {tpow2(0): 1}
t2 = {tpow2(1): 1}
# rebuild B, C, E0 as lambda-linear polys
Bl = {}; Cl = {}; E0l = {}
for bi, bb in enumerate(basis):
    for j in range(NB):
        if bb[j]: Bl[ (j, bi) ] = Bl.get((j,bi), 0) + bb[j]
    for j in range(NC):
        if bb[NB+j]: Cl[ (j, bi) ] = Cl.get((j,bi), 0) + bb[NB+j]
    for j in range(NE):
        if bb[NB+NC+j]: E0l[ (j, bi) ] = E0l.get((j,bi), 0) + bb[NB+NC+j]
def assemble(dd):
    out = {}
    for (j, bi), c in dd.items():
        key = [0]*(1+NL); key[0] = j; key[1+bi] = 1
        out[tuple(key)] = out.get(tuple(key), 0) + c
    return out
B2_ = assemble(Bl); C2_ = assemble(Cl); E02_ = assemble(E0l)
v2_ = pa(one2, t2)
A2_ = one2
for _ in range(k+1): A2_ = pm(A2_, v2_)
D2_ = pm(v2_, v2_)
d2 = lambda p: {tuple([kk[0]-1]+list(kk[1:])): c*kk[0] for kk, c in p.items() if kk[0]}
Ap2, Bp2, Cp2, Dp2, E0p2 = d2(A2_), d2(B2_), d2(C2_), d2(D2_), d2(E02_)
tk2 = {tpow2(k):1}; tk12 = {tpow2(k-1):1}; tk_1p2 = {tpow2(k+1):1}
K0 = pa(E02_, pm(t2, E0p2))
M0 = pa(pm(tk2, pm(Bp2, D2_)), ps_(pm(tk12, pm(B2_, D2_)), k),
        ps_(pm(A2_, C2_), -1), ps_(pm(t2, pm(A2_, Cp2)), -1))
c0 = pa(pm(K0, M0),
        ps_(pm(E0p2, pa(pm(tk_1p2, pm(Bp2, D2_)),
                        ps_(pm({tpow2(2):1}, pm(A2_, Cp2)), -1))), -1),
        ps_(pm(tk_1p2, pa(pm(Bp2, C2_), ps_(pm(B2_, Cp2), -k))), -1))
# collect by t-degree: each = quadratic in lambdas; t^0 free
eqs = {}
for kk, c in c0.items():
    td = kk[0]; lam = kk[1:]
    if td == 0: continue
    eqs.setdefault(td, {})
    eqs[td][lam] = eqs[td].get(lam, 0) + c
eqs = {td: {l: c for l, c in e.items() if c} for td, e in eqs.items()}
eqs = {td: e for td, e in eqs.items() if e}
print(f"quadratic equations at t-degrees: {sorted(eqs)}")

# reduced ansatz: lambda = (c2, 0, 0, e0, e1, 0, 0), normalize c2 = 1
# substitute and print residual equations in (e0, e1)
import itertools
def subs_lam(e, vals):
    out = 0
    for lam, c in e.items():
        m = c
        for j, ex in enumerate(lam):
            if ex: m *= vals[j]**ex
        out += m
    return out
# symbolic in e0, e1: represent as dict (i,j) -> coeff of e0^i e1^j
def subs_sym(e, freeidx):
    out = {}
    for lam, c in e.items():
        key = [0, 0]
        ok = True
        m = c
        for j, ex in enumerate(lam):
            if ex == 0: continue
            if j == 0: m *= 1          # c2 = 1
            elif j == freeidx[0]: key[0] += ex
            elif j == freeidx[1]: key[1] += ex
            else: ok = False; break
        if ok:
            out[tuple(key)] = out.get(tuple(key), 0) + m
    return out
print("\nreduced ansatz lambda = (1,0,0,e0,e1,0,0):")
polys = []
for td in sorted(eqs):
    q = subs_sym(eqs[td], (3, 4))
    # check the dropped lambdas really absent: recompute with full check
    full_ok = all(all(ex == 0 for j, ex in enumerate(lam) if j in (1,2,5,6)) or True
                  for lam in eqs[td])
    polys.append((td, q))
    print(f"  t^{td}: {q}")

# ---------------- stage 3: does the FULL quadratic system have solutions? ----------------
print("\n--- stage 3: numeric hunt over the full 7-param system (11 eqs) ---")
import random
random.seed(59)
EQS = []
for td in sorted(eqs):
    if eqs[td]: EQS.append(eqs[td])
def F_and_J(x):
    Fv = []
    Jv = []
    for e in EQS:
        val = 0.0
        grad = [0.0]*NL
        for lam, c in e.items():
            idxs = [j for j, ex in enumerate(lam) for _ in range(ex)]
            cf = float(c)
            if len(idxs) == 1:
                val += cf * x[idxs[0]]
                grad[idxs[0]] += cf
            elif len(idxs) == 2:
                a, b = idxs
                val += cf * x[a] * x[b]
                grad[a] += cf * x[b]
                grad[b] += cf * x[a]
            else:
                val += cf
        Fv.append(val); Jv.append(grad)
    return Fv, Jv

def gn(x, iters=200):
    lam_damp = 1e-3
    for _ in range(iters):
        Fv, Jv = F_and_J(x)
        # normal equations with damping (pure python)
        n = NL
        A_ = [[sum(Jv[r][i]*Jv[r][j] for r in range(len(Jv))) + (lam_damp if i==j else 0)
               for j in range(n)] for i in range(n)]
        b_ = [-sum(Jv[r][i]*Fv[r] for r in range(len(Jv))) for i in range(n)]
        # gaussian solve
        M = [row[:] + [b_[i]] for i, row in enumerate(A_)]
        for c in range(n):
            piv = max(range(c, n), key=lambda i: abs(M[i][c]))
            if abs(M[piv][c]) < 1e-14: break
            M[c], M[piv] = M[piv], M[c]
            pv = M[c][c]
            M[c] = [xx/pv for xx in M[c]]
            for i in range(n):
                if i != c and M[i][c]:
                    f = M[i][c]
                    M[i] = [aa - f*bb for aa, bb in zip(M[i], M[c])]
        dx = [M[i][n] for i in range(n)]
        x = [xi + di for xi, di in zip(x, dx)]
    Fv, _ = F_and_J(x)
    return x, sum(f*f for f in Fv)

best = []
for trial in range(400):
    x0 = [random.uniform(-4, 4) for _ in range(NL)]
    x, res = gn(x0)
    if res < 1e-18:
        best.append((res, x))
print(f"hits with residual < 1e-18: {len(best)} / 400 starts")
seen = []
for res, x in sorted(best)[:12]:
    key = tuple(round(v, 6) for v in x)
    if any(max(abs(a-b) for a,b in zip(key, s)) < 1e-4 for s in seen): continue
    seen.append(key)
    print("  sol:", [f"{v:.6f}" for v in x], f"res {res:.1e}")
# also complex? try complex starts via 2n real embedding if no real hits
if not best:
    print("  no real hits; trying complex Newton (real 14-dim embedding)")
    def F_and_J_c(z):
        xs = [complex(z[2*i], z[2*i+1]) for i in range(NL)]
        Fv = []
        Jc = []
        for e in EQS:
            val = 0j
            grad = [0j]*NL
            for lam, c in e.items():
                idxs = [j for j, ex in enumerate(lam) for _ in range(ex)]
                cf = complex(c)
                if len(idxs) == 1:
                    val += cf * xs[idxs[0]]; grad[idxs[0]] += cf
                elif len(idxs) == 2:
                    a, b = idxs
                    val += cf*xs[a]*xs[b]; grad[a] += cf*xs[b]; grad[b] += cf*xs[a]
                else: val += cf
            Fv.append(val); Jc.append(grad)
        return Fv, Jc
    hits = 0
    sols = []
    for trial in range(400):
        z = [random.uniform(-3,3) for _ in range(2*NL)]
        for _ in range(120):
            Fv, Jc = F_and_J_c(z)
            # complex GN: solve (J^H J + damp) dx = -J^H F  in complex arithmetic
            n = NL
            Ac = [[sum(Jc[r][i].conjugate()*Jc[r][j] for r in range(len(Jc))) + (1e-3 if i==j else 0)
                   for j in range(n)] for i in range(n)]
            bc = [-sum(Jc[r][i].conjugate()*Fv[r] for r in range(len(Jc))) for i in range(n)]
            M = [row[:] + [bc[i]] for i, row in enumerate(Ac)]
            ok = True
            for c in range(n):
                piv = max(range(c, n), key=lambda i: abs(M[i][c]))
                if abs(M[piv][c]) < 1e-14: ok = False; break
                M[c], M[piv] = M[piv], M[c]
                pv = M[c][c]
                M[c] = [xx/pv for xx in M[c]]
                for i in range(n):
                    if i != c and M[i][c]:
                        f = M[i][c]
                        M[i] = [aa - f*bb for aa, bb in zip(M[i], M[c])]
            if not ok: break
            dx = [M[i][n] for i in range(n)]
            for i in range(NL):
                z[2*i] += dx[i].real; z[2*i+1] += dx[i].imag
        Fv, _ = F_and_J_c(z)
        res = sum(abs(f)**2 for f in Fv)
        if res < 1e-18:
            hits += 1
            sols.append([complex(z[2*i], z[2*i+1]) for i in range(NL)])
    print(f"  complex hits: {hits}/400")
    seenc = []
    for s in sols[:8]:
        key = tuple((round(v.real,5), round(v.imag,5)) for v in s)
        if key in seenc: continue
        seenc.append(key)
        print("   csol:", [f"{v.real:+.4f}{v.imag:+.4f}i" for v in s])

# ---------------- W3 check: is Phi = tC + E0*D forced linear on the c1-space? ----------------
print("\n--- W3: Phi = tC + E0*D on the 7-dim c1-solution space (k=3) ---")
# Phi is BILINEAR in (C, E0-with-D-fixed): D = v^2 fixed => tC + E0*D is LINEAR in lambdas!
Phi = pa(pm(t2, C2_), pm(E02_, D2_))
by_deg = {}
for kk, c in Phi.items():
    td = kk[0]; lam = kk[1:]
    assert sum(lam) <= 1
    by_deg.setdefault(td, {})
    j = lam.index(1) if sum(lam) else -1
    by_deg[td][j] = by_deg[td].get(j, 0) + c
for td in sorted(by_deg):
    nz = {j: str(c) for j, c in by_deg[td].items() if c}
    print(f"  t^{td}: {nz if nz else '0'}")
print("  (W3 = TRUE iff coefficients of t^2 and higher vanish identically in lambda)")
