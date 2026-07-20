#!/usr/bin/env python3
"""
jacobian_lifting_ladder_boxeph_S142.py  (HYP-8105)

THE LIFTING PROGRAM EXECUTED: hunt for new Keller kernels at other defect classes.

FRAME (from S141): every C*-equivariant z-linear map with weights (-1,1,2)->(2,1,-1)
is  F = ( z A(s) + y^2 B(s),  y g1(s) + xz g2(s),  x (h0(s) - r) ),
s = xy, r = x^2 z, five univariate polynomials.  The known kernel decodes as
A = u^3 (u = 1+s), g2 = A', B = u + A', g1 = 1+12s+9s^2, h0 = 2-3s, and its
collision resolvent  Phi(s) = s g1 + h0 g2 = 2(2s+3)  has the single root -3/2
(=> axis fiber 1 + 2 = 3, matching the S3 triple cover).

THE HUNT: for each d, set A = u^d, g2 = A', B = u + A' (the decoded generative
identities), leave g1 (deg d-1) and h0 (deg d-2) UNKNOWN, impose det JF == const
SYMBOLICALLY (unknown coefficients as extra formal variables), solve the system
by triangular elimination, and for every solution:
  * verify det == const exactly;
  * compute Phi_d and its rational roots -> explicit collisions, verified;
  * mod-p census (p = 11, 13): fiber histogram (degree/monodromy fingerprint)
    and non-bijectivity certificate (Keller + injective => automorphism =>
    mod-p bijective; so non-bijectivity certifies a genuine counterexample).

d = 3 must reproduce the known kernel (control).  boxeph-2026-07-19-S142.
"""

from fractions import Fraction as Fr
from itertools import product

# ---------- n-variable exact polynomial arithmetic (tuple-keyed) ----------------
def mul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            k = tuple(i + j for i, j in zip(e, f))
            out[k] = out.get(k, 0) + c * d
    return {k: v for k, v in out.items() if v}

def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}

def sc(p, c):
    return {k: v * c for k, v in p.items() if v * c}

def diff(p, var):
    out = {}
    for e, c in p.items():
        if e[var]:
            k = list(e); k[var] -= 1
            out[tuple(k)] = out.get(tuple(k), 0) + c * e[var]
    return out

def NV_var(nv, i):
    e = [0] * nv; e[i] = 1
    return {tuple(e): Fr(1)}

def det3(M):
    return add(mul(M[0][0], add(mul(M[1][1], M[2][2]), sc(mul(M[1][2], M[2][1]), -1))),
               sc(mul(M[0][1], add(mul(M[1][0], M[2][2]), sc(mul(M[1][2], M[2][0]), -1))), -1),
               mul(M[0][2], add(mul(M[1][0], M[2][1]), sc(mul(M[1][1], M[2][0]), -1))))

def hunt(d, verbose=True):
    ng1, nh0 = d, d - 1                      # #coeffs: g1 deg d-1, h0 deg d-2
    nv = 3 + ng1 + nh0                       # x,y,z + gammas + etas
    X, Y, Z = NV_var(nv, 0), NV_var(nv, 1), NV_var(nv, 2)
    ONE = {tuple([0] * nv): Fr(1)}
    S = mul(X, Y)
    U = add(ONE, S)
    def upow(k):
        r = ONE
        for _ in range(k): r = mul(r, U)
        return r
    A = upow(d)
    A2 = sc(upow(d - 1), d)                  # g2 = A' (d/ds)
    B = add(U, A2)                           # B = u + A'
    g1 = add(*[mul(NV_var(nv, 3 + i), mulpow(S, i, nv)) for i in range(ng1)]) if False else None
    def spow(k):
        r = ONE
        for _ in range(k): r = mul(r, S)
        return r
    g1 = add(*[mul(NV_var(nv, 3 + i), spow(i)) for i in range(ng1)])
    h0 = add(*[mul(NV_var(nv, 3 + ng1 + i), spow(i)) for i in range(nh0)])
    f1 = add(mul(Z, A), mul(mul(Y, Y), B))
    f2 = add(mul(Y, g1), mul(mul(X, Z), A2))
    f3 = add(mul(X, h0), sc(mul(mul(mul(X, X), X), Z), -1))
    F = [f1, f2, f3]
    J = [[diff(f, v) for v in range(3)] for f in F]
    D = det3(J)
    # constraint system: for each xyz-monomial != 0, coefficient-poly in unknowns = 0
    sysm = {}
    const_terms = {}
    for e, c in D.items():
        xyz, rest = e[:3], e[3:]
        if xyz == (0, 0, 0):
            const_terms[rest] = const_terms.get(rest, 0) + c
        else:
            sysm.setdefault(xyz, {})[rest] = sysm.get(xyz, {}).get(rest, 0) + c
    eqs = [v for v in sysm.values()]
    # triangular solve: repeatedly find an equation linear in exactly one unfixed unknown
    vals = {}
    nunk = ng1 + nh0
    def reduce_eq(eq):
        out = {}
        for rest, c in eq.items():
            coef = c
            newrest = list(rest)
            for i, ex in enumerate(rest):
                if i in vals and ex:
                    coef *= vals[i] ** ex
                    newrest[i] = 0
            k = tuple(newrest)
            out[k] = out.get(k, 0) + coef
        return {k: v for k, v in out.items() if v}
    progress = True
    inconsistent = False
    while progress and not inconsistent:
        progress = False
        for eq in eqs:
            r = reduce_eq(eq)
            if not r: continue
            unk = sorted({i for k in r for i, ex in enumerate(k) if ex})
            if not unk:
                inconsistent = True; break
            if len(unk) == 1:
                i = unk[0]
                lin = all(k[i] <= 1 for k in r)
                if lin:
                    a = sum(v for k, v in r.items() if k[i] == 1)
                    b = sum(v for k, v in r.items() if k[i] == 0)
                    if a != 0:
                        vals[i] = Fr(-b, a) if not isinstance(b, Fr) else -b / a
                        progress = True
    solved = len(vals) == nunk and not inconsistent
    # check all equations vanish
    ok = solved and all(not reduce_eq(eq) for eq in eqs)
    return solved, ok, vals, ng1, nh0, inconsistent

def mulpow(S, i, nv):
    pass

results = {}
print("=" * 96)
print("THE LADDER HUNT: A = u^d, g2 = A', B = u + A'; solve det==const for (g1, h0)")
print("=" * 96)
for d in range(2, 8):
    solved, ok, vals, ng1, nh0, inc = hunt(d)
    if ok:
        g1c = [vals.get(i, None) for i in range(ng1)]
        h0c = [vals.get(ng1 + i, None) for i in range(nh0)]
        print("d=%d: SOLVED  g1 = %s  h0 = %s" % (d, g1c, h0c))
        results[d] = (g1c, h0c)
    else:
        print("d=%d: %s (fixed %d/%d unknowns)" % (d, "INCONSISTENT" if inc else "not fully solved by triangular pass", len(vals), ng1 + nh0))
print()

# ---------------- assemble + certify each solved rung ---------------------------
def build_F(d, g1c, h0c):
    def P3(mon): return {mon: Fr(1)}
    X, Y, Z, ONE = P3((1,0,0)), P3((0,1,0)), P3((0,0,1)), P3((0,0,0))
    S = mul(X, Y); U = add(ONE, S)
    def pw(p, k):
        r = ONE
        for _ in range(k): r = mul(r, p)
        return r
    A = pw(U, d); A2 = sc(pw(U, d - 1), d); B = add(U, A2)
    g1 = add(*[sc(pw(S, i), c) for i, c in enumerate(g1c)])
    h0 = add(*[sc(pw(S, i), c) for i, c in enumerate(h0c)])
    f1 = add(mul(Z, A), mul(mul(Y, Y), B))
    f2 = add(mul(Y, g1), mul(mul(X, Z), A2))
    f3 = add(mul(X, h0), sc(mul(pw(X, 3), Z), -1))
    return [f1, f2, f3], g1, h0, A2

def pev(p, x, y, z):
    return sum(c * x**e[0] * y**e[1] * z**e[2] for e, c in p.items())

for d, (g1c, h0c) in sorted(results.items()):
    Fd, g1, h0, g2 = build_F(d, g1c, h0c)
    J = [[diff(f, v) for v in range(3)] for f in Fd]
    D = det3(J)
    const = list(D.items())
    assert len(const) == 1 and const[0][0] == (0, 0, 0), ("det not const?!", d)
    detval = const[0][1]
    # Phi(s) = s g1 + h0 g2 in one variable
    Phi = {}
    for e, c in add(mul({(1,1,0):Fr(1)}, g1), mul(h0, g2)).items():
        Phi[e[0]] = Phi.get(e[0], 0) + c    # power of x == power of y == power of s
    Phi = {k: v for k, v in Phi.items() if v}
    # rational roots of Phi
    import math
    roots = []
    if Phi:
        deg = max(Phi)
        lead = Phi[deg]; c0 = Phi.get(0, 0)
        cands = set()
        def divisors(n):
            n = abs(int(n)); out = set()
            for i in range(1, min(n, 2000) + 1):
                if n % i == 0: out.add(i); out.add(n // i)
            return out or {1}
        num0 = c0.numerator if isinstance(c0, Fr) else int(c0)
        den0 = lead.numerator if isinstance(lead, Fr) else int(lead)
        for pn in divisors(num0 if num0 else 1) | {0}:
            for qd in divisors(den0):
                for sgn in (1, -1):
                    cands.add(Fr(sgn * pn, qd))
        for cand in cands:
            if sum(v * cand**k for k, v in Phi.items()) == 0:
                roots.append(cand)
    print("d=%d: det = %s ; Phi(s) coeffs %s ; rational roots %s" % (d, detval, dict(sorted(Phi.items())), sorted(roots)))
    # explicit collisions: x=0 sheet point + the s0-branch points at x = +-1
    for s0 in roots:
        z0 = pev(h0, 1, s0, 0)      # h0(s0): evaluate with x=1,y=s0 making s=s0
        # branch point (x, y, z) = (1, s0, h0(s0)) since r = x^2 z = h0(s0) ✓
        p1 = (Fr(1), s0, z0)
        p2 = (Fr(-1), -s0, z0)
        w = tuple(pev(f, *p1) for f in Fd)
        w2 = tuple(pev(f, *p2) for f in Fd)
        p0 = (Fr(0), Fr(0), Fr(w[0]) / pev({(0,0,0):Fr(1)}, 0, 0, 0) if False else None)
        # x=0 sheet: F(0,0,z) = (z*A(0), 0, 0); A(0) = 1 => z0ax = w[0]
        p0 = (Fr(0), Fr(0), w[0])
        w0 = tuple(pev(f, *p0) for f in Fd)
        print("   collision @ s0=%s: F%s = %s ; F%s = %s ; F%s = %s  EQUAL: %s"
              % (s0, p1, w, p2, w2, p0, w0, w == w2 == w0))
    # mod-p census
    for p in (11, 13):
        fs = [{k: int(v) % p for k, v in f.items()} for f in Fd]
        cnt = {}
        for x, y, z in product(range(p), repeat=3):
            img = tuple(sum(c * pow(x, i, p) * pow(y, j, p) * pow(z, k, p)
                            for (i, j, k), c in f.items()) % p for f in fs)
            cnt[img] = cnt.get(img, 0) + 1
        hist = {}
        for v in cnt.values(): hist[v] = hist.get(v, 0) + 1
        defc = p**3 - len(cnt)
        print("   mod %d: fibers %s deficiency %d %s" % (p, dict(sorted(hist.items())), defc,
              "NON-BIJECTIVE => genuine counterexample" if defc else "bijective?!"))
print("DONE.")
