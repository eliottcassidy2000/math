#!/usr/bin/env python3
"""
jacobian_port_weight3_boxeph_S143.py  (HYP-8110)

THE METHOD PORTED to weights (-1,1,3) -> (3,1,-1), z-linear:
   F = ( z A(s) + y^3 B(s),  y g1(s) + x^2 z g2(s),  x h0(s) - x^4 z ),
   s = xy, r = x^3 z, w = h0 - r.
Derived by hand (verified below): det JF * w^3 = -Jac_{s,r}(P,Q) with
P = F2 F3 = (s g1 + r g2) w,  Q = F1 F3^3 = (r A + s^3 B) w^3;
Jac = w^2 (b1 + b2 w + b3 w^2), b0 = 0 automatic, and Keller <=> bracket = -c w:
   b1 = PHI PSI' - 3 PHI' PSI  (= -c),   PHI = s g1 + h0 g2,  PSI = h0 A + s^3 B
   b2 = 4 A PHI' - A' PHI + 3 g2' PSI - 2 g2 PSI'  (= 0)
   b3 = 2 g2 A' - 4 A g2'  (= 0)   <=>  A = alpha v^2, g2 = gamma v (SQUARE LAW).
Log obstruction (same proof shape): PHI constant (automorphism-type) or LINEAR;
axis fiber = 1 + 3 * #roots(PHI): candidate kernels are DEGREE-4 covers with a
Z/3 torus element replacing sigma.  b2's leading coefficient ~ (2-k):
k = deg v = 2 distinguished; k = 1 forces K = 0 in PSI = K PHI^3 + c/(3 phi1).

THE HUNT: solve both branches exactly, impose the reconstruction congruences
   g1 = (PHI - h0 g2)/s  poly  <=>  PHI(0) = h0(0) g2(0)
   B  = (PSI - h0 A)/s^3 poly  <=>  PSI - h0 A == 0 mod s^3   (3 conditions),
assemble any solution into F, verify det symbolically, find explicit collisions,
mod-p census (p = 7, 11, 13 -- cube-root Chebotarev fingerprint), certify
non-bijectivity.  boxeph-2026-07-19-S143.
"""

import random
from fractions import Fraction as Fr
from itertools import product

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
        for k, v in p.items(): out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}
def sc(p, c): return {k: v * c for k, v in p.items() if v * c}
def diff(p, var):
    out = {}
    for e, c in p.items():
        if e[var]:
            k = list(e); k[var] -= 1
            out[tuple(k)] = out.get(tuple(k), 0) + c * e[var]
    return out
def det3(M):
    return add(mul(M[0][0], add(mul(M[1][1],M[2][2]), sc(mul(M[1][2],M[2][1]),-1))),
               sc(mul(M[0][1], add(mul(M[1][0],M[2][2]), sc(mul(M[1][2],M[2][0]),-1))),-1),
               mul(M[0][2], add(mul(M[1][0],M[2][1]), sc(mul(M[1][1],M[2][0]),-1))))

rng = random.Random(143)
def rpoly(deg):
    out = {}
    for k in range(deg + 1):
        c = rng.randint(-3, 3)
        if c: out[(k,)] = Fr(c)
    out[(deg,)] = Fr(rng.randint(1, 3))
    return out

print("=" * 96)
print("(M') master factorization + b-forms for the (-1,1,3) class: 6 random instances")
print("=" * 96)
def build_xyz(A1, B1, G11, G21, H01):
    def sub_s(p1):
        return {(k, k, 0): c for (k,), c in p1.items()}
    X, Y, Z = {(1,0,0):Fr(1)}, {(0,1,0):Fr(1)}, {(0,0,1):Fr(1)}
    A, B, g1, g2, h0 = map(sub_s, (A1, B1, G11, G21, H01))
    f1 = add(mul(Z, A), mul(mul(mul(Y, Y), Y), B))
    f2 = add(mul(Y, g1), mul(mul(mul(X, X), Z), g2))
    f3 = add(mul(X, h0), sc(mul(mul(mul(mul(X,X),X),X), Z), -1))
    return [f1, f2, f3]
for trial in range(6):
    A1, B1 = rpoly(rng.randint(1,2)), rpoly(rng.randint(0,2))
    G11, G21, H01 = rpoly(rng.randint(0,2)), rpoly(rng.randint(0,2)), rpoly(rng.randint(0,1))
    Fs = build_xyz(A1, B1, G11, G21, H01)
    D = det3([[diff(f, v) for v in range(3)] for f in Fs])
    def c2(p1): return {(k, 0): c for (k,), c in p1.items()}
    Sv, Rv = {(1,0):Fr(1)}, {(0,1):Fr(1)}
    A2, B2, g12, g22, h02 = map(c2, (A1, B1, G11, G21, H01))
    w = add(h02, sc(Rv, -1))
    P = mul(add(mul(Sv, g12), mul(Rv, g22)), w)
    Q = mul(add(mul(Rv, A2), mul(mul(mul(Sv,Sv),Sv), B2)), mul(mul(w, w), w))
    Jpq = add(mul(diff(P,0), diff(Q,1)), sc(mul(diff(P,1), diff(Q,0)), -1))
    def embed(p2):
        return {(i + 3*j, i, j): c for (i, j), c in p2.items()}
    okM = (mul(D, embed(mul(mul(w,w),w))) == sc(embed(Jpq), -1))
    PHI = add(mul(Sv, g12), mul(h02, g22))
    PSI = add(mul(h02, A2), mul(mul(mul(Sv,Sv),Sv), B2))
    b1 = add(mul(PHI, diff(PSI,0)), sc(mul(diff(PHI,0), PSI), -3))
    b2 = add(sc(mul(A2, diff(PHI,0)), 4), sc(mul(diff(A2,0), PHI), -1),
             sc(mul(PSI, diff(g22,0)), 3), sc(mul(g22, diff(PSI,0)), -2))
    b3 = add(sc(mul(g22, diff(A2,0)), 2), sc(mul(A2, diff(g22,0)), -4))
    recon = mul(mul(mul(w, w), w), add(b1, mul(b2, w), mul(b3, mul(w, w))))
    okB = (Jpq == recon)
    print("trial %d: master %s ; bracket %s" % (trial, okM, okB))
    assert okM and okB
print("=> ported identities CONFIRMED.")

print("\n" + "=" * 96)
print("(HUNT) assemble kernels: gauge v = 1+s (k=1, K=0) and v = quadratics (k=2)")
print("=" * 96)

def try_branch(vco, Kfix, label, h0deg=2, gridvals=(Fr(0),Fr(1),Fr(-1),Fr(2),Fr(-2),Fr(3),Fr(-3),Fr(1,2),Fr(-1,2),Fr(4),Fr(6),Fr(-6),Fr(2,3),Fr(-3,2),Fr(9),Fr(1,3))):
    """v fixed; unknowns alpha,gamma,phi0,phi1,K(,eta*). Solve:
       b2 == 0, PHI(0)=h0(0)g2(0), PSI-h0*A == 0 mod s^3, all exact.
       Strategy: small grid on (alpha,gamma,phi1,K) with phi0,eta* solved linearly."""
    v = {(k,): c for k, c in enumerate(vco) if c}
    found = []
    for alpha in (Fr(1),):
        for gamma in (Fr(1), Fr(2), Fr(3), Fr(1,2)):
            A1 = sc(mul(v, v), alpha); g21 = sc(v, gamma)
            for phi1 in (Fr(1), Fr(2), Fr(3), Fr(4), Fr(1,2)):
                for c in (Fr(-3), Fr(-2), Fr(-1), Fr(1), Fr(3)):
                    E = c / (3 * phi1)
                    for K in ((Fr(0),) if Kfix else (Fr(0), Fr(1), Fr(1,2), Fr(1,27), Fr(-1), Fr(1,8), Fr(1,64))):
                        for phi0 in gridvals:
                            PHI = {(0,): phi0, (1,): phi1}
                            PHI3 = mul(PHI, mul(PHI, PHI))
                            PSI = add(sc(PHI3, K), {(0,): E})
                            b2 = add(sc(mul(A1, diff(PHI,0)), 4), sc(mul(diff(A1,0), PHI), -1),
                                     sc(mul(PSI, diff(g21,0)), 3), sc(mul(g21, diff(PSI,0)), -2))
                            if b2 != {}: continue
                            # reconstruction: h0 (deg <= h0deg) with PSI - h0 A == 0 mod s^3
                            # solve linearly for eta coefficients from the s^0,s^1,s^2 conditions
                            import itertools
                            nA = {k[0]: cc for k, cc in A1.items()}
                            nP = {k[0]: cc for k, cc in PSI.items()}
                            # (h0 A)_m = sum eta_i nA_{m-i}
                            # unknowns eta_0..eta_{h0deg}; 3 linear equations
                            rowsM = []
                            for m in range(3):
                                rowsM.append([nA.get(m - i, Fr(0)) for i in range(h0deg + 1)] + [nP.get(m, Fr(0))])
                            # gaussian solve (underdetermined ok -> parameterize free vars = 0)
                            M = [r[:] for r in rowsM]; piv = []
                            rr = 0
                            for ccol in range(h0deg + 1):
                                p = None
                                for r2 in range(rr, 3):
                                    if M[r2][ccol] != 0: p = r2; break
                                if p is None: continue
                                M[rr], M[p] = M[p], M[rr]
                                pr = M[rr]; inv = 1 / pr[ccol]
                                for r2 in range(3):
                                    if r2 != rr and M[r2][ccol] != 0:
                                        f2_ = M[r2][ccol] * inv
                                        for c2_ in range(h0deg + 2): M[r2][c2_] -= f2_ * pr[c2_]
                                piv.append(ccol); rr += 1
                            if any(all(M[r2][c2_] == 0 for c2_ in range(h0deg + 1)) and M[r2][h0deg + 1] != 0 for r2 in range(3)):
                                continue
                            eta = [Fr(0)] * (h0deg + 1)
                            for idx, ccol in enumerate(piv):
                                eta[ccol] = M[idx][h0deg + 1] / M[idx][ccol]
                            h0 = {(k,): e for k, e in enumerate(eta) if e}
                            # g1 polynomial: PHI(0) = h0(0) g2(0)
                            if PHI.get((0,), Fr(0)) != (h0.get((0,), Fr(0)) * g21.get((0,), Fr(0))):
                                continue
                            # assemble and verify
                            g1num = add(PHI, sc(mul(h0, g21), -1))
                            if any(k[0] == 0 for k in g1num): continue
                            G11 = {(k[0]-1,): cc for k, cc in g1num.items()}
                            Bnum = add(PSI, sc(mul(h0, A1), -1))
                            if any(k[0] < 3 for k in Bnum): continue
                            B1 = {(k[0]-3,): cc for k, cc in Bnum.items()}
                            Fs = build_xyz(A1, B1, G11, g21, h0)
                            D = det3([[diff(f, vv) for vv in range(3)] for f in Fs])
                            if list(D.keys()) == [(0,0,0)]:
                                found.append((dict(A1), dict(B1), dict(G11), dict(g21), dict(h0), D[(0,0,0)], PHI, PSI))
    print("branch %s: Keller assemblies found: %d" % (label, len(found)))
    return found

sols = try_branch([Fr(1), Fr(1)], Kfix=True, label="k=1 (v=1+s, K=0)")
sols += try_branch([Fr(1), Fr(1), Fr(1)], Kfix=False, label="k=2 (v=1+s+s^2)")
sols += try_branch([Fr(1), Fr(2), Fr(1)], Kfix=False, label="k=2 (v=(1+s)^2)")
sols += try_branch([Fr(1), Fr(0), Fr(1)], Kfix=False, label="k=2 (v=1+s^2)")

def pev(p, x, y, z):
    return sum(c * x**e[0] * y**e[1] * z**e[2] for e, c in p.items())
seen_shapes = set()
for (A1, B1, G11, g21, h0, detv, PHI, PSI) in sols[:6]:
    key = (tuple(sorted(A1.items())), tuple(sorted(h0.items())))
    if key in seen_shapes: continue
    seen_shapes.add(key)
    Fs = build_xyz(A1, B1, G11, g21, h0)
    print("\nKERNEL CANDIDATE: det =", detv)
    print("  A =", A1, " B =", B1, " g1 =", G11, " g2 =", g21, " h0 =", h0)
    # resolvent root(s)
    phi0, phi1 = PHI.get((0,), Fr(0)), PHI.get((1,), Fr(0))
    rho = -phi0 / phi1
    psi_rho = sum(c * rho**k[0] for k, c in PSI.items())
    print("  resolvent root s0 = %s ; PSI(s0) = %s" % (rho, psi_rho))
    # collision test: x=1 branch point (1, rho, h0(rho)) vs axis (0,0, z*) if images match
    h0r = sum(c * rho**k[0] for k, c in h0.items())
    p1 = (Fr(1), rho, h0r)
    img1 = tuple(pev(f, *p1) for f in Fs)
    A0 = A1.get((0,), Fr(0))
    p0 = (Fr(0), Fr(0), img1[0] / A0)
    img0 = tuple(pev(f, *p0) for f in Fs)
    print("  collision: F%s = %s ; F%s = %s ; EQUAL: %s" % (p1, img1, p0, img0, img1 == img0))
    for p in (7, 11, 13):
        fs = [{k: int(vv) % p for k, vv in f.items()} for f in Fs]
        cnt = {}
        for x, y, z in product(range(p), repeat=3):
            im = tuple(sum(cc * pow(x, i, p) * pow(y, j, p) * pow(z, k, p)
                           for (i, j, k), cc in f.items()) % p for f in fs)
            cnt[im] = cnt.get(im, 0) + 1
        hist = {}
        for vv in cnt.values(): hist[vv] = hist.get(vv, 0) + 1
        defc = p**3 - len(cnt)
        print("  mod %2d: fibers %s deficiency %d (p%%3=%d) %s" % (p, dict(sorted(hist.items())), defc, p % 3,
              "NON-BIJECTIVE" if defc else "bijective"))
print("\nDONE.")
