#!/usr/bin/env python3
r"""
.scratch/twojet_g1_gb_wildcert_kps.py  (kind-pasteur 2026-07-26, G1 probe)

GB EMPTINESS with WILDNESS CERTIFICATE for the 2-jet Keller system of
THM-2446 on small coefficient boxes, protocol of
04-computation/jacobian_k45_groebner_verdict_kps_S128c99.py
(Rabinowitsch over GF(32003), GB == [1] test, grevlex).

Map shape:  F = A(x,y) z^2 + B(x,y) z + C(x,y),  A,B,C : C^2 -> C^3.
Keller <=> D5=D4=D3=D2=D1=0 in C[x,y], D0 = nonzero constant
(THM-2446 eq (3)); Rabinowitsch d00*t - 1 forces D0(0,0) != 0 and the
nonconstant coefficients of D0 are equated to zero.

WILDNESS CERTIFICATES (THM-2446 hunt-design warning: the box contains
tame maps, plain consistency means NOTHING):

(a) alpha-conic nondegeneracy (box (2,2,2) only): the 3x3 matrix M of
    A's degree-2 part (rows = components, cols = x^2, xy, y^2 coeffs)
    has rank 3, imposed by Rabinowitsch det(M)*s - 1 = 0.
(b) fiber-degree->=4 certificate (boxes with dA = 1): auxiliary fiber
    points.  Gauge normalization (all in-box: source translation +
    constant z-translation + target translation): fiber point 1 =
    (0,0,0), target P = F(0,0,0) = C(0,0) = 0 (so c*_00 = 0 is fixed).
    Aux vars (x_i,y_i,z_i), i=2,3,4 with F(x_i,y_i,z_i) = 0 and
    Rabinowitsch  z2*z3*z4*(z2-z3)*(z2-z4)*(z3-z4)*s - 1 = 0
    (four fiber points with pairwise-distinct z; for a quasi-finite
    dominant map 4 distinct fiber points => field degree >= 4).

Gauge fact used everywhere: C enters every D_k only through C_x, C_y
(third bracket slots are always A or B), so the constant part of C is
pure target-translation gauge; we fix C(0,0) = 0 in every box run.

Modes: control | box111 | box122 | box222 | box111nocert
"""
import sys
import time
import sympy as sp

x, y, z, v = sp.symbols("x y z v")
t_, s_ = sp.symbols("t_ s_")
PRIME = 32003


def bracket(u, vv, w):
    return sp.Matrix.hstack(sp.Matrix(u), sp.Matrix(vv), sp.Matrix(w)).det()


def d_coeffs(A, B, C):
    """THM-2446 eq (2); verbatim convention of the S131 companion."""
    Ax = [sp.diff(a, x) for a in A]; Ay = [sp.diff(a, y) for a in A]
    Bx = [sp.diff(b, x) for b in B]; By = [sp.diff(b, y) for b in B]
    Cx = [sp.diff(c, x) for c in C]; Cy = [sp.diff(c, y) for c in C]
    D5 = 2 * bracket(Ax, Ay, A)
    D4 = bracket(Ax, Ay, B) + 2 * bracket(Ax, By, A) + 2 * bracket(Bx, Ay, A)
    D3 = (bracket(Ax, By, B) + bracket(Bx, Ay, B)
          + 2 * (bracket(Ax, Cy, A) + bracket(Cx, Ay, A) + bracket(Bx, By, A)))
    D2 = (bracket(Bx, By, B) + bracket(Ax, Cy, B) + bracket(Cx, Ay, B)
          + 2 * (bracket(Bx, Cy, A) + bracket(Cx, By, A)))
    D1 = bracket(Bx, Cy, B) + bracket(Cx, By, B) + 2 * bracket(Cx, Cy, A)
    D0 = bracket(Cx, Cy, B)
    return [D0, D1, D2, D3, D4, D5]


def generic_vec(name, deg, drop00=False):
    coeffs, comp = [], []
    for i in range(3):
        e = sp.Integer(0)
        for dx in range(deg + 1):
            for dy in range(deg + 1 - dx):
                if drop00 and dx == 0 and dy == 0:
                    continue  # gauge: C(0,0) = 0 fixed by target translation
                sname = sp.Symbol(f"{name}{i}_{dx}{dy}")
                coeffs.append(sname)
                e += sname * x ** dx * y ** dy
        comp.append(e)
    return comp, coeffs


def graded_equations(A, B, C):
    Ds = d_coeffs(A, B, C)
    eqs = []
    for k in range(1, 6):
        Dk = sp.expand(Ds[k])
        if Dk == 0:
            continue
        eqs += [c for c in sp.Poly(Dk, x, y).coeffs() if c != 0]
    D0 = sp.expand(Ds[0])
    d00 = D0.subs({x: 0, y: 0})
    rest = sp.expand(D0 - d00)
    if rest != 0:
        eqs += [c for c in sp.Poly(rest, x, y).coeffs() if c != 0]
    eqs = list(dict.fromkeys(eqs))
    return eqs, d00


def fiber_cert(A, B, C):
    """Certificate (b): 3 extra fiber points over P = C(0,0) = 0 with
    pairwise-distinct nonzero z (point 1 = origin, z1 = 0)."""
    P = [c.subs({x: 0, y: 0}) for c in C]
    fib, eqs, zs = [], [], []
    for i in (2, 3, 4):
        xi, yi, zi = sp.symbols(f"xf{i} yf{i} zf{i}")
        fib += [xi, yi, zi]
        zs.append(zi)
        for j in range(3):
            eqs.append(sp.expand(
                A[j].subs({x: xi, y: yi}) * zi ** 2
                + B[j].subs({x: xi, y: yi}) * zi
                + C[j].subs({x: xi, y: yi}) - P[j]))
    dist = zs[0] * zs[1] * zs[2] * (zs[0] - zs[1]) * (zs[0] - zs[2]) * (zs[1] - zs[2])
    return fib, eqs, sp.expand(dist)


def alpha_matrix_det(ca):
    """Certificate (a): det of the 3x3 quadratic-form coefficient matrix
    of A's degree-2 part (cols x^2, xy, y^2)."""
    d = {str(sname): sname for sname in ca}
    M = sp.Matrix(3, 3, lambda i, j:
                  d[f"a{i}_{['20','11','02'][j]}"])
    return sp.expand(M.det())


def run_gb(eqs, gens, label, t0):
    print(f"  system: {len(eqs)} eqs, {len(gens)} gens  [{time.time()-t0:.0f}s]",
          flush=True)

    def clearden(e):
        _, pe = sp.Poly(sp.expand(e), *gens, domain="QQ").clear_denoms()
        return pe.as_expr()

    eqs = [clearden(e) for e in eqs]
    G = sp.groebner(eqs, *gens, modulus=PRIME, order="grevlex")
    exprs = list(G.exprs)
    if exprs == [sp.Integer(1)] or exprs == [1]:
        print(f"  GB = [1]  =>  EMPTY over GF({PRIME})-bar "
              f"[{label}] [{time.time()-t0:.0f}s]", flush=True)
        return "EMPTY"
    print(f"  GB size {len(exprs)} != [1]  =>  CONSISTENT mod {PRIME} "
          f"[{label}] [{time.time()-t0:.0f}s]", flush=True)
    return "CONSISTENT"


def box_run(dA, dB, dC, cert, label):
    t0 = time.time()
    print(f"\n==== {label}: box (dA,dB,dC)=({dA},{dB},{dC}), cert={cert} ====",
          flush=True)
    A, ca = generic_vec("a", dA)
    B, cb = generic_vec("b", dB)
    C, cc = generic_vec("c", dC, drop00=True)
    unknowns = ca + cb + cc
    eqs, d00 = graded_equations(A, B, C)
    ngraded = len(eqs)
    eqs = eqs + [sp.expand(d00 * t_ - 1)]
    gens = list(unknowns)
    if cert == "fiber4":
        fib, feqs, dist = fiber_cert(A, B, C)
        eqs += feqs + [sp.expand(dist * s_ - 1)]
        gens += fib + [s_, t_]
        print(f"  box coeffs {len(unknowns)} (c*_00 gauge-fixed to 0), "
              f"fiber vars {len(fib)}, +s,t; graded eqs {ngraded}, "
              f"fiber eqs {len(feqs)}, dist-Rabinowitsch 1, d00-Rabinowitsch 1",
              flush=True)
    elif cert == "alpha3":
        detM = alpha_matrix_det(ca)
        eqs += [sp.expand(detM * s_ - 1)]
        gens += [s_, t_]
        print(f"  box coeffs {len(unknowns)} (c*_00 gauge-fixed to 0), +s,t; "
              f"graded eqs {ngraded}, alpha-rank3-Rabinowitsch 1, "
              f"d00-Rabinowitsch 1", flush=True)
    else:
        gens += [t_]
        print(f"  box coeffs {len(unknowns)} (c*_00 gauge-fixed to 0), +t; "
              f"graded eqs {ngraded}, d00-Rabinowitsch 1; NO WILDNESS CERT "
              f"(tame witness (x+z^2,y,z) is in-box => consistency expected, "
              f"meaningless for G1)", flush=True)
    return run_gb(eqs, gens, label, t0)


# ---------------- controls ----------------

def control():
    u = 1 + x * y
    Bc = [u ** 3, 3 * x * u ** 2, -x ** 3]
    Cc = [y ** 2 * u * (4 + 3 * x * y),
          y + 3 * x * y ** 2 * (4 + 3 * x * y),
          2 * x - 3 * x ** 2 * y]
    Z3 = [sp.Integer(0)] * 3

    # C1: THM-1310 z-affine embed; graded system must hold exactly,
    # GB with d00-Rabinowitsch, NO alpha certificate => CONSISTENT.
    t0 = time.time()
    print("\n==== C1: THM-1310 embed (A=0), alpha-Rabinowitsch DROPPED ====",
          flush=True)
    Ds = d_coeffs(Z3, Bc, Cc)
    assert all(sp.expand(Ds[k]) == 0 for k in range(1, 6)), "graded system fails"
    assert sp.expand(Ds[0]) == -2, "D0 != -2"
    print("  exact: D5=D4=D3=D2=D1=0, D0=-2 (S131 cancellation reproduced)",
          flush=True)
    r = run_gb([sp.expand(Ds[0] * t_ - 1)], [t_], "C1", t0)
    assert r == "CONSISTENT", "C1 must be CONSISTENT"

    # C2: same + alpha-rank-3 certificate => must flip to EMPTY
    # (field degree 3 map correctly OUTSIDE the certified stratum).
    t0 = time.time()
    print("\n==== C2: THM-1310 embed + alpha-rank3 Rabinowitsch ====", flush=True)
    detM = sp.Integer(0)  # A = 0 => quadratic-form matrix = 0
    r = run_gb([sp.expand(Ds[0] * t_ - 1), sp.expand(detM * s_ - 1)],
               [s_, t_], "C2", t0)
    assert r == "EMPTY", "C2 must be EMPTY"

    # C3: tame witness (x+z^2, y, z) in box (1,1,1) + fiber4 certificate
    # => must be EMPTY (field degree 1: no 4-point fiber exists).
    t0 = time.time()
    print("\n==== C3: tame witness (x+z^2,y,z) + fiber4 certificate ====",
          flush=True)
    At = [sp.Integer(1), sp.Integer(0), sp.Integer(0)]
    Bt = [sp.Integer(0), sp.Integer(0), sp.Integer(1)]
    Ct = [x, y, sp.Integer(0)]
    Dst = d_coeffs(At, Bt, Ct)
    assert all(sp.expand(Dst[k]) == 0 for k in range(1, 6)) and Dst[0] == 1
    fib, feqs, dist = fiber_cert(At, Bt, Ct)
    r = run_gb(feqs + [sp.expand(dist * s_ - 1), sp.expand(Dst[0] * t_ - 1)],
               fib + [s_, t_], "C3", t0)
    assert r == "EMPTY", "C3 must be EMPTY"

    # C4: k=3-style family control (port of S128c99 verdict(3,4,2,3,3),
    # determinant routed through the 2-jet graded brackets with A=0).
    t0 = time.time()
    print("\n==== C4: k=3-style family (S128c99 box, graded route, "
          "no alpha cert) ====", flush=True)
    kap = sp.Symbol("kap")
    syms = []
    W = sp.Integer(0)
    for i in range(5):
        for j in range(5 - i):
            cs = sp.Symbol(f"w{i}_{j}"); syms.append(cs); W += cs * x**i * y**j
    S = sp.Integer(0)
    for i in range(3):
        for j in range(3 - i):
            cs = sp.Symbol(f"sc{i}_{j}"); syms.append(cs); S += cs * x**i * y**j
    R = sp.Integer(0)
    for m_ in range(4):
        for n_ in range(4):
            cs = sp.Symbol(f"r{m_}_{n_}"); syms.append(cs); R += cs * x**m_ * u**n_
    B1 = u * W; B2 = S + kap * x * W
    P = sp.expand(R - x**3 * B1)
    Pv = sp.together(P.subs(y, (v - 1) / x))
    num, _ = sp.fraction(sp.cancel(Pv)); num = sp.expand(num)
    pn = sp.Poly(num, v)
    lin = []
    for j in range(min(3, pn.degree() + 1)):
        for c_ in sp.Poly(pn.nth(j), x).all_coeffs():
            if c_ != 0:
                lin.append(c_)
    lin = list(dict.fromkeys(lin))
    sol1 = sp.solve(lin, syms, dict=True)
    sub1 = sol1[0] if sol1 else {}
    W1, S1, R1, B11, B21 = [e.subs(sub1) for e in (W, S, R, B1, B2)]
    B31 = sp.cancel((R1 - x**3 * B11) / u**3)
    Bfam = [u**3, kap * u**2 * x, -x**3]     # 2-jet B (z-coefficient)
    Cfam = [sp.expand(B11), sp.expand(B21), sp.expand(B31)]  # 2-jet C
    Dsf = d_coeffs(Z3, Bfam, Cfam)
    assert all(sp.expand(Dsf[k]) == 0 for k in (3, 4, 5)), "A=0 => D3,D4,D5 == 0"
    eqs = []
    for k in (1, 2):
        Dk = sp.expand(Dsf[k])
        if Dk != 0:
            eqs += [c for c in sp.Poly(Dk, x, y).coeffs() if c != 0]
    D0f = sp.expand(Dsf[0])
    d00 = D0f.subs({x: 0, y: 0})
    rest = sp.expand(D0f - d00)
    if rest != 0:
        eqs += [c for c in sp.Poly(rest, x, y).coeffs() if c != 0]
    eqs = list(dict.fromkeys(eqs))
    rem = sorted(set().union(*[e.free_symbols for e in eqs + [d00]]) - {x, y},
                 key=str)
    r = run_gb(eqs + [sp.expand(d00 * t_ - 1)], rem + [t_], "C4", t0)
    assert r == "CONSISTENT", "C4 must be CONSISTENT (THM-1310 family lives here)"
    print("\nCONTROLS: ALL PASS", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "control"
    if mode == "control":
        control()
    elif mode == "box111":
        box_run(1, 1, 1, "fiber4", "G1 box (1,1,1) + fiber-degree>=4 cert")
    elif mode == "box122":
        box_run(1, 2, 2, "fiber4", "G1 box (1,2,2) + fiber-degree>=4 cert")
    elif mode == "box222":
        box_run(2, 2, 2, "alpha3", "G1 box (2,2,2) + alpha-rank-3 cert")
    elif mode == "box111nocert":
        box_run(1, 1, 1, None, "box (1,1,1) NO certificate (tame-inhabited)")
    else:
        raise SystemExit("unknown mode " + mode)
    print("DONE " + mode, flush=True)
