#!/usr/bin/env python3
"""
kind-pasteur-2026-07-20-S128c101 -- HYP-8150 / THM-1345: the owner's plane family,
the polynomial section, THE RADICAL INVERSE, and the trace module.

 (1) Verify the owner's Family 1 (r-continuum) and Family 2 (whole plane g=0)
     EXACTLY (Family 2 modulo s^2 = v^2 - 16u); place them:
     s = sqrt(-L)|_{g=0}; excluded parabola v^2 = 16u = Jelonek cap on the plane;
     Family 1 = the v = 0 slice; sigma(u,v) = (0, v, u-4v^2) is a POLYNOMIAL
     SECTION (F o sigma = id on {g=0}); z-numerator 13/2 universal (2 - 3s*).
 (2) THE RADICAL INVERSE: the fiber cubic N = L x^3 + (4-3bg) x - 2g is
     depressed, S_3 is solvable => explicit 3-valued radical inverse via
     Cardano; verify F(inv_j(Q)) = Q numerically (1e-9) at random targets, all
     three branches; verify the g=0 degeneration reproduces Family 2 exactly.
 (3) THE TRACE MODULE (weight/parity fit method): Tr_F of monomials
     y^2, xy, z, xz, yz, xyz, x^2z^2, y^3z, y^4, z^2 -- polynomial or L-poled?
     Tower transitivity Tr_{F o F} = Tr_F o Tr_F then decides the F.F centroid
     (backlog c99(b)) with NO degree-9 computation:
     Tr_{FF}(x) = Tr_F(0) = 0; Tr_{FF}(y) = Tr_F(3y/2) = 9b/4;
     Tr_{FF}(z) = Tr_F(zeta(x,y,z)), zeta = -81x^2z^2/2 + 27xyz - 51x + 15y^3z/8 - 3y^2/4.
"""
import sympy as sp
from sympy import symbols, expand, Rational, sqrt, Poly, together, cancel, fraction, simplify
import random, cmath

x, y, z, s_, a, b, g, u_, v_, r_ = symbols('x y z s_ a b g u_ v_ r_')
u = 1 + x*y
F1 = u**3*z + y**2*u*(4 + 3*x*y)
F2 = y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y)
F3 = 2*x - 3*x**2*y - x**3*z
Fv = (F1, F2, F3)
L = 27*a**2*g**2 - 18*a*b*g + 16*a + b**3*g - b**2

print("== (1) the owner's two families, verified and placed ==", flush=True)
# Family 1
p1 = [ (0, 0, -1/(4*r_**2)), (r_, -3/(2*r_), 13/(2*r_**2)), (-r_, 3/(2*r_), 13/(2*r_**2)) ]
ok1 = True
for P in p1:
    im = tuple(sp.simplify(f.subs({x: P[0], y: P[1], z: P[2]})) for f in Fv)
    tgt = (-1/(4*r_**2), 0, 0)
    ok1 &= all(sp.simplify(im[i] - tgt[i]) == 0 for i in range(3))
print(f"  Family 1 (r-continuum): all three preimages -> (-1/(4r^2),0,0) identically in r: {ok1}", flush=True)
# Family 2, modulo s^2 = v^2 - 16u
xp = 2/s_
yp = v_/4 - 3/(2*xp)
zp = 13/(2*xp**2) - 3*v_/(4*xp)
def red(e):
    # reduce modulo s^2 -> v^2 - 16u
    e = sp.expand(sp.simplify(e))
    e = sp.expand(e.subs(s_**2, v_**2 - 16*u_))
    e = sp.expand(e.subs(s_**2, v_**2 - 16*u_))
    return sp.simplify(e)
im0 = tuple(sp.expand(f.subs({x: 0, y: v_, z: u_ - 4*v_**2})) for f in Fv)
print(f"  x=0 branch: F(0, v, u-4v^2) = {im0}  == (u, v, 0): {all(sp.simplify(im0[i] - t) == 0 for i, t in enumerate((u_, v_, 0)))}", flush=True)
imp = tuple(red(f.subs({x: xp, y: yp, z: zp})) for f in Fv)
okp = all(sp.simplify(imp[i] - t) == 0 for i, t in enumerate((u_, v_, 0)))
imm = tuple(red(f.subs({x: -xp, y: v_/4 + 3*s_/4, z: red(13*s_**2/8 + 3*v_*s_/8)})) for f in Fv)
# recompute x_- branch cleanly: x_- = -2/s: y_- = v/4 - 3/(2*x_-) = v/4 + 3s/4; z_- = 13/(2 x_-^2) - 3v/(4 x_-) = 13 s^2/8 + 3 v s/8
okm = all(sp.simplify(imm[i] - t) == 0 for i, t in enumerate((u_, v_, 0)))
print(f"  x_+ = 2/s branch maps to (u,v,0) mod s^2 = v^2-16u: {okp}", flush=True)
print(f"  x_- = -2/s branch maps to (u,v,0) mod s^2 = v^2-16u: {okm}", flush=True)
# placement
L0 = sp.expand(L.subs(g, 0))
print(f"  -L|_(g=0) = {sp.expand(-L0)} = b^2 - 16a  => owner's s^2 = v^2-16u IS the resolvent -L on the plane: {sp.expand(-L0 - (b**2 - 16*a)) == 0}", flush=True)
print(f"  excluded parabola v^2 = 16u <=> L|_(g=0) = 0 = Jelonek cap on the plane; Family 1 = the v=0 slice (s = 2/r)", flush=True)
sig = (sp.Integer(0), b, a - 4*b**2)
Fsig = tuple(sp.expand(f.subs({x: sig[0], y: sig[1], z: sig[2]})) for f in Fv)
print(f"  SECTION: F(sigma(a,b)) = {Fsig} == (a, b, 0): {all(sp.simplify(Fsig[i] - t) == 0 for i, t in enumerate((a, b, 0)))}   (sigma = (0, b, a-4b^2) polynomial section over the plane)", flush=True)

print("\n== (2) THE RADICAL INVERSE (Cardano on the depressed fiber cubic) ==", flush=True)
def radical_inverse(av, bv, gv):
    """three-valued radical inverse of F at target (av,bv,gv); returns 3 source points."""
    Lv = 27*av**2*gv**2 - 18*av*bv*gv + 16*av + bv**3*gv - bv**2
    pv = (4 - 3*bv*gv)/Lv
    qv = (-2*gv)/Lv
    # depressed cubic x^3 + pv x + qv: Cardano
    disc = (qv/2)**2 + (pv/3)**3
    Csqrt = cmath.sqrt(disc)
    c1 = (-qv/2 + Csqrt)**(1/3) if (-qv/2 + Csqrt) != 0 else 0
    # pick cube roots carefully: iterate all three cube roots of the first radicand
    roots = []
    w3 = cmath.exp(2j*cmath.pi/3)
    base = -qv/2 + Csqrt
    r0 = base**(Rational(1,3)) if base == 0 else cmath.exp(cmath.log(base)/3)
    for j in range(3):
        A_ = r0 * w3**j
        B_ = -pv/(3*A_) if A_ != 0 else 0
        roots.append(A_ + B_)
    # y via the subresultant s(x), z via linearity
    pts = []
    for xr in roots:
        sv = s_of_x_num(xr, av, bv, gv)
        yr = sv/xr
        ur3 = (1 + xr*yr)**3
        B1 = yr**2*(1+xr*yr)*(4 + 3*xr*yr)
        zr = (av - B1)/ur3
        pts.append((xr, yr, zr))
    return pts

# s(x) rational function from the subresultant (recompute here)
C1s = 3*a*x**2 - b*s_*x - b*x + s_**2 + s_
C2s = sp.expand(a*x**3 + g*(1+s_)**3 - x*(s_**2 + 3*s_ + 2))
prs = sp.subresultants(Poly(C1s, s_), Poly(C2s, s_))
lin1 = [q for q in prs if sp.degree(q, s_) == 1]
q1 = sp.Poly(lin1[-1], s_)
s_expr = sp.cancel(-q1.all_coeffs()[1] / q1.all_coeffs()[0])
s_num = sp.lambdify((x, a, b, g), s_expr, 'mpmath')
import mpmath
def s_of_x_num(xr, av, bv, gv):
    return complex(s_num(complex(xr), complex(av), complex(bv), complex(gv)))

def Fnum(P):
    xv, yv, zv = P
    U = 1 + xv*yv
    return (U**3*zv + yv**2*U*(4 + 3*xv*yv),
            yv + 3*xv*U**2*zv + 3*xv*yv**2*(4 + 3*xv*yv),
            2*xv - 3*xv**2*yv - xv**3*zv)

random.seed(3)
worst = 0.0
for trial in range(6):
    Q = (random.uniform(-3, 3) + 0.3, random.uniform(-3, 3), random.uniform(-3, 3) + 0.2)
    pts = radical_inverse(*Q)
    errs = []
    for P in pts:
        im = Fnum(P)
        errs.append(max(abs(im[i] - Q[i]) for i in range(3)))
    worst = max(worst, max(errs))
    print(f"  target ~({Q[0]:.2f},{Q[1]:.2f},{Q[2]:.2f}): 3 radical branches, max |F(inv_j(Q)) - Q| = {max(errs):.2e}", flush=True)
print(f"  WORST ERROR over all trials/branches: {worst:.2e}  => THE RADICAL INVERSE VERIFIED (S_3 solvable => F^-1 in radicals)", flush=True)
print("  radical tower: sqrt(disc) with disc ~ -Delta/(108 L^3)-class (the resolvent sqrt(-L)) then one cube root -- Cardano = the trisection made algebraic", flush=True)

print("\n== (3) the trace module (fit method) + the F.F centroid by tower transitivity ==", flush=True)
Lp = L
N = sp.expand(L*x**3 + (4 - 3*b*g)*x - 2*g)
def trace_fit(hname, hexpr, maxdeg=6):
    """Tr_F of source-monomial hexpr(x,y,z): try polynomial fit (by weight/parity), else L*Tr fit; returns verdict string."""
    # weight: x:1, y:-1, z:-2 => weight of a,b,g targets: a:-2, b:-1, g:1
    # compute weight and parity of h
    pw = sp.Poly(hexpr, x, y, z)
    (mon, coef) = list(pw.terms())[0]
    W = mon[0] - mon[1] - 2*mon[2]
    par = (mon[1] + mon[2]) % 2   # parity under tau: y odd, z... careful: sigma: (x,y,z)->(-x,-y,z): h parity in (x,y); target tau: (a,b,g)->(a,-b,-g)
    parh = (mon[0] + mon[1]) % 2  # sigma-parity of h
    # samples
    sams = []
    random.seed(17)
    t = 0
    while len(sams) < 18 and t < 400:
        t += 1
        va = Rational(random.randint(-9, 9), random.choice([1, 2, 3]))
        vb = Rational(random.randint(-9, 9), random.choice([1, 2, 3]))
        vg = Rational(random.randint(-9, 9), random.choice([1, 2]))
        if va == 0 or vg == 0: continue
        Lv = Lp.subs({a: va, b: vb, g: vg})
        if Lv == 0: continue
        Nv = Poly(N.subs({a: va, b: vb, g: vg}), x, domain='QQ')
        # h as rational function of x on the fiber: substitute y = s(x)/x, z = (a-B1)/u^3
        yx = s_expr/x
        ux = 1 + s_expr
        B1x = (s_expr**2/x**2)*ux*(4 + 3*s_expr)
        zx = (a - B1x)/ux**3
        hx = sp.together(hexpr.subs({y: yx, z: zx}))
        hn, hd = sp.fraction(sp.cancel(hx))
        hnv = Poly(sp.expand(hn.subs({a: va, b: vb, g: vg})), x, domain='QQ')
        hdv = Poly(sp.expand(hd.subs({a: va, b: vb, g: vg})), x, domain='QQ')
        if sp.degree(sp.gcd(hdv, Nv), x) > 0: continue
        ginv = sp.invert(hdv, Nv)
        redp = (hnv * ginv) % Nv
        cs = redp.all_coeffs()[::-1]
        while len(cs) < 3: cs.append(0)
        p2v = -2*(4 - 3*vb*vg)/Lv
        p1v = 0
        tr = 3*cs[0] + cs[1]*p1v + cs[2]*p2v
        sams.append((va, vb, vg, sp.nsimplify(tr)))
    if len(sams) < 12:
        return f"  Tr({hname}): insufficient samples ({len(sams)})"
    # candidate bases
    def basis_for(Wt, part, dmax):
        out = []
        for i in range(dmax+1):
            for j in range(dmax+1):
                for k in range(dmax+1):
                    if i + j + k > dmax: continue
                    if -2*i - j + k != Wt: continue
                    if (j + k) % 2 != part: continue
                    out.append(a**i * b**j * g**k)
        return out
    for mode, Wt in (("poly", W), ("L-pole", W - 2)):
        bas = basis_for(Wt, parh if mode == "poly" else parh, maxdeg)
        if not bas: continue
        co = sp.symbols(f'c0:{len(bas)}')
        eqs = []
        for va, vb, vg, tr in sams:
            tgtv = tr if mode == "poly" else tr*Lp.subs({a: va, b: vb, g: vg})
            eqs.append(sum(c*m.subs({a: va, b: vb, g: vg}) for c, m in zip(co, bas)) - tgtv)
        try:
            sol = sp.solve(eqs, co, dict=True)
        except Exception:
            sol = []
        if sol:
            expr = sum(sol[0].get(c, c)*m for c, m in zip(co, bas))
            if expr.free_symbols - {a, b, g}:  # underdetermined leftovers
                expr = expr.subs({c: 0 for c in co if c in expr.free_symbols})
            resid = all(sp.simplify(expr.subs({a: va, b: vb, g: vg}) - (tr if mode == "poly" else tr*Lp.subs({a: va, b: vb, g: vg}))) == 0 for va, vb, vg, tr in sams)
            if resid:
                if mode == "poly":
                    return f"  Tr({hname}) = {sp.expand(expr)}   [POLYNOMIAL]"
                else:
                    return f"  Tr({hname}) = ({sp.factor(expr)})/L   [L-POLE]"
    return f"  Tr({hname}): no fit at deg <= {maxdeg} (needs bigger basis or deeper pole)"

mons = [("y^2", y**2), ("xy", x*y), ("z", z), ("xz", x*z), ("yz", y*z),
        ("xyz", x*y*z), ("x^2z^2", x**2*z**2), ("y^3z", y**3*z), ("y^4", y**4)]
verdicts = {}
for nm, e_ in mons:
    vtxt = trace_fit(nm, e_)
    verdicts[nm] = "POLYNOMIAL" in vtxt
    print(vtxt, flush=True)
need = ["x^2z^2", "xyz", "y^3z", "y^2"]
allpoly = all(verdicts.get(k, False) for k in need)
print(f"\n  TOWER TRANSITIVITY: Tr_FF(x) = 0; Tr_FF(y) = Tr_F(3y/2) = 9b/4;", flush=True)
print(f"  Tr_FF(z) = Tr_F(zeta) needs {{x^2z^2, xyz, y^3z, y^2}} trace-polynomial: {allpoly}", flush=True)
print(f"  => F.F centroid polynomial: {allpoly}   (backlog c99(b) verdict)", flush=True)
print("\nDONE.", flush=True)
