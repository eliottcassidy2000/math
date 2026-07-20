"""
opus-2026-07-19-S414 part B: DEFORMATION THEORY of the counterexample (the
robust replacement for the full coefficient solve).

Equivariant ansatz G = (zP + y^2 Q, yR + xzS, xT), P..T in C[s,w] as in part A
(15 coefficients).  Constraint: det J(G) == const.  At the known point F:
  (i)  verify det = -2 and all non-constant coefficient equations vanish;
  (ii) LINEARIZE: the Jacobian of the constraint system at F; kernel dimension
       = infinitesimal moduli (including trivial directions);
  (iii) identify the TRIVIAL directions (reparametrizations): weighted torus
       rescalings of source/target that preserve the ansatz shape — count them
       and subtract;
  (iv) attempt to INTEGRATE one nontrivial kernel direction: line-search a
       finite deformation with det still constant; if found, re-derive its
       exotic orbit ({T=0} cap {sR+wS=0}) => a NEW explicit counterexample.
"""
import sympy as sp

x, y, z, s, w, t = sp.symbols('x y z s w t')
a = sp.symbols('a0:4'); b = sp.symbols('b0:3'); c = sp.symbols('c0:3')
d = sp.symbols('d0:3'); e = sp.symbols('e0:2')
coeffs = [*a, *b, *c, *d, *e]
sub_sw = {s: x*y, w: x**2*z}
Pg = a[0]+a[1]*s+a[2]*s**2+a[3]*s**3
Qg = b[0]+b[1]*s+b[2]*s**2
Rg = c[0]+c[1]*s+c[2]*s**2
Sg = d[0]+d[1]*s+d[2]*s**2
Tg = e[0]+e[1]*s-w
G = sp.Matrix([z*Pg.subs(sub_sw)+y**2*Qg.subs(sub_sw),
               y*Rg.subs(sub_sw)+x*z*Sg.subs(sub_sw),
               x*Tg.subs(sub_sw)])
detG = sp.expand(G.jacobian([x, y, z]).det())
poly = sp.Poly(detG, x, y, z)
eqs, const = [], sp.Integer(0)
for mono, coeff in poly.terms():
    if mono == (0, 0, 0): const = coeff
    else: eqs.append(sp.expand(coeff))
eqs = sorted(set(eqs), key=sp.default_sort_key)
known = {a[0]:1,a[1]:3,a[2]:3,a[3]:1, b[0]:4,b[1]:7,b[2]:3,
         c[0]:1,c[1]:12,c[2]:9, d[0]:3,d[1]:6,d[2]:3, e[0]:2,e[1]:-3}
print(f"(i) #distinct constraint equations: {len(eqs)}; "
      f"known point satisfies all: {all(q.subs(known) == 0 for q in eqs)}; "
      f"det at known point: {const.subs(known)}")

Jc = sp.Matrix([[sp.diff(q, v) for v in coeffs] for q in eqs]).subs(known)
ns = Jc.nullspace()
print(f"(ii) linearized constraint rank: {Jc.rank()} of {len(coeffs)} "
      f"coefficients -> tangent/kernel dimension: {len(ns)}")

print("(iii) kernel directions (as coefficient vectors):")
for i, v in enumerate(ns):
    print(f"   dir{i}: {list(v.T)}")

# (iv) integrate: try known + t*dir for each kernel direction; find t-family
print("(iv) integration attempts:")
for i, v in enumerate(ns):
    trial = {coeffs[j]: known[coeffs[j]] + t*v[j] for j in range(len(coeffs))}
    residuals = [sp.expand(q.subs(trial)) for q in eqs]
    if all(r == 0 for r in residuals):
        dt = sp.expand(const.subs(trial))
        print(f"   dir{i}: EXACT 1-parameter family! det(t) = {dt}")
        Rt = Rg.subs(trial); St = Sg.subs(trial); Tt = Tg.subs(trial)
        line = sp.expand((s*Rt + w*St).subs(w, Tt.subs(w, 0)))
        print(f"     exotic-orbit equation on T=0: {sp.collect(line, s)} = 0")
        for tv in (1, -1, 2):
            inst = {k: sp.nsimplify(val.subs(t, tv)) for k, val in trial.items()}
            roots = sp.solve(line.subs(t, tv), s)
            print(f"     t={tv}: coeffs {[inst[k] for k in coeffs]}; "
                  f"det = {const.subs(inst)}; s-roots {roots}")
    else:
        nz = [r for r in residuals if r != 0]
        deg1 = all(sp.degree(sp.Poly(r, t)) >= 2 or r == 0 for r in residuals)
        print(f"   dir{i}: obstructed at order >= 2: {deg1} "
              f"({len(nz)} nonzero residuals)")
