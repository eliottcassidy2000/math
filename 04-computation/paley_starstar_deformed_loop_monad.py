#!/usr/bin/env python3
"""
paley_starstar_deformed_loop_monad.py
monad-explorer-2026-06-07 (deep-research, 11th session)

Search for a DEFORMED loop equation for V(s,y) = U(x,y), s=x/(1-x):
   V(s,y) = sum_{m>=1} Q_m(s) y^m
At y=-1: F=1+V(s,-1) solves xF^2+F-1=0, i.e. (in s,V):  s V^2 + (1+3s) V + s = 0.

We look for a(s,y) V^2 + b(s,y) V + c(s,y) = 0 with a,b,c polynomials of bounded
degree, matching the known columns Q_1..Q_4 (EXACT, all s) and partial Q_5,Q_6.
If found, analyzing at s=-1 should yield BOTH handoffs:
   V(-1,y) = -y ;  V_s(-1,y) = y/((1-y)(1-2y)).
"""
import sympy as sp

s, y = sp.symbols('s y')

# exact Q_m (from P_m via Q_m = s^m (1+s)^{m-1} P_m(s/(1+s)) )
Q = {}
Q[1] = s
Q[2] = 3*s**2 + 3*s**3
Q[3] = 13*s**3 + 33*s**4 + 20*s**5
Q[4] = 69*s**4 + 304*s**5 + 416*s**6 + 181*s**7
# partial Q5 (with c2,c3,c4): from run:
c2,c3,c4 = sp.symbols('c2 c3 c4')
Q[5] = 421*s**5 + 2740*s**6 + (5694+c2)*s**7 + (4852+2*c2+c3)*s**8 + (1477+c2+c3+c4)*s**9

NY = 5  # use y up to y^4 fully (Q1..Q4 exact); y^5 partial as a check
Vfull = sum(Q[m]*y**m for m in range(1, 6))

def fit_quadratic(degA_s, degA_y, label):
    """Try a(s,y)V^2+b(s,y)V+c(s,y)=0 with each of a,b,c poly deg<=degA_s in s, <=degA_y in y.
       Solve linear system from y^1..y^4 coefficients (s-exact) = 0. c2,c3,c4 left symbolic
       but we DON'T use Q5 here (only Q1..Q4), so result is independent of them."""
    V = sum(Q[m]*y**m for m in range(1, 5))  # exact part only
    # unknown coeff symbols
    A = {}; B = {}; C = {}
    syms = []
    for (D, name) in [(A,'a'),(B,'b'),(C,'c')]:
        for i in range(degA_s+1):
            for j in range(degA_y+1):
                u = sp.Symbol(f'{name}_{i}_{j}')
                D[(i,j)] = u; syms.append(u)
    a = sum(A[(i,j)]*s**i*y**j for i in range(degA_s+1) for j in range(degA_y+1))
    b = sum(B[(i,j)]*s**i*y**j for i in range(degA_s+1) for j in range(degA_y+1))
    cc= sum(C[(i,j)]*s**i*y**j for i in range(degA_s+1) for j in range(degA_y+1))
    expr = sp.expand(a*V**2 + b*V + cc)
    # collect powers of y up to NY-1; each coeff is a poly in s -> all its s-coeffs = 0
    eqs = []
    p = sp.Poly(expr, y)
    for j in range(0, 5):  # y^0..y^4 (V exact to y^4, V^2 exact to y^4 since V=O(y))
        cj = p.coeff_monomial(y**j) if j>0 else expr.coeff(y,0)
        cj = sp.expand(cj)
        if cj == 0:
            continue
        ps = sp.Poly(cj, s)
        for co in ps.all_coeffs():
            eqs.append(co)
    sol = sp.solve(eqs, syms, dict=True)
    print(f"--- {label}: degA_s={degA_s},degA_y={degA_y}, #unknowns={len(syms)}, #eqs={len(eqs)} ---")
    if not sol:
        print("   NO solution"); return None
    # find a nontrivial solution
    sol0 = sol[0]
    nontrivial = any(sol0.get(u,u)!=0 for u in syms)
    # substitute, normalize
    aa = sp.expand(a.subs(sol0)); bb = sp.expand(b.subs(sol0)); ccc = sp.expand(cc.subs(sol0))
    free = (aa.free_symbols|bb.free_symbols|ccc.free_symbols) - {s,y}
    if free:
        # set first free param to 1, rest 0
        subm = {f:0 for f in free}
        fl = sorted(free, key=str)[0]; subm[fl]=1
        aa=sp.expand(aa.subs(subm)); bb=sp.expand(bb.subs(subm)); ccc=sp.expand(ccc.subs(subm))
    if aa==0 and bb==0 and ccc==0:
        print("   only trivial"); return None
    print("   a =", aa); print("   b =", bb); print("   c =", ccc)
    # verify at y=-1 matches loop eqn s V^2+(1+3s)V+s (up to scalar)
    a1=aa.subs(y,-1); b1=bb.subs(y,-1); c1=ccc.subs(y,-1)
    print("   at y=-1: a,b,c =", sp.simplify(a1), sp.simplify(b1), sp.simplify(c1),
          " (loop wants prop to s,1+3s,s)")
    # cross-check on Q5 (partial): plug Vfull mod y^6 and see y^5 coeff constraints on c2,c3,c4
    chk = sp.expand(aa*Vfull**2 + bb*Vfull + ccc)
    pc = sp.Poly(chk, y).coeff_monomial(y**5)
    print("   [y^5] residual (should be 0; constrains c2,c3,c4):", sp.expand(pc))
    return (aa,bb,ccc)

for (ds,dy) in [(2,2),(3,2),(3,3),(4,3),(4,4),(5,4)]:
    r = fit_quadratic(ds,dy,f"quad s<= {ds}, y<= {dy}")
    if r:
        break
