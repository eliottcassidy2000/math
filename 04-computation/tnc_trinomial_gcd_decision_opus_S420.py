import sympy as sp
u,a=sp.symbols('u a')
def CT(Rexpr,N,mm):
    return sp.Poly(sp.expand(Rexpr**mm),u).coeff_monomial(u**(N*mm))
print("THE FINISHING COMPUTATION.  Gauge-fix r0=r_d=1, free middle coeff = a.")
print("Each CT(m) is a POLYNOMIAL in a.  A trinomial nullcone violator = a common NONZERO")
print("root of ALL CT(m).  If gcd of a FINITE set of them has no nonzero root, NONE exists.")
print()
print("   trinomial (N; charges)        nonzero CT(m) polys in a        gcd nonzero-root?")
allclosed=True
for N,j,d in [(2,3,6),(2,1,4),(2,1,5),(2,3,5),(3,1,5),(3,2,7),(2,1,7),(3,1,4),(3,4,5),(2,5,6)]:
    R=1+a*u**j+u**d
    polys=[]
    for m in range(1,16):
        c=sp.Poly(sp.expand(CT(R,N,m)),a)
        if c.total_degree()>0 or (c.as_expr()!=0 and c.as_expr().free_symbols):
            if c.as_expr()!=0: polys.append(c.as_expr())
        if len(polys)>=4: break
    if not polys:
        print(f"   N={N} charges(-{N},{j-N},{d-N}): all CT constant in a -- degenerate gauge, skip")
        continue
    g=polys[0]
    for p in polys[1:]: g=sp.gcd(sp.Poly(g,a),sp.Poly(p,a)).as_expr()
    # nonzero roots of the gcd
    if g.free_symbols:
        rts=[r for r in sp.solve(g,a) if r!=0]
    else:
        rts=[]   # gcd is a nonzero constant -> no common root at all
    closed=(len(rts)==0)
    allclosed=allclosed and closed
    print(f"   N={N} (-{N},{j-N},{d-N}):   {len(polys)} polys, gcd={sp.factor(g)}   "
          f"nonzero common root: {rts if rts else 'NONE -> TNC'}")
print()
print(f"ALL trinomial patterns closed by finite gcd (no nonzero common root): {allclosed}")
print()
print("DETAIL on the witness {-2,1,4}: show CT(3) and CT(6) are COPRIME in a.")
R=1+a*u**3+u**6
c3=sp.factor(CT(R,2,3)); c6=sp.factor(CT(R,2,6))
g=sp.gcd(sp.Poly(CT(R,2,3),a),sp.Poly(CT(R,2,6),a)).as_expr()
print(f"   CT(3) = {c3}")
print(f"   CT(6) = {c6}")
print(f"   gcd(CT(3),CT(6)) in a = {sp.factor(g)}   -> {'COPRIME (only a=0)' if sp.Poly(g,a).total_degree()==0 or all(r==0 for r in sp.solve(g,a)) else 'shares a root'}")
print("   CT(3)=0 at a=+-i; is CT(6)=0 there?", [sp.simplify(CT(R,2,6).subs(a,r)) for r in [sp.I,-sp.I]])
