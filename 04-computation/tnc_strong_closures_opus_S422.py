import sympy as sp
u=sp.symbols('u')
def CT(R,N,mm):
    return sp.Poly(sp.expand(R**mm),u).coeff_monomial(u**(N*mm))
print("STRONG CLOSURE 1: common-ray coefficients.  If every nonzero r_k = rho_k e^{i phi}")
print("(SAME argument phi), then every rep of 0 at level m uses EXACTLY m factors, so")
print("CT(m) = e^{i m phi} * (positive sum) != 0.  Hence R is TNC-clear.")
print("  Verify with a complex common phase phi (all coeffs = positive * e^{i phi}):")
phi=sp.pi/5
for coeffs,exps,N in [([1,1,1],[0,3,6],2),([2,3,1],[0,1,5],3),([1,1,1,1],[0,3,4,5],2)]:
    R=sum(c*sp.exp(sp.I*phi)*u**e for c,e in zip(coeffs,exps))
    m0=next((m for m in range(1,10) if sp.simplify(CT(R,N,m))!=0),None)
    val=sp.simplify(CT(R,N,m0))
    print(f"   coeffs {coeffs} e^(i pi/5), exps {exps}, N={N}: CT(m0={m0}) = {val}  (nonzero: {val!=0})")
print("  => the nullcone MISSES the entire common-ray locus (real-codim structure). PROVED.")
print()
print("STRONG CLOSURE 2: the (k-1)-level bound -- do {m0,2m0,...,(k-1)m0} saturate?")
a,b,c,w=sp.symbols('a b c w')
def empty_via_levels(R,N,params,mults,m0):
    polys=[sp.expand(CT(R,N,mm)) for mm in [k*m0 for k in mults]]
    polys=[p for p in polys if p!=0 and not p.is_number]
    prod=1
    for pr in params: prod*=pr
    G=sp.groebner(polys+[1-w*prod],*params,w,order='grevlex')
    return list(G)==[sp.Integer(1)], len(polys)
print("  5-NOMIAL (3 params): levels {m0,2m0,3m0,4m0} -- do 4 = k-1 levels suffice?")
for N,exps in [(2,[3,4,5,6]),(3,[1,2,4,5])]:
    j,l,p,d=exps
    R=1+a*u**j+b*u**l+c*u**p+u**d
    m0=next((m for m in range(1,10) if sp.expand(CT(R,N,m))!=0 and not sp.expand(CT(R,N,m)).is_number),None)
    for K in [3,4,5]:
        emp,ne=empty_via_levels(R,N,[a,b,c],list(range(1,K+1)),m0)
        print(f"   N={N} charges({','.join(str(e-N) for e in [0]+exps)}): m0={m0}, "
              f"levels 1..{K}*m0 ({ne} eqs) -> empty: {emp}")
        if emp: break
