import sympy as sp
from sympy import factorint
u,a=sp.symbols('u a')
def CT(R,N,mm):
    return sp.Poly(sp.expand(R**mm),u).coeff_monomial(u**(N*mm))
def acontent(c):
    q=sp.Poly(c,a)
    while q.total_degree()>0 and q.eval(0)==0: q=sp.Poly(sp.quo(q.as_expr(),a),a)
    return q
print("ANGULAR UNIFORM (HYP-8535): does the largest bad prime of Res(CT(m0),CT(2m0)) stay")
print("BOUNDED as the SPAN grows (fixed charge-count 3)?  If yes -> a uniform good prime.")
print("The tunable pattern {-N, 1, M} at growing M (N=2, middle charge 1):")
print()
print("   pattern (N=2)     m0   |Res|-largest-prime   smallest good prime")
for d in [4,6,8,10,12,14]:   # R=1+a u^{d-1}... use charges {-2,1,d-2}: j so charge 1 -> j=3
    # {-2,1,M}: j-N=1 => j=3; d-N=M => d=M+2, span M+... use R=1+a u^3+u^{M+2}? keep charge 1 fixed
    R=1+a*u**3+u**d   # charges -2,1,d-2
    N=2
    m0=next((m for m in range(1,20) if sp.expand(CT(R,N,m))!=0 and not sp.expand(CT(R,N,m)).is_number),None)
    if m0 is None: 
        print(f"   d={d}: no tunable level"); continue
    c0=acontent(sp.expand(CT(R,N,m0))); c1=acontent(sp.expand(CT(R,N,2*m0)))
    if c0.total_degree()<1: 
        print(f"   d={d} (charges -2,1,{d-2}): m0={m0} CT(m0) monomial (unique-min) -> TNC, no Res needed"); continue
    res=sp.resultant(c0.as_expr(),c1.as_expr(),a)
    if res==0: print(f"   d={d}: RES=0!"); continue
    fac=factorint(int(abs(res))); maxp=max(fac); 
    from sympy import primerange
    good=next(p for p in primerange(2,500) if int(abs(res))%p!=0)
    print(f"   d={d} (charges -2,1,{d-2}): m0={m0}  max bad prime={maxp}  good p={good}")
print()
print("If the largest bad prime GROWS with span, no single uniform prime -- but a good prime")
print("<= (next prime after max multinomial) always exists by Kummer height. Check growth rate.")
