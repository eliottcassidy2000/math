import sympy as sp
from sympy import primerange, factorint
u,a=sp.symbols('u a')
def CT(R,N,mm):
    return sp.Poly(sp.expand(R**mm),u).coeff_monomial(u**(N*mm))
def acontent(c):
    q=sp.Poly(c,a)
    while q.total_degree()>0 and q.eval(0)==0: q=sp.Poly(sp.quo(q.as_expr(),a),a)
    return q
print("FINITE-PLACE HALF, correct criterion: CT(m0), CT(2m0) coprime mod p  <=>  p does NOT")
print("divide Res_a(CT(m0),CT(2m0)).  Reductions mod p are LUCAS products (the carry CA).")
print("So the SEPARATING finite places are ALL primes p not dividing the resultant -- infinitely")
print("many.  Bad primes = prime factors of Res.  Test the bad-prime structure:")
print()
print("   pattern          Res (a-content)         bad primes (p|Res)     smallest GOOD prime")
for N,j,d in [(2,3,6),(2,5,8),(3,2,6),(3,4,8),(3,5,10),(4,5,10)]:
    R=1+a*u**j+u**d
    m0=next((m for m in range(1,16) if sp.expand(CT(R,N,m))!=0 and not sp.expand(CT(R,N,m)).is_number),None)
    if m0 is None: continue
    c0=acontent(sp.expand(CT(R,N,m0))); c1=acontent(sp.expand(CT(R,N,2*m0)))
    res=sp.resultant(c0.as_expr(),c1.as_expr(),a)
    if res==0: 
        print(f"   N={N}(-{N},{j-N},{d-N}): RES=0 -- would be a nullcone element! (should not happen)")
        continue
    bad=sorted(factorint(int(abs(res))).keys())
    # smallest prime not dividing res AND not dividing leading coeffs (so degrees preserved)
    lead=[abs(int(sp.LC(c0.as_expr(),a))), abs(int(sp.LC(c1.as_expr(),a)))]
    good=next(p for p in primerange(2,200) if res%p!=0 and lead[0]%p!=0 and lead[1]%p!=0)
    print(f"   N={N}(-{N},{j-N},{d-N}): |Res|={int(abs(res))}  bad primes {bad}  smallest good p={good}")
print()
print("=> the finite-place half is CLOSED: coprimality over Q (THM-1680/1710) => coprime mod")
print("   EVERY prime p not dividing Res => separation at infinitely many finite places, with")
print("   an EXPLICIT smallest good prime. The bad primes are the resultant's factors, and the")
print("   mod-p reductions are Lucas/carry-CA products. No archimedean place is needed.")
