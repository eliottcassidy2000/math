"""Independent all-A2 carrier audit and Hermite signs, no producer imports."""
from math import comb,gcd
import sympy as S
import json,hashlib,sys
sys.stdout.reconfigure(newline='\n')
t,u=S.symbols('t u');GATES=0
def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)
def choose(n,k):return comb(n,k) if 0<=k<=n else 0
def coefficients(B,x,y):
    return {j:comb(x+y-2*j,x+(B-2)*j) for j in range(-(x//(B-2)),y//B+1)}
def square(b):
    return {k:sum(v*b.get(k-j,0) for j,v in b.items()) for k in range(2*min(b),2*max(b)+1)}
def hermite_negative(p,response):
    n=p.degree();basis=[]
    for j in range(n):
        r=S.Poly(S.rem(t**(j+1),p.as_expr(),t),t)
        basis.append([r.nth(i) for i in range(n)])
    T=S.Matrix(basis).T
    inv=S.invert(t,p.as_expr(),t)
    rr=S.Poly(S.rem(S.expand(t*response)*inv,p.as_expr(),t),t)
    op=sum((rr.nth(i)*T**i for i in range(n)),S.zeros(n))
    traces=[S.trace(op*T**k) for k in range(2*n-1)]
    H=S.Matrix(n,n,lambda i,j:traces[i+j])
    for k in range(1,n+1):need((-1)**k*H[:k,:k].det()>0,'independent negative-definite trace form')
def maps():
    records=[]
    for B in (3,5,8,9):
      for h in (1,2,3):
       for x in (0,3):
        for z in (0,1):
            r=1;y=B*h+r;m=x+y+z;L=x//(B-2);q=L+h;k=2*h+z
            b=coefficients(B,x,y)
            G=S.Poly(sum(v*t**(j+L) for j,v in b.items()),t)
            p=S.Poly(sum(v*choose(m,z+2*j)*t**j for j,v in b.items() if j>=0),t)
            w=sum(v*choose(2*m,2*z+2*j)*t**j for j,v in square(b).items() if choose(2*m,2*z+2*j))
            need(G.degree()==q and G.nth(0)>0 and G.nth(q)>0,'complete beta endpoints')
            need(G.count_roots(-S.oo,0)==q,'full beta root reality at fresh parameters')
            need(p.degree()==h and p.count_roots(-S.oo,0)==h,'complete original cancellation roots')
            rho=-S.Rational(3,2)
            # The carrier is built by substituting in the entire polynomial,
            # and exact ordinary multiplication computes its square.
            H=S.Poly(S.expand((1+u)**m*sum(v*rho**(j+L)*u**(2*(q-j-L)) for j,v in b.items())),u)
            need(H.degree()==m+2*q and H.nth(0)!=0 and 0<k<H.degree(),'actual interior and endpoints')
            need(H.nth(k)==rho**L*p.eval(rho),'original first coefficient map')
            need((H*H).nth(2*k)==rho**(2*L)*w.subs(t,rho),'complete alpha square map')
            hermite_negative(p,w)
            records.append([B,h,x,r,z])
    print('Fresh general A2 rows',len(records),'cover B3,5,8,9; both residues; x0 included')
    return records
def coupled():
    records=[]
    for h,x in ((1,2),(2,3),(3,2),(4,6),(5,2),(6,6),(7,2),(8,6)):
        q=x+h;g=x+3*h+1
        fs=[]
        for n in (3*q,3*q-1,3*q-2):
            fs.append(S.Poly(sum(comb(n-2*j,j)*t**j for j in range(n//3+1)),t))
        KB=sum(v*t**j*u**(2*q-2*j) for (j,),v in fs[0].terms())
        KC=sum(v*t**j*u**(2*q-2-2*j) for (j,),v in fs[1].terms())
        KD=sum(v*t**j*u**(2*q-2-2*j) for (j,),v in fs[2].terms())
        need(S.expand(S.diff(KB,u)-S.Rational(2,3)*u*((q+2)*KC+u*S.diff(KC,u)))==0,'first coupled Euler identity')
        need(S.expand((q+1)*KD+u*S.diff(KD,u)-2*KC-S.Rational(3,2)*u*S.diff(KC,u))==0,'second coupled Euler identity')
        # Direct composition polynomial product, then binomial convolution,
        # keeps the mixed coefficient before converting to the Laurent row.
        CD=fs[1]*fs[2]
        mixed=sum(v*choose(2*g,4*h-4*q+4+2*j)*t**j for (j,),v in CD.terms())
        skip=S.expand(2*t**(1-2*x)*mixed)
        b=coefficients(3,x,3*h)
        p=S.Poly(sum(v*choose(g,1+2*j)*t**j for j,v in b.items() if j>=0),t)
        need(p.count_roots(-S.oo,0)==h,'complete fresh coupled roots')
        hermite_negative(p,skip)
        records.append([h,x,g,gcd(g,6*h+3)])
    print('Fresh mixed-product rows',json.dumps(records),'all rootwise skip signs negative, FINITE-EXACT only')
    return records
if __name__=='__main__':
    records=[maps(),coupled()]
    need(S.Poly(u*u+1,u).count_roots(-S.oo,S.oo)==0,'positive phase cannot preserve real roots')
    need(S.Poly(u**3-1,u).count_roots(-S.oo,S.oo)==1,'higher pullback exponent hostile')
    print('RECORD_SHA256',hashlib.sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)
