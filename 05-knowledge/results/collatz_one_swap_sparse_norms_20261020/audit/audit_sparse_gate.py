from math import gcd
from collections import Counter
import sympy as s

def ceildiv(a,b): return -(-a//b)
def carry(w):
    b=0
    for i,e in enumerate(w):
        if e:b=3*b+2**i
    return b
z=s.Symbol('z')
checks=Counter()
maxratio=(0,None)
for p in range(2,24):
    norms={}
    for a in range(p):
        for eps in (-1,1):
            Q=z**(a+1)-z**a+1 if eps==1 else z**(p-a)-z+1
            N=abs(int(s.resultant(z**p-2,Q,z)))
            M=N
            for l in (2,3):
                while M and M%l==0:M//=l
            assert 0<M<=3**p-2,(p,a,eps,N,M)
            norms[a,eps]=M
            checks['norms']+=1
    for r in range(1,p):
        if gcd(p,r)>1:continue
        for k in range(1,16):
            q=p*k+r
            c=[((i+1)*p)//q-(i*p)//q for i in range(q)]
            positions=[i for i,e in enumerate(c) if e]
            D=2**q-3**p
            for i in range(q):
                j=(i+1)%q
                if c[i]==c[j]:continue
                w=c.copy(); w[i],w[j]=w[j],w[i]
                moved=i if c[i] else j
                eps=1 if c[i] else -1
                jone=positions.index(moved)+1
                t=jone-1
                eta=lambda m:ceildiv((t+m)*r,p)-ceildiv(t*r,p)
                a=(-(t+1)*r)%p
                h=max(0,-eta(1)-eps)
                P=sum(3**(p-1-m)*2**(h+eta(m)+eps*(m==1)+k*m) for m in range(p))
                g=gcd(abs(D),carry(w))
                assert gcd(abs(D),P)==g,('P',p,r,k,i,eps)
                M=norms[a,eps]
                assert M%g==0,('divis',p,r,k,i,g,M)
                assert gcd(gcd(abs(D),P),M)==g
                if eps==-1 and a==0:
                    assert g==gcd(3**p-2,2**k*3**(r-1)-1),('endminus',p,r,k,i,g)
                    checks['endpoint_minus']+=1
                if eps==1 and a==p-1:
                    assert g==gcd(3**p-2**(p-1),2**(k+r)-3**(r+1)),('endplus',p,r,k,i,g)
                    checks['endpoint_plus']+=1
                checks['swaps']+=1
                if k==1:checks['k1_swaps']+=1
print(dict(checks))
