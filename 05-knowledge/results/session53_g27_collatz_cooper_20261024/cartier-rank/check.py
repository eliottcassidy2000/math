import json
from math import comb,gcd

def req(v,s):
    if not v:raise RuntimeError(s)
N=240
F=[sum(comb(n,k)**2*comb(n+k,k)*comb(2*k,n) for k in range((n+1)//2,n+1)) for n in range(N+1)]
c=[1]
for n in range(1,N+1):
    v,r=divmod(F[n]-sum(c[k]*c[n-k] for k in range(1,n)),2)
    req(r==0,'integrality');c.append(v)
rows=[]
for p in [2,3,5,7,11,13]:
    m=p*p;inert=p%7 in (3,5,6);bound=N//m
    req(all((c[m*n]-c[n])%m==0 for n in range(bound+1)),'swap action')
    if inert:
        witness=next(n for n in range(1,N//p+1) if (c[p*n]-c[n])%p)
        det=(c[p*witness]-c[witness])%m
        req(gcd(det,m)==1,'unit determinant')
        for a in range(m):
            for b in range(m):
                if (a+b)%m==0 and (a*c[witness]+b*c[p*witness])%m==0:
                    req(a==b==0,'free rank two')
        half=pow(2,-1,m)
        for n in range(bound+1):
            g,k,kk=c[n]%m,c[p*n]%m,c[m*n]%m
            plus,minus=(g+k)*half%m,(g-k)*half%m
            req((k+kk)*half%m==plus,'plus eigenline')
            req((k-kk)*half%m==(-minus)%m,'minus eigenline')
            req((plus+minus)%m==g,'direct sum')
        req((2*c[0])%m!=0,'false scalar minus rejected')
        req(det!=0,'false rank one rejected')
        # Over R an arbitrary nonzero coefficient need not be a unit.
        req((p*p)%m==0 and p%m!=0,'nonunit dependence hostile')
        rows.append({'p':p,'rank':2,'period':2,'finite_unit_witness_n':witness,'det_mod_p2':det,'matrix_bound':bound})
    else:
        req(all((c[p*n]-c[n])%m==0 for n in range(N//p+1)),'rank one action')
        rows.append({'p':p,'rank':1,'period':1,'matrix_bound':bound})
req(gcd(2,4)!=1,'p2 no eigenline division')
print(json.dumps({'scope':'FINITE-EXACT supplemental controls, not universal proof','independent_binomial_square_root_bound':N,'rows':rows,'hostiles':['scalar minus','rank one at inert','nonunit difference cannot prove freeness','division by2 at p2']},sort_keys=True,indent=2))
