import json, hashlib
from pathlib import Path

def require(x,msg):
    if not x: raise RuntimeError(msg)
N=10000
c=[1,2]
for n in range(1,N):
    num=(26*n*n+13*n+2)*c[n]+3*(3*n-1)*(3*n-2)*c[n-1]
    q,r=divmod(num,(n+1)**2)
    require(r==0,f'integrality {n+1}')
    c.append(q)
primes=[p for p in range(2,101) if all(p%d for d in range(2,int(p**.5)+1))]
rows=[]
for p in primes:
    failures=[n for n in range(1,N//(p*p)+1) if (c[p*p*n]-c[n])%(p*p)]
    require(not failures,f'two-step {p}: {failures}')
    bad3=next((n for n in range(1,N//(p*p)+1) if (c[p*p*n]-c[n])%(p**3)),None)
    rows.append({'p':p,'chi':0 if p==7 else (1 if p%7 in (1,2,4) else -1),'positive_n_checked':N//(p*p),'first_mod_p3_failure_n_in_range':bad3})
# Independent exact q-series check of f*=theta^3(1-3t/4), multiplied by 8.
M=80
chi=lambda d: 0 if d%7==0 else (1 if d%7 in (1,2,4) else -1)
A1=[0]+[sum(chi(d)*d*d for d in range(1,n+1) if n%d==0) for n in range(1,M+1)]
A2=[0]+[sum(chi(n//d)*d*d for d in range(1,n+1) if n%d==0) for n in range(1,M+1)]
th=[1]+[2*sum(chi(d) for d in range(1,n+1) if n%d==0) for n in range(1,M+1)]
def mul(a,b):
    return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(M+1)]
th3=mul(mul(th,th),th)
eta=[1]+[0]*M
for k in (1,7):
    for n in range(k,M+1,k):
        for _ in range(3):
            for j in range(M,n-1,-1): eta[j]-=eta[j-n]
E=[0]+eta[:-1]
for n in range(M+1):
    require((8 if n==0 else -7*A1[n]+49*A2[n])==8*th3[n]-6*E[n],f'fstar identity q^{n}')
# Divisor-sum U_p recurrences are exact, including the p^2 correction.
for p in primes:
    if p in (2,7): continue
    for n in range(1,M//p+1):
        prev1=A1[n//p] if n%p==0 else 0
        prev2=A2[n//p] if n%p==0 else 0
        require(A1[p*n]==(1+chi(p)*p*p)*A1[n]-chi(p)*p*p*prev1,'A1 Hecke')
        require(A2[p*n]==(chi(p)+p*p)*A2[n]-chi(p)*p*p*prev2,'A2 Hecke')
require((c[9]-c[1])%9==0 and (c[9]-c[1])%27!=0,'sharpness p3')
report={'status':'PASS','coefficient_bound':N,'identity_q_bound':M,'two_step_tests':rows,'p3_sharpness':{'c9':c[9],'c1':c[1],'difference_mod27':(c[9]-c[1])%27},'checks_survive_optimization':True}
print(json.dumps(report,indent=2))
