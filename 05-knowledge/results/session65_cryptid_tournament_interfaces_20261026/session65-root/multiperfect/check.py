from math import isqrt
from fractions import Fraction
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)
def factor(n):
    out=[];p=2
    while p*p<=n:
        a=0
        while n%p==0:n//=p;a+=1
        if a:out.append((p,a))
        p+=1
    if n>1:out.append((n,1))
    return out
def sigma(n):
    v=1
    for p,a in factor(n):v*=(p**(a+1)-1)//(p-1)
    return v
def direct(n):
    return sum(d+(n//d if d*d!=n else 0) for d in range(1,isqrt(n)+1) if n%d==0)
def v2(n):
    return (n&-n).bit_length()-1
rows=[]
for r in range(2,25):
    apex=(1<<r)-1;dim=(1<<(r-1))*apex
    sd=sigma(dim);sm=sigma(apex)
    need(sd==apex*sm,'multiplicative arc dimension')
    need(Fraction(sd,dim)==Fraction(sm,1<<(r-1)),'exact abundancy transport')
    need((sd%dim==0)==(v2(sm)>=r-1),'exact multiperfect valuation test')
    prime=len(factor(apex))==1 and factor(apex)[0][1]==1
    need((sd==2*dim)==prime,'Euclid-Euler family control')
    rows.append({'exponent':r,'order':1<<r,'arc_dimension':dim,'apex':apex,'sigma_apex':sm,'binary_valuation':v2(sm),'abundancy':str(Fraction(sd,dim)),'multiperfect':sd%dim==0})
for n in range(1,10001):
    need(sigma(n)==direct(n),'independent divisor enumeration')
    if n%2:
        budget=sum(v2(p+1)+v2(a+1)-1 for p,a in factor(n) if a%2)
        need(v2(sigma(n))==budget,'odd Euler valuation budget')
        need((sigma(n)%2==1)==(isqrt(n)**2==n),'odd sigma parity iff square')
need(sigma(120)==360 and sigma(2016)!=2*2016,'apex prime/multiperfect distinction')
need(sigma(21)==32 and sigma(21)%21!=0,'valuation alone is not a generic odd witness')
print(json.dumps({'status':'PASS','family':rows,'independent_sigma_values':10000,'odd_budget_values':5000,'hostiles':{'120':'three-perfect, not perfect','21':'sigma=32 but21 does not divide32; power-of-two cancellation test only applies to the specified arc family'}},sort_keys=True,indent=2))
