"""Target-specific progression check of the older forty-template fixed bank.

This is an exploratory sidecar, not part of the main central-word audit.
The constant matrix file is copied from the reviewed session198 family.
"""
from functools import lru_cache
from pathlib import Path
from math import gcd,isqrt
import json

F=json.loads((Path(__file__).parent/'forms.json').read_text())
@lru_cache(None)
def g(n):
    if n==1:return 1
    if n%2==0:return g(n//2)
    h=(n+1)//2
    return g(h)+h-h%2

def check(n,expected):
    tested=0
    for u in range(2,isqrt((n-1)//3)+1):
        matches=[]
        s=(g(u),u,u%2,1)
        for name,A in F.items():
            alpha,beta,gamma,delta=[sum(s[i]*A[i][j] for i in range(4)) for j in range(4)]
            if alpha<=0 or beta<=0:raise RuntimeError('positive inverse coefficients')
            bound=(n-1)//(3*u)
            for q in (0,1):
                lo=(u-q+1)//2;hi=(bound-q)//2;h=gcd(2*beta,alpha)
                rhs=n-delta-(beta+gamma)*q
                if rhs%h:continue
                mod=alpha//h
                r=0 if mod==1 else rhs//h*pow(2*beta//h,-1,mod)%mod
                for t in range(lo+(r-lo)%mod,hi+1,mod):
                    v=2*t+q;tested+=1
                    if alpha*g(v)+beta*v+gamma*q+delta==n:
                        matches.append({'form':name,'left_index':u,'right_index':v,
                                        'g_left':g(u),'g_right':g(v),
                                        'row':[alpha,beta,gamma,delta]})
        if matches:
            if u!=expected:raise RuntimeError('minimum smaller index mismatch')
            return {'target':n,'minimum_smaller_index':u,'matches':matches,'progression_tests':tested}
    raise RuntimeError('expected representation absent')

if __name__=='__main__':
    print(json.dumps({'status':'AUTHOR CHECK PASS','scope':'two fixed-target minima only',
                      'targets':[check(1348681,41),check(6888203,86)]},indent=2,sort_keys=True))
