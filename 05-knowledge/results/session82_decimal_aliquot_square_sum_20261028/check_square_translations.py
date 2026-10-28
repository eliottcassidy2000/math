#!/usr/bin/env python3
"""Exact translation classification; stdlib only, no assert-dependent checks."""
import json
from math import isqrt
from itertools import combinations

def check(p,msg):
    if not p: raise RuntimeError(msg)
def square(x):return x>=0 and isqrt(x)**2==x
def factor(x):
    out={};p=2
    while p*p<=x:
        while x%p==0:out[p]=out.get(p,0)+1;x//=p
        p=3 if p==2 else p+2
    if x>1:out[x]=out.get(x,0)+1
    return out

def aliquot(x):
    f=factor(x);sigma=1
    for p,e in f.items():sigma*=(p**(e+1)-1)//(p-1)
    return sigma-x

def triple(n):return {'a':n*n+1,'b':n+1,'c':n*n+n}

def candidates(n):
    values=triple(n);out=[]
    for center,m in values.items():
        (kx,x),(ky,y)=sorted(((k,v) for k,v in values.items() if k!=center),key=lambda z:z[1])
        d=y-x
        for u in range(1,isqrt(d)+1):
            if d%u:continue
            v=d//u
            if (u+v)%2:continue
            r=(v-u)//2;s=(v+u)//2
            if (r*r-m-x)%2:continue
            t=(r*r-m-x)//2
            if min(values.values())+t<=0:continue
            path=[x+t,m+t,y+t]
            check(path[0]+path[1]==r*r and path[1]+path[2]==s*s,'edge identity')
            out.append(dict(center=center,offset=t,path=path,roots=[r,s],factor_pair=[u,v]))
    return sorted(out,key=lambda z:(z['offset'],z['center']))

def main():
    # Independent direct square-root enumeration, exhaustive finite upper bound.
    for n in range(2,101):
        expected=set()
        vals=triple(n)
        for mid,m in vals.items():
            x,y=sorted(v for k,v in vals.items() if k!=mid);d=y-x
            for r in range((d-1)//2+1):
                if square(r*r+d) and (r*r-m-x)%2==0:
                    t=(r*r-m-x)//2
                    if min(vals.values())+t>0:expected.add((mid,t))
        actual={(r['center'],r['offset']) for r in candidates(n)}
        check(actual==expected,'independent enumeration')
        check(all(mid!='b' for mid,t in actual),'forbidden center')
    n=1000; vals=triple(n);rows=candidates(n)
    for row in rows:
        row['aliquot']=[aliquot(x) for x in row['path']]
        row['factors']=[factor(x) for x in row['path']]
    check(len(rows)==24,'complete count')
    check(all(row['offset']>0 for row in rows),'no nonpositive offsets')
    for row in rows:
        p,h,q=row['path'];disc=1+4*(q-p)
        check(square(disc),'recovery discriminant')
        recovered_n=(1+isqrt(disc))//2;recovered_t=p-recovered_n-1
        check(recovered_n==1000 and recovered_t==row['offset'],'source recovery')
        check(h==recovered_n**2+recovered_n+recovered_t,'center recovery')
    check(rows[0]['offset']==177612,'least offset')
    check(rows[0]['aliquot']==[1,981406,1],'prime endpoint fiber')
    check([r['offset'] for r in rows if r['aliquot'][0]==r['aliquot'][2]==1]==[177612], 'unique prime endpoint offset')
    # Explicit unbounded family n=8k, k>=1, with center c.
    for n in range(8,8001,8):
        a,b,c=triple(n).values();d=n*(n-1);r=d//4-1;s=d//4+1
        t=(r*r-b-c)//2
        check((r*r-b-c)%2==0 and t>=0,'family positivity')
        check(b+c+2*t==r*r and a+c+2*t==s*s,'family edges')
    original_edges=[(x,y) for x,y in combinations(vals.values(),2) if square(x+y)]
    print(json.dumps(dict(status='FINITE-EXACT with proved enumeration bound',original=vals,
      original_edges=original_edges,original_aliquot={k:aliquot(v) for k,v in vals.items()},
      translation_count=len(rows),prime_endpoint_offsets=[r['offset'] for r in rows if r['aliquot'][0]==r['aliquot'][2]==1],
      translations=rows,independent_crosscheck='all n=2..100',family_check='n=8,16,...,8000'),indent=2))
if __name__=='__main__':main()
