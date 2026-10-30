#!/usr/bin/env python3
"""Ramified metallic valuation layers; exact deterministic independent tests."""
import json
from pathlib import Path

def check(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def mul(x,y,a,n):
    b,c=x;d,e=y
    return ((b*d+c*e)%n,(b*e+c*d+a*c*e)%n)

def trace(x,a,n):return (2*x[0]+a*x[1])%n

def power(x,m,a,n):
    y=(1,0)
    while m:
        if m%2:y=mul(y,x,a,n)
        x=mul(x,x,a,n);m//=2
    return y

rows=[]
for p,a in [(5,1),(13,3),(17,8),(29,5),(41,18)]:
    D=a*a+4
    check(D%p==0 and D%(p*p)!=0,(p,a,D))
    for k in [1,2,3]:
        n=p**k; delta=(-a,2);g=power((0,1),4,a,n)
        cases=0
        for val in range(2*k):
            # Multiplication by rational units leaves the prescribed delta depth.
            for unit in [1,2,p-1,p+1]:
                x=mul((unit,0),power(delta,val,a,n),a,n)
                r=2*k-val;full=p**(r//2);obs=p**((r-1)//2)
                check(mul(power(g,full,a,n),x,a,n)==x,('full upper',p,a,k,val))
                if full>1:
                    check(mul(power(g,full//p,a,n),x,a,n)!=x,('full minimal',p,a,k,val))
                def same_pair(d):
                    y=mul(power(g,d,a,n),x,a,n)
                    return trace(x,a,n)==trace(y,a,n) and trace(mul(g,x,a,n),a,n)==trace(mul(g,y,a,n),a,n)
                check(same_pair(obs),('trace upper',p,a,k,val))
                if obs>1:check(not same_pair(obs//p),('trace minimal',p,a,k,val))
                cases+=1
        rows.append(dict(prime=p,metallic_parameter=a,k=k,cases=cases))
result={'status':'PASS','cases':sum(r['cases'] for r in rows),'rows':rows}
print(json.dumps(result,indent=2))
Path(__file__).with_name('metallic_results.json').write_text(json.dumps(result,indent=2)+'\n')
