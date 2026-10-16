#!/usr/bin/env python3
"""Full-ring orbit walks versus valuation/Burnside counts; elementary prime controls."""
import json
from collections import Counter
from pathlib import Path
import sys

def need(b,s):
    if not b: raise RuntimeError(s)

def mul(x,y,a,M):
    return ((x[0]*y[0]+x[1]*y[1])%M,
            (x[0]*y[1]+x[1]*y[0]+a*x[1]*y[1])%M)

def power(x,n,a,M):
    y=(1,0)
    while n:
        if n%2: y=mul(x,y,a,M)
        x=mul(x,x,a,M);n//=2
    return y

def order(x,a,M,cap):
    y=(1,0)
    for j in range(1,cap+1):
        y=mul(y,x,a,M)
        if y==(1,0): return j
    raise RuntimeError('order cap')

def valuation(x,p,k):
    if x==0: return k
    h=0
    while x%p==0: x//=p;h+=1
    return h

def census(a,p,k):
    need(p>2 and (a*a+4)%p!=0,'unramified odd prime')
    pi=order((0,1),a,p,2*(p+1));m=pi
    while m%2==0: m//=2
    q=pi//m;M=p**k;g=power((0,1),q,a,M)
    lift=power((0,1),pi,a,p**8)
    h=min(valuation((lift[0]-1)%p**8,p,8),valuation(lift[1],p,8))
    need(h<8,'valuation precision insufficient')
    expected=Counter()
    for l in range(1,k+1):
        d=m*p**max(0,l-h)
        need((p**(2*l)-p**(2*l-2))%d==0,'integral shell orbits')
        expected[d]+=(p**(2*l)-p**(2*l-2))//d
    visited=bytearray(M*M);visited[0]=1;actual=Counter();shells=Counter()
    for idx in range(1,M*M):
        if visited[idx]: continue
        x=(idx%M,idx//M);y=x;length=0;orbit=set()
        while True:
            key=y[0]+M*y[1]
            need(not visited[key],'orbit intersects previous cycle')
            visited[key]=1;orbit.add(key);length+=1;y=mul(g,y,a,M)
            if y==x: break
        need((-x[0])%M+M*((-x[1])%M) not in orbit,'negation fixed orbit')
        l=k-min(valuation(x[0],p,k),valuation(x[1],p,k))
        need(length==m*p**max(0,l-h),'content-shell length')
        actual[length]+=1;shells[l]+=1
    need(actual==expected,'full orbit histogram')
    d=order(g,a,M,m*p**max(0,k-h))
    fixed=0;y=(1,0)
    for j in range(d):
        v=min(valuation((y[0]-1)%M,p,k),valuation(y[1],p,k))
        fixed+=p**(2*v);y=mul(y,g,a,M)
    need(fixed%d==0 and fixed//d-1==sum(actual.values()),'Burnside fixed kernels')
    total=sum(actual.values());need(total%2==0,'paired orbits')
    return {'a':a,'p':p,'k':k,'split':pow((a*a+4)%p,(p-1)//2,p)==1,
            'pi':pi,'odd_order':m,'lift_delay':h,'nonzero_orbits':total,
            'tournament_exponent':total//2,'length_histogram':dict(actual)}

def isprime(n):
    if n<2: return False
    d=2
    while d*d<=n:
        if n%d==0: return False
        d+=1
    return True

def factors(n):
    result=[];d=2
    while d*d<=n:
        if n%d==0:
            result.append(d)
            while n%d==0:n//=d
        d+=1
    if n>1:result.append(n)
    return result

def prime_controls():
    p=40487;need(isprime(p),'non-generous prime')
    fs=factors(p-1)
    primitive=lambda r: all(pow(r,(p-1)//q,p)!=1 for q in fs)
    need(all(not primitive(r) for r in range(1,5)) and primitive(5),'least primitive root')
    need(pow(5,p-1,p*p)==1 and pow(5,p-1,p**3)!=1,'exact non-generous delay two')
    wag=[]
    for q in (3,5,7,11,13,17,19):
        w=(2**q+1)//3
        if not isprime(w):continue
        x=1;d=0
        while True:
            x=x*2%w;d+=1
            if x==1:break
        need(d==(2 if q==3 else 2*q),'Wagstaff order including exceptional q=3')
        if q>3:need(w%4==3 and (w-1)%(2*q)==0,'odd square subgroup')
        wag.append({'q':q,'prime':w,'order_2':d,'invariant_tournament_exponent':(w-1)//(2*q) if q>3 else None})
    return {'non_generous':{'p':p,'least_primitive_root':5,'delay':2},'wagstaff':wag}

def main():
    cases=[]
    for a in range(1,5):
        for p in (3,5,7,11,13,17,19):
            if (a*a+4)%p==0:continue
            for k in (1,2):cases.append(census(a,p,k))
            if p<=7:cases.append(census(a,p,3))
    cases.append(census(2,31,2))
    data={'status':'FINITE-EXACT checks of an all-depth theorem proved separately',
          'cases':cases,'prime_controls':prime_controls(),
          'total_ring_elements':sum(r['p']**(2*r['k']) for r in cases)}
    Path(sys.argv[1]).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'cases':len(cases),'ring_elements':data['total_ring_elements'],
                      'exceptional':[(r['a'],r['p'],r['k'],r['nonzero_orbits']) for r in cases if r['lift_delay']>1],
                      'prime_controls':data['prime_controls']},indent=2))

if __name__=='__main__':main()
