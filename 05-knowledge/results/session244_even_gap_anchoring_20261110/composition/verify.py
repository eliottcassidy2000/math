#!/usr/bin/env python3
import importlib.util,json,sys
from pathlib import Path
from itertools import product
from math import gcd,isqrt
spec=importlib.util.spec_from_file_location('physical',Path(__file__).with_name('physical.py'))
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)

def act(z,w):
    D,E,U,d=z
    for c in w:
        if c=='O':D,E,U,d=D,0,2*U,d
        elif c=='Z':D,E,U,d=D+U-E,d,2*U-d,d
        elif c=='J':D,E,U,d=D+U-E,D-E,D+U,D
        else:raise ValueError(c)
    return D,E,U,d

def literal(a,r,w):
    for c in w:
        n=len(a)
        a=p.join(a,[0]) if c=='O' else p.graft(a,r,[0])
        if c=='J':r=n
    return a,r

def matrix(w):
    cols=[act(tuple(int(i==j) for i in range(4)),w) for j in range(4)]
    return [list(x) for x in zip(*cols)]

def endpoint(N):
    n=max(4,N.bit_length()+2);a=p.endpoint(n-2)
    selected={i for i in range(1,n-1) if N>>(i-1)&1}
    a=[x|((1<<n) if i in selected else 0) for i,x in enumerate(a)]
    a.append(sum(1<<i for i in range(n) if i not in selected))
    return a,n

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def main():
    if '--fail' in sys.argv:p.need(False,'intentional failure')
    words=['JJJJZ','JOZJZ']; expected=[[[9,-3,11,0],[3,-1,3,0],[9,-3,11,0],[3,-1,3,0]],[[9,-1,11,0],[3,-1,3,0],[9,-1,11,0],[3,-1,3,0]]]
    p.need([matrix(w) for w in words]==expected,'full matrix identity')
    checks=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            a=p.labelled(n,code)
            for r in range(n):
                z=p.state(a,r)
                for c in 'OZJ':
                    b,s=literal(a,r,c)
                    p.need(p.state(b,s)==act(z,c),'one letter physical law');checks+=1
    examples=[]
    for N in range(1,25):
        a,r=endpoint(N);z=p.state(a,r);p.need(z[1]==N,'endpoint target')
        pair=[]
        for w in words:
            b,s=literal(a,r,w);h,end=p.tables(b);v=act(z,w)
            p.need(h[-1]==v[0] and end[-1][s]==v[1] and h[((1<<len(b))-1)^(1<<s)]==v[3],'pair graph count and endpoint and deletion')
            p.need(p.strong(b),'strong output');pair.append(h[-1])
        p.need(pair[1]-pair[0]==2*N,'exact gap')
        examples.append(dict(N=N,source=z,counts=pair,order=len(b)))
    primes=[]
    for a,r in [([20,1,2,7,14],2),([24,1,3,6,14],0),([24,1,3,6,14],3)]:
        z=p.state(a,r);pair=[]
        for w in words:
            b,s=literal(a,r,w);v=p.state(b,s)
            p.need(v==act(z,w) and p.strong(b),'prime physical full state');p.need(prime(v[0]),'trial division prime');pair.append(v[0])
        primes.append(dict(rows=a,root=r,state=z,counts=pair))
    primitive=0
    for h in range(1,501,2):
        for e in range(1,h+1,2):
            lo,hi,k=20*h-3*e,20*h-e,6*h-e
            p.need((lo-3*k)//2==h and 3*lo-10*k==e and (lo-3*k)%2==0,'exact Fminus inverse')
            p.need(gcd(lo,k)==gcd(h,e),'odd determinant two gcd')
            p.need(gcd(lo,hi)==gcd(5*h,e),'pair gcd identity')
            primitive+=1
    # Exact row classes at length five, requiring terminal Z.
    rows={tuple(matrix(''.join(w))[0]) for w in product('OZJ',repeat=5) if w[-1]=='Z'}
    pairs=sorted((r,(r[0],r[1]+2,r[2],r[3])) for r in rows if (r[0],r[1]+2,r[2],r[3]) in rows)
    p.need(len(pairs)==2,'length five row classes')
    print(json.dumps(dict(status='PASS',matrices=dict(zip(words,expected)),one_letter_checks=checks,endpoint_examples=examples,prime_examples=primes,gcd_integer_checks=primitive,length_five_terminal_Z_row_pairs=pairs,scope='Universal endpoint gaps and invariant face only; no arbitrary placement, completeness, or prime infinitude.'),sort_keys=True,indent=2))
if __name__=='__main__':main()
