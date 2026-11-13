#!/usr/bin/env python3
"""Independent audit: no imports from the packet under review.

Eta path: E(q)^r / E(q^N)^r using formal division, cross-checked against a
fresh filtered-divisor recurrence. Signed actions: literal invariance and
pointwise set stabilizers, no author orbital transport or linear elimination.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path
import json


def need(ok, detail):
    if not ok:
        raise RuntimeError(detail)


def mul(a,b,D):
    out=[0]*(D+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:D+1-i]):
            out[i+j]+=x*y
    return out


def power(a,r,D):
    out=[1]+[0]*D
    for _ in range(r):
        out=mul(out,a,D)
    return out


def divide(a,b,D):
    need(b[0]==1, "formal unit divisor")
    q=[0]*(D+1)
    for n in range(D+1):
        q[n]=a[n]-sum(b[j]*q[n-j] for j in range(1,n+1))
    need(mul(b,q,D)==a, "formal division identity")
    return q


def ordinary_euler(D):
    out=[1]+[0]*D
    for d in range(1,D+1):
        for n in range(D,d-1,-1):
            out[n]-=out[n-d]
    return out


def v2(n):
    if n==0:
        return None
    return (abs(n)&-abs(n)).bit_length()-1


def eta_tests():
    D=80
    E=ordinary_euler(D)
    levels=[N for N in range(2,26) if 24%(N-1)==0]
    result=[]
    for N in levels:
        r=24//(N-1)
        EN=[E[n//N] if n%N==0 else 0 for n in range(D+1)]
        # Full Euler numerator/denominator, not a forbidden-d factor product.
        A=divide(power(E,r,D),power(EN,r,D),D)
        divisor=[0]+[sum(d for d in range(1,n+1) if n%d==0 and d%N!=0)
                     for n in range(1,D+1)]
        rec=[1]
        for n in range(1,D+1):
            numerator=-r*sum(divisor[j]*rec[n-j] for j in range(1,n+1))
            need(numerator%n==0,("integer recursion",N,n))
            rec.append(numerator//n)
        need(A==rec,("independent eta paths",N))
        depth=v2(r)
        for s in range(depth+1):
            root=divide(power(E,r//(1<<s),D),power(EN,r//(1<<s),D),D)
            need(power(root,1<<s,D)==A,("integral root",N,s))
        need(A[1]%(1<<(depth+1))!=0,("first impossible root",N))
        odd=[n for n in range(1,D+1) if A[n]%2]
        need(odd[0]==1<<depth,("first odd exponent",N))
        need(all(n%(1<<depth)==0 for n in odd),("Frobenius support",N))
        for n in range(1,D+1):
            need(A[n]%(1<<max(0,depth-v2(n)))==0,("valuation bound",N,n))
        # Exact divisor cusp orders (including width) in the trivial-character chart.
        cusp={}
        for c in range(1,N+1):
            if N%c==0:
                cusp[c]=Fraction(r*(N-c*c),24*gcd(c*c,N))
        need(cusp[N]==-1 and cusp[1]==1,("infinity/origin cusp",N))
        need(all(v>=0 for c,v in cusp.items() if c!=N),("no other pole",N))
        need(all(v==0 for c,v in cusp.items() if c not in (1,N)),("square cusp zero",N))
        need(isqrt(N**r)**2==N**r,("trivial character square",N))
        # Mobius inverse on a rational real/imaginary coordinate, exactly.
        y=Fraction(1,100*N*N)
        a,b=Fraction(1,N),y
        u,v=1-N*a,-N*b
        den=u*u+v*v
        re,im=(a*u+b*v)/den,(b*u-a*v)/den
        need(re==-Fraction(1,N) and im==1/(N*N*y),("cusp inverse",N))
        result.append({"N":N,"r":r,"degree":D,"root_depth":depth,
                       "first_odd_degree":odd[0],"first13":A[:13],
                       "cusp_orders":{str(c):str(v) for c,v in cusp.items()}})
    return result


def induced(p):
    pts=list(combinations(range(len(p)),2))
    loc={x:i for i,x in enumerate(pts)}
    return tuple(loc[tuple(sorted((p[a],p[b])))] for a,b in pts)


def stabilizer_test(n,action):
    pairs=list(combinations(range(n),2))
    for a,b in pairs:
        for p,negative in action:
            if {p[a],p[b]}=={a,b}:
                swap=int(p[a]==b)
                if swap!=negative:
                    return {"pair":[a,b],"permutation":list(p),"negative":negative,"swap":swap}
    return None


def obeys(mask,n,action):
    pairs=list(combinations(range(n),2))
    loc={x:i for i,x in enumerate(pairs)}
    for p,neg in action:
        for i,(a,b) in enumerate(pairs):
            u,v=p[a],p[b]
            j=loc[tuple(sorted((u,v)))]
            if ((mask>>i)&1)^((mask>>j)&1) != neg^int(u>v):
                return False
    return True


def unsigned_masks(n, generator):
    pairs=list(combinations(range(n),2))
    loc={x:i for i,x in enumerate(pairs)}
    todo=set(pairs)
    bits=[]
    while todo:
        a,b=min(todo)
        cycle=[]
        u,v=a,b
        while True:
            key=tuple(sorted((u,v)))
            cycle.append((loc[key],int(u<v)))
            todo.discard(key)
            u,v=generator[u],generator[v]
            if (u,v)==(a,b):
                break
        bits.append(cycle)
    for select in range(1<<len(bits)):
        mask=0
        for k,cycle in enumerate(bits):
            for i,bit in cycle:
                mask|=(bit^((select>>k)&1))<<i
        yield mask


def action_tests():
    result=[]
    for m in (3,5,7,9):
        plus=[(tuple((i+a)%m for i in range(m)),0) for a in range(m)]
        minus=[(tuple((a-i)%m for i in range(m)),1) for a in range(m)]
        full=plus+minus
        need(stabilizer_test(m,full) is None,("natural signed criterion",m))
        lifted=[(induced(p),s) for p,s in full]
        n=m*(m-1)//2
        obstruction=stabilizer_test(n,lifted)
        need((obstruction is None)==(m==3),("pair signed obstruction",m))
        need(stabilizer_test(n,lifted[:m]) is None,("cyclic survivor",m))
        reflection=lifted[m][0]
        fixed=[i for i in range(n) if reflection[i]==i]
        need(len(fixed)==(m-1)//2,("reflected pair fixed points",m))
        row={"m":m,"pair_vertices":n,"negative_reflection_fixed_vertices":fixed,
             "signed_obstruction":obstruction}
        if m<=5:
            masks=[mask for mask in range(1<<(m*(m-1)//2)) if obeys(mask,m,full)]
            need(len(masks)==1<<((m-1)//2),("literal signed census",m))
            row["literal_signed_tournaments"]=len(masks)
        if m==3:
            count=sum(obeys(mask,n,lifted) for mask in range(1<<(n*(n-1)//2)))
            need(count==2,"m3 signed pair positive")
            row["literal_signed_pair_tournaments"]=count
        if m==5:
            states=list(unsigned_masks(n,lifted[1][0]))
            need(len(states)==len(set(states))==512,"m5 cyclic orientation count")
            need(all(obeys(mask,n,lifted[:m]) for mask in states),"literal C5 invariance")
            need(not any(obeys(mask,n,lifted) for mask in states),"all C5 states fail signed lift")
            row["literal_cyclic_pair_tournaments"]=len(states)
        result.append(row)
    # Negative involution on three points may fix one vertex, but not two.
    positive=[((0,1,2),0),((1,0,2),1)]
    negative=[((0,1,2,3),0),((1,0,2,3),1)]
    need(stabilizer_test(3,positive) is None,"one fixed point allowed")
    need(stabilizer_test(4,negative) is not None,"two fixed points forbidden")
    return result


if __name__=="__main__":
    here=Path(__file__).resolve().parent
    author=here.parent/'sporadic'
    eta=eta_tests()
    original=json.loads((author/'check.json').read_text())
    for row,expected in zip(eta,original['eta']):
        need(row['N']==expected['level'] and row['first13']==expected['first_13_coefficients'],
             ('author receipt coefficients',row['N']))
    print(json.dumps({'status':'PASS: independent finite controls and source-input audit',
                      'eta':eta,'signed':action_tests(),
                      'source_hashes':{name:sha256((author/name).read_bytes()).hexdigest()
                                       for name in ['BRIDGES.md','check.py','check.json']},
                      'audit_code_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},
                     indent=2,sort_keys=True))
