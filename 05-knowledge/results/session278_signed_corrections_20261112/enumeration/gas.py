"""Exact rational divisibility-defect corrections for ordinary tournament counts."""
from fractions import Fraction as F
from math import factorial,gcd,comb
from collections import Counter
import json,time

def parts(n,lo=1):
    if not n: yield ()
    for a in range(lo,n+1,2):
        for p in parts(n-a,a): yield (a,)+p

def weight(p, gas=False):
    e=sum((a-1)//2 for a in p)
    defects=[]
    for i,a in enumerate(p):
        for b in p[i+1:]:
            baseline=gcd(a,b) if gas and a==3 else min(a,b)
            e+=baseline
            d=baseline-gcd(a,b)
            if d: defects.append(F((1<<d)-1,1<<d))
    z=1
    for a,c in Counter(p).items(): z*=a**c*factorial(c)
    return F(1<<e,z),defects

def elementary(q,k=3):
    e=[F(1)]+[F(0)]*k
    for x in q:
        for j in range(k,0,-1): e[j]+=x*e[j-1]
    return e

def envelope_dp(n):
    # All coefficient arithmetic rational; size is a sufficient response state.
    d=[F(0)]*(n+1);d[0]=F(1)
    for a in range(1,n+1,2):
        z=[F(0)]*(n+1)
        for m,v in enumerate(d):
            if not v:continue
            for j in range((n-m)//a+1):
                e=j*m+j*(a-1)//2+a*comb(j,2)
                z[m+a*j]+=v*F(1<<e,a**j*factorial(j))
        d=z
    return d[n]

def census(n):
    sums=[[F(0)]*5 for _ in range(2)]; lows=[[F(0)]*4 for _ in range(2)]; exact=F(0);count=0;max_edges=0
    for p in parts(n):
        count+=1
        for g in (0,1):
            w,q=weight(p,bool(g)); es=elementary(q,4)
            for k in range(5): sums[g][k]+=w*es[k]
            for k in range(4):
                if len(q)>k:lows[g][k]+=w*F(k+1,len(q))*es[k+1]
            if not g:
                exact+=w*__import__('functools').reduce(lambda a,b:a*(1-b),q,F(1))
                max_edges=max(max_edges,len(q))
    assert exact.denominator==1
    assert sums[0][0]==envelope_dp(n)
    out={'n':n,'odd_partitions':count,'exact':str(exact.numerator),'max_incompatibility_edges':max_edges}
    for g,label in [(0,'plain'),(1,'three_resummed')]:
        cs=sums[g]; bs=[sum(((-1)**j*cs[j] for j in range(k+1)),F(0)) for k in range(4)]
        assert bs[1]<=exact and bs[3]<=exact and exact<=bs[0] and exact<=bs[2]
        for k in range(4): assert lows[g][k]<=abs(bs[k]-exact)<=cs[k+1]
        lo=max(bs[k]+lows[g][k] if k%2 else bs[k]-cs[k+1] for k in range(4))
        hi=min(bs[k]+cs[k+1] if k%2 else bs[k]-lows[g][k] for k in range(4))
        ilo=-(-lo.numerator//lo.denominator); ihi=hi.numerator//hi.denominator
        assert ilo<=exact<=ihi
        bits=(ihi-ilo).bit_length()
        modulus=1<<bits; residue=int(exact)%modulus
        recovered=ilo+(residue-ilo)%modulus
        assert recovered==exact and recovered+modulus>ihi
        out[label]={'interval_low':str(ilo),'interval_high':str(ihi),'residue_bits_sufficient':bits,'remainder_lower':[str(v) for v in lows[g]],'remainder_upper':[str(v) for v in cs[1:]],'signed_errors':[str(b-exact) for b in bs], 'error_bits':[round(__import__('math').log2(abs((b-exact).numerator))-__import__('math').log2((b-exact).denominator),4) if b!=exact else None for b in bs]}
    return out

if __name__=='__main__':
    st=time.time(); rows=[census(n) for n in list(range(25))+[28,32,36,40,48,56,64]]
    # Independent labelled Burnside at n<=7: enumerate permutations and edge orbits.
    import itertools
    for n in range(8):
        total=0
        for perm in itertools.permutations(range(n)):
            seen=set();cycles=[]
            for v in range(n):
                if v in seen:continue
                c=[];u=v
                while u not in seen:seen.add(u);c.append(u);u=perm[u]
                cycles.append(len(c))
            if any(a%2==0 for a in cycles):continue
            edges={(i,j) for i in range(n) for j in range(i+1,n)};o=0
            while edges:
                u,v=next(iter(edges));start=(u,v);o+=1
                while True:
                    edges.discard(tuple(sorted((u,v))));u,v=perm[u],perm[v]
                    if (u,v)==start:break
            total+=1<<o
        assert str(total//factorial(n))==rows[n]['exact'] and total%factorial(n)==0
    for r in range(1,8):
        for f in range(8):
            for a,b,g,q in [(3,5,False,F(3,4)),(5,7,True,F(15,16))]:
                w,qs=weight((1,)*f+(a,)*r+(b,),g)
                assert qs==[q]*r
    assert all(x['three_resummed']['signed_errors'][3]=='0' for x in rows if x['n']<=23)
    assert rows[24]['three_resummed']['signed_errors'][3]!='0'
    print(json.dumps({'status':'PASS','seconds':time.time()-st,'rows':rows},indent=2))
