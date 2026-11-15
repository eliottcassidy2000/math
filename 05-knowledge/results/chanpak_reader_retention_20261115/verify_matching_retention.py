#!/usr/bin/env python3
"""Independent visit optimizer and literal slot-profile controls; stdlib only."""
from itertools import combinations, product
from functools import lru_cache
from random import Random
from fractions import Fraction
from math import factorial, comb

def graph(n,code):
    a=[0]*n
    for i,(x,y) in enumerate(combinations(range(n),2)):
        if code>>i&1: a[x]|=1<<y
        else: a[y]|=1<<x
    return a

def matching(a):
    n=len(a)
    O=[v for v in range(1,n) if a[0]>>v&1]
    I=[v for v in range(1,n) if a[v]&1]
    best=0; options={0:1}
    for o in O:
        new=options.copy()
        for mask,mult in options.items():
            for i in I:
                w=1<<i
                if a[o]&w and not mask&w:
                    new[mask|w]=new.get(mask|w,0)+mult
        options=new
    best=max(x.bit_count() for x in options)
    return best,sum(c for s,c in options.items() if s.bit_count()==best)

def visit_counts(a,maxr):
    n=len(a); full=(1<<n)-2
    # Enumerate one path by adding unused non-port vertices, or another
    # copy of the port when the preceding letter is non-port.
    paths=[{} for _ in range(1<<n)]
    @lru_cache(None)
    def trail(mask,last):
        ans={(mask,0):1}
        for v in range(n):
            if last!=-1 and not (a[last]>>v&1):continue
            if v and mask>>v&1:continue
            for (s,k),ways in trail(mask|(1<<v if v else 0),v).items():
                key=s,k+(v==0)
                ans[key]=ans.get(key,0)+ways
        return ans
    for (s,k),ways in trail(0,-1).items():paths[s][k]=ways
    R=[{} for _ in range(maxr+1)]
    D={(0,0):1}
    for r in range(1,maxr+1):
        B={}
        for (s,k),prev in D.items():
            m=full^s; sub=m
            while True:
                for j,c in paths[sub].items():
                    key=s|sub,k+j
                    B[key]=B.get(key,0)+prev*c
                if not sub:break
                sub=(sub-1)&m
        D=B
        R[r]={k:c for (s,k),c in B.items() if s==full and k}
    return R

def path_table(a):
    n=len(a); F=(1<<n)-1
    d=[[0]*n for _ in range(F+1)]
    for i in range(n):d[1<<i][i]=1
    for s in range(F+1):
        for i,val in enumerate(d[s]):
            if val:
                m=a[i]&~s
                while m:
                    bit=m&-m;m^=bit
                    d[s|bit][bit.bit_length()-1]+=val
    return [1]+[sum(d[s]) for s in range(1,F+1)]

def pi(a,r):
    h=path_table(a); F=len(h)-1
    f=[1]+[0]*F
    for _ in range(r):
        out=[]
        for s in range(F+1):
            score=0;sub=s
            while True:
                score+=f[sub]*h[s^sub]
                if not sub:break
                sub=(sub-1)&s
            out.append(score)
        f=out
    return f[-1]

def plug(a,b):
    m=len(b);groups=[list(range(m))]+[[m+i-1]for i in range(1,len(a))]
    x=b+[0]*(len(a)-1)
    for i,s in enumerate(groups):
        for j,t in enumerate(groups):
            if a[i]>>j&1:
                for z in s:
                    x[z]|=sum(1<<w for w in t)
    return x

def main():
    count=sat=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=graph(n,code); v,b=matching(a);R=visit_counts(a,3)
            for r in (1,2,3):
                assert set(R[r])==set(range(1,r+v+1))
                if n-1==2*v:
                    assert R[r][r+v]==b*factorial(v)*comb(v+r-1,r-1)
                    sat+=1
                count+=1
    tri=[2,4,1]
    assert visit_counts(tri,2)[2]=={1:12,2:10,3:2}
    rng=Random(1001)
    for _ in range(80):
        n=rng.randrange(1,6);a=graph(n,rng.randrange(1<<(n*(n-1)//2)))
        b=plug(tri,a)
        for r in range(1,5):
            assert pi(b,r)==r*(pi(a,r+1)+r*pi(a,r)-pi(a,r-1))
    a=graph(8,87558995); b=graph(8,249732159)
    assert [pi(a,r) for r in (1,2,3)]==[345,5640,45225]
    assert [pi(b,r) for r in (1,2,3)]==[345,5640,44505]
    x,y=plug(tri,a),plug(tri,b)
    assert pi(x,1)==pi(y,1)==5985 and pi(x,2)-pi(y,2)==1440
    assert all(pi(a,r)-pi(b,r)==2*r*r*(r*r-1)*(r*r-4) for r in range(1,6))
    g=[Fraction(0),Fraction(1),Fraction(0)]
    for r in range(2,31):
        g.append(g[r-1]-(r-Fraction(1,r))*g[r])
    assert g[3:6]==[1,Fraction(-8,3),11]
    for r in range(3,len(g)-1):
        assert g[r]*g[r+1]<0 and abs(g[r])<abs(g[r+1])
    print(count,'visit demands,',sat,'saturated formulas, 320 literal triangle reads passed')

if __name__=='__main__':main()
