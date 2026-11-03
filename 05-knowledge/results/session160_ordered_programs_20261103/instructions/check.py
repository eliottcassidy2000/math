#!/usr/bin/env python3
"""Independent literal-graph controls for the decoded-state constructions."""
from fractions import Fraction

def require(condition, *args):
    if not condition: raise AssertionError(args)

def trans(n):
    return [[int(i<j) for j in range(n)] for i in range(n)]

def ordinal(a,b):
    n,m=len(a),len(b)
    return [[(a[i][j] if i<n and j<n else
              b[i-n][j-n] if i>=n and j>=n else
              int(i<n and j>=n)) for j in range(n+m)] for i in range(n+m)]

def wrap(a):
    n=len(a); r=[row+[1,0] for row in a]
    r.append([0]*n+[0,1]); r.append([1]*n+[0,0])
    return r

def strong(a):
    n=len(a)
    for v in range(n):
        seen={v}; todo=[v]
        while todo:
            i=todo.pop()
            for j in range(n):
                if a[i][j] and j not in seen:
                    seen.add(j);todo.append(j)
        if len(seen)!=n:return False
    return True

def subset_h(a):
    n=len(a); size=1<<n
    dp=[[0]*n for _ in range(size)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,size):
        for i in range(n):
            if dp[s][i]:
                for j in range(n):
                    if not (s>>j)&1 and a[i][j]:dp[s|1<<j][j]+=dp[s][i]
    return [1]+[sum(dp[s]) for s in range(1,size)]

def cover(a,r):
    hs=subset_h(a); full=len(hs)-1
    c=[0]*(full+1);c[0]=1
    for _ in range(r):
        d=[0]*(full+1)
        for s in range(full+1):
            t=s
            while True:
                d[s]+=c[s^t]*hs[t]
                if t==0:break
                t=(t-1)&s
        c=d
    return c[full]

def apex(n):
    bits=bin(n)[2:]+'0';a=trans(len(bits));m=len(bits)
    for row in a:row.append(0)
    a.append([0]*(m+1))
    for i,b in enumerate(bits):a[m][i]=int(b=='1');a[i][m]=int(b=='0')
    return a

def root_int(x,b):
    lo,hi=0,1
    while hi**b<x:hi*=2
    while lo<hi:
        m=(lo+hi)//2
        if m**b<x:lo=m+1
        else:hi=m
    return lo

for n in range(1,8):
    e=wrap(trans(n)); hs=subset_h(e)
    require(strong(e) and hs[-1]==1+2**n,'E',n)
    require((hs[-1]-1 & -(hs[-1]-1)).bit_length()-1==n,'decode',n)
    if n<=5:
        for r in range(5):
            expected=r**(n+2)+r*((r+1)**n-(r-1)**n)
            require(cover(e,r)==expected,'profile',n,r)
        for k in range(3):
            w=wrap(ordinal(e,trans(k)));h=subset_h(w)[-1]
            c=2**(n+2)+2*(3**n-1)
            require(strong(w) and h==1+2**n+2**k*c,'W',n,k,h)
    if n<=3:
        x=hs[-1]
        for bit in (0,1):
            y=subset_h(wrap(trans(2*n+bit)))[-1]
            require(y==1+(2 if bit else 1)*(x-1)**2,'digit',n,bit)
print('PASS strong E(n), profiles, both digit moves and W(n,k) by literal subset DP')

A=[1,2]
for n in range(2,128):
    q=n//2
    A.append(A[q]+q if n%2==0 else A[q]+1-(q%2))
for n in range(1,128):
    h=n+A[n]; t=n%2
    require(2*n+A[2*n] == h+2*n if 2*n<128 else True,'zero',n)
    require(2*n+1+A[2*n+1] == h+n+2-t if 2*n+1<128 else True,'one',n)
    if n<64:
        a=apex(n);require(strong(a) and subset_h(a)[-1]==h,'apex',n)
require([n+A[n] for n in (2,3,4,6)]==[5,5,9,11], 'collision')
print('PASS one-apex automaton and same-count, different-successor collision')

for n in range(1,51):
    for b in range(1,10):
        z=2**n; q=root_int(z,b)
        require((q**b==z)==(n%b==0),'root guard',n,b)
program=[Fraction(3,2),Fraction(5,3)];n=12;trace=[n]
while True:
    allowed=[f for f in program if n*f.denominator and n%f.denominator==0]
    if not allowed:break
    f=allowed[0];n=n*f.numerator//f.denominator;trace.append(n)
require(trace==[12,18,27,45,75,125],'fractran',trace)
for i in range(len(trace)-1):
    n=trace[i];f=next(f for f in program if n%f.denominator==0)
    require(1+root_int(2**n,f.denominator)**f.numerator==1+2**trace[i+1],'count fractran',n)
for m in range(1,512):
    y=2*m+1
    terms=[(1+2**(i+1))-1 for i in range(m.bit_length()) if m>>i&1]
    require(1+sum(terms)==y,'tile',m)
print('PASS first-applicable FRACTRAN, exact root guards and dyadic sums')
