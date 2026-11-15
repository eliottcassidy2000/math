"""Replay the optimal-carry injection and the exact carry generating series."""
from collections import Counter

def p(x):
    v=0
    while x:
        u=x//2
        if x%2==0: v+=u+u%2
        x=u
    return v

def F(x,y):
    return 100*x*y+89*x+11*(p(x)+p(y)+1+x%2)

def weight(w,l):
    a,b=1,0
    for i in range(l):
        if w>>i&1: a,b=max(2*a,b),b
        else: a,b=a,max(a,2*b)
    return a

def carry(a,b,l):
    c=0
    for i in range(l):
        if ((a>>i&1)+(b>>i&1)+c)&2: c=1
        else: c=0
        yield c

def need(caption,test):
    if not test: raise RuntimeError(caption)

expected=[1,3,9,27,83,259,819,2611,8371,26931,86835,280371,906035]
for l,a in enumerate(expected):
    s=sum(weight(w,l) for w in range(1<<l))
    need((l,s,a),s==a)
    # (1-z)(1-2z)(1-2z-4z^2)=1-5z+4z^2+8z^3-8z^4.
    v=sum(c*expected[l-i] for i,c in enumerate([1,-5,4,8,-8]) if i<=l)
    need(('gf',l,v),v==([1,-2,-2,2][l] if l<4 else 0))
for r in range(3,10):
    V=1<<r
    ys=set()
    for w in range(V//8):
        groups={}
        for b in range(w+1):
            groups.setdefault(tuple(carry(b,w-b,r-3)),[]).append(b)
        selected=max(groups.values(),key=len)
        need((r,w),len(selected)==weight(w,r-3))
        for u in range(V,2*V):
            for b in selected:
                value=F(V*u+b,V*u+V//2+w-b)
                need((r,u,w,b,value), value not in ys)
                need((r,value),100*V**4<=value<=500*V**4)
                ys.add(value)
    need(('count',r),len(ys)==V*expected[r-3])
    print(r,len(ys))
