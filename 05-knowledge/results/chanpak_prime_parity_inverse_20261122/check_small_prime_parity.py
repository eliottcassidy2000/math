#!/usr/bin/env python3
"""Literal grafts and local controls for the three truncated signs."""
from itertools import product
from fractions import Fraction
from math import gcd, isqrt

ROWS=((31035,48343),(31045,48353),(1,1))
BLOCKS=((220,117,24,18,192,29,172,46),(158,180,40,66,76,25,39,124))
MATS=(((3,3),(6,8)),((9,9),(30,32)),((3,3),(12,16)))


def need(ok, message):
    if not ok:
        raise ArithmeticError(message)


def dot(row,x):
    return sum(a*b for a,b in zip(row,x))


def add(T,k,u):
    n=len(T)
    out=list(T)
    out[-1]|=1<<(n+k)
    bits=f'{u:0{k}b}'
    for j,b in enumerate(bits):
        r=((1<<n)-1)|(((1<<k)-(1<<(j+1)))<<n)
        if b=='0':
            r|=1<<(n+k)
        out.append(r)
    out.append(((1<<(n-1))-1)|(int(bits[::-1],2)<<n))
    return tuple(out)


def assemble(mid, right):
    n=10+len(right)
    full=(1<<n)-1
    # Left is the converse of the marked forward pair. Marked i=1.
    arcs=[row|(full^3) for row in (0,1)]
    arcs += [(row<<2)|(full^1023) for row in mid]
    arcs += [row<<10 for row in right]
    for i,j in ((1,2),(3,n-1)):
        need((arcs[i]>>j)&1,'original forward seam')
        arcs[i]^=1<<j
        arcs[j]^=1<<i
    for i in range(n):
        reached={i}
        while (expanded := reached | {j for a in reached for j in range(n)
                                     if arcs[a]>>j&1}) != reached:
            reached=expanded
        need(len(reached)==n,'strong graft')
    return tuple(arcs)


def hp(T):
    n=len(T)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        dp[1<<i][i]=1
    for S,row in enumerate(dp):
        for i,c in enumerate(row):
            if c:
                todo=T[i]&~S
                while todo:
                    b=todo&-todo
                    todo^=b
                    dp[S|b][b.bit_length()-1]+=c
    return sum(dp[-1])


def main():
    for mid,row in zip(BLOCKS,ROWS):
        need(hp(assemble(mid,add((2,0),4,6)))==dot(row,(1,2)), 'literal graft')
    need(tuple(ROWS[1][i]-ROWS[0][i] for i in range(2))==(10,10),'scaled gap')
    need(all(gcd(*s)==1 for s in ROWS),'primitive rows')
    det=lambda x,y:x[0]*y[1]-x[1]*y[0]
    need([abs(det(ROWS[i],ROWS[j])) for i,j in ((0,1),(0,2),(1,2))]
         ==[173080,17308,17308],'exceptional primes')
    need(all(4327%d for d in range(2,isqrt(4327)+1)),'4327 prime')
    for p in (5,7,11,17):
        q=p*p
        primitive=[x for x in product(range(q),repeat=2) if gcd(*x,p)==1]
        for row in ROWS:
            for j in (1,2):
                m=p**j
                prob=Fraction(sum(dot(row,x)%m==0 for x in primitive),len(primitive))
                need(prob==Fraction(1,(p+1)*p**(j-1)),'prime-power tail')
        union=sum(any(dot(r,x)%p==0 for r in ROWS) for x in primitive)
        need(Fraction(union,len(primitive))==Fraction(2 if p==5 else 3,p+1),
             'number of kernel lines')
    print("two literal strong grafts; primitive rows; exceptional determinants; 24 prime-power tail controls")


if __name__=='__main__':
    main()
