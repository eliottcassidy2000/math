#!/usr/bin/env python3
"""Replay the physical contact matrices, binary inverse, and odd CRT generators."""
from collections import deque
from itertools import product
from math import gcd
from functools import cache

@cache
def F(u):
    if u == 0:
        return 1
    v = u//2
    return F(v) if u%2 else F(v) + v + v%2

def matrix(k,u):
    p = 1<<k
    b = int(f"{u:0{k}b}"[::-1],2)
    return F(u),u+1,p-b,p

def act(A,v):
    a,b,c,d = A
    x,y = v
    return a*x+b*y,c*x+d*y

def mul(A,B):
    a,b,c,d=A
    e,f,g,h=B
    return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h

def determinant(A):
    a,b,c,d=A
    return a*d-b*c

def inverse_at(A,v):
    a,b,c,d=A; det=determinant(A)
    x,y=act((d,-b,-c,a),v)
    if x%det or y%det:
        raise ArithmeticError("not exact")
    q,r=x//det,y//det
    if min(q,r)<=0: raise ArithmeticError("not positive")
    return q,r

def require(value,message):
    if not value: raise ArithmeticError(message)

A=matrix(3,2)
B=matrix(4,4)
W=(1,0,0,1)
for k,u in [(5,10),(3,2),(5,14),(5,9),(4,4)]:
    M=matrix(k,u)
    require(determinant(M)!=0,"rank-one actual reset step")
    W=mul(M,W)
seed=act(B,(1,2))

def decode(g,gg):
    require((2*g-gg)%15==0 and (gg-g)%19==0,"inexact reading")
    v=((2*g-gg)//15,(gg-g)//19)
    out=[]
    while True:
        tag={2:"A",4:"B"}[ (15*v[0]+19*v[1])%15 ]
        ready=inverse_at(A if tag=="A" else B,v)
        require(tuple(x%15 for x in ready)==(0,1),"not ready")
        out.append(tag)
        if ready==seed:
            return "".join(reversed(out))
        v=inverse_at(W,ready)

def branches():
    require(tuple(x%15 for x in W)==(0,0,1,11),"reset matrix")
    require(seed==(15,46),"seed")
    ma=matrix(5,18); mb=matrix(6,37)
    require(ma[:2]==(15,19) and mb[:2]==(15,38),"readouts")
    for bits in range(1,11):
        for word in product("AB",repeat=bits):
            v=seed
            for j,tag in enumerate(word):
                if j: v=act(W,v)
                require(tuple(x%15 for x in v)==(0,1),"online phase")
                v=act(A if tag=="A" else B,v)
            h,R=v
            require((h%9!=0 and h%3==0 and h%5!=0) if word[-1]=="A"
                    else (h%25!=0 and h%5==0 and h%3!=0),"exact tag")
            g,gg=15*h+19*R,15*h+38*R
            require(gcd(g,15)==gcd(gg,15)==1,"forced divisor")
            require(g%15==(2 if word[-1]=="A" else 4),"visible tag")
            require(decode(g,gg)=="".join(word),"inverse")
    print("word inverse: 2046 binary words through depth ten")

def grow(k,u,L):
    return k+2*L,(((u<<L)+1)<<L)-1

def crt():
    for m in range(3,62,2):
        L=next(j for j in range(1,m+1) if pow(2,j,m)==1)
        k,u=grow(2,1,L)
        P1=matrix(k,u)
        P2=matrix(*grow(k,u,L))
        Q=matrix(k+1,u)
        T=matrix(*grow(1,1,L))
        for actual,residue in ((P1,(2,2,2,4)),(P2,(3,2,2,4)),(Q,(2,2,4,8)),(T,(2,2,1,2))):
            require(all((a-b)%m==0 for a,b in zip(actual,residue)),"contact loop")
        def reach(all_residues=False):
            bank=[P1,P2,Q]
            if all_residues: bank.append(matrix(1,0))
            seen={(1,2)}
            todo=deque(seen)
            while todo:
                v=todo.popleft()
                for g in bank:
                    w=tuple(x%m for x in act(g,v))
                    if w not in seen:
                        seen.add(w)
                        todo.append(w)
            return {tuple(x%m for x in act(T,v)) for v in seen}
        expected={(a,b) for a in range(m) for b in range(m) if gcd(gcd(a,b),m)==1}
        require(reach()==expected,"primitive profile")
        require(len(reach(True))==m*m,"nonprimitive profile")
    print("joint CRT: every primitive and nonprimitive vector modulo every odd m<=61")

if __name__=="__main__":
    branches()
    crt()
