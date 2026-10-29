#!/usr/bin/env python3
"""Scalar plethystic jet inversion; root CHILD slots, not vertex ports."""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb
import json
from check_marked_components import channel, mul, inv, add, require, relabel_mask, union_scc

def power(a,j):
    v=[1]+[0]*(len(a)-1)
    for _ in range(j): v=mul(v,a)
    return v

def jets(N,J):
    A=[[2**comb(j,2)*v for v in channel(N,1,j)] for j in range(J+1)]
    ai=inv(A[0]); ai2=mul(ai,ai)
    S=[add([1]+[0]*N,ai,-1)]
    for j in range(1,J+1):
        v=mul(A[j],ai2)
        for i in range(1,j): v=add(v,mul(ai,mul(A[i],S[j-i])),-comb(j,i))
        S.append(v)
    K=[v[:] for v in S]
    if N: K[0][1]-=1
    if J: K[1][0]-=1
    # Exponential partial Bell polynomial B[j,h](A1,A2,...).
    bell=[[[0]*(N+1) for _ in range(J+1)] for _ in range(J+1)]
    bell[0][0]=[1]+[0]*N
    for j in range(1,J+1):
        for h in range(1,j+1):
            for i in range(1,j-h+2):
                bell[j][h]=add(bell[j][h],mul(A[i],bell[j-i][h-1]),comb(j-1,i-1))
    H=[K[0]]; slots=[K[0]]
    F=A[0][:]; F[0]-=1
    for j in range(1,J+1):
        v=K[j][:]
        for h in range(1,j): v=add(v,mul(H[h],bell[j][h]),-1)
        # A1 has constant1, so integral formal division is possible.
        H.append(mul(v,power(inv(A[1]),j)))
        slots.append(mul(power(F,j),H[j]))
    require(all(c>=0 for row in slots for c in row),'negative slot inventory')
    return slots,H

def strong_modules(n,mask,edges):
    out=[set() for _ in range(n)]
    for k,(a,b) in enumerate(edges):
        if mask>>k&1: out[a].add(b)
        else: out[b].add(a)
    V=set(range(n)); modules=[]
    for m in range(1,1<<n):
        M={i for i in range(n) if m>>i&1}
        if all(len(out[v]&M) in (0,len(M)) for v in V-M): modules.append(M)
    strong=[M for M in modules if all(not(M&L) or M<=L or L<=M for L in modules)]
    return [M for M in strong if M!=V and not any(M<L and L!=V for L in strong)]

def literal(N,J):
    vals=[[0]*(N+1) for _ in range(J+1)]
    checks=0
    for n in range(3,N+1):
        edges=list(combinations(range(n),2)); pp=list(permutations(range(n)))
        seen=set()
        for m in range(1<<len(edges)):
            if m in seen: continue
            images=[relabel_mask(m,p,edges) for p in pp]; seen.update(images)
            if len(union_scc(n,(m,),edges))!=1: continue
            aut=[p for p,im in zip(pp,images) if im==m]
            blocks=strong_modules(n,m,edges)
            require(set().union(*blocks)==set(range(n)) and sum(map(len,blocks))==n,'bad strong partition')
            for j in range(J+1):
                z=Fraction()
                for p in aut:
                    fix=sum({p[v] for v in B}==B for B in blocks)
                    fall=1
                    for i in range(j): fall*=fix-i
                    z+=Fraction(fall,len(aut))
                require(z.denominator==1,'bad slot Burnside')
                vals[j][n]+=int(z); checks+=1
    return vals,checks

def main():
    vals,H=jets(20,4); lit,checks=literal(6,4)
    for j in range(5):
        for n in range(7): require(vals[j][n]==lit[j][n],f'slot mismatch {n,j}')
    require([vals[j][3] for j in range(4)]==[1,1,2,2],'C3 control')
    require(vals[1][4]==3,'strong order4 has three root slots, but four vertex ports')
    print(json.dumps({'scope':'ordinary prime-root child slots, ordered distinct named slots; not actual vertex ports',
        'literal_checks':checks, 'rows':{str(n):[vals[j][n] for j in range(5)] for n in range(3,21)},
        'all_H_integral':all(isinstance(c,int) for row in H for c in row)},indent=2))
if __name__=='__main__': main()
