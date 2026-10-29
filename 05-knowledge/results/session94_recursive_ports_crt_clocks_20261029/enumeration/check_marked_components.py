#!/usr/bin/env python3
"""Exact marked-SCC recurrences; no third-party packages and no class input."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial, gcd
import json
from pathlib import Path

def require(p, msg):
    if not p: raise RuntimeError(msg)

@lru_cache(None)
def oddparts(n, low=1):
    if n == 0: return ((),)
    return tuple((d,)+p for d in range(low,n+1,2) for p in oddparts(n-d,d))

def channel(N, layers, bits):
    require(layers>=1 and bits>=0, 'illegal channel')
    out=[]
    for n in range(N+1):
        ans=Fraction()
        for p in oddparts(n):
            e=sum((d-1)//2 for d in p)+sum(gcd(a,b) for a,b in combinations(p,2))
            z=1
            for d,m in Counter(p).items(): z*=d**m*factorial(m)
            ans+=Fraction(2**(layers*e+bits*len(p)),z)
        require(ans.denominator==1,'nonintegral Burnside')
        out.append(int(ans))
    return out

def mul(a,b):
    N=len(a)-1
    return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(N+1)]

def inv(a):
    require(a[0]==1,'unit required')
    b=[1]+[0]*(len(a)-1)
    for n in range(1,len(a)): b[n]=-sum(a[i]*b[n-i] for i in range(1,n+1))
    return b

def add(a,b,c=1): return [x+c*y for x,y in zip(a,b)]

def marked(N,J,layers=1,bits=0):
    # One independently specified Burnside channel per distinct-port count.
    C=[channel(N-i,layers,bits+layers*i) for i in range(min(J,N)+1)]
    st=[[0]*(J+1) for _ in range(J+1)]; st[0][0]=1
    for j in range(1,J+1):
        for i in range(1,j+1): st[j][i]=st[j-1][i-1]+i*st[j-1][i]
    R=[]
    for j in range(J+1):
        row=[0]*(N+1)
        for i in range(min(j,N)+1):
            scale=st[j][i]*2**(layers*comb(i,2)+bits*i)
            for n in range(i,N+1): row[n]+=scale*C[i][n-i]
        R.append(row)
    A=R[0]; Ai=inv(A); Ai2=mul(Ai,Ai)
    I=[add([1]+[0]*N,Ai,-1)]
    for j in range(1,J+1):
        val=mul(R[j],Ai2)
        for i in range(1,j): val=add(val,mul(Ai,mul(R[i],I[j-i])),-comb(j,i))
        I.append(val)
    # G[b][j] marks exactly b occupied SCCs among j named ports (collisions allowed).
    G=[[[0]*(N+1) for j in range(J+1)] for b in range(J+1)]
    G[0][0]=A
    for b in range(1,J+1):
        for j in range(b,J+1):
            row=[0]*(N+1)
            for i in range(1,j+1): row=add(row,mul(I[i],G[b-1][j-i]),comb(j,i))
            G[b][j]=mul(A,row)
    for j in range(J+1):
        for n in range(N+1):
            require(sum(G[b][j][n] for b in range(J+1))==R[j][n], 'occupied block conservation')
            require(all(G[b][j][n]>=0 for b in range(J+1)), 'negative count')
            require(I[j][n]>=0,'negative strong count')
    return R,I,G

def relabel_mask(mask,p,edges):
    pos={e:i for i,e in enumerate(edges)}; out=0
    for i,(a,b) in enumerate(edges):
        x,y=p[a],p[b]; orient=(mask>>i)&1
        if x>y: x,y=y,x; orient^=1
        out|=orient<<pos[x,y]
    return out

def union_scc(n,masks,edges):
    reach=[1<<i for i in range(n)]
    for mask in masks:
        for k,(a,b) in enumerate(edges):
            if (mask>>k)&1: reach[a]|=1<<b
            else: reach[b]|=1<<a
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1: reach[i]|=reach[k]
    comps=[]; seen=set()
    for i in range(n):
        if i not in seen:
            c=tuple(j for j in range(n) if reach[i]>>j&1 and reach[j]>>i&1)
            comps.append(c); seen.update(c)
    return comps

def literal(N,J,layers,bits):
    tables={}
    for n in range(N+1):
        edges=list(combinations(range(n),2)); pp=list(permutations(range(n)))
        maps=[[relabel_mask(m,p,edges) for m in range(1<<len(edges))] for p in pp]
        cmaps=[[sum(((m>>i)&1)<<p[i] for i in range(n)) for m in range(1<<n)] for p in pp]
        seen=set(); rr=[0]*(J+1); ii=[0]*(J+1); gg=[[Fraction() for j in range(J+1)] for b in range(J+1)]
        for masks in product(range(1<<len(edges)),repeat=layers):
            for cols in product(range(1<<n),repeat=bits):
                state=masks+cols
                if state in seen: continue
                orbit={tuple(maps[k][m] for m in masks)+tuple(cmaps[k][m] for m in cols) for k in range(len(pp))}
                seen.update(orbit)
                aut=[p for k,p in enumerate(pp) if tuple(maps[k][m] for m in masks)+tuple(cmaps[k][m] for m in cols)==state]
                comps=union_scc(n,masks,edges); compid={v:i for i,c in enumerate(comps) for v in c}
                for j in range(J+1):
                    count=Fraction(sum(sum(p[v]==v for v in range(n))**j for p in aut),len(aut))
                    rr[j]+=count
                    if len(comps)==1: ii[j]+=count
                    for p in aut:
                        fixed=[v for v in range(n) if p[v]==v]
                        for ports in product(fixed,repeat=j):
                            b=len({compid[v] for v in ports})
                            gg[b][j]+=Fraction(1,len(aut))
        tables[n]=(rr,ii,gg)
    return tables

def main():
    checks=0
    receipts={}
    for layers,bits,N in [(1,0,5),(2,0,4),(2,1,3),(1,1,4)]:
        R,I,G=marked(20,4,layers,bits)
        lit=literal(N,3,layers,bits)
        for n,(rr,ii,gg) in lit.items():
            for j in range(4):
                require(R[j][n]==rr[j],f'literal full mismatch {layers,bits,n,j}'); checks+=1
                require(I[j][n]==ii[j],f'literal strong mismatch {layers,bits,n,j}'); checks+=1
                for b in range(4):
                    require(G[b][j][n]==gg[b][j],f'literal occupied mismatch {layers,bits,n,j,b}'); checks+=1
        receipts[f'layers{layers}_bits{bits}']={str(n):{'classes':R[0][n],'union_strong':I[0][n],
            'ports1':R[1][n], 'ports2':R[2][n], 'ports3':R[3][n],
            'two_same_component':G[1][2][n], 'two_different_components':G[2][2][n],
            'three_component_occupancy':[G[b][3][n] for b in range(1,4)]} for n in (2,3,11,12,20)}
    R,I,G=marked(12,3)
    require(R[0][:13]==[1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168],'A000568 control')
    require(G[1][2][2]==2 and G[2][2][2]==2,'TT2 hostile')
    R2,I2,G2=marked(2,2,2,0)
    require(I2[0][2]==1 and I[0][2]==0,'union strong is not individual strong')
    receipts['checks']=checks; receipts['range']='formal recurrences through n20 and 4 ports; literal class/automorphism checks to stated small bounds'
    print(json.dumps(receipts,indent=2))
if __name__=='__main__': main()
