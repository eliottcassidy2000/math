#!/usr/bin/env python3
"""Standard-library exact checks. No assert-dependent validation."""
from itertools import combinations, permutations, product
from collections import Counter
from fractions import Fraction
import json

def require(c,m):
    if not c: raise RuntimeError(m)
def induced(A,vs): return tuple(tuple(A[i][j] for j in vs) for i in vs)
def isos(A,B):
    n=len(A)
    if n!=len(B): return []
    da=list(map(sum,A)); db=list(map(sum,B)); p=[-1]*n; used=set(); out=[]
    order=sorted(range(n),key=lambda i:db.count(da[i]))
    def rec(k):
        if k==n: out.append(tuple(p)); return
        i=order[k]
        for j in range(n):
            if j in used or da[i]!=db[j]: continue
            if any(A[i][u]!=B[j][p[u]] for u in order[:k]): continue
            p[i]=j; used.add(j); rec(k+1); used.remove(j)
    rec(0); return out

def root(A):
    n=len(A); universe=set(range(n)); mods=[]
    for bits in range(1,1<<n):
        M={i for i in range(n) if bits>>i&1}
        if all(len({A[v][u] for u in M})==1 for v in universe-M): mods.append(M)
    strong=[M for M in mods if all(not(M&N and M-N and N-M) for N in mods)]
    proper=[M for M in strong if len(M)<n]
    blocks=sorted([sorted(M) for M in proper if not any(M<N for N in proper)])
    Q=tuple(tuple(A[B[0]][C[0]] if i!=j else 0 for j,C in enumerate(blocks)) for i,B in enumerate(blocks))
    return blocks,Q

def cyc(p):
    seen=set(); out=[]
    for i in range(len(p)):
        if i in seen: continue
        j=i; size=0
        while j not in seen: seen.add(j);size+=1;j=p[j]
        out.append(size)
    return tuple(sorted(out))
def poly(G): return {k:Fraction(v,len(G)) for k,v in Counter(map(cyc,G)).items()}
def mul(P,Q):
    R={}
    for a,x in P.items():
        for b,y in Q.items():
            c=tuple(sorted(a+b));R[c]=R.get(c,0)+x*y
    return R

def predicted(A):
    if len(A)==1:return [(0,)],{(1,):Fraction(1)}
    blocks,Q=root(A); children=[induced(A,B) for B in blocks]
    childgroups=[predicted(C)[0] for C in children]
    maps={(i,j):isos(C,D) for i,C in enumerate(children) for j,D in enumerate(children)}
    H=[h for h in isos(Q,Q) if all(maps[i,h[i]] for i in range(len(blocks)))]
    group=[]
    for h in H:
        for lifts in product(*(maps[i,h[i]] for i in range(len(blocks)))):
            p=[-1]*len(A)
            for i,B in enumerate(blocks):
                for a,v in enumerate(B):p[v]=blocks[h[i]][lifts[i][a]]
            group.append(tuple(p))
    Z={}
    for h in H:
        seen=set(); term={():Fraction(1)}
        for i in range(len(h)):
            if i in seen:continue
            j=i;l=0
            while j not in seen:seen.add(j);l+=1;j=h[j]
            factor={tuple(l*d for d in c):v for c,v in poly(childgroups[i]).items()}
            term=mul(term,factor)
        for c,v in term.items():Z[c]=Z.get(c,0)+v/len(H)
    require(Z==poly(group),'cycle recursion mismatch')
    return group,Z

def tournament(n,bits):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j]=(bits>>k)&1;A[j][i]=1-A[i][j]
    return tuple(map(tuple,A))
def inflate(Q,cs):
    vertices=[(i,a) for i,C in enumerate(cs) for a in range(len(C))]
    return tuple(tuple(cs[i][a][b] if i==j else Q[i][j] for j,b in vertices) for i,a in vertices)
def triples(A):return [list(t) for t in combinations(range(len(A)),3) if all(sum(A[i][j] for j in t)==1 for i in t)]
def serial(Z):return {','.join(map(str,k)):str(v) for k,v in sorted(Z.items())}
def main():
    cases=0
    for n in range(1,6):
        for bits in range(1<<(n*(n-1)//2)):
            A=tournament(n,bits);G,Z=predicted(A)
            require(set(G)==set(isos(A,A)),f'full group mismatch {n},{bits}');cases+=1
    C=((0,1,0),(0,0,1),(1,0,0));one=((0,),);chain=((0,1),(0,0))
    A=inflate(chain,[one,C]);dual=tuple(zip(*A));T=inflate(C,[A,A,A]);U=inflate(C,[dual,A,A])
    require(triples(T)==triples(U),'hostile cyclic triples differ')
    data={}
    for name,X,size in [('T',T,81),('U',U,27)]:
        G,Z=predicted(X); direct=isos(X,X)
        require(set(G)==set(direct) and len(G)==size,'hostile group mismatch')
        data[name]={'order':len(G),'cycle_index':serial(Z),'adjacency':X}
    require(len(triples(T))==67,'cyclic triple count')
    print(json.dumps({'status':'PASS','all_labelled_orders_1_through_5':cases,'hostile_same_C3':True,'hostile_C3_edges':67,'hostiles':data},indent=2))
if __name__=='__main__':main()
