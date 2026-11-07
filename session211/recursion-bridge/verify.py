#!/usr/bin/env python3
"""Exact self-contained tests. No assertions: checks survive python -O."""
from itertools import combinations, permutations
import json
from pathlib import Path

def require(ok, label):
    if not ok: raise RuntimeError(label)

def tournament(n, code):
    A=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j]=(code>>k)&1; A[j][i]=1-A[i][j]
    return A

def opposite(A): return [list(r) for r in zip(*A)]
def induced(A,S): return [[A[i][j] for j in S] for i in S]
def paths(A):
    n=len(A)
    if not n:return (1,[1],[])
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    h=[0]*(1<<n);h[0]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask>>v)&1:continue
            old=mask^(1<<v)
            if old:dp[mask][v]=sum(dp[old][u] for u in range(n) if A[u][v])
        h[mask]=sum(dp[mask])
    return h[-1],h,dp[-1]

def padding(A, gate=False):
    n=len(A);B=[r+[0,1] for r in A]
    B += [[1]*n+[0,0 if gate else 1],[0]*n+[1 if gate else 0,0]]
    return B

def double(A):
    n=len(A);B=[[0]*(2*n) for _ in range(2*n)]
    for i in range(n):
        for j in range(n):
            B[i][j]=A[i][j]
            B[n+i][n+j]=A[j][i]
            B[i][n+j]=1 if i==j else A[i][j]
            B[n+i][j]=0 if i==j else A[i][j]
    return B

def ispath(A,p):return all(A[u][v] for u,v in zip(p,p[1:]))
def energy(A,p):return sum(A[p[j]][p[i]] for i in range(len(p)) for j in range(i+1,len(p)))
def gate_count(A):
    n=len(A);out=[sum(r) for r in A];inn=[n-1-d for d in out]
    return sum(A[u][v] and out[u]==inn[v]==1 for u in range(n) for v in range(n))
def triangles(A):
    return sum(A[a][b]*A[b][c]*A[c][a]+A[a][c]*A[c][b]*A[b][a] for a,b,c in combinations(range(len(A)),3))
def join(A,B):
    return [r+[1]*len(B) for r in A]+[[0]*len(A)+r for r in B]

def main():
    summaries=[]; total=0; reflected=0; sector_direct=0; gates=0
    for n in range(1,6):
        gate_values=set(); cyclic_gate_values=set(); doubled_values={}
        for code in range(1<<(n*(n-1)//2)):
            A=tournament(n,code); h,hs,e=paths(A); s=paths(opposite(A))[2]
            sigma=padding(A); G=padding(A,True); D=double(A)
            p2=sum(hs[S]*hs[((1<<n)-1)^S] for S in range(1<<n))
            hg=paths(G)[0]; hd=paths(D)[0]
            require(paths(sigma)[0]==h,'padding')
            require(hg==h+p2,'apex gate increment')
            require(hg>=1+2**n+triangles(A)*(2+2**(n-1)),'triangle tax')
            require(hg>=1+2**n and ((hg==1+2**n)==(triangles(A)==0)),'gate minimum')
            if triangles(A):
                require(hg>=3+3*2**(n-1),'gate second minimum')
                require((hg==3+3*2**(n-1))==(triangles(A)==1),'gate second equality')
                cyclic_gate_values.add(hg)
            gate_values.add(hg)
            e2=sum(x*x for x in e);s2=sum(x*x for x in s)
            b01=(h*h+e2)//2;b10=(h*h-s2)//2
            require(hd>=b01+b10,'two-block lower bound')
            require(hd>=h*(h+1)//2,'endpoint-free lower bound')
            # Separate explicit list products verify sector counts, not DP states.
            if n<=4:
                words=[p for p in permutations(range(n)) if ispath(A,p)]
                invwords=[p for p in permutations(range(n)) if ispath(opposite(A),p)]
                c01=sum(ispath(D,p+tuple(v+n for v in q)) for p in words for q in invwords)
                c10=sum(ispath(D,tuple(v+n for v in q)+p) for p in words for q in invwords)
                require((c01,c10)==(b01,b10),'independent two-block enumeration');sector_direct+=1
            for p in permutations(range(n)):
                q=p+tuple(v+n for v in reversed(p))
                require(energy(D,q)==n*(n-1)//2+2*energy(A,p),'folded energy')
                require(ispath(D,q)==ispath(A,p),'folded sinks')
                reflected+=1
            require(min(map(sum,double(sigma)))>0,'D Sigma has no sink')
            require(min(map(sum,padding(padding(D))))==0,'Sigma Sigma D has a sink')
            require(gate_count(double(G))==0,'D Gamma is ungated')
            require(gate_count(padding(padding(D,True),True))==1,'Gamma Gamma D is gated');gates+=1
            doubled_values.setdefault(h,set()).add(hd)
            total+=1
        summaries.append({'n':n,'tournaments':1<<(n*(n-1)//2),'gate_minimum':min(gate_values),'gate_second_minimum':min(cyclic_gate_values) if cyclic_gate_values else None,'double_values_by_H':{str(k):sorted(v) for k,v in sorted(doubled_values.items())}})
    C=tournament(3,2);require(triangles(C)==1,'C3 encoding')
    X=join(C,[[0]]);Y=join([[0]],C)
    hx,_,ex=paths(X);hy,_,ey=paths(Y)
    hostile={'C3_then_1':{'H':hx,'H_double':paths(double(X))[0],'end_norm_square':sum(v*v for v in ex),'start_norm_square':sum(v*v for v in paths(opposite(X))[2])},'1_then_C3':{'H':hy,'H_double':paths(double(Y))[0],'end_norm_square':sum(v*v for v in ey),'start_norm_square':sum(v*v for v in paths(opposite(Y))[2])}}
    result={'status':'PASS','universe':'all labelled tournaments of orders 1 through 5; independent explicit two-block enumeration through order 4','tournaments':total,'folded_permutation_words':reflected,'two_block_direct_cases':sector_direct,'gate_noncommutation_cases':gates,'by_order':summaries,'scalar_hostile':hostile,'claims_not_made':['No scalar H+2 or 2H recursion','No proof of unrestricted two-gap spectrum','No classification of all doubled Hamilton paths by two-block sectors']}
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
