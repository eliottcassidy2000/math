#!/usr/bin/env python3
"""Independent subset-DP audit of cycle-reversal score-fiber gate blindness."""
import itertools, json

def gate(A):
    n=len(A); B=[[0]*(n+2) for _ in range(n+2)]
    for i in range(n):
        for j in range(n): B[i][j]=A[i][j]
        B[n+1][i]=1; B[i][n]=1
    B[n][n+1]=1
    return B

def gates(A):
    n=len(A)
    return [(a,b) for a in range(n) for b in range(n) if a!=b and A[a][b] and all(A[b][v] and A[v][a] for v in range(n) if v not in (a,b))]

def count(A,cyclic=False):
    n=len(A)
    if n==0: return int(not cyclic)
    dp=[[0]*n for _ in range(1<<n)]
    for i in (range(1) if cyclic else range(n)): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if not dp[mask][i]: continue
            for j in range(n):
                if not mask>>j&1 and A[i][j]: dp[mask|1<<j][j]+=dp[mask][i]
    return sum(dp[-1][i]*(A[i][0] if cyclic else 1) for i in range(n))

def strong(A):
    n=len(A)
    for st in range(n):
        seen={st}; stack=[st]
        while stack:
            i=stack.pop()
            for j in range(n):
                if A[i][j] and j not in seen: seen.add(j);stack.append(j)
        if len(seen)!=n: return False
    return True

def check(cond,msg):
    if not cond: raise RuntimeError(msg)

cases=0; graphs=0; strong_ungated_min=None; witness=None
for n in range(1,6):
    pairs=list(itertools.combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs): A[i][j]=(mask>>k)&1; A[j][i]=1-A[i][j]
        B=gate(A); h=count(A)
        check(count(B,True)==h,('gate HC',n,mask))
        if n>=3 and strong(A) and not gates(A) and (strong_ungated_min is None or n<strong_ungated_min): strong_ungated_min=n
        for v in range(n):
            U=[row[:] for row in B]
            for i,j in [(n,n+1),(n+1,v),(v,n)]: U[i][j],U[j][i]=U[j][i],U[i][j]
            D=[[A[i][j] for j in range(n) if j!=v] for i in range(n) if i!=v]
            check(list(map(sum,B))==list(map(sum,U)),('score',n,mask,v))
            check(strong(U),('strong',n,mask,v))
            check(count(U,True)==count(D),('cycle deletion',n,mask,v))
            if n>=3: check(bool(gates(U))==(sum(A[v]) in (0,n-1)),('gate iff',n,mask,v))
            if n==3 and h==3 and witness is None:
                witness={'core':A,'vertex':v,'gate':B,'reversed':U,'scores':list(map(sum,B)),'gate_HC':count(B,True),'reversed_HC':count(U,True),'gates_before':gates(B),'gates_after':gates(U)}
            cases+=1
        graphs+=1
check(strong_ungated_min==5,('minimal strong ungated',strong_ungated_min))
print(json.dumps({'status':'PASS','method':'direct adjacency matrices, independent subset-DP path/cycle counts, BFS strongness, explicit gate search','core_universe':'all labeled tournaments on 1..5 vertices','cores':graphs,'pointed_reversals':cases,'minimum_strong_ungated_order_at_least_3':strong_ungated_min,'witness':witness},indent=2))
