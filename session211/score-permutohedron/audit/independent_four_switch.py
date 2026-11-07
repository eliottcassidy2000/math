#!/usr/bin/env python3
"""Independent direct-matrix checks of two score-fiber cycle-count witnesses."""
from itertools import permutations
import json

def check(x,label):
    if not x:raise RuntimeError(label)

def matrix(rows):return [[row>>j&1 for j in range(len(rows))] for row in rows]
def count(A,cycle=False):
    n=len(A)
    if n==0:return int(not cycle)
    dp=[[0]*n for _ in range(1<<n)]
    for i in (range(1) if cycle else range(n)):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            if dp[mask][i]:
                for j in range(n):
                    if not mask>>j&1 and A[i][j]:dp[mask|1<<j][j]+=dp[mask][i]
    return sum(dp[-1][i]*(A[i][0] if cycle else 1) for i in range(n))
def gates(A):
    n=len(A)
    return [(a,b) for a in range(n) for b in range(n) if a!=b and A[a][b] and all(A[b][v] and A[v][a] for v in range(n) if v not in (a,b))]
def connected(A):
    for root in range(len(A)):
        reached={root};todo=[root]
        while todo:
            i=todo.pop()
            for j in range(len(A)):
                if A[i][j] and j not in reached:reached.add(j);todo.append(j)
        if len(reached)!=len(A):return False
    return True

cases=[([24,1,3,6,14],0,1,[56,33,35,38,46,64,31],[26,96,35,38,46,1,61],7),([30,44,56,48,34,1],1,4,[94,108,120,112,98,65,128,63],[94,60,120,112,224,65,2,111],21)]
out=[]
for T,u,v,G,U,expected in cases:
    A=matrix(T);B=matrix(G);C=matrix(U);n=len(T)
    for i in range(n):
        check(B[i][:n]==A[i],'core restriction')
        check(B[n+1][i] and B[i][n],'cross gates')
    check(B[n][n+1],'gate arc')
    R=[row[:] for row in B]
    for x,y in [(n,n+1),(n+1,v),(v,u),(u,n)]:
        check(R[x][y]==1,'directed four-cycle')
        R[x][y],R[y][x]=R[y][x],R[x][y]
    check(R==C and list(map(sum,B))==list(map(sum,C)),'score-preserving switch')
    A[u][v],A[v][u]=A[v][u],A[u][v]
    z=0
    for mask in range(1<<n):
        S=[x for x in range(n) if mask>>x&1]
        if u not in S or v not in S:continue
        middle=[x for x in S if x not in (u,v)]
        xcount=sum(all(A[p[i]][p[i+1]] for i in range(len(p)-1)) for mid in permutations(middle) for p in [(u,)+mid+(v,)])
        rem=[[A[i][j] for j in range(n) if j not in S] for i in range(n) if i not in S]
        z+=xcount*count(rem)
    h=count(C);hc=count(C,True)
    direct_h=sum(all(C[p[i]][p[i+1]] for i in range(len(C)-1)) for p in permutations(range(len(C))))
    direct_hc=sum(all(C[p[i]][p[(i+1)%len(C)]] for i in range(len(C))) for q in permutations(range(1,len(C))) for p in [(0,)+q])
    check(z==hc==direct_hc==expected,'three cycle count methods')
    check(h==direct_h,'two path count methods')
    check(connected(C) and not gates(C),'strong and ungated')
    out.append({'core_rows':T,'u':u,'v':v,'gate_rows':G,'switched_rows':U,'scores':list(map(sum,B)),'H_gate':count(B),'H_switched':h,'HC_gate':count(B,True),'HC_switched':hc,'forest_count':z,'strong':True,'gates':[]})
print(json.dumps({'status':'PASS','method':'independent adjacency, subset DP, full permutations and endpoint-permutation/forest convolution','examples':out},indent=2))
