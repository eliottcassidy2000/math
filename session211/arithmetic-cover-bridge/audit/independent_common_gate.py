#!/usr/bin/env python3
"""Standalone permutation audit of common digraph cover and gate identities."""
import itertools, json

def choose(x,k):
    r=1
    for j in range(k): r=r*(x-j)//(j+1)
    return r

def poly(A,q):
    n=len(A)
    return sum(choose(q+sum(A[p[i]][p[i+1]] for i in range(n-1)),n) for p in itertools.permutations(range(n)))

def gate(A):
    n=len(A); B=[[0]*(n+2) for _ in range(n+2)]
    for i in range(n):
        for j in range(n): B[i][j]=A[i][j]
        B[n+1][i]=1; B[i][n]=1
    B[n][n+1]=1
    return B

def hc(A):
    n=len(A)
    if n<2: return 0
    return sum(all(A[p[i]][p[(i+1)%n]] for i in range(n)) for p0 in itertools.permutations(range(1,n)) for p in [(0,)+p0])

def check(a,b,context):
    if a!=b: raise RuntimeError((context,a,b))

checked=0
for n in range(1,4):
    edges=[(i,j) for i in range(n) for j in range(n) if i!=j]
    for m in range(1<<len(edges)):
        A=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(edges): A[i][j]=(m>>k)&1
        C=[[int(i!=j and not A[i][j]) for j in range(n)] for i in range(n)]
        B=gate(A)
        for q in range(-n-3,n+4):
            check(poly(A,-q),(-1)**n*poly(C,q),('reciprocity',n,m,q))
            check(poly(B,q),q*q*poly(A,q)+q*(poly(A,q+1)-poly(A,q-1)),('gate',n,m,q))
        check(hc(B),poly(A,1),('cycle bijection',n,m))
        checked+=1
empty_gap={q: [poly(gate([]),q),q*q] for q in (1,2,3)}
check(empty_gap,{1:[1,1],2:[4,4],3:[9,9]},'empty core polynomial boundary')
check(hc(gate([])),0,'empty gate has no directed Hamiltonian cycle')
check(poly([],1),1,'empty-core path convention')
print(json.dumps({'status':'PASS','method':'direct permutation word enumeration and rooted cycle enumeration; no imported producer functions','digraphs_checked':checked,'universes':'all labeled loopless digraphs on 1..3 vertices, allowing opposite arcs and missing arcs','identities':['word polynomial complement reciprocity','general nonempty digraph gate polynomial','HC(gate D)=H(D)'],'empty_core_gate_operator_valid':empty_gap,'empty_core_cycle_bijection_excluded':{'HC_gate_empty':0,'P_empty_at_1':1}},indent=2))
