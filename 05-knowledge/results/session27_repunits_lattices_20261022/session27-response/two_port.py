"""Audit exact recovery of endpoint and ordered-two-path-cover data.
Includes all labelled cores through k4, deterministic k5/k6 samples,
and cyclic vertex-transitive controls. Full tournament endpoint DP is
independent of the core-permutation path-cover computation.
"""
from common import *

def recover(Q):
    k=len(Q);h=response(Q,());single=[response(Q,(1<<u,)) for u in range(k)]
    L=[x-h for x in single]
    A=[response(Q,(1<<u,0))-response(Q,(0,1<<u)) for u in range(k)]
    J=[[0]*k for _ in range(k)];C=[[0]*k for _ in range(k)];D=[[0]*k for _ in range(k)]
    for u in range(k):
        for v in range(u+1,k):
            D[u][v]=response(Q,(1<<u,1<<v))-response(Q,(1<<v,1<<u));D[v][u]=-D[u][v]
            J[u][v]=A[u]-A[v]-D[u][v];J[v][u]=-J[u][v]
            symmetric=h+L[u]+L[v]-response(Q,((1<<u)|(1<<v),))
            if (symmetric+J[u][v])%2:raise RuntimeError('nonintegral C')
            C[u][v]=(symmetric+J[u][v])//2;C[v][u]=(symmetric-J[u][v])//2
    if any((A[u]-L[u])%2 for u in range(k)):raise RuntimeError('nonintegral S')
    S=[(A[u]-L[u])//2 for u in range(k)]
    E=[L[u]+S[u]-sum(C[u]) for u in range(k)]
    return h,S,E,C,A,D

def audit(Q):
    k=len(Q);h,S,E,C,A,D=recover(Q);actual_h,l,c=coefficients(Q);s=[0]*k;e=[0]*k
    for p in permutations(range(k)):
        if all(Q[p[j-1]][p[j]] for j in range(1,k)):s[p[0]]+=1;e[p[-1]]+=1
    if (h,S,E,C)!=(actual_h,s,e,c):raise RuntimeError(('reconstruction failure',Q,(h,S,E,C),(actual_h,s,e,c)))
    for u in range(k):
        for v in range(k):
            for w in range(k):
                if D[u][v]+D[v][w]+D[w][u] != -(C[u][v]-C[v][u]+C[v][w]-C[w][v]+C[w][u]-C[u][w]):raise RuntimeError('curl mismatch')
    return {'k':k,'H':h,'S':S,'E':E,'A':A,'D01':D[0][1] if k>1 else None}

if __name__=='__main__':
    n=0
    for k in range(1,5):
        for mask in range(1<<(k*(k-1)//2)):audit(core(k,mask));n+=1
    for k in [5,6]:
        maximum=(1<<(k*(k-1)//2))-1
        for m in [0,1,2,3,5,17,41,maximum//3,maximum//2,maximum-1,maximum]:audit(core(k,m));n+=1
    print('Exact sidecar recovery and triangle curl checks PASS:',n,'cores',flush=True)
    for k in [3,5,7]:
        Q=[[int(1<=(v-u)%k<=k//2) for v in range(k)] for u in range(k)]
        print('Cyclic vertex-transitive control:',audit(Q),flush=True)
    Q=core(2,1);D=recover(Q)[-1][0][1]
    if D!=0:raise RuntimeError('two-vertex hostile')
    print('Hostile to D_uv=±2H(Q−{u,v}): k2 arc0→1, D01=0 but H(empty)=1, so RHS=±2.')

    Q=core(3,7);h,S,E,C,A,D=recover(Q)
    if C!=[[0,1,0],[1,0,1],[2,1,0]] or C[0][2]-C[2][0]!=-2:raise RuntimeError('orientation hostile')
    print('Orientation hostile: transitive3 arc0→2 has J02=-2, opposite to its arc; adjacent arc currents vanish.')
