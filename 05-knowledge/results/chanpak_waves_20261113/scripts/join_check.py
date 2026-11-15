"""Wave 25: brute-force check of the cyclic-join identity (chanpak_absolute_reflection_20261114):
T = A -> B -> c -> A (all arcs A->B, B->c, c->A),  H(T) = (f_A - h_A)(f_B - h_B) + 2 h_A h_B,
with h = H, f = F_2 = sum over subsets S of H(S) H(complement) (H(empty) = 1)."""
import random, itertools
def hsub(adj,n):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            if not dp[S][v]: continue
            for u in range(n):
                if not (S>>u)&1 and adj[v][u]: dp[S|1<<u][u]+=dp[S][v]
    H=[sum(dp[S]) for S in range(1<<n)]; H[0]=1
    return H
def prof(adj,n):
    H=hsub(adj,n); full=(1<<n)-1
    return H[full], sum(H[S]*H[full^S] for S in range(1<<n)), H
def rand_t(n):
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if random.random()<0.5: adj[i][j]=1
            else: adj[j][i]=1
    return adj
random.seed(1); bad=0
for trial in range(300):
    a=random.randint(1,5); b=random.randint(1,5)
    A=rand_t(a); B=rand_t(b); n=a+b+1
    adj=[[0]*n for _ in range(n)]
    for i in range(a):
        for j in range(a): adj[i][j]=A[i][j]
    for i in range(b):
        for j in range(b): adj[a+i][a+j]=B[i][j]
    c=a+b
    for i in range(a):
        for j in range(b): adj[i][a+j]=1
    for j in range(b): adj[a+j][c]=1
    for i in range(a): adj[c][i]=1
    hA,fA,_=prof(A,a); hB,fB,_=prof(B,b); hT,fT,_=prof(adj,n)
    if hT!=(fA-hA)*(fB-hB)+2*hA*hB: bad+=1
print("cyclic join identity: 300 random trials, mismatches:",bad)
