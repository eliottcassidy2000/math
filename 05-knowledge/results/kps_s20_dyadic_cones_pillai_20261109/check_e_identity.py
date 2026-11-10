# check e_v(T) (# Hamiltonian paths ending at v) == O - E, where O/E = sum over paths P from v with |V(P)| odd/even
# of H(T - V(P)), on random tournaments of order <= 9 (brute-force subset DP in Python)
import random, itertools
random.seed(11)
def counts(n,adj):
    N=1<<n
    end=[[0]*n for _ in range(N)]; st=[[0]*n for _ in range(N)]; H=[0]*N; H[0]=1
    for S in range(1,N):
        for v in range(n):
            if not S>>v&1: continue
            P=S&~(1<<v)
            if P==0: end[S][v]=1; st[S][v]=1; continue
            end[S][v]=sum(end[P][w] for w in range(n) if P>>w&1 and adj[w][v])
            st[S][v]=sum(st[P][w] for w in range(n) if P>>w&1 and adj[v][w])
        H[S]=sum(end[S][v] for v in range(n) if S>>v&1)
    return end,st,H
bad=tot=0
for n in range(1,10):
    for _ in range(12 if n<9 else 4):
        adj=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                if random.random()<0.5: adj[i][j]=1
                else: adj[j][i]=1
        end,st,H=counts(n,adj); F=(1<<n)-1
        for v in range(n):
            O=E=0
            for S in range(1,1<<n):
                if S>>v&1:
                    c=st[S][v]*H[F&~S]
                    if bin(S).count('1')%2: O+=c
                    else: E+=c
            tot+=1
            if end[F][v]!=O-E: bad+=1
print('ports checked',tot,'violations of e = O - E:',bad)
