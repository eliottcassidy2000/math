"""Wave 25 certificate: the prime target N = 100003 chooses its host.
First core A = TT_2, (h, f) = (1, 4): the line 3 f_B - h_B = N. The census supplies an order-9 core B on it,
(h_B, f_B) = (1997, 34000): 3*34000 - 1997 = 100003. Build T = A -> B -> c -> A and count Hamiltonian paths."""
B_UPPER="111100011111010111111111011110111101"   # gentourng order-9 line: bit = 1 means i -> j (i < j)
def tri(line,n):
    adj=[[0]*n for _ in range(n)]; k=0
    for i in range(n):
        for j in range(i+1,n):
            if line[k]=='1': adj[i][j]=1
            else: adj[j][i]=1
            k+=1
    return adj
def ham(adj,n):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            x=dp[S][v]
            if x:
                for u in range(n):
                    if not (S>>u)&1 and adj[v][u]: dp[S|1<<u][u]+=x
    return dp
def profile(adj,n):
    dp=ham(adj,n); H=[sum(dp[S]) for S in range(1<<n)]; H[0]=1
    return H[-1], sum(H[S]*H[((1<<n)-1)^S] for S in range(1<<n))
A=[[0,1],[0,0]]; B=tri(B_UPPER,9)
print("A = TT_2 profile", profile(A,2), "; B profile", profile(B,9))
na,nb=2,9; n=na+nb+1; adj=[[0]*n for _ in range(n)]
for i in range(na):
    for j in range(na): adj[i][j]=A[i][j]
for i in range(nb):
    for j in range(nb): adj[na+i][na+j]=B[i][j]
for i in range(na):
    for j in range(nb): adj[i][na+j]=1
for j in range(nb): adj[na+j][n-1]=1
for i in range(na): adj[n-1][i]=1
H=sum(ham(adj,n)[(1<<n)-1])
print("H(TT_2 -> B -> c -> TT_2) =",H,"(prime target 100003)")
print("outgoing masks:",[sum(1<<j for j in range(n) if adj[i][j]) for i in range(n)])
