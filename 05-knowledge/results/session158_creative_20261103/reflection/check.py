from itertools import permutations
from collections import Counter

def tiling(n, bits):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    T=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            reverse=(j>=i+2 and ((bits>>free.index((i,j)))&1))
            T[i][j]=not reverse; T[j][i]=bool(reverse)
    return T

def chord(T):
    n=len(T)
    return [[False if i==j else (T[i][j] if abs(i-j)==1 else not T[i][j]) for j in range(n)] for i in range(n)]

def subset_tables(T):
    n=len(T); N=1<<n
    dp=[[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v]=1
    for m in range(1,N):
        for v in range(n):
            if not(m>>v&1): continue
            p=m^(1<<v)
            if p: dp[m][v]=sum(dp[p][u] for u in range(n) if p>>u&1 and T[u][v])
    H=[sum(row) for row in dp]; H[0]=1
    return H,dp

def sandwich(A):
    k=len(A); n=2*k+1; T=[[False]*n for _ in range(n)]
    def f(i,j):
        if i<k and j<k: return A[i][j]
        if i>k and j>k: return not A[n-1-i][n-1-j] if i!=j else False
        return i<j
    for i in range(n):
        for j in range(n):
            if i!=j: T[i][j]=f(i,j)
    return T

def brute(T):
    n=len(T)
    return sum(all(T[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))

def main():
    for k in range(1,6):
        m=(k-1)*(k-2)//2 if k>=2 else 0
        pairs=Counter(); records=[]
        for x in range(1<<m):
            A=tiling(k,x); X=chord(A); S=sandwich(A); K=chord(S); n=2*k+1
            ha,_=subset_tables(A); hx,end=subset_tables(X); hs,_=subset_tables(S); hk,_=subset_tables(K)
            h=hx[-1]; a=k-1; starts,_=subset_tables([[X[j][i] for j in range(k)] for i in range(k)])
            # Number starting at a in X equals number ending at a in X opposite.
            _,op_end=subset_tables([[X[j][i] for j in range(k)] for i in range(k)])
            s=op_end[-1][a]
            full=(1<<k)-1
            r=sum(end[u][a]*hx[full^u] for u in range(1,full+1) if u>>a&1)
            if (r-s)%2: raise RuntimeError(("portal parity",k,x,r,s))
            expect=r*r+2*h*r+(h-s)**2
            if (expect-1-2*(s%2))%4: raise RuntimeError("mod 4")
            if hs[-1]!=ha[-1]**2 or hk[-1]!=expect: raise RuntimeError((k,x,ha[-1],hs[-1],h,s,r,hk[-1],expect))
            if not all(S[i][i+1] and K[i][i+1] for i in range(n-1)) or not K[n-1][0]: raise RuntimeError('cycle')
            if any(S[i][j]!=S[n-1-j][n-1-i] or K[i][j]!=K[n-1-j][n-1-i] for i in range(n) for j in range(n)): raise RuntimeError('reflection')
            if k<=3 and (brute(S)!=hs[-1] or brute(K)!=hk[-1]): raise RuntimeError('brute')
            pairs[(ha[-1],hk[-1])]+=1;records.append((x,ha[-1],h,s,r,hk[-1]))
        print('k=',k,'cases=',1<<m,'ordinary and optimized arithmetic check: OK', 'pairs=', sorted(pairs.items()))
        if k==4:
            for rec in records: print('  mask, H(A), H(X), start, r, H(K) =',*rec)
if __name__ == "__main__":
    main()
