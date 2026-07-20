# Arborescence count A(T)=sum_r a_r (out-arborescences rooted at r, Matrix-Tree) vs Ham-path H(T).
# a_r = det((D_in - A) with row/col r deleted).  Contrast the {7,21} rule.
from itertools import product
def bareiss_det(M):  # exact integer determinant
    M=[row[:] for row in M]; n=len(M)
    if n==0: return 1
    sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            sw=next((i for i in range(k+1,n) if M[i][k]!=0),None)
            if sw is None: return 0
            M[k],M[sw]=M[sw],M[k]; sign=-sign
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[n-1][n-1]
def ham(n,arc):
    N=1<<n; dp=[0]*(N*n)
    for v in range(n): dp[(1<<v)*n+v]=1
    for S in range(N):
        b=S*n
        for v in range(n):
            c=dp[b+v]
            if c:
                for w in range(n):
                    if arc[v][w] and not (S>>w)&1: dp[(S|(1<<w))*n+w]+=c
    f=(N-1)*n; return sum(dp[f:f+n])
def arb_total(n,arc):  # sum_r out-arborescences rooted at r
    indeg=[sum(arc[j][i] for j in range(n)) for i in range(n)]
    L=[[ (indeg[i] if i==j else 0) - arc[j][i] for j in range(n)] for i in range(n)]  # D_in - A^T? check below
    tot=0
    for r in range(n):
        idx=[i for i in range(n) if i!=r]
        sub=[[L[i][j] for j in idx] for i in idx]
        tot+=bareiss_det(sub)
    return tot
def tour(n,bits):
    arc=[[0]*n for _ in range(n)]; idx=0
    for i in range(n):
        for j in range(i+1,n):
            if (bits>>idx)&1: arc[i][j]=1
            else: arc[j][i]=1
            idx+=1
    return arc
# sanity: transitive n=3 -> A=2 (=(n-1)!); 3-cycle -> A=3
tr=[[0,1,1],[0,0,1],[0,0,0]]; cyc=[[0,1,0],[0,0,1],[1,0,0]]
print("sanity: transitive n=3 A=",arb_total(3,tr),"(want 2); 3-cycle A=",arb_total(3,cyc),"(want 3)")
Aach={}; Hach={}
for n in range(2,7):
    A=set(); H=set(); parities=[0,0]
    for bits in range(1<<(n*(n-1)//2)):
        arc=tour(n,bits); a=arb_total(n,arc); A.add(a); parities[a%2]+=1; H.add(ham(n,arc))
    Aach[n]=A; Hach[n]=H
    print(f"n={n}: ARB values={sorted(A)}")
    print(f"       (even count {parities[0]}, odd count {parities[1]}); HAM values={sorted(H)}")
UA=set().union(*Aach.values()); print("\nARB achieved (n<=6):",sorted(UA))
print("ARB forbidden ints <= max:", [k for k in range(1,max(UA)) if k not in UA])
print("HAM forbidden odds:", [k for k in range(1,15,2) if k not in set().union(*Hach.values())])
