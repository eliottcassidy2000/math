from itertools import product,permutations

def adj(b):
    m=len(b)
    return [[int(i<j) if i<m and j<m else b[i] if i<m else 1-b[j] if j<m else 0 for j in range(m+1)] for i in range(m+1)]

def nullity(B):
    n=len(B); A=[r[:] for r in B]
    for i in range(n): A[i][i]=sum(A[i])%2
    P=[sum(((sum(A[i][k]*A[k][j] for k in range(n))+A[i][j]+(i==j))%2)<<j for j in range(n)) for i in range(n)]
    lead={}
    for r in P:
        while r:
            p=r.bit_length()
            if p not in lead:
                lead[p]=r
                break
            r^=lead[p]
    return n-len(lead)

tested=0
for m in range(1,11):
    strong=total=0
    for b in product((0,1), repeat=m):
        nu=nullity(adj(b))
        D=sum((1 if (m-i-1)%2==0 else -1)*bit for i,bit in enumerate(b))
        expected=2 if D%3==1+m%2 else 0
        if nu!=expected: raise ValueError(b,nu,expected)
        total+=nu==2
        strong+=b[0]==0 and b[-1]==1 and nu==2
        tested+=1
    if total != (2**m-1-m%2)//3: raise ValueError(m,total)
    if m>=2 and strong != (2**(m-2)+(2 if m%2==0 else -2))//3: raise ValueError(m,strong)
for b,nu_expected in [((0,0,0,1),2),((0,1,0,1),0)]:
    B=adj(b)
    H=sum(all(B[i][j] for i,j in zip(p,p[1:])) for p in permutations(range(len(B))))
    if (H,nullity(B))!=(9,nu_expected): raise ValueError(b,H,nullity(B))
print(f'{tested} binary rank checks; two direct permutation checks; all pass')
