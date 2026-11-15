from itertools import combinations

def need(t,*a):
    if not t:raise AssertionError(a)
def orient(n,back):
    rows=[0]*n
    for i,j in combinations(range(n),2):
        a,b=(j,i) if (i,j) in back else (i,j)
        rows[a]|=1<<b
    return rows
def DP(rows):
    n=len(rows)
    D=[[0]*n for s in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for s in range(1,1<<n):
        for i,v in enumerate(D[s]):
            if not v:continue
            a=rows[i]&~s
            while a:
                u=a&-a;a-=u;D[s|u][u.bit_length()-1]+=v
    return D
def profile(R):
    n=len(R);q=(1<<n)-1
    D=DP(R)
    V=DP([sum(1<<j for j in range(n) if R[j]>>i&1)for i in range(n)])
    H=[1]+[sum(v) for v in D[1:]]
    return H[-1],[sum(V[s][i]*H[q^s] for s in range(1,q+1)
                           if s.bit_count()%2==0) for i in range(n)]
for j in range(3,13):
    H,E=profile(orient(j+1,{(0,j)}))
    need(H==1+(1<<(j-1)),j)
    A=[(1<<(j-1))-1,(1<<(j-2))+3**(j-2)-1]
    for i in range(2,j):
        d=j-i-1
        A.append((1<<d)+(1<<(i-2))*(2*3**d-1))
    A.append(1<<(j-2))
    need(E==A,j,E,A)
out=[]
for B in ({(0,1),(1,2),(2,3)},{(0,2),(1,2),(1,3)}):
    h,E=profile(orient(5,{(i,j+1) for i,j in B}))
    out.append((h,h+2*E[0]))
need(out==[(11,37),(9,31)],out)
print('one-edge port spectra through order 13; ordered P4 loads:',out)
