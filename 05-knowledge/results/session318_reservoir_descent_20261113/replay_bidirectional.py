from itertools import product
def subtable(A):
    n=len(A);D=[[0]*n for _ in range(1<<n)]
    for j in range(n):D[1<<j][j]=1
    for S in range(1,1<<n):
        for j,c in enumerate(D[S]):
            if c:
                N=A[j]&~S
                while N:
                    b=N&-N;N^=b;D[S|b][b.bit_length()-1]+=c
    return D
for m in range(1,9):
    n=m+2; N=(1<<n)-1
    A=[2, N&~3]+[1+sum(1<<j for j in range(i+1,n)) for i in range(2,n)]
    D=subtable(A)
    R=[sum(1<<j for j,r in enumerate(A) if r>>i&1) for i in range(n)]
    starts=subtable(R)[N]
    h=[1]+list(map(sum,D[1:]))
    if h[-1]!=2**m+1 or starts[:2]!=[1,1]:raise ValueError((m,h[-1],starts))
    for j in range(1,m+1):
        E=sum(D[S][j+1]*h[N^S] for S in range(1,N+1))
        want=2**(m-j+1)*3**(j-1)+2**j-(j==m)
        if E!=want or starts[j+1]!=2**(m-j):raise ValueError((m,j,E,want))
print('36 marked starts and split-end weights checked by independent subset DP')
