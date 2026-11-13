from itertools import combinations
from replay_bsd_parity import phi
def rowzero(a):return tuple(x^((x.bit_count()%2)<<i)for i,x in enumerate(a))
def mul(a,b):
    d=[]
    for t in a:
        s=0
        for i,r in enumerate(b):
            if t>>i&1:s^=r
        d.append(s)
    return d
for n in range(2,7):
    E=[(i,j)for i,j in combinations(range(n),2)if j>i+1]
    L=rowzero(tuple(sum(1<<j for j in range(n) if abs(j-i)==1)for i in range(n)))
    for code in range(1<<len(E)):
        T=[1<<(i+1) if i+1<n else 0 for i in range(n)]
        U=T[:]
        for k,(i,j) in enumerate(E):
            if code>>k&1:T[i]|=1<<j;U[j]|=1<<i
            else:T[j]|=1<<i;U[i]|=1<<j
        A=rowzero(T);p,q=phi(T),phi(U)
        D=tuple(x^y for x,y in zip(p,q))
        want=tuple(a^b^c^d for a,b,c,d in zip(mul(A,L),mul(L,A),mul(L,L),L))
        if D!=want or any(d>>i&1 or d.bit_count()%2 for i,d in enumerate(D)):
            raise ValueError((n,code))
        for i,j in combinations(range(1,n-1),2):
            if j-i>=3 and ((D[i]>>j)^(A[i]>>(j-1))^(A[i]>>(j+1))^(A[i-1]>>j)^(A[i+1]>>j))&1:
                raise ValueError('stencil')
print('line stencil checked on every path-pinned tournament through order six')
