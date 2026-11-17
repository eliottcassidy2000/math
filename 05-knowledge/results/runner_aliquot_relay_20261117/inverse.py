from math import isqrt

A,g,C=13190,37,-5884

def f(x):
    y=1
    while x:
        q=x//2
        if not x&1:y+=q+(q&1)
        x=q
    return y

def read(z):
    return C+A*(z+1)+g*f(z)

def ceil(a,b):return -((-a)//b)

def inverse(N):
    Y=N-C
    if Y<A:return []
    L=2
    while A*L*L<Y:L*=2
    k=L.bit_length()-1
    Q=Y//(A*L)
    F=[1]*max(L,Q+1)
    rev=[0]*L
    for x in range(1,len(F)):
        q=x//2;F[x]=F[q]+(0 if x&1 else q+(q&1))
        if x<L:rev[x]=(rev[q]>>1)|((x&1)<<(k-1))
    ans=[]
    for r in range(L):
        B=A*L+g*(L-1-rev[r])
        C0=Y-A*(r+1)-g*(F[r]-1)
        q=(C0-g)//B
        if q<0 or q<ceil(C0-2*g,B+g):continue
        if (C0-B*q)%g:continue
        if C0==B*q+g*(F[q]+(r<L//2)*(q&1)):
            z=q*L+r
            ans.append(z)
    return sorted(ans)

def s(n):
    if n<1:raise ValueError
    S=0
    for x in range(1,isqrt(n)+1):
        if n%x==0:
            S+=x
            if x*x<n:S+=n//x
    return S-n

if __name__=='__main__':
    from collections import defaultdict
    data=defaultdict(list)
    for z in range(2048):
        N=read(z)
        data[N].append(z)
    for N,ss in data.items():
        if N>A*2048+C:continue
        if inverse(N)!=ss:raise AssertionError((N,ss))
    for p in [11,13,179,181,1049,1051,3889007,5417681,30011,30013,1000000007,1000000009]:
        print(p,inverse(7*p))
    for q in [11,13,1049,1051,25,35,121,127,129]:
        if q%7:
            if s(7*q)!=q+8*s(q):raise AssertionError(q)
            print('aliquot',q,s(7*q),s(q)-1)

    def full(rows):
        n=len(rows);F=(1<<n)-1
        if not all(not(r>>i&1) and all(((r>>j&1)+(rows[j]>>i&1))==1 for j in range(i)) for i,r in enumerate(rows)):
            raise AssertionError
        def dp(R):
            D=[[0]*n for _ in range(F+1)]
            h=[1]+[0]*F
            for i in range(n):D[1<<i][i]=1
            for M in range(1,F+1):
                for i,d in enumerate(D[M]):
                    h[M]+=d
                    if d:
                        S=R[i]&~M
                        while S:
                            a=S&-S;S^=a;D[M|a][a.bit_length()-1]+=d
            return D,h
        D,h=dp(rows)
        P,_=dp([F^(1<<i)^r for i,r in enumerate(rows)])
        return h[-1],sum(P[M][-1]*h[F^M] for M in range(1,F+1) if not M.bit_count()&1),h[F^(1<<(n-1))]
    Arows=(160,13,73,97,15,86,147,318,127)
    n=9
    B=tuple(q+((i!=8)<<n) for i,q in enumerate(Arows))+(1<<8,)
    if full(Arows)!=(1059,3149,533) or full(B)[0]!=7357 or s(7357)!=1059:raise AssertionError
    print('literal reverse aliquot lift 1059 → 7357')
    Crows=(2,36,313,83,355,457,135,287,75)
    Drows=tuple(q+(int(not(226>>i&1))<<9) for i,q in enumerate(Crows))+(226,)
    if full(Crows)[0]!=1057 or full(Drows)[0]!=7343 or s(7343)!=1057:raise AssertionError
    for z,q in ((2063,3889007),(2870,5417681)):
        if read(z)!=7*q or any(any(p%d==0 for d in range(2,isqrt(p)+1)) for p in (q,q+2)):
            raise AssertionError
    print('literal reverse aliquot lift 1057 → 7343')
