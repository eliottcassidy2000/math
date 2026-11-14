from itertools import permutations,product,groupby
from functools import lru_cache
@lru_cache(None)
def compositions(n,prev=0):
    if n==0:return [()]
    return [(i,)+r for i in range(1,n+1) if not (i==prev==1) for r in compositions(n-i,i)]
def flip_full(L):
    n=sum(L);R=[0]*n
    ends=tuple(accumulate((0,)+L));blocks=[tuple(range(a,b)) for a,b in zip(ends,ends[1:])]
    for i in range(n):R[i]=(1<<i)-1
    for a,b in zip(blocks,blocks[1:]):
        for i,j in ((a[-1],b[0]),(a[0],b[-1])):
            R[i]^=1<<j;R[j]^=1<<i
    return R
from itertools import accumulate
def dp(R):
    n=len(R);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for m,S in enumerate(D):
        for i,c in enumerate(S):
            r=R[i]&~m
            while c and r:
                j=r&-r;r-=j;D[m+j][j.bit_length()-1]+=c
    return sum(D[-1])
def canonical(R):
    n=len(R);order=sorted(range(n),key=lambda i:R[i].bit_count())
    P=product(*(tuple(permutations(g)) for _,g in groupby(order,lambda i:R[i].bit_count())))
    return min(tuple(sum((R[u]>>v&1)<<j for j,v in enumerate(p)) for u in p) for p0 in P for p in [sum(p0,())])
def triangles(R):
    return sum((R[i]>>j&1)==(R[j]>>k&1)==(R[k]>>i&1) for i in range(len(R)) for j in range(i+1,len(R)) for k in range(j+1,len(R)))
if __name__=='__main__':
    for n in range(4,10):
        P=[p for p in compositions(n) if p==p[::-1]]
        RP=[(p,flip_full(p)) for p in P]
        d={}
        for p,R in RP:
            b=0 if len(p)==1 else 2*n-p[0]-p[-1]-2-(p[0]==1)-(p[-1]==1)
            assert b==triangles(R)
            d.setdefault((b,dp(R)),[]).append((p,canonical(R),tuple(sorted(x.bit_count() for x in R))))
        print(n,len(P),len({u for a in d.values() for p,u,s in a}),sorted(d),flush=True)
        for h,v in d.items():
            if len({u for p,u,s in v})>1:print('collision',h,[(p,s) for p,u,s in v],flush=True)

    for a in range(2,7):
        assert dp(flip_full((a,a)))==1+2**(2*a-2)
        k=a-2
        assert dp(flip_full((a,1,a)))==1+6*2**k+4*4**k+2*(k>=1)
    P,Q=(1,2,2,2,1),(2,2,2,2)
    def triangle_degrees(R):
        ins=[0]*len(R)
        for i in range(len(R)):
            for j in range(i+1,len(R)):
                for k in range(j+1,len(R)):
                    if (R[i]>>j&1)==(R[j]>>k&1)==(R[k]>>i&1):
                        for v in (i,j,k):ins[v]+=1
        return sorted(ins)
    assert triangle_degrees(flip_full(P))==[2,2,3,3,4,4,6,6]
    assert triangle_degrees(flip_full(Q))==[2,2,3,3,5,5,5,5]
    assert dp(flip_full(P))==dp(flip_full(Q))==173
