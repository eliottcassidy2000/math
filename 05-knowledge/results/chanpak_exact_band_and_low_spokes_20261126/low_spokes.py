"""Finite exact controls; the general statements are proved in README.md."""
from functools import lru_cache

@lru_cache(None)
def F(v):
    if v <= 1:
        return 1
    q, bit = divmod(v, 2)
    return F(q) + (1-bit)*(q+(q&1))

def H(b,e):
    return b*e+F(b-1)*F(e-1)-(b&1)*(e&1)

def cycle_dp(rows):
    n=len(rows); d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for S in range(1,1<<n):
        for i,count in enumerate(d[S]):
            B=rows[i]&~S
            while count and B:
                bit=B&-B; B-=bit
                d[S|bit][bit.bit_length()-1]+=count
    return sum(d[-1])

def core(u,n):
    return [sum(1<<k for k in range(j+1,n))+((u>>j&1)<<n) for j in range(n)]+[((1<<n)-1)^u]

def seam(b,e,n,m):
    A,B=core(b-1,n),core(e-1,m); a=n+1
    A[-1]|=1<<(a+m)
    for j in range(m+1):
        A.append(((1<<a)-1)-((j==m)<<n)+
            sum((B[i]>>j&1)<<(a+i) for i in range(m+1)))
    return A

for b in range(1,18):
    for e in range(1,18):
        assert cycle_dp(seam(b,e,max(1,(b-1).bit_length()),max(1,(e-1).bit_length())))==H(b,e)
fork=[]
for p in range(256):
    for q in range(256):
        for x in (0,1):
            for y in (0,1):
                b,e=2*p+1+x,2*q+1+y
                z=H(b,e)
                if {H(2*p+2-x,e)-z,H(b,2*q+2-y)-z}=={-2,2}:
                    fork.append((b,e,z))
assert fork==[(1,2,3),(2,1,3)]
for s in range(3,81):
    b=1<<s
    for q,W,sign in [(b-2,2*b*b-3*b//2-1,1),(b-3,2*b*b-15*b//4-1,1),
                (b+2,2*b*b+27*b//4+5,-1),(b+1,2*b*b+9*b//2+3,-1)]:
        assert H(b,2*q+1)==W and H(b,2*q+2)==W+2*sign
print('PASS 289 literal seams, 262144 fork settings, and four rays for s=3..80')
