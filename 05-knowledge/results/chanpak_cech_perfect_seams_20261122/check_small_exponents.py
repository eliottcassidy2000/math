"""Replay literal p=4 and p=6 seams; exclude order 7 at p=4."""
from itertools import combinations, permutations
from functools import lru_cache

@lru_cache(None)
def H(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i,num in enumerate(dp[s]):
            if not num:
                continue
            r=rows[i]&~s
            while r:
                b=r&-r
                r^=b
                dp[s|b][b.bit_length()-1]+=num
    return sum(dp[-1])

def flip(rows,i,j):
    rows=list(rows)
    rows[i]^=1<<j
    rows[j]^=1<<i
    return tuple(rows)

def strong(rows):
    for i in range(len(rows)):
        s=1<<i
        while True:
            t=s
            for j,r in enumerate(rows):
                if s>>j&1:
                    t|=r
            if t==s:
                break
            s=t
        if s.bit_count()!=len(rows):
            return False
    return True

def reps(n):
    E=list(combinations(range(n),2))
    seen=set()
    for m in range(1<<len(E)):
        if m in seen:
            continue
        r=[0]*n
        for b,(i,j) in enumerate(E):
            if m>>b&1:r[i]|=1<<j
            else:r[j]|=1<<i
        seen.update(sum(((r[p[i]]>>p[j])&1)<<b for b,(i,j) in enumerate(E))
                    for p in permutations(range(n)))
        yield r

if __name__=='__main__':
    records=[
        (4,(192,13,1,5,47,207,30,94),(0,1)),
        (6,(266,436,121,194,201,345,3,101,220),(4,8)),
    ]
    for p,r,ij in records:
        E=(1<<(p-1))*((1<<p)-1)
        s=flip(r,*ij)
        assert strong(r) and strong(s)
        assert (H(r),H(s))==(E-1,E+1)
        print('p',p,'order',len(r),'flip',ij,'counts',H(r),H(s))
    r6=list(reps(6))
    assert len(r6)==56
    hits=0; lower=0
    for r in r6:
        for a in range(64):
            rr=tuple(v|(((a>>i&1)^1)<<6) for i,v in enumerate(r))+(a,)
            if H(rr)==119:
                lower+=1
                hits+=sum(H(flip(rr,i,j))==121
                          for i,j in combinations(range(7),2))
    assert hits==0
    print('order-seven lower presentations',lower,'and bridge presentations',hits)
