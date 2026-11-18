"""Exact replay for the partial uniform parallel-orbit bound and Q6 certificates."""
from fractions import Fraction as F
from math import comb,prod
import sys

def upper_sq(d,h):
    n=d-1
    if h==n:
        A=d*(1<<(d-2))
        return F(A*A,prod(6*comb(n,p) for p in range((n+1)//2)))
    g=n-h;c=comb(h,(h-1)//2)
    A=d*(1<<(d-2))*comb(n,h)
    return F(A*A,prod(6*c*comb(g,r) for r in range(g+1)))

def hist_cert(mask,d=6):
    S=[n for n in range(1<<d) if mask>>n&1]
    G=set()
    for a in range(1<<d):
        for i in range(d):
            if a>>i&1:continue
            b=a^(1<<i);H=[0]*d
            for v in S:H[min((a^v).bit_count(),(b^v).bit_count())]+=1
            G.add(tuple(H))
    assert len(G)==d*(1<<(d-1))
    return len(S)

def main():
    caps=tuple(F(a,b) for a,b in ((1,100000),(8,10000),(13,10000),(44,1000),(132,1000),(19,1000)))
    assert sum(caps)<F(1,5)
    for h,c in enumerate(caps[:5],1):assert upper_sq(11,h)<c*c
    assert all(upper_sq(d,d-1)<caps[-1]**2 for d in (11,12))
    assert hist_cert(0x02283022a042a00a)==15
    if len(sys.argv)>1:
        with open(sys.argv[1]) as f:rows=[int(s.strip(),16) for s in f if s.strip() and not s.startswith('#')]
        assert len(rows)==len(set(rows))==229
        assert {hist_cert(x) for x in rows}=={15}
        print('recounted certificates 229')
    print('strict uniform partial union bound 1/5')
if __name__=='__main__':main()
