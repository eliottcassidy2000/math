"""Finite controls for physical norm profiles; proofs are all-order."""
from itertools import combinations
from math import gcd

def hp_subsets(n, bits):
    arcs = [0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if bits >> k & 1:
            arcs[i] |= 1 << j
        else:
            arcs[j] |= 1 << i
    end = [[0]*n for s in range(1<<n)]
    h = [1] + [0]*((1<<n)-1)
    for mask in range(1,1<<n):
        for j in range(n):
            if mask & (1<<j):
                pre = mask ^ (1<<j)
                end[mask][j] = 1 if not pre else sum(
                    end[pre][i] for i in range(n) if arcs[i] >> j & 1
                )
        h[mask] = sum(end[mask])
    return h

def check():
    tally = 0
    h3 = set()
    for n in range(1,6):
        full=(1<<n)-1
        for bits in range(1 << (n*(n-1)//2)):
            hp=hp_subsets(n,bits)
            h=hp[-1]
            q=sum(hp[i]*hp[full^i] for i in range(full+1))
            z=q-h
            if not (z>=n*h and (q-(-1)**n*h)%3==0):
                raise AssertionError((n,bits,h,q))
            if h==3:
                h3.add((n,q))
            if z>h and gcd(z,h)==1:
                for a,b in ((z+2*h,h),(2*z+h,z),(2*z-h,z)):
                    if a>=3*b and (a+b)%4==0:
                        raise AssertionError((z,h,a,b))
            tally+=1
    if h3 != {(3,12),(4,24),(5,48)}:
        raise AssertionError(h3)
    print(f"{tally} labeled profile controls; no primitive physical edge; fixed-H=3 rays at n≤5")
    hp=hp_subsets(5,12)
    if hp[-1]!=11 or sum(hp[i]*hp[31^i] for i in range(32))!=100:
        raise AssertionError("block profile")
    count=0
    for L in range(3,9):
        for a in range(1,5):
            M=2**L*3**a
            P=3**(a-1)
            actual=set()
            for k in range(2**(L-3)*P):
                h=pow(11,k,M)
                q=pow(100,k,M)*pow(2,L,M)%M
                for _ in range(2*P):
                    actual.add(((q-h)**2+2*h*h)%M)
                    q=2*q%M
            target={v for v in range(3,M,8) if v%3!=1}
            if actual!=target:
                raise AssertionError((L,a))
            count+=len(target)
    print(f"{count} primitive joint 2,3-adic residues across 24 moduli")

if __name__=='__main__':
    check()
