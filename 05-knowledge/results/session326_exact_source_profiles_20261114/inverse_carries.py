"""Exact joint one-apex-over-chain inverse; arbitrary optional-two-slot count."""
from itertools import product

def logical(ok):
    if not ok:
        raise RuntimeError('certificate mismatch')

def pair(bits):
    n=len(bits)
    m=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if bits[i] and not bits[j]: m[j-i]+=1
    H=1+sum(m[d]*((1<<(d-1))+(d==1)) for d in range(1,n))
    Q=(1<<(n+1))+sum(m[d]*(1<<(n-d))*(3**(d-1)+1) for d in range(1,n))
    return H,Q

def v2(a):
    return (a&-a).bit_length()-1

def q_candidates(n,Q):
    c=[(1<<(k+1))*(3**(n-2-k)+1) for k in range(n-1)]
    R=Q-(1<<(n+1))
    if R<c[0]:
        return
    half=n//2
    def go(a,b,t,P):
        k=len(a)
        if k< len(c) and (R-P)%(1<<v2(c[k])):
            return
        if k==half:
            for mid in ((),) if n%2==0 else ((0,),(1,)):
                w=tuple(a)+mid+tuple(1-z for z in reversed(b))
                if pair(w)[1]==Q:
                    yield w
            return
        for x,y in product((0,1),repeat=2):
            z=x+y+sum(a[i]*b[k-i] for i in range(1,k))
            yield from go(a+(x,), b+(y,), t+1, P+c[k]*z)
    yield from go((1,), (1,), 0, c[0])

def inverse(H,Q):
    if H<3 or H%2==0 or Q<=0: return ()
    lower=H.bit_length()
    upper=(H-1).bit_length()+1
    return tuple(w for n in range(lower,upper+1) for w in q_candidates(n,Q) if pair(w)[0]==H)

def main():
    for n in range(2,12):
        buckets={}
        for tail in product((0,1),repeat=n-2):
            w=(1,)+tail+(0,)
            buckets.setdefault(pair(w)[1],set()).add(w)
        for Q,W in buckets.items():
            logical(set(q_candidates(n,Q))==W)
        profiles={pair(w) for W in buckets.values() for w in W}
        for H,Q in profiles:
            found=inverse(H,Q)
            logical(bool(found) and all(pair(w)==(H,Q) for w in found))
    print('arbitrary-Q carry inverse: all active cuts through chain length 11')
if __name__=='__main__': main()
