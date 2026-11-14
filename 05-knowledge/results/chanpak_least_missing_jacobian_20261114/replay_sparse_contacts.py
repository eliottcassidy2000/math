from itertools import product
from functools import cache

def hp(T):
    n=len(T); m=(1<<n)-1
    @cache
    def dfs(s,i):
        R=s^(1<<i)
        return sum(dfs(R,j) for j in range(n) if ((T[i]&R)>>j)&1) if R else 1
    return sum(dfs(m,i) for i in range(n))

def exact(base,S,edge):
    n=len(base); mask=(1<<n)-1; a=S&mask; b=S>>n
    T=[x + (1<<n if not (a>>i)&1 else 0) + (1<<(n+1) if not (b>>i)&1 else 0) for i,x in enumerate(base)]
    T.extend((a+(1<<(n+1) if edge else 0), b+(1<<n if not edge else 0)))
    return hp(T)

def check(x):
    if not x:raise AssertionError('sparse probe degree verification')

for B,edge in product([(2,4,1),(6,4,0),(6,12,8,1)],(0,1)):
    n=len(B); lo=(1<<n)-1
    probes={S:exact(B,S,edge) for S in range(1<<(2*n))
            if (S&lo).bit_count()<=2 and (S>>n).bit_count()<=2}
    coeff={}
    for S in probes:
        v=0; r=S
        while True:
            v+=(-1)**((S^r).bit_count())*probes[r]
            if not r:break
            r=(r-1)&S
        coeff[S]=v
    for S in range(1<<(2*n)):
        check(exact(B,S,edge)==sum(v for r,v in coeff.items() if r&S==r))
print('768 two-apex exact outputs reconstructed from their separate size-two contact probes.')
