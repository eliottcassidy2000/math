from itertools import combinations
from functools import cache

def make(W,back,Z=None):
    # x,y,block,z followed by A,B. Block Z wins against both apices.
    n=len(W)
    out=[0]*(n+2)
    for i in range(n):
        for j in range(i+1,n):
            edge = (Z[i-2]>>(j-2)&1) if Z is not None and 2<=i<j<n-1 else 1
            a,b=(i,j) if edge else (j,i)
            out[a]|=1<<b
        for j,v in enumerate([2,1]):
            if W[i]&v:
                out[n+j]|=1<<i
            else:
                out[i]|=1<<(n+j)
    out[n+back]|=1<<(n+1-back)
    return out

def count_all(out):
    @cache
    def f(m,i):
        prev=m^(1<<i)
        if not prev: return 1
        return sum(f(prev,j) for j in range(len(out)) if (prev&out[i])>>j&1)
    H=[1]+[sum(f(m,j) for j in range(len(out)) if m>>j&1) for m in range(1,1<<len(out))]
    return H

def check():
    cases=0
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for code in range(1<<len(pairs)):
            Z=[0]*n
            for i,(a,b) in enumerate(pairs):
                if code>>i&1:a,b=b,a
                Z[a]|=1<<b
            H=count_all(Z); F=(1<<n)-1
            @cache
            def P(m,k):
                if k==1:return H[m]
                return sum(H[s]*P(m^s,k-1) for s in range(m+1) if s&m==s)
            for back in (0,1):
                out=make([1,2]+[0]*n+[2],back,Z)
                assert count_all(out)[-1]==4*P(F,3)+9*P(F,2)+(2*back-3)*H[F]
            cases+=1
    print('Both literal shield orientations checked for all',cases,'labelled blocks on 1–4 vertices.')
if __name__=='__main__':check()
