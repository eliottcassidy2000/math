"""Odd-cycle packing polynomial, counted by literal subset DP."""
from functools import cache
def layers(rows):
    n=len(rows);N=1<<n
    cyc=[0]*N
    for i in range(n):
        B=1<<i; d={B:{i:1}}
        for s in range(B,N,B*2):
            if s not in d:continue
            for u,h in d[s].items():
                if s.bit_count()%2 and s.bit_count()>=3 and rows[u]>>i&1:cyc[s]+=h
                r=rows[u]&~s&~(B-1)
                while r:
                    b=r&-r;r-=b;j=b.bit_length()-1
                    a=d.setdefault(s|b,{})
                    a[j]=a.get(j,0)+h
    @cache
    def pack(s):
        if not s:return (1,)
        r=s&-s
        out=list(pack(s^r));c=s
        while c:
            if c&r and cyc[c]:
                a=pack(s^c)
                out.extend([0]*max(0,len(a)+1-len(out)))
                for j,k in enumerate(a,1):out[j]+=cyc[c]*k
            c=(c-1)&s
        return tuple(out)
    return pack(N-1)

def host(rows,x,a,flip=False):
    n=len(rows)
    add=[x,a]
    out=[r+sum((m>>i&1)<<(n+j) for j,m in enumerate(add)) for i,r in enumerate(rows)]
    out += [((1<<n)-1)^m for m in add]
    out[n+flip]|=1<<(n+1-flip)
    return out

if __name__=='__main__':
    rows=(0,53,57,83,225,201,135,15)
    for x,a,h in ((34,68,2325),(35,69,12263),(42,76,2805),(43,77,11623)):
        p,q=(layers(host(rows,x,a,b)) for b in (False,True))
        assert tuple(v-u for u,v in zip(p,q))==(0,5,-2,0)
        assert sum(k*2**i for i,k in enumerate(p))==h
        print(x,a,p,q)
