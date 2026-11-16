from functools import lru_cache
from itertools import combinations
from math import isqrt

def counts(rows):
    n=len(rows)
    @lru_cache(None)
    def path_end(s,v):
        s ^= 1<<v
        if not s:return 1
        return sum(path_end(s,u) for u in range(n) if (s>>u)&(rows[u]>>v)&1)
    def H(s):
        return sum(path_end(s,v) for v in range(n) if s>>v&1) if s else 1
    return H

def run():
    hist={}
    for n in range(1,7):
        arcs=list(combinations(range(n),2))
        tri=list(combinations(range(n),3))
        f=(1<<n)-1
        for code in range(1<<len(arcs)):
            rows=[0]*n
            for j,(a,b) in enumerate(arcs):
                if code>>j&1: rows[a]|=1<<b
                else: rows[b]|=1<<a
            cycles=[(1<<a)|(1<<b)|(1<<c) for a,b,c in tri
                    if (rows[a]>>b&1)==(rows[b]>>c&1)!=(rows[a]>>c&1)]
            H=counts(tuple(rows)); hp=H(f)
            for x in range(n):
                base=f^(1<<x); h=H(base)
                r,odd=divmod(hp-h,2)
                assert not odd and r>=0
                if r:
                    tau=sum(c&base==c for c in cycles)
                    qmax=max(H(f^c) for c in cycles if c>>x&1)
                    assert tau+(qmax-1)//2<=r
                    if r<=5:hist.setdefault(r,set()).add((h,hp))
    print('Order at most 6, gap pairs:',{k:sorted(v) for k,v in sorted(hist.items())})
    active={}
    for a in range(1,7825):
        for b in range(a+1,isqrt(7825*7825-a*a)+1):
            c=isqrt(a*a+b*b)
            if c*c==a*a+b*b:
                active.setdefault(c,[]).append((a,b))
    ac={v for c,ab in active.items() for a,b in ab for v in (a,b,c)}
    assert len(ac)==6494 and len(active[7825])==7
    print('7825 triples:',active[7825])
if __name__=='__main__':run()

def check_cubes():
    for target,row,a,b in [(7,(22,28,24,17,0),2,4),
                           (21,(22,60,56,49,32,1),0,4),
                           (35,(14,52,56,50,33,1),1,4)]:
        row=list(row)
        f=(1<<len(row))-1
        before=counts(tuple(row))(f)
        row[a]^=1<<b
        row[b]^=1<<a
        after=counts(tuple(row))(f)
        assert {before,after}=={target-2,target+2}
    print('Partial template endpoint counts verified.')

if __name__=='__main__':check_cubes()
