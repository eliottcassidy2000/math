"""Pure helpers vendored from the frozen first Session58 packet."""
import itertools
def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def pairs(n): return list(itertools.combinations(range(n),2))
def relabel(x,n,p,tournament=False):
    ps=pairs(n); idx={e:k for k,e in enumerate(ps)}; y=0
    for k,(i,j) in enumerate(ps):
        a,b=p[i],p[j]; bit=(x>>k)&1
        if tournament and a>b: bit^=1
        if bit:y|=1<<idx[tuple(sorted((a,b)))]
    return y

