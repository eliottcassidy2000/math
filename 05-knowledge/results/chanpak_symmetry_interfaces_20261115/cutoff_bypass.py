from math import comb
from verify_modular import step
def partner(r,op):
    if op==5:return r**4
    return r*sum(comb(r,k)*b for k,b in enumerate((0,5,22,36,24)) if k<=r)
def comp(v,M,op):
    if op in (2,3):return step(v,M,op)
    x=(0,)+v+(0,); b=[partner(r,op)%M for r in range(M+1)]
    return tuple((r*((b[r+1]-b[r])*x[r+1]+(3*b[r]-b[r+1]-b[r-1])*x[r]-(b[r]-b[r-1])*x[r-1]))%M for r in range(1,M))
def cert(ops):
    M=(3,5); w=max(ops); R={i:set() for i in range(2-w,1)}; R[1]={((1,2),(1,2,3,4))}
    full={1:False}; seen={}
    for n in range(2,401):
        R[n]={tuple(comp(v,p,op) for v,p in zip(V,M)) for op in ops for V in R[n-op]}
        full[n]=len({tuple(v[0] for v in V) for V in R[n]})==15
        key=tuple(frozenset(R[j]) for j in range(n-w+1,n+1))
        if key in seen:
            m=seen[key]
            assert all(full[k] for k in range(m,n)),(m,n,full)
            return (m,n,max((k for k in range(1,m) if not full[k]),default=0)+1)
        seen[key]=n
    raise RuntimeError(ops)
if __name__=='__main__':
    for ops,want in (((2,3,6),(83,85,42)),((2,3,5,6),(68,70,30))):
        got=cert(ops)
        assert got==want,(ops,got,want)
        print(ops,got)
