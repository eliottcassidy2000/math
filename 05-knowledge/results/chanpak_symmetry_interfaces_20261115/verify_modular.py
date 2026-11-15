from math import gcd

def step(state,M,mode):
    s=(0,)+state+(0,)
    if mode==2:
        return tuple((r*(s[r+1]+r*s[r]-s[r-1]))%M for r in range(1,M))
    return tuple((r*((2*r+1)*s[r+1]+(r*r-2)*s[r]-(2*r-1)*s[r-1]))%M for r in range(1,M))

def cert(M,max_n=100):
    R={-1:set(),0:set(),1:{tuple(range(1,M))}}
    seen={}; full={1:False}
    for n in range(2,max_n+1):
        R[n]={step(v,M,t) for t in (2,3) for v in R[n-t]}
        full[n]=len({v[0] for v in R[n]})==M
        key=tuple(frozenset(R[j]) for j in range(n-2,n+1))
        if key in seen:
            m=seen[key]
            assert all(full[k] for k in range(m,n))
            threshold=max((k for k in range(1,m) if not full[k]),default=0)+1
            return (m,n,threshold)
        seen[key]=n
    raise RuntimeError((M,max_n))

def joint15():
    mods=(3,5); R={-1:set(),0:set(),1:{tuple(tuple(range(1,m)) for m in mods)}}
    for n in range(2,149):
        R[n]={tuple(step(v,m,k) for v,m in zip(V,mods)) for k in (2,3) for V in R[n-k]}
    assert tuple(R[j] for j in (144,145,146))==tuple(R[j] for j in (146,147,148))
    for n in (146,147):
        project={(V[0][0],V[1][0]) for V in R[n]}
        missing=[a for a in range(15) if (a%3,a%5) not in project]
        assert missing==([5] if n%2 else [10])
    return (146,148,10,5)

if __name__=='__main__':
    wanted={3:(14,16),5:(27,29),7:(42,44),9:(52,54)}
    for M,w in wanted.items():
        ans=cert(M)
        assert ans[:2]==w,(M,ans,w)
        print(M,*ans)
    print('joint mod15 repeated endpoints and missing even-/odd-order residues:',joint15())
