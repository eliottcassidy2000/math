import itertools as it, json, collections as co, fractions as fr
from functools import lru_cache

def tiles(n): return [(x,y) for x in range(n) for y in range(x-1)]
def adj(n,t):
    a=[0]*n
    for x in range(1,n): a[x]|=1<<(x-1)
    for k,(x,y) in enumerate(tiles(n)):
        a[x if (t>>k)&1 else y]|=1<<(y if (t>>k)&1 else x)
    return a

def canonical(a):
    n=len(a); sc=[r.bit_count() for r in a]; gr=co.defaultdict(list)
    for v,s in enumerate(sc):gr[s].append(v)
    best=None
    for pp in it.product(*(it.permutations(gr[s]) for s in sorted(gr))):
        p=sum(pp,()); code=0
        for x in range(1,n):
            for y in range(x):code=(code<<1)|((a[p[x]]>>p[y])&1)
        if best is None or code<best:best=code
    return best

def converse(a):
    n=len(a);return [((1<<n)-1)^r^(1<<v) for v,r in enumerate(a)]
def c3(a):return len(a)*(len(a)-1)*(len(a)-2)//6-sum(r.bit_count()*(r.bit_count()-1)//2 for r in a)
def make(n):
    m=len(tiles(n)); N=1<<m; keys={}; pi=[]; reps=[]; cs=[]; scores=[]; scflag=[]
    for t in range(N):
        a=adj(n,t); cc=canonical(a); oc=canonical(converse(a)); key=min(cc,oc)
        if key not in keys:
            keys[key]=len(keys); reps.append(t);cs.append(c3(a)); scores.append(sorted(r.bit_count() for r in a));scflag.append(cc==oc)
        pi.append(keys[key])
    q=len(keys); A=[[[0]*q for _ in range(q)] for c in range(2)]; rows=[co.Counter() for _ in range(q)]; refl=[]
    tt=tiles(n); perm=[tt.index((n-1-y,n-1-x)) for x,y in tt]
    for t in range(N):
        r=sum(((t>>k)&1)<<v for k,v in enumerate(perm));refl.append(r)
        u,v=pi[t],pi[t^(N-1)];blue=int(r==t)
        A[blue][u][v]+=1;rows[u][cs[v]-cs[u]]+=1
        a=adj(n,t); assert cs[v]-cs[u]==a[-1].bit_count()-a[0].bit_count()-1
    return dict(n=n,m=m,N=N,q=q,pi=pi,reps=reps,c3=cs,scores=scores,selfconverse=scflag,Ablack=A[0],Ablue=A[1],D=list(co.Counter(pi)[u] for u in range(q)),reflection=refl)
if __name__=='__main__':
    import pathlib
    for n in range(3,8):
        d=make(n); pathlib.Path(f'n{n}.json').write_text(json.dumps(d,indent=2));print(n,d['q'],flush=True)
