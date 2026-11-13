from itertools import permutations, combinations
from collections import Counter
from math import lcm

def require(ok):
    if not ok: raise ValueError('control failed')

def parity(x): return x.bit_count()&1

def run(n):
    es=list(combinations(range(n),2))
    zs=[(i,j) for i,j in es if j-i>1]
    m=len(zs)
    perm=tuple(zs.index((n-1-j,n-1-i)) for i,j in zs)
    def sigma(x): return sum(((x>>i)&1)<<j for i,j in enumerate(perm))
    X=list(range(1<<m)); F=[x for x in X if x==sigma(x)]
    reps=[x for x in X if x==min(x^f for f in F)]
    def canon(x):
        rows={}
        for i,j in es:
            a=(j==i+1) or ((x>>zs.index((i,j)))&1)
            rows[i,j]=bool(a); rows[j,i]=not a
        return min(sum(rows[p[i],p[j]]<<a for a,(i,j) in enumerate(es)) for p in permutations(range(n)))
    c=[canon(x) for x in X]; ds=Counter(c)
    cs=sorted(ds); N=len(cs); idx={v:i for i,v in enumerate(cs)}
    cls=[idx[i] for i in c]
    rho={}
    for x in X:
        a,b=cls[x],cls[sigma(x)]
        require(a not in rho or rho[a]==b)
        rho[a]=b
    def zero(): return [[0]*N for _ in cs]
    As={}
    for f in F:
        M=zero()
        for x in X: M[cls[x]][cls[x^f]]+=1
        As[f]=M
        require(all(M[a][b]==M[b][a]==M[rho[a]][rho[b]] for a in range(N) for b in range(N)))
    any_noncommute=False
    L=lcm(*ds.values())
    mass=[L//ds[a] for a in cs]
    for f in F:
        for g in F:
            a,b=As[f],As[g]
            if any(sum((a[i][k]*b[k][j]-b[i][k]*a[k][j])*mass[k] for k in range(N)) for i in range(N) for j in range(N)):
                if not any_noncommute: print('noncommuting translation pair',n,zs,f,g)
                any_noncommute=True
                break
        if any_noncommute: break
    for k in {tuple(parity(k&f) for f in F):k for k in X}.values():
        G=zero()
        pos=neg=0
        for x in reps:
            d=x^sigma(x)
            sign=(-1)**parity(k&d)
            pos+=sign==1; neg+=sign==-1
            b=[0]*N
            for f in F: b[cls[x^f]]+=(-1)**parity(k&f)
            require(all(b[rho[i]]==sign*b[i] for i in range(N)))
            for i in range(N):
                for j in range(N): G[i][j]+=b[i]*b[j]
        require((pos,neg) in [(len(reps),0),(len(reps)//2,len(reps)//2)])
        for i in range(N):
            for j in range(N):
                require(G[i][j]==sum((-1)**parity(k&f)*As[f][i][j] for f in F))
    print(n,m,len(F),len(reps),N,'noncommute',any_noncommute)

if __name__=='__main__':
    for n in range(3,7): run(n)
