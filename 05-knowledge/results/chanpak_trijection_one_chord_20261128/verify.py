"""Independent literal checks of the marked-tree and one-chord constructions."""
from itertools import combinations, product
from pathlib import Path
from collections import Counter, deque
from sys import argv

def make(w):
    n=len(w)+1
    T=[[False]*n for _ in range(n)]
    F=[]; root=0
    for b,letter in enumerate(w,1):
        for a in range(b):
            T[a][b]=(a!=root or letter=="P")
            T[b][a]=not T[a][b]
        F.append((root,b))
        if letter=="M":root=b
    return T,F

def path(n,F,a,b):
    G=[[] for _ in range(n)]
    for i,(x,y) in enumerate(F):
        G[x].append((y,i));G[y].append((x,i))
    q=deque([a]);seen={a:None}
    while q:
        x=q.popleft()
        if x==b:break
        for y,k in G[x]:
            if y not in seen:seen[y]=(x,k);q.append(y)
    result=[]
    while b!=a:
        b,k=seen[b];result.append(k)
    return result[::-1]

def reduced(T,F):
    L={}
    for a,b in combinations(range(len(T)),2):
        p=path(len(T),F,a,b)
        if len(p)==1:continue
        x,y=p[0],p[-1]
        assert (min(x,y),max(x,y)) not in L
        L[min(x,y),max(x,y)]=T[a][b] if x<y else T[b][a]
    assert len(L)==len(F)*(len(F)-1)//2
    return L

def count(T):
    n=len(T); inc=[sum(1<<i for i in range(n) if T[i][j]) for j in range(n)]
    b=[[0]*n for _ in range(1<<n)]
    for s in range(1,1<<n):
        for j in range(n):
            if s>>j&1:
                r=s^(1<<j)
                b[s][j]=sum(b[r][i] for i in range(n) if (r&inc[j])>>i&1) if r else 1
    return sum(b[-1])

def degrees(n,F):
    r=Counter()
    for a,b in F:r[a]+=1;r[b]+=1
    return r

if __name__=='__main__':
    checks=0
    for l in range(7):
        for w in product("PKM",repeat=l):
            T,F=make(w)
            assert all(reduced(T,F).values())
            checks+=1
    print("transitive compressed controls",checks)
    fin=Path(argv[1]) if len(argv)>1 else Path(__file__).with_name('witnesses.txt')
    types=Counter(); cycle_lengths=Counter(); reductions=Counter(); results=set()
    for line in fin.read_text().splitlines():
        a,w,i,j,old=line.split()
        h,i,j,old=map(int,(a,i,j,old)); T,F=make(w)
        n=len(T)
        p=path(n,F,i,j); assert len(p)>=2 and T[i][j]
        assert count(T)==old
        R=reduced(T,F); assert all(R.values())
        T[i][j],T[j][i]=False,True
        assert count(T)==h
        diff={edge for edge,value in reduced(T,F).items() if value!=R[edge]}
        expected={tuple(sorted((p[0],p[-1])))}
        assert diff==expected
        (r,s),=expected
        reductions['adjacent' if s==r+1 else 'nonadjacent']+=1
        d=degrees(n,F)
        b=next(k for k,(a,b) in enumerate(F) if k not in p and (d[a]==1 or d[b]==1))
        F2=[e for k,e in enumerate(F) if k!=b]+[(i,j)]
        D=degrees(n,F2)
        assert len(F2)==n-1 and len(D)==n-1   # exactly one isolate
        active=list(D)
        for x in active:path(n,F2,active[0],x)
        cycle_lengths[len(p)+1]+=1
        types[n]+=1
        results.add(h)
    assert len(results)==27
    print("Independently verified repaired heights",len(results))
    print("Orders",sorted(types.items()),"cycle lengths",sorted(cycle_lengths.items()),
          "one compressed flip",sorted(reductions.items()))
    print("Every witness gives unicyclic component plus isolate, with no padding.")
