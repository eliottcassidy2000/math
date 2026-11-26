from pathlib import Path
from extra_routes import load
import json,random,sys

def solve(grid,dom,seed,budget):
    m,D=load(grid)
    R=[list(map(int,l.split())) for l in Path(dom).read_text().splitlines()]
    x,n=R[0];R=R[1:];assert len(R)==n
    inv={r[0]:i for i,r in enumerate(R)}
    slots={}
    for i,r in enumerate(R):
        for k in range(1,5):
            v=(k,r[k]);slots[v]=slots.get(v,0)|(1<<i)
    for y in range(m*m):
        if y//m in (0,28) or len(D[y])!=m-1 or y==x:continue
        want=int(y//m!=x//m and y%m!=x%m)-int(y in D[x])-len(D[y]&D[x])
        if want:
            assert want==1
            z=sum(1<<inv[t] for t in D[y] if t in inv)
            slots[(5,y)]=z
    cover=[0]*n; active=[]
    for j,(k,bits) in enumerate(slots.items()):
        active.append(bits)
        z=bits
        while z:
            t=z&-z;z^=t;cover[t.bit_length()-1]|=1<<j
    N=len(active)
    rejects=[]
    for i,(y,*_) in enumerate(R):
        reach=set(D[y])
        for z in D[y]:reach.update(D[z])
        z=sum(1<<inv[t] for t in reach if t in inv)|(1<<i)
        fg=cover[i]
        while fg:
            t=fg&-fg;fg^=t
            z |= active[t.bit_length()-1]
        rejects.append(z)
    rng=random.Random(seed)
    seen=0; best=0; partial=[]
    def walk(U,W,run):
        nonlocal seen,best,partial
        if seen==budget:return None
        seen+=1
        if len(run)>best:best=len(run);partial=run.copy()
        if not U:return run
        K=U; MM=0; d=1<<30
        while K:
            t=K&-K;K^=t
            S=active[t.bit_length()-1]&W
            k=S.bit_count()
            if not k:return None
            if k<d:d=k;MM=S
        options=[]
        while MM:
            t=MM&-MM;MM^=t
            options.append(t.bit_length()-1)
        rng.shuffle(options)
        for i in options:
            Z=walk(U&~cover[i], W&~rejects[i],run+[i])
            if Z is not None:return Z
            if seen==budget:return None
        return None
    solution=walk((1<<N)-1,(1<<n)-1,[])
    out={"vertex":x,"conditions":N,"choices":n,"nodes":seen,"selected":best,"witness":None}
    if solution is not None:
        T=[R[i][0] for i in solution];assert len(T)+len(D[x])==55
        for y in T:
            assert y not in D[x]; D[x].add(y);D[y].add(x)
        for y in D[x]:
            assert not (D[x]&D[y])
            for z in D[y]-{x}:
                assert len(D[x]&D[z])==1
        out["witness"]=T
    print(json.dumps(out))
if __name__=="__main__":solve(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]))
