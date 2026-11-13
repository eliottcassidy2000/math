from itertools import permutations, combinations
from collections import Counter

def run(n):
    es=list(combinations(range(n),2))
    zs=[(i,j) for i,j in es if j-i>1]
    m=len(zs); X=range(1<<m)
    perm=tuple(zs.index((n-1-j,n-1-i)) for i,j in zs)
    def sigma(x): return sum(((x>>i)&1)<<j for i,j in enumerate(perm))
    def canon(x):
        rows={}
        for i,j in es:
            a=(j==i+1) or ((x>>zs.index((i,j)))&1)
            rows[i,j]=bool(a); rows[j,i]=not a
        return min(sum(rows[p[i],p[j]]<<a for a,(i,j) in enumerate(es)) for p in permutations(range(n)))
    cls=[canon(x) for x in X]
    A=Counter((cls[x],cls[x^((1<<m)-1)]) for x in X)
    rho={cls[x]:cls[sigma(x)] for x in X}
    ys=sorted({min(a,b) for a,b in rho.items()})
    B=Counter((min(a,rho[a]),min(b,rho[b])) for x in X
              for a,b in [(cls[x],cls[x^((1<<m)-1)])])
    def refine(key):
        d={v:i for i,v in enumerate(dict.fromkeys(key.values()))}
        return {v:d[key[v]] for v in ys}
    colors=refine({v:(sum(B[v,w] for w in ys),) for v in ys})
    rounds=0
    while True:
        new=refine({v:(colors[v],tuple(sorted((colors[w],B[v,w]) for w in ys if B[v,w])))
                    for v in ys})
        rounds+=1
        if len(set(new.values()))==len(set(colors.values())): break
        colors=new
    ch={i for i in ys if rho[i]!=i}
    seen=set(); sizes=[]
    for i in sorted(ch):
        if i in seen:continue
        comp={i}; stack=[i]
        while stack:
            k=stack.pop()
            for l in ch-comp:
                if A[k,l]!=A[k,rho[l]]:comp.add(l);stack.append(l)
        seen|=comp; sizes.append(len(comp))
    out=(len(ys),len(set(colors.values())),len(ch),tuple(sizes))
    want={4:(3,3,1,(1,)),5:(10,10,2,(1,1)),6:(34,34,22,(22,))}
    if out!=want[n]:raise ValueError((n,out,want[n]))
    print(n,out,"refinement rounds",rounds)
for n in (4,5,6):run(n)
