def variables(n,q,anti=True):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    info={};orbit=[]
    for p in pairs:
        if p in info:continue
        comp={p:0};todo=[p]
        while todo:
            u,v=todo.pop();s=comp[u,v]
            trans=[((u+q)%n,(v+q)%n,False)]
            if anti:trans.append(((n-1-u)%n,(n-1-v)%n,True))
            for a,b,rev in trans:
                if a>b:a,b=b,a;rev=not rev
                t=s^rev
                if (a,b) in comp:
                    assert comp[a,b]==t,('bad',n,q,p)
                else:comp[a,b]=t;todo.append((a,b))
        k=len(orbit);orbit.append(comp)
        for p,s in comp.items():info[p]=(k,s)
    forced={}
    for i in range(n):
        j=(i+1)%n;p=tuple(sorted((i,j)));k,s=info[p]
        t=(i<j)^s
        assert k not in forced or forced[k]==t
        forced[k]=t
    return orbit,forced
if __name__=='__main__':
    for n,q in [(10,2),(14,2),(15,3),(21,3),(12,4),(20,4),(28,4)]:
        for anti in (True,False):
            orb,f=variables(n,q,anti)
            print(n,q,'anti',anti,'orbits',len(orb),'forced',len(f),'free',len(orb)-len(f))
