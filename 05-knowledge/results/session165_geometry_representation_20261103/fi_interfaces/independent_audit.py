from itertools import combinations,permutations

def need(ok,label):
    if not ok:raise RuntimeError(label)

def comps(n,r):
    if r==1:
        yield(n,);return
    for k in range(n+1):
        for c in comps(n-k,r-1):yield(k,)+c

def covers(adj,vs,r):
    out=set()
    for p in permutations(vs):
        for c in comps(len(vs),r):
            pos=0;paths=[]
            for k in c:
                paths.append(p[pos:pos+k]);pos+=k
            if all(all(b in adj[a] for a,b in zip(path,path[1:])) for path in paths):out.add(tuple(paths))
    return out

def slot(c,q):return next(t for t,path in enumerate(c) if q in path)
def starts(c,q):return c[slot(c,q)][0]==q
def before(c,t,a):return c[:t]+((a,)+c[t],)+c[t+1:]
def after(c,q,a):
    t=slot(c,q);path=c[t];j=path.index(q)+1
    return c[:t]+(path[:j]+(a,)+path[j:],)+c[t+1:]
def exceptional(c,rho,q,a):
    t=slot(c,rho);return c[:t]+((q,a)+c[t],)+c[t+1:]
def check_image(images,expected,label):
    need(len(images)==len(set(images)),label+' injectivity')
    need(set(images)==expected,label+' surjectivity')

configurations=0;instances=0
for n in range(2,5):
    es=list(combinations(range(n),2))
    for mask in range(1<<len(es)):
        adj=[set() for _ in range(n)]
        for k,(a,b) in enumerate(es):
            if mask>>k&1:adj[a].add(b)
            else:adj[b].add(a)
        for q in range(n):
            inc=[v for v in range(n) if q in adj[v]]
            if len(inc)!=1:continue
            rho=inc[0];a=n;out=[set(x) for x in adj]+[set(range(n))-{q}];out[q].add(a)
            configurations+=1
            for r in range(1,4):
                old=covers(adj,tuple(range(n)),r)
                S={c for c in old if starts(c,q)};U=old-S
                d={c for c in covers(adj,tuple(v for v in range(n) if v!=q),r) if starts(c,rho)}
                new=covers(out,tuple(range(n+1)),r)
                cs=[before(c,t,a) for c in S for t in range(r) if t!=slot(c,q)]+[before(c,t,a) for c in U for t in range(r)]
                cu=[after(c,q,a) for c in old]+[exceptional(c,rho,q,a) for c in d]
                check_image(cs,{c for c in new if starts(c,a)},'C S')
                check_image(cu,{c for c in new if not starts(c,a)},'C U')
                check_image(list(S),{c for c in covers(out,tuple(range(n)),r) if starts(c,q)},'C d')
                ds=[after(c,q,a) if t==slot(c,q) else before(c,t,a) for c in S for t in range(r)]+[exceptional(c,rho,q,a) for c in d]
                du=[before(c,t,a) for c in U for t in range(r)]+[after(c,q,a) for c in U]
                dd=[before(c,t,a) for c in d for t in range(r) if t!=slot(c,rho)]
                check_image(ds,{c for c in new if starts(c,q)},'D S')
                check_image(du,{c for c in new if not starts(c,q)},'D U')
                check_image(dd,{c for c in covers(out,tuple(v for v in range(n+1) if v!=q),r) if starts(c,rho)},'D d')
                instances+=6
print({'guarded_configurations':configurations,'six_coordinate_bijection_checks':instances,'slot_sizes':[1,2,3],'status':'PASS'})
