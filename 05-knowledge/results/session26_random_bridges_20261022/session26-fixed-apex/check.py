import itertools, math, json
from fractions import Fraction
from functools import lru_cache
out=[]
for n in range(1,6):
    pairs=list(itertools.combinations(range(n),2)); ids={e:i for i,e in enumerate(pairs)}
    full=(1<<len(pairs))-1
    perms=list(itertools.permutations(range(n)))
    def arc(t,u,v):
        return bool(t>>ids[min(u,v),max(u,v)]&1)==(u<v)
    def trans(t,vs):
        return all(not(arc(t,a,b)==arc(t,b,c)==arc(t,c,a)) for a,b,c in itertools.combinations(vs,3))
    @lru_cache(None)
    def tau(t):
        return next(n-s for s in range(n,-1,-1) if any(trans(t,vs) for vs in itertools.combinations(range(n),s)))
    def strong(t):
        reach=[{i}|{j for j in range(n) if i!=j and arc(t,i,j)} for i in range(n)]
        for k in range(n):
            for i in range(n):
                if k in reach[i]: reach[i]|=reach[k]
        return all(len(r)==n for r in reach)
    def relabel(t,p):
        return sum((1<<i) for i,(u,v) in enumerate(pairs) if arc(t,p[u],p[v]))
    classes={}; checked=0
    for t in range(full+1):
        h=0
        for p in perms:
            if all(arc(t,p[i],p[i+1]) for i in range(n-1)):
                h+=1; checked+=1
                pathmask=sum(1<<ids[min(p[i],p[i+1]),max(p[i],p[i+1])] for i in range(n-1))
                tp=t^full^pathmask
                assert tau(t)+tau(tp)>=n//3
                assert tau(tp)<=(n+tau(t))//2
        c=min(relabel(t,p) for p in perms)
        if c==t:
            aut=sum(relabel(t,p)==t for p in perms)
            classes[t]={'tau':tau(t),'H':h,'aut':aut,'strong':strong(t),'sc':min(relabel(full^t,p) for p in perms)==t}
    row={'n':n,'labelled_tournaments':full+1,'ordinary_classes':len(classes),'path_transforms_checked':checked,'strata':[]}
    for k in range(1,n):
        m=n-k
        subset={t:d for t,d in classes.items() if d['tau']<=k}
        weighted=sum((Fraction(d['H'],d['aut']) for d in subset.values()),Fraction())
        fh=sum((Fraction(sum(trans(t,b) for b in itertools.combinations(range(n),m))*d['H'],d['aut']) for t,d in subset.items()),Fraction())
        # Presentation rows use fixed ascending base 0,...,m-1; core m,...,n-1.
        forced=sum(1<<ids[u,v] for u,v in itertools.combinations(range(m),2))
        free=[i for i,(u,v) in enumerate(pairs) if v>=m]
        L=1<<len(free); sumh=bad=collision=sccbad=strongcount=0
        for bits in range(L):
            t=forced|sum(1<<i for j,i in enumerate(free) if bits>>j&1)
            h=sum(all(arc(t,p[i],p[i+1]) for i in range(n-1)) for p in perms)
            sumh+=h
            f=sum(trans(t,b) for b in itertools.combinations(range(n),m))
            bad+=f>1
            profiles=[tuple(arc(t,u,b) for b in range(m)) for u in range(m,n)]
            collision+=len(set(profiles))<k
            # Reachability criterion for all core vertices lying in one SCC.
            reach=[{i}|{j for j in range(n) if i!=j and arc(t,i,j)} for i in range(n)]
            for v in range(n):
                for u in range(n):
                    if v in reach[u]:reach[u]|=reach[v]
            coregood=all(v in reach[u] for u in range(m,n) for v in range(m,n))
            sccbad+=not coregood
            if m>=2 and coregood:
                endpoint=(any(arc(t,u,0) for u in range(m,n)) and any(arc(t,m-1,u) for u in range(m,n)))
                assert strong(t)==endpoint
            strongcount+=strong(t)
        assert Fraction(sumh,math.factorial(k))==fh
        bound=sum(math.comb(k,r)*math.comb(m,r)*Fraction(m-r+1,2**(m-r))**r for r in range(1,min(k,m)+1))
        assert Fraction(bad,L)<=bound
        assert Fraction(collision,L)<=Fraction(math.comb(k,2),2**m)
        assert Fraction(sccbad,L)<=k*(k-1)*Fraction(3,4)**m
        row['strata'].append({'k':k,'m':m,'classes':len(subset),'strong_classes':sum(d['strong'] for d in subset.values()),'self_converse_classes':sum(d['sc'] for d in subset.values()),'weighted_H':str(weighted),'weighted_fH':str(fh),'presentations':L,'presentation_H_sum':sumh,'nonunique_presentations':bad,'collision_presentations':collision,'core_SCC_bad_presentations':sccbad,'strong_presentations':strongcount})
    out.append(row)
print(json.dumps({'status':'all assertions passed','scope':'All labelled tournaments n=1..5 and every directed Hamilton path; presentation strata 1<=k<n.','results':out},indent=2))
