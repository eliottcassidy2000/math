"""Independent checks of primitive seams and the conflict-atom recursion."""
from collections import Counter
from functools import cache
from itertools import combinations, permutations, product
from math import factorial


@cache
def pairs(n):
    return tuple(combinations(range(n), 2))


@cache
def paths(n):
    return tuple(permutations(range(n)))


def permutation_action(p):
    """Orientation signs of the old arc coordinates under relabeling."""
    ps = pairs(len(p))
    return tuple((ps.index(tuple(sorted((p[a],p[b])))), p[a]>p[b]) for a,b in ps)


def transform(x, action):
    return sum((((x>>i)&1)^rev)<<j for i,(j,rev) in enumerate(action))


@cache
def ham(n, x):
    mat = [0]*n
    for j,(a,b) in enumerate(pairs(n)):
        if (x>>j)&1:
            a,b=b,a
        mat[a] |= 1<<b
    d = [[0]*n for _ in range(1<<n)]
    for a in range(n):
        d[1<<a][a]=1
    for S in range(1,1<<n):
        for a,v in enumerate(d[S]):
            if v:
                A=mat[a]&~S
                while A:
                    low=A&-A
                    d[S|low][low.bit_length()-1]+=v
                    A^=low
    return sum(d[-1])


def cuts(n,x):
    p=pairs(n)
    for S in range(1,(1<<n)-1):
        if all(bool(x>>j&1)==bool(S>>b&1) for j,(a,b) in enumerate(p)
               if (S>>a^S>>b)&1):
            yield S


def induced(n,x,S):
    s=[i for i in range(n) if S>>i&1]
    ps=pairs(n)
    return len(s),sum(((x>>ps.index((s[a],s[b])))&1)<<j
                       for j,(a,b) in enumerate(pairs(len(s))))


def test_seams():
    prim=Counter()
    all_unit=Counter()
    for n in range(2,7):
        acts=[permutation_action(p) for p in paths(n)]
        remaining=set(range(1<<len(pairs(n))))
        while remaining:
            x=min(remaining)
            remaining.difference_update(transform(x,a) for a in acts)
            aut=[a for a in acts if transform(x,a)==x]
            cx=set(cuts(n,x))
            seen=set()
            for k in range(len(pairs(n))):
                if k in seen:
                    continue
                seen.update(a[k][0] for a in aut)
                y=x^(1<<k)
                h,j=ham(n,x),ham(n,y)
                if j!=h+2:
                    continue
                assert all(a[k][0]!=k or a==acts[0] for a in aut)
                cy=set(cuts(n,y))
                common=sorted(cx&cy)
                all_unit[n,h]+=1
                if not common:
                    prim[n,h]+=1
                    if h>=5:
                        assert not cx and not cy
                    continue
                u,v=pairs(n)[k]
                bounds=[0]+common+[(1<<n)-1]
                blocks=[b^a for a,b in zip(bounds,bounds[1:])]
                middle=next(S for S in blocks if S>>u&1)
                assert middle>>v&1
                assert all(S.bit_count()==1 for S in blocks if S!=middle)
                assert ham(*induced(n,x,middle))==h
                assert ham(*induced(n,y,middle))==j
        assert all_unit[n,1]==max(0,n-2)
        assert all_unit[n,3]==2*max(0,n-3)
    assert {k:v for k,v in prim.items() if k[1]<5}=={(3,1):1,(4,3):2}
    for (n,h),v in all_unit.items():
        assert v==sum((n-m+1)*prim[m,h] for m in range(2,n+1))
    print("Primitive seam tests through n=6:", sum(all_unit.values()),
          "unit classes; weak primitive counts", {k:v for k,v in prim.items() if k[1]<5})
    print("First nontrivial primitive n=5 fibers:", sorted((h,v) for (n,h),v in prim.items() if n==5))


def path_atoms(n,g,marked):
    """Return (number of free variables, (prescribed bits, u-weight, v-weight))."""
    coords={}
    for a,b in pairs(n):
        if marked and (a,b)==(0,1):
            continue
        if (a,b) in coords:
            continue
        i=0 if not coords else 1+max(j for j,s in coords.values())
        u,v=a,b
        while tuple(sorted((u,v))) not in coords:
            coords[tuple(sorted((u,v)))]=(i,u>v)
            u,v=g[u],g[v]
        if (u,v)!=(a,b):
            return None
    ps=set(paths(n))
    weights={}
    while ps:
        P=min(ps)
        members=[P]
        Q=tuple(map(g.__getitem__,P))
        while Q!=P:
            members.append(Q);Q=tuple(map(g.__getitem__,Q))
        ps.difference_update(members)
        vals={}
        mu=mv=len(members)
        compatible=True
        for a,b in zip(P,P[1:]):
            q=tuple(sorted((a,b)))
            if marked and q==(0,1):
                if a==0:mv=0
                else:mu=0
                continue
            j,s=coords[q]
            val=(a>b)^s
            if j in vals and vals[j]!=val:
                compatible=False
                break
            vals[j]=val
        if compatible:
            k=tuple(sorted(vals.items()))
            pu,pv=weights.get(k,(0,0))
            weights[k]=(pu+mu,pv+(mv if marked else 0))
    m=0 if not coords else 1+max(j for j,s in coords.values())
    return m,tuple(weights.items())


def polynomial(n,g,marked):
    A=path_atoms(n,g,marked)
    if A is None:return {}
    m,atoms=A

    @cache
    def go(i,D):
        if i==len(atoms):
            return {(0,0):1<<(m-len(D))}
        k,(u,v)=atoms[i]
        drop=go(i+1,D)
        left,right=dict(D),dict(k)
        if any(a in left and left[a]!=b for a,b in k):
            return drop
        keep=go(i+1,tuple(sorted((left|right).items())))
        res=Counter(drop)
        for (a,b),c in keep.items():
            res[a+u,b+v]+=c
            res[a,b]-=c
        return {t:c for t,c in res.items() if c}

    return go(0,())


def test_conflict():
    for n in range(2,5):
        group=tuple(permutations(range(n)))
        for marked in (False,True):
            gp=[g for g in group if not marked or g[:2]==(0,1)]
            total=Counter()
            for g in gp:
                direct=Counter()
                act=permutation_action(g)
                for x in range(1<<len(pairs(n))):
                    if (marked and x&1) or transform(x,act)!=x:
                        continue
                    direct[ham(n,x),ham(n,x^1) if marked else 0]+=1
                value=polynomial(n,g,marked)
                assert dict(direct)==value, (n,marked,g,direct,value)
                total.update(value)
            assert all(v%len(gp)==0 for v in total.values())
    print("Both path-conflict transforms agree permutation-by-permutation through n=4.")


if __name__=='__main__':
    test_seams()
    test_conflict()
