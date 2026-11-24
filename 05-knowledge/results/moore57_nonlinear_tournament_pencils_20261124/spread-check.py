"""Exact PG(3,2) pointed-spread CSP and PG(3,7) local regular-spread catalog."""
from itertools import combinations, product
from collections import Counter

def line_catalog_2():
    return sorted({frozenset((a,b,a^b)) for a,b in combinations(range(1,16),2)}, key=lambda l:sorted(l))

def solve_2():
    L=line_catalog_2()
    spreads=[frozenset(s) for s in combinations(range(35),5)
             if len(set().union(*(L[j] for j in s)))==15]
    domains=[tuple(s for s in spreads if i in s) for i in range(35)]
    assert len(spreads)==56 and set(map(len,domains))=={8}
    full=(1<<8)-1
    pair={}
    for i,j in product(range(35),repeat=2):
        if i==j: continue
        masks=[]
        for S in domains[i]:
            m=0
            for b,T in enumerate(domains[j]):
                cross=j in S
                if cross==(i in T) and ((S&T)=={i,j} if cross else len(S&T)<=1):
                    m|=1<<b
            masks.append(m)
        pair[i,j]=masks
    def propagate(D,queue):
        while queue:
            i=queue.pop()
            for j in range(35):
                if i==j: continue
                m=0
                for a in range(8):
                    if D[i]>>a&1: m|=pair[i,j][a]
                d=m&D[j]
                if not d: return False
                if d!=D[j]:
                    D[j]=d
                    queue.append(j)
        return True
    calls=0
    def search(D):
        nonlocal calls
        calls+=1
        if max(map(int.bit_count,D))==1: return D
        i=min((i for i in range(35) if D[i].bit_count()>1),key=lambda i:D[i].bit_count())
        for a in range(8):
            if not D[i]>>a&1: continue
            E=D.copy()
            E[i]=1<<a
            if propagate(E,[i]):
                ans=search(E)
                if ans: return ans
        return None
    ans=search([full]*35)
    assert ans
    selected=[domains[i][v.bit_length()-1] for i,v in enumerate(ans)]
    assert len(set(selected))==35
    H=[P-{i} for i,P in enumerate(selected)]
    G=[{j+15 for j,l in enumerate(L) if i in l} for i in range(1,16)]
    G += [{p-1 for p in l}|{j+15 for j in H[i]} for i,l in enumerate(L)]
    assert all(len(x)==7 for x in G)
    assert all(len(G[i]&G[j])==int(j not in G[i]) and (j in G[i])==(i in G[j])
               for i,j in combinations(range(50),2))
    print("PG(3,2):",len(spreads),"spreads; 8 per line; CSP search nodes",calls,
          "; selected 35 distinct; Moore certificate passed")
    profiles=Counter()
    binary=[]
    for S in spreads:
        d=[len(T&S) for T in H]
        e=sum(d[z] for z in S)//2
        profiles[e,max(d)]+=1
        if e==0 and max(d)<=2: binary.append(S)
    def resolutions(left):
        if not left:
            yield ()
            return
        for S in spreads:
            if min(left) in S and S<=left:
                for more in resolutions(left-S):
                    yield (S,)+more
    res=list(resolutions(set(range(35))))
    assert len(res)==240
    assert Counter(sum(S in selected for S in r) for r in res)=={4:210,7:30}
    covers=sum(all(S in binary for S in r) for r in res)
    print("HS resolutions and selected spreads per resolution:",len(res),dict(Counter(sum(S in selected for S in r) for r in res)))
    print("Frame profiles (internal edges, maximum contacts):",sorted(profiles.items()))
    print("Independent binary frames and full resolutions:",len(binary),covers)
    direction_profiles=Counter()
    for i,U in enumerate(L):
        # Coordinates for a line skew to U: binary graph of a map from a
        # complement W into U, packed as two 2-bit columns.
        a,b=sorted(U)[:2]
        c=next(x for x in range(1,16) if x not in U)
        d=next(x for x in range(1,16) if x not in ({0}|U|{c}|{c^y for y in U}))
        dec={((a if v&1 else 0)^(b if v&2 else 0)^(c if v&4 else 0)^(d if v&8 else 0)):v for v in range(16)}
        def coord(j):
            l={dec[x]//4:dec[x]%4 for x in L[j]}
            assert set(l)=={1,2,3}
            return l[1]+4*l[2]
        flats=[{coord(k) for k in selected[j] if k!=i} for j in H[i]]
        assert len(set().union(*flats))==16
        dirs={frozenset(x^min(F) for x in F) for F in flats}
        common=set.intersection(*(set(R) for R in dirs))
        generated={0}
        for R in dirs:
            for v in R:
                generated |= {w^v for w in generated}
        direction_profiles[len(dirs),len(common),len(generated)]+=1
    print("HS local direction profiles (#directions, common intersection size, span size):",dict(direction_profiles))
    anchors=Counter()
    for i,j in combinations(range(35),2):
        if j not in H[i]: continue
        A=H[i]-{j}; B=H[j]-{i}; rest=set(range(35))-selected[i]-selected[j]
        cross=[(a,b) for a in A for b in B if L[a].isdisjoint(L[b])]
        profiles2=Counter((len(H[x]&A),len(H[x]&B)) for x in rest)
        assert len(cross)==profiles2[1,1]==0
        grid={x for x in rest if H[x].isdisjoint(A|B)}
        assert len(grid)==9 and all(not L[x].isdisjoint(L[i]) and not L[x].isdisjoint(L[j]) for x in grid)
        anchors[sum(len(H[x]&grid) for x in grid)//2]+=1
    print("HS grid edge counts across 70 anchors:",dict(anchors))
    return len(binary),covers

def regular_spreads_pair(q=7):
    I=(1,0,0,1)
    Z=(0,0,0,0)
    def mul(A,B):
        a,b,c,d=A
        e,f,g,h=B
        return ((a*e+b*g)%q,(a*f+b*h)%q,(c*e+d*g)%q,(c*f+d*h)%q)
    def det(A):
        a,b,c,d=A
        return (a*d-b*c)%q
    GL={A for A in product(range(q),repeat=4) if det(A)}
    fields=set()
    for A in GL:
        a,b,c,d=A
        if any((x*x-(a+d)*x+det(A))%q==0 for x in range(q)): continue
        F=frozenset(tuple((s*i+t*v)%q for i,v in zip(I,A)) for s,t in product(range(q),repeat=2))
        fields.add(F)
    chart=set()
    for F in fields:
        rem=GL.copy()
        while rem:
            A=min(rem)
            AF=frozenset(mul(A,B) for B in F)
            chart.add(AF)
            rem-=AF
    print(f"q={q}: |GL2|={len(GL)}, field embeddings={len(fields)}, regular spreads through a fixed skew pair={len(chart)}")
    for S in chart:
        assert len(S)==q*q and Z in S
        assert all(det(tuple((a-b)%q for a,b in zip(A,B))) for A,B in combinations(S,2))
    n=Counter(A for S in chart for A in S if A!=Z)
    print("Third-line incidences on the chart:",Counter(n.values()))
    if q==7:
        assert len(fields)==21 and len(chart)==882 and Counter(n.values())=={21:2016}
        # All nonzero elements invertible and chart planes 2-dimensional.
        S=min(chart,key=lambda T:sorted(T))
        T=sorted((U for U in chart if U&S=={Z}),key=lambda U:sorted(U))
        num=[sum(len(U&V)==q for V in T) for U in T]
        assert len(T)==721
        pairs=sum(num)//2
        print("Transverse directions and compatible pairs:",len(T),pairs,"partner histogram",sorted(Counter(num).items()))
        print("Exactly two-direction partitions for a fixed central chart:",pairs*(2**q-2))
        U,V=next((U,V) for U,V in combinations(T,2) if len(U&V)==q)
        def add(x,y):return tuple((a+b)%q for a,b in zip(x,y))
        W={add(u,v) for u in U for v in V}
        assert len(W)==q**3
        # Use U in the zero slice, V elsewhere. Each point of S
        # points to one of the resulting disjoint affine planes.
        flats=[{add(y,x) for x in (U if y in W else V)} for y in S]
        assert len(flats)==q*q and len(set().union(*flats))==q**4
        assert all(len(F&S)==1 for F in flats)
        # In M_2(F_7)/<I>, normalize projective trace-zero coordinates.
        P=[v for v in product(range(q),repeat=3) if next((x for x in v if x),0)==1]
        square={x*x%q for x in range(q)}
        inner=[(a,b,c) for a,b,c in P if (a*a+b*c)%q not in square]
        inside=Counter(sum((s*a+t*b+u*c)%q==0 for a,b,c in inner) for s,t,u in P)
        assert len(P)==57 and len(inner)==21 and max(inside)==4
        print("Regular directions in a line pencil, and projective line loads:",len(inner),sorted(inside.items()))



if __name__=="__main__":
    solve_2()
    regular_spreads_pair(2)
    regular_spreads_pair(7)
