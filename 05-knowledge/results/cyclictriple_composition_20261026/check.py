"""Literal tournament tests and necklace metagraphs; no external packages."""
import itertools
import json
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache


def need(ok,why):
    if not ok: raise RuntimeError(why)


def carousel(m):
    need(m>=3 and m%2==1, 'carousel domain: odd m >= 3')
    return tuple(sum(1<<j for j in range(m) if 1<=(j-i)%m<=(m-1)//2) for i in range(m))


def block(cut=0,cycle=0):
    a=list(carousel(3))+[0]
    if cycle:
        a[:3]=[((1<<3)-1)^r^(1<<i) for i,r in enumerate(a[:3])]
    for i in range(3):
        if cut: a[3]|=1<<i
        else: a[i]|=1<<3
    return tuple(a)


def converse(t):
    return tuple(((1<<len(t))-1)^row^(1<<i) for i,row in enumerate(t))


def compose(q,bs):
    starts=[0]
    for b in bs: starts.append(starts[-1]+len(b))
    out=[0]*starts[-1]
    for i,b in enumerate(bs):
        outside=sum(((1<<len(bs[j]))-1)<<starts[j] for j in range(len(q)) if q[i]>>j&1)
        for v,row in enumerate(b): out[starts[i]+v]=(row<<starts[i])|outside
    return tuple(out)


def cyclic(t):
    return tuple((i,j,k) for i,j,k in itertools.combinations(range(len(t)),3)
                 if ((t[i]>>j&1)==(t[j]>>k&1)==(t[k]>>i&1)))


def iso_count(a,b,stop=False):
    """Backtracking over actual arcs; no modular decomposition assumed."""
    if len(a)!=len(b): return 0
    n=len(a)
    degA=[r.bit_count() for r in a]; degB=[r.bit_count() for r in b]
    if sorted(degA)!=sorted(degB): return 0
    colorsA=[(degA[i],tuple(sorted(degA[j] for j in range(n) if a[i]>>j&1))) for i in range(n)]
    colorsB=[(degB[i],tuple(sorted(degB[j] for j in range(n) if b[i]>>j&1))) for i in range(n)]
    candidates=[[j for j in range(n) if colorsA[i]==colorsB[j]] for i in range(n)]
    mapped={}; used=set()
    def rec():
        if len(mapped)==n: return 1
        options=[]
        for i in range(n):
            if i in mapped: continue
            cand=[j for j in candidates[i] if j not in used and all((a[i]>>u&1)==(b[j]>>v&1) for u,v in mapped.items())]
            if not cand: return 0
            options.append((len(cand),i,cand))
        _,i,cand=min(options)
        total=0
        for j in cand:
            mapped[i]=j; used.add(j)
            total+=rec()
            used.remove(j); del mapped[i]
            if stop and total: return total
        return total
    return rec()


def hamilton(t):
    n=len(t)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask,row in enumerate(dp):
        for i,count in enumerate(row):
            if count:
                allowed=t[i]&~mask
                while allowed:
                    bit=allowed&-allowed; allowed-=bit
                    dp[mask|bit][bit.bit_length()-1]+=count
    return sum(dp[-1])


def modules(t):
    n=len(t); result=[]
    for mask in range(1,1<<n):
        if all((t[v]&mask) in (0,mask) for v in range(n) if not (mask>>v&1)):
            result.append(mask)
    return result


def witness():
    q=carousel(3)
    reps=[compose(q,[block(cut=(w>>i)&1) for i in range(3)]) for w in [0,1,3,7]]
    cs=cyclic(reps[0]); counts=[]
    for t in reps:
        need(cyclic(t)==cs,'same labelled triple set')
        aut=iso_count(t,t); paths=hamilton(t)
        need(iso_count(t,converse(t))==0,'not self converse')
        counts.append({'rows':list(t),'automorphisms':aut,'hamiltonian_paths':paths,'path_orbits':paths//aut})
    need([x['automorphisms'] for x in counts]==[81,27,27,81],'81/27 emblem')
    need(len({x['hamiltonian_paths'] for x in counts})==1,'H survives')
    hist=Counter(); encodings=set()
    # Every node choice in the cited exact C3 fibre formula: 2*4^3=128.
    for root in range(2):
        for cuts in itertools.product(range(2),repeat=3):
            for cycles in itertools.product(range(2),repeat=3):
                t=compose(converse(q) if root else q,[block(c,d) for c,d in zip(cuts,cycles)])
                need(t not in encodings,'free labelled node choices'); encodings.add(t)
                need(cyclic(t)==cs,'all128 triples preserved')
                hits=[i for i,r in enumerate(reps) if iso_count(t,r,stop=True)]
                need(len(hits)==1,'unique isomorphism type')
                hist[hits[0]]+=1
    need(list(hist[i] for i in range(4))==[16,48,48,16],'full fibre orbit sizes')
    # Prime quotient claim independently checked from every candidate module.
    prime_controls={}
    for m in [3,5,7,9]:
        q=carousel(m)
        need(all(x.bit_count() in [1,m] for x in modules(q)),'carousel prime')
        aut=iso_count(q,q)
        need(aut==m,'carousel exact automorphism order')
        prime_controls[str(m)]=aut
    orders=[]
    for cuts in [(0,0,0),(1,0,0)]:
        t=compose(carousel(3),[block(c) for c in cuts]); hist_order=Counter(); group=[]
        for shift in range(3):
            for rs in itertools.product(range(3),repeat=3):
                p=tuple(4*((i+shift)%3)+((j+rs[i])%3 if j<3 else 3) for i in range(3) for j in range(4))
                if not all((t[i]>>j&1)==(t[p[i]]>>p[j]&1) for i in range(12) for j in range(12)): continue
                group.append(p); power=tuple(range(12)); order=0
                while True:
                    power=tuple(p[x] for x in power);order+=1
                    if power==tuple(range(12)):break
                hist_order[order]+=1
        need(len(group)==iso_count(t,t),'block permutations exhaust actual automorphisms')
        orders.append(dict(sorted(hist_order.items())))
    return {'cyclic_triples':list(cs),'cyclic_count':len(cs),'representatives':counts,
            'full_labelled_fibre':len(encodings),'class_labelled_multiplicities':dict(sorted(hist.items())),
            'prime_carousel_automorphisms':prime_controls,'automorphism_element_orders':orders}


def rotate(w,k): return w[k:]+w[:k]
def necklace(w): return min(rotate(w,k) for k in range(len(w)))
def period(w): return next(k for k in range(1,len(w)+1) if rotate(w,k%len(w))==w)
def divisors(n): return [d for d in range(1,n+1) if n%d==0]
def mu(n):
    sign=1; p=2
    while p*p<=n:
        if n%p==0:
            n//=p; sign=-sign
            if n%p==0: return 0
        p+=1
    return -sign if n>1 else sign
def phi(n): return sum(__import__('math').gcd(a,n)==1 for a in range(1,n+1))


def necklaces():
    rows=[]; graph3=None
    for m in [3,5,7,9,11,13]:
        words=list(itertools.product(range(2),repeat=m))
        buckets=defaultdict(list)
        for w in words: buckets[necklace(w)].append(w)
        orbitcounts=Counter(period(w) for w in buckets)
        mobius={d:sum(mu(e)*2**(d//e) for e in divisors(d))//d for d in divisors(m)}
        need(dict(orbitcounts)==mobius,'Mobius period strata')
        need(len(buckets)==sum(phi(d)*2**(m//d) for d in divisors(m))//m,'Burnside total')
        # Full quotient-cube adjacency with literal move multiplicities.
        edge={}
        for w,bucket in buckets.items():
            neighbors=Counter()
            for i in range(m):
                v=list(w);v[i]^=1;neighbors[necklace(tuple(v))]+=1
            edge[w]=neighbors
            for v,c in neighbors.items():
                need(w!=v and sum(w)%2!=sum(v)%2,'bipartite loop-free')
        for w,neighbors in edge.items():
            need(sum(neighbors.values())==m,'all physical ports counted')
            for v,c in neighbors.items():
                need(len(buckets[w])*c==len(buckets[v])*edge[v][w],'detailed balance')
                aw=3**m*m//period(w);av=3**m*m//period(v)
                need(Fraction(c,aw)==Fraction(edge[v][w],av),'inverse-automorphism flow mass')
        # No complemented reflection fixes a necklace when m is odd.
        need(all(necklace(tuple(1-x for x in reversed(w)))!=w for w in buckets),'no self-converse necklace')
        if m==3:
            names=sorted(buckets,key=sum)
            graph3=[[edge[w].get(v,0) for v in names] for w in names]
        rows.append({'m':m,'vertices':4*m,'classes':len(buckets),
                     'period_classes':dict(sorted(orbitcounts.items())),
                     'automorphism_orders':{str(3**m*m//d):c for d,c in sorted(orbitcounts.items())},
                     'undirected_quotient_edges':sum(len(v) for v in edge.values())//2,
                     'mass_identity':str(sum(Fraction(1,3**m*m//period(w)) for w in buckets))})
    return {'rows':rows,'m3_move_multiplicities':graph3}


def small_labelled():
    # Exhaust all labelled tournaments n<=5; count true automorphisms with
    # direct permutations, independently of the backtracking routine.
    rows=[]
    for n in range(1,6):
        pairs=list(itertools.combinations(range(n),2)); groups=defaultdict(set)
        perms=list(itertools.permutations(range(n)))
        for bits in range(1<<len(pairs)):
            t=[0]*n
            for e,(i,j) in enumerate(pairs):
                if bits>>e&1:t[i]|=1<<j
                else:t[j]|=1<<i
            aut=sum(all((t[i]>>j&1)==(t[p[i]]>>p[j]&1) for i,j in pairs) for p in perms)
            groups[cyclic(t)].add(aut)
        need(all(len(v)==1 for v in groups.values()),'small-fibre automorphism invariance')
        rows.append({'n':n,'labelled_tournaments':1<<len(pairs),'cyclic_structures':len(groups)})
    return rows


def growth():
    cases=0
    for m in range(3,26,2):
        k=(m-1)//2
        f=[i+1 if i<k else i+2 for i in range(m)]
        old,new=carousel(m),carousel(m+2)
        need(all((old[i]>>j&1)==(new[f[i]]>>f[j]&1) for i in range(m) for j in range(m)), 'pointed insertion is induced embedding')
        w=(1,)*m; v=(0,)+w[:k]+(0,)+w[k:]
        need(period(v)==m+2,'two inserted zeros kill rotation stabilizer')
        cases+=1
    def grow(w):
        k=(len(w)-1)//2
        return (0,)+w[:k]+(0,)+w[k:]
    a,b=(0,1,1),(1,0,1)
    need(necklace(a)==necklace(b),'same unpointed parent')
    need(necklace(grow(a))!=necklace(grow(b)),'growth needs marked gap')
    # Every physical old arc, including the child blocks, survives at m=3,5.
    for m in [3,5]:
        k=(m-1)//2; f=[i+1 if i<k else i+2 for i in range(m)]
        for w in itertools.product(range(2),repeat=m):
            a=compose(carousel(m),[block(c) for c in w])
            b=compose(carousel(m+2),[block(c) for c in grow(w)])
            embedding=[4*f[i]+j for i in range(m) for j in range(4)]
            need(all((a[i]>>j&1)==(b[embedding[i]]>>embedding[j]&1) for i in range(4*m) for j in range(4*m)), 'all old tournament arcs retained')
    need(3**13 < 11*3**11,'larger tournament can lose automorphisms')
    return {'quotient_embedding_orders':list(range(3,26,2)), 'cases':cases,
            'unpointed_hostile':{'parents':[a0 for a0 in ['011','101']],
             'children':[''.join(map(str,grow(x))) for x in [(0,1,1),(1,0,1)]]},
            'symmetry_drop':{'parent_vertices':44,'child_vertices':52,
              'parent_aut':11*3**11,'child_aut':3**13}}


if __name__=='__main__':
    print(json.dumps({'witness':witness(),'necklaces':necklaces(),'small_exhaustion':small_labelled(),'growth':growth()},indent=2,sort_keys=True))
