"""Finite exact controls for the affine compiler and overlapping odd prefixes.

No external imports, isomorphism/SCC filters, random selection or disabled
assertions. The producer parity graph is checked against direct matrices.
"""
from itertools import combinations, permutations, product
from collections import Counter
import json
from explore import Parity, edge, edges, source, words, feasible

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def literal(n,bits):
    T=[[False]*n for _ in range(n)]
    for (u,v),bit in zip(combinations(range(n),2),bits):T[u][v]=bool(bit);T[v][u]=not bit
    return T

def anti(T,q):
    r=dict(zip(q,q[::-1]));n=len(q)
    return all(T[u][v]==T[r[v]][r[u]] for u in range(n) for v in range(n) if u!=v)

def directed(T,q):return all(T[u][v] for u,v in zip(q,q[1:]))

def flip(T,q):
    E=edges(q);n=len(q)
    return [[False if u==v else T[u][v]^(edge(u,v) not in E) for v in range(n)] for u in range(n)]

def cycles(g):
    seen=set();out=[]
    for v in range(len(g)):
        if v in seen:continue
        C=[];u=v
        while u not in seen:seen.add(u);C.append(u);u=g[u]
        out.append(C)
    return out

def forced(d,u,v):
    r,x=d.find(d.ids[edge(u,v)]);s,y=d.find(d.c)
    return (x^y^(u>v)) if r==s else None

def dimension(d):
    root=d.find(d.c)[0]
    return len({d.find(i)[0] for i in range(d.c)}-{root})

def affine_move_blue(d,q,flips):d.anti(q,flips);return d.ok

def matchings(xs):
    if not xs:yield ();return
    a=xs[0]
    for b in xs[1:]:
        for rest in matchings(tuple(x for x in xs[1:] if x!=b)):yield ((a,b),)+rest

def reversal_words(n):
    for fixed in range(n) if n%2 else (None,):
        for pairs in matchings(tuple(x for x in range(n) if x!=fixed)):
            yield tuple(a for a,b in pairs)+(() if fixed is None else (fixed,))+tuple(b for a,b in pairs[::-1])

def main():
    counts=Counter();small_hist={}
    # Independent literal tournament path/reversal census versus parity feasibility.
    for n in range(2,6):
        P=tuple(range(n));qs=list(permutations(range(n)));N=n*(n-1)//2
        actual=Counter()
        for bits in product((0,1),repeat=N):
            T=literal(n,bits)
            if not(anti(T,P) and directed(T,P)):continue
            for q in qs:
                if anti(T,q) and directed(T,q):actual[q]+=1
                counts['literal_pair_comparisons']+=1
        for q in qs:
            d=Parity(n);d.anti(P);d.anti(q)
            for u,v in zip(P,P[1:]):d.force(u,v)
            for u,v in zip(q,q[1:]):d.force(u,v)
            expected=(1<<dimension(d)) if d.ok else 0
            need(expected==actual[q],f'literal pair count n{n} q{q}')
            counts['exact_two_blue_families']+=1
        small_hist[str(n)]=dict(sorted(Counter(actual.values()).items()))

    # Exact target compiler and legality over the complete tiny ledger universe.
    n=3;P=tuple(range(n));qs=list(permutations(range(n)));pairs=list(combinations(range(n),2));allpairs=set(pairs)
    legal_ledgers=0
    for bits in product((0,1),repeat=3):
        T=literal(n,bits)
        for ledger in product(qs,repeat=3):
            now=T;flips=set();legal=True
            for q in ledger:
                bybits=all(bool(bits[pairs.index(edge(u,v))]^(u>v)^(edge(u,v) in flips)) for u,v in zip(q,q[1:]))
                need(bybits==directed(now,q),'intermediate affine legality')
                legal&=bybits
                now=flip(now,q);flips^=allpairs^edges(q)
            expected=flips==(allpairs^edges(P))
            need((now==flip(T,P))==expected,'literal target vs affine target')
            legal_ledgers+=legal;counts['three_move_ledgers']+=1

    # All affine families in a bounded guarded theorem universe, not one sample/family.
    population=0;new_cells=0
    for n in range(11,37):
        P=tuple(range(n));allpairs=set(combinations(range(n),2))
        for ell in range(1,n-3,2):
            if n<2*ell+3:continue
            Q=words(n,ell)
            need(edges(Q[0])^edges(Q[1])^edges(Q[2])==edges(P),'ledger XOR')
            for h in range(2,n-3):
                if h in {2,ell+2} or (ell>=3 and h in {4,ell-1}):continue
                d=source(n,h,ell)
                need(d.ok,'theorem source family nonempty')
                population+=1<<dimension(d)
                flips=set()
                for i,q in enumerate(Q):
                    for u,v in zip(q,q[1:]):
                        need(forced(d,u,v)==(1^(edge(u,v) in flips)),'universal join implication')
                    f=source(n,h,ell);f.anti(q,flips)
                    need(not f.ok,'universal blackness through affine contradiction')
                    flips^=allpairs^edges(q)
                counts['guarded_affine_cells']+=1
                if n<2*ell+6:new_cells+=1

    # The overlap cycle table, using direct positional reversals, through L=300.
    for L in range(8,301,2):
        for j in (1,2,3):
            n=2*L-j;q=words(n,L-3)[2];r=dict(zip(q,q[::-1]));g=[r[n-1-v] for v in range(n)]
            cs=cycles(g);oa=next(c for c in cs if 0 in c);ports={1,n-3,n-2,n-1}
            if L%3 and j in (1,3):need(len(cs)==1,'single overlap cycle')
            elif L%3:
                need(len(oa)==L-1 and len(set(oa)&ports)<=2,'two-overlap orbit capacity')
                tail=2 if L%6==2 else 1
                top=2*L-4 if tail==2 else 2*L-5
                need(oa==list(range(0,top+1,3))+list(range(tail,L-2,3)),'explicit a orbit')
            elif j in (1,2):need(any(len(c)==L//3 for c in cs),'unchanged even residue cycle')
            else:need(len(cs)==3 and all(len(c)==2*L//3-1 for c in cs) and set(oa)&ports=={n-3},'three residue cycles')
            counts['overlap_cycle_tables']+=1

    # All long odd cuts in this finite box. These are evidence, not an all-order theorem.
    failures=[]
    for n in range(15,41):
        for ell in range(1,n-3,2):
            for h in range(2,n-3):
                if h in {2,ell+2}:continue
                for j in (1,2):
                    f=feasible(n,h,ell,j)
                    if f:failures.append([n,ell,h,j+1])
                    counts['long_cut_blue_queries']+=1
    need(failures==[[15,9,3,3],[15,9,6,3],[18,9,3,3],[18,9,6,3]],'retained genuine blue final hostiles')

    # Complete positional-reversal universe through ten for the iff criterion.
    for n in range(2,11):
        P=tuple(range(n))
        for q in reversal_words(n):
            r=dict(zip(q,q[::-1]));g=[r[n-1-v] for v in range(n)]
            odd=all(len(c)%2 for c in cycles(g))
            d=Parity(n);d.anti(P);d.anti(q)
            need(d.ok==odd,'sharp common anti-reflection criterion')
            counts['common_reflection_families']+=1

    # Positive/hostile signed cycle, and path-edge XOR without legality.
    d=Parity(3);d.join(0,1,0);d.join(1,2,0);need(d.ok,'positive cycle prefix');d.join(2,0,1);need(not d.ok,'odd signed cycle')
    P=(0,1,2);T=literal(3,(1,1,1));bad=(P,P[::-1],P)
    need(edges(bad[0])^edges(bad[1])^edges(bad[2])==edges(P),'hostile XOR correct')
    need(not directed(flip(T,P),bad[1]),'hostile illegal intermediate join')
    print(json.dumps({'status':'PASS','scope':'finite exact; unrestricted theorem is in README',
      'counts':dict(counts),'new_overlap_cells':new_cells,'source_population_with_cell_multiplicity':str(population),
      'literal_small_family_histogram':small_hist,'legal_tiny_ledgers':legal_ledgers,
      'long_cut_blue_failures':failures,'controls':['odd signed cycle','literal reversed-join hostile','genuine blue final hostiles']},sort_keys=True,indent=2))

if __name__=='__main__':main()
