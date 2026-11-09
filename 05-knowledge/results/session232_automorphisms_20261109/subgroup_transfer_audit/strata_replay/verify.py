#!/usr/bin/env python3
"""Exact, self-contained subgroup orbit census through order 7.
No assertion statement is used; python3 -O retains every gate.
"""
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path
import hashlib
import json

OUT=Path(__file__).resolve().parent

def require(ok, message):
    if not ok: raise RuntimeError(message)

def degree_orders(rows):
    bydegree=defaultdict(list)
    for v,row in enumerate(rows): bydegree[row.bit_count()].append(v)
    return product(*(permutations(bydegree[d]) for d in sorted(bydegree)))

def relabel(rows, order):
    return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(order)) for v in order)

def canonical(rows):
    return min(relabel(rows,tuple(v for block in blocks for v in block)) for blocks in degree_orders(rows))

def automorphisms(rows):
    n=len(rows); groups=defaultdict(list)
    for i,row in enumerate(rows): groups[row.bit_count()].append(i)
    keys=sorted(groups)
    ans=[]
    for images in product(*(permutations(groups[d]) for d in keys)):
        p=[0]*n
        for d,block in zip(keys,images):
            for i,j in zip(groups[d],block): p[i]=j
        if all(((rows[u]>>v)&1)==((rows[p[u]]>>p[v])&1) for u in range(n) for v in range(u+1,n)):
            ans.append(tuple(p))
    return ans

def hpaths(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            if not c: continue
            rem=rows[v]&~mask
            while rem:
                bit=rem&-rem; rem-=bit; w=bit.bit_length()-1
                dp[mask|bit][w]+=c
    return sum(dp[-1])

def literal_paths(rows):
    n=len(rows)
    return [p for p in permutations(range(n)) if all((rows[p[i]]>>p[i+1])&1 for i in range(n-1))]

def subgroups(group):
    n=len(group[0]); ident=group.index(tuple(range(n)))
    ix={g:i for i,g in enumerate(group)}
    mul=[[ix[tuple(a[b[v]] for v in range(n))] for b in group] for a in group]
    def generated(h,x):
        generators=list(h)+[x]; seen={ident}; todo=[ident]
        while todo:
            a=todo.pop()
            for g in generators:
                b=mul[a][g]
                if b not in seen: seen.add(b); todo.append(b)
        return frozenset(seen)
    found={frozenset([ident])}; todo=list(found)
    while todo:
        h=todo.pop()
        for x in range(len(group)):
            if x not in h:
                k=generated(h,x)
                if k not in found: found.add(k); todo.append(k)
    return sorted(found,key=lambda h:(len(h),tuple(sorted(h))))

def normalized_census(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    fibres=Counter()
    for mask in range(1<<len(free)):
        rows=[0]*n
        for i in range(n-1): rows[i]|=1<<(i+1)
        for k,(i,j) in enumerate(free):
            if mask>>k&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
        fibres[canonical(rows)]+=1
    return fibres

def cayley(n,steps):
    return tuple(sum(1<<((v+s)%n) for s in steps) for v in range(n))

def encode_graph(rows):
    n=len(rows); require(n>=3,'graph encoding scope')
    adj=[set() for _ in rows]
    def vertex(): adj.append(set()); return len(adj)-1
    def edge(u,v): adj[u].add(v); adj[v].add(u)
    def pendant(v,length):
        for _ in range(length):
            w=vertex(); edge(v,w); v=w
    arcs={}
    for u in range(n):
        for v in range(n):
            if rows[u]>>v&1:
                a,b=vertex(),vertex(); edge(u,a); edge(a,b); edge(b,v)
                pendant(a,1); pendant(b,2); arcs[u,v]=(a,b)
    return adj,arcs

def decode_graph(adj):
    # Recover the 2-core by leaf peeling, independent of supplied vertex roles.
    alive=set(range(len(adj))); degree=[len(nb) for nb in adj]
    todo=[v for v in alive if degree[v]<2]
    while todo:
        v=todo.pop()
        if v not in alive: continue
        alive.remove(v)
        for w in adj[v]:
            if w in alive:
                degree[w]-=1
                if degree[w]<2: todo.append(w)
    marker={}
    for v in alive:
        outside=adj[v]-alive
        if not outside: marker[v]=0
        else:
            require(len(outside)==1,'unique marker arm')
            w=next(iter(outside)); prev=v; length=1
            while len(adj[w])!=1:
                candidates=adj[w]-{prev}
                require(len(candidates)==1,'marker is simple arm')
                prev,w=w,next(iter(candidates)); length+=1
            marker[v]=length
    originals=sorted(v for v in alive if marker[v]==0)
    mapping={v:i for i,v in enumerate(originals)}; rows=[0]*len(originals)
    for u in originals:
        for a in adj[u]:
            if marker[a]!=1: continue
            b=next(x for x in adj[a]&alive if marker[x]==2)
            v=next(x for x in adj[b]&alive if marker[x]==0)
            rows[mapping[u]]|=1<<mapping[v]
    return tuple(rows)

def graph_automorphisms(adj):
    n=len(adj)
    def refine(c,d):
        while True:
            signatures=[]
            for colors in (c,d):
                signatures.append([(colors[v],tuple(sorted(colors[w] for w in adj[v]))) for v in range(n)])
            palette={sig:i for i,sig in enumerate(sorted(set(signatures[0]+signatures[1])))}
            cc,dd=([palette[sig] for sig in sigs] for sigs in signatures)
            if Counter(cc)!=Counter(dd): return None
            if len(set(cc))==len(set(c)): return cc,dd
            c,d=cc,dd
    result=[]
    def visit(c,d):
        pair=refine(c,d)
        if pair is None:return
        c,d=pair; cells=defaultdict(list); targets=defaultdict(list)
        for v in range(n): cells[c[v]].append(v); targets[d[v]].append(v)
        nontrivial=[(len(vs),color) for color,vs in cells.items() if len(vs)>1]
        if not nontrivial:
            perm=tuple(targets[c[v]][0] for v in range(n))
            require(all({perm[w] for w in adj[v]}==adj[perm[v]] for v in range(n)),'full graph aut validation')
            result.append(perm);return
        _,color=min(nontrivial); v=cells[color][0]; fresh=max(c)+1
        for w in targets[color]:
            cc=c.copy();dd=d.copy();cc[v]=fresh;dd[w]=fresh;visit(cc,dd)
    colors=[len(a) for a in adj];visit(colors,colors)
    return result

def graph_checks(rows):
    adj,arcs=encode_graph(rows); n=len(rows)
    require(len(adj)==n+5*n*(n-1)//2,'graph order')
    require(decode_graph(adj)==rows,'graph decode')
    edges=[(u,v) for u,nb in enumerate(adj) for v in nb if u<v]
    # Independent full uncoloured graph aut enumeration, with no gadget roles.
    autos=graph_automorphisms(adj)
    tgroup=set(automorphisms(rows))
    require({tuple(p[v] for v in range(n)) for p in autos}==tgroup,'graph automorphism restrictions')
    require(len(autos)==len(tgroup),'graph has no hidden automorphism kernel')
    signs=Counter(sum(1 for u,v in edges if p[u]>p[v])%2 for p in autos)
    require(not signs[1],'graph sign-even')
    paths=literal_paths(rows)
    encoded=[]
    for path in paths:
        word=[path[0]]
        for u,v in zip(path,path[1:]): word.extend((*arcs[u,v],v))
        require(len(set(word))==len(word),'encoded path simple')
        require(all(v in adj[u] for u,v in zip(word,word[1:])),'encoded graph traversal')
        encoded.append(tuple(word))
    require(len(set(encoded))==hpaths(rows),'path encoding injective')
    leaves=sum(len(x)==1 for x in adj)
    require(leaves>2,'ordinary graph Hamilton paths must be absent')
    return {'vertices':len(adj),'edges':len(edges),'leaves':leaves,
        'aut_order':len(autos),'edge_reversal_sign_counts':dict(signs),
        'admissible_paths':len(encoded),'ordinary_graph_hamilton_paths':0}

def valuation(n,p):
    k=0
    while n%p==0: n//=p;k+=1
    return k

def ordinal_sum(*blocks):
    n=sum(map(len,blocks)); rows=[]; offset=0
    for block in blocks:
        above=((1<<n)-1)^((1<<(offset+len(block)))-1)
        rows.extend((r<<offset)|above for r in block);offset+=len(block)
    return tuple(rows)

def sccs(rows):
    reach=[rows[v]|(1<<v) for v in range(len(rows))]
    for k in range(len(rows)):
        for v in range(len(rows)):
            if reach[v]>>k&1: reach[v]|=reach[k]
    left=set(range(len(rows))); result=[]
    while left:
        v=min(left); block=sorted(w for w in left if (reach[v]>>w&1) and (reach[w]>>v&1))
        left.difference_update(block)
        result.append(tuple(sum(((rows[u]>>w)&1)<<j for j,w in enumerate(block)) for u in block))
    return result

def equality_stack(rows,alpha,beta):
    counts=Counter(canonical(x) for x in sccs(rows))
    allowed={canonical((0,)),canonical(cayley(3,[1])),canonical(cayley(5,[1,2])),canonical(cayley(7,[1,2,4]))}
    return (set(counts)<=allowed and counts[canonical(cayley(5,[1,2]))]==alpha
            and counts[canonical(cayley(7,[1,2,4]))]==beta)

def main():
    expected={1:1,2:1,3:2,4:4,5:12,6:56,7:456}; summary=[]; witnesses=[]
    for n in range(1,8):
        fibres=normalized_census(n); require(len(fibres)==expected[n],f'class census n={n}')
        quotient_values=set(); pairs=0; group_hist=Counter(); low_prime_gate=0; prime_cost_tests=0; equality_tests=0
        for rows,mass in sorted(fibres.items()):
            group=automorphisms(rows); h=hpaths(rows); require(h%len(group)==0,'free divisibility')
            q=h//len(group); require(q==mass,'normalized path mass')
            group_hist[len(group)]+=1
            if q<9:
                z=len(group)
                for prime in (3,5):
                    while z%prime==0: z//=prime
                require(z==1,'low-quotient prime restriction'); low_prime_gate+=1
            paths=literal_paths(rows)
            require(len(paths)==h,'independent literal Hamilton count')
            for subgroup in subgroups(group):
                pairs+=1
                reps={min(tuple(group[g][v] for v in path) for g in subgroup) for path in paths}
                count=len(reps); require(count*len(subgroup)==h,'literal subgroup free orbit count')
                require(count not in (7,21),'subgroup gap')
                require(count==len(group)//len(subgroup)*q,'restriction index law')
                if count<23:
                    prime_cost_tests+=1
                    z=len(subgroup)
                    for prime in (3,5,7):
                        while z%prime==0:z//=prime
                    require(z==1,'small quotient prime support')
                    alpha,beta=valuation(len(subgroup),5),valuation(len(subgroup),7)
                    price=3**(alpha+2*beta)
                    require(count>=price,'small quotient symmetry price')
                    structure=(len(subgroup)==len(group) and equality_stack(rows,alpha,beta))
                    require((count==price)==structure,'equality iff stack')
                    if count==price:equality_tests+=1
                quotient_values.add(count)
        require(sum(fibres.values())==1<<((n-1)*(n-2)//2),'pointed mass')
        summary.append({'n':n,'classes':len(fibres),'pointed_mass':sum(fibres.values()),
            'literal_subgroups_over_classes':pairs,'prime_cost_tests':prime_cost_tests,'equality_cases':equality_tests,'automorphism_order_histogram':dict(sorted(group_hist.items())),
            'subgroup_orbit_values':sorted(quotient_values),'q_less_than_9_classes':low_prime_gate})
    for name,steps in [('carousel7',[1,2,3]),('paley7',[1,2,4])]:
        rows=cayley(7,steps); group=automorphisms(rows); h=hpaths(rows)
        subgroup_counts=Counter(len(k) for k in subgroups(group))
        witnesses.append({'name':name,'rows':rows,'H':h,'aut_order':len(group),'q':h//len(group),
             'subgroup_order_counts':dict(sorted(subgroup_counts.items())),
             'orbit_counts':sorted({h//s for s in subgroup_counts}),
             'graph':graph_checks(rows)})
    sharp=[]
    for name,rows,expected in [
        ('R5',cayley(5,[1,2]),(15,5,3)),
        ('R5 then R5',ordinal_sum(cayley(5,[1,2]),cayley(5,[1,2])),(225,25,9)),
        ('P7 then C3',ordinal_sum(cayley(7,[1,2,4]),cayley(3,[1])),(567,63,9))]:
        h=hpaths(rows);a=len(automorphisms(rows));q=h//a
        require((h,a,q)==expected,'sharp control')
        require(equality_stack(rows,valuation(a,5),valuation(a,7)),'sharp equality structure')
        sharp.append({'name':name,'rows':rows,'H':h,'aut_order':a,'q':q})
    data={'status':'FINITE-EXACT','universe':'All ordinary tournament isomorphism classes, n=1..7; every literal automorphism subgroup of each representative',
          'checks':summary,'explicit_graph_controls':witnesses,'sharp_controls':sharp,
          'external_software':'none; Python standard library only',
          'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out=json.dumps(data,indent=2,sort_keys=True)+'\n'; (OUT/'receipt.json').write_text(out); print(out,end='')
if __name__=='__main__':main()
