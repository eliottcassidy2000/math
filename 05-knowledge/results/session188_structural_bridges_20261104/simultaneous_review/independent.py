#!/usr/bin/env python3
"""Literal-word full-response banks and original-module XOR audit.
No imports from the producer's response or module code.
"""
from itertools import combinations,permutations,product
from collections import defaultdict,Counter
from pathlib import Path
import json


def need(ok,why):
    if not ok:raise RuntimeError(why)


def basis(n):
    masks=[s for s in range(1<<n) if s.bit_count()<=2]
    return sorted(a|(c<<n) for a in masks for c in masks)


def edge_table(n):return list(combinations(range(n),2))


def rows(n,code):
    a=[0]*n
    for b,(u,v) in enumerate(edge_table(n)):
        if code>>b&1:a[u]|=1<<v
        else:a[v]|=1<<u
    return a


def word_skeletons(n):
    index={p:i for i,p in enumerate(edge_table(n))}
    result=Counter()
    for word in permutations(range(n+2)):
        mask=value=ones=zeros=0
        for u,v in zip(word,word[1:]):
            if u<n and v<n:
                bit=1<<index[min(u,v),max(u,v)]
                mask|=bit
                if u<v:value|=bit
            elif u>=n and v>=n:
                if (u,v)!=(n,n+1):break
            else:
                old,port=(u,v) if u<n else (v,u)
                bit=1<<(old+n*(port-n))
                if u<n:ones|=bit
                else:zeros|=bit
        else:
            need(not ones&zeros,'conflicting word contact')
            result[mask,value,ones,zeros]+=1
    return result


def polynomial_bank(n):
    B=basis(n);index={s:i for i,s in enumerate(B)}
    e=n*(n-1)//2;all_edges=(1<<e)-1
    bank=[[0]*len(B) for _ in range(1<<e)]
    for (mask,value,ones,zeros),weight in word_skeletons(n).items():
        terms=[];s=zeros
        while True:
            terms.append((index[ones|s],(-1 if s.bit_count()%2 else 1)*weight))
            if not s:break
            s=(s-1)&zeros
        free=all_edges^mask;s=free
        while True:
            vector=bank[value|s]
            for j,c in terms:vector[j]+=c
            if not s:break
            s=(s-1)&free
    return B,[tuple(v) for v in bank]


def modules(a):
    n=len(a);answer=[]
    for M in range(1,1<<n):
        if M.bit_count()<2:continue
        good=True
        for v in range(n):
            if not M>>v&1:
                edges=[bool(a[v]>>u&1) for u in range(n) if M>>u&1]
                if min(edges)!=max(edges):good=False;break
        if good:answer.append(M)
    return answer


def hull(a,u,v):
    M=(1<<u)|(1<<v)
    while True:
        old=M
        for x in range(len(a)):
            if not M>>x&1 and 0<(a[x]&M)<M:M|=1<<x
        if M==old:return M


def flip(n,M):return sum(1<<i for i,(u,v) in enumerate(edge_table(n)) if M>>u&1 and M>>v&1)


def category(ms):
    if len(ms)==1:return 'single'
    intersections=[(a,b,a&b) for a,b in combinations(ms,2)]
    if all(not i for a,b,i in intersections):return 'pairwise_disjoint'
    if all(not i or i==a or i==b for a,b,i in intersections):return 'nested_laminar'
    return 'overlapping'


def literal_full(a):
    n=len(a);patterns=Counter()
    for word in permutations(range(n+2)):
        ones=zeros=0
        for u,v in zip(word,word[1:]):
            if u<n and v<n:
                if not a[u]>>v&1:break
            elif u>=n and v>=n:
                if (u,v)!=(n,n+1):break
            else:
                old,port=(u,v) if u<n else (v,u)
                bit=1<<(old+n*(port-n))
                if u<n:ones|=bit
                else:zeros|=bit
        else:patterns[ones,zeros]+=1
    polynomial=Counter();truth=[0]*(1<<(2*n));full=(1<<(2*n))-1
    for (ones,zeros),weight in patterns.items():
        s=zeros
        while True:
            polynomial[ones|s]+=weight*(-1 if s.bit_count()%2 else 1)
            if not s:break
            s=(s-1)&zeros
        free=full^(ones|zeros);s=free
        while True:
            truth[ones|s]+=weight
            if not s:break
            s=(s-1)&free
    return {s:c for s,c in polynomial.items() if c},truth


def main():
    root=Path(__file__).resolve().parent.parent
    class_counts=Counter();tournaments=families=coordinates=hulls=module_count=0
    census=[]
    for n in range(1,6):
        B,bank=polynomial_bank(n)
        published=json.loads((root/'simultaneous_response'/f'full_response_n{n}.json').read_text())
        need(published['basis']==[[s&((1<<n)-1),s>>n] for s in B],'basis')
        need(len(set(bank))==len(bank),'full-response injectivity in finite bank')
        for rec in published['records']:
            code=rec['code']
            need(bank[code]==tuple(rec['vector']) and rows(n,code)==rec['rows'],'full exact bank')
        local=Counter();pairs=Counter()
        for code in range(len(bank)):
            a=rows(n,code);mods=modules(a);module_count+=len(mods)
            generated={}
            for u,v in edge_table(n):
                H=hull(a,u,v);hulls+=1;generated.setdefault(H,(u,v))
            need(set(generated)==set(mods),'pair hull generation theorem')
            ordered=sorted(mods,key=lambda M:(M.bit_count(),M))
            for M in ordered:
                u,v=generated[M]
                for N in ordered:
                    need(bool(N>>u&1 and N>>v&1)==(M&N==M),'inclusion zeta minor')
            need(len(mods)<=n*(n-1)//2,'module count bound')
            # A tournament is transitive exactly when some order realizes all arcs.
            transitive=any(all(a[u]>>v&1 for i,u in enumerate(w) for v in w[i+1:]) for w in permutations(range(n)))
            need((len(mods)==n*(n-1)//2)==transitive,'sharp module count characterization')
            targets=[0]*(1<<len(mods));edgeflips=[flip(n,M) for M in mods]
            for f in range(1,1<<len(mods)):
                bit=f&-f;j=bit.bit_length()-1;targets[f]=targets[f^bit]^edgeflips[j]
                selected=[M for j,M in enumerate(mods) if f>>j&1]
                kind=category(selected);local[kind]+=1;families+=1
                if len(selected)==2:pairs[kind]+=1
                target=code^targets[f]
                need(target!=code,'nonempty identity family')
                need(bank[target]!=bank[code],'full exact response collision')
            tournaments+=1
        coordinates+=len(B)*len(bank);class_counts.update(local)
        census.append(dict(n=n,families_by_type=dict(local),two_member_families=dict(pairs)))
    receipt=json.loads((root/'simultaneous_response'/'receipt.json').read_text())
    need(families==130942 and tournaments==1099,'family universe')
    need(dict(class_counts)==receipt['families_by_type'],'family classification counts')
    examples=[];cache={}
    for ex in receipt['examples']:
        a,b=ex['T'],ex['U'];n=len(a);B=basis(n)
        for adj in (a,b):
            if tuple(adj) not in cache:cache[tuple(adj)]=literal_full(adj)
        p,fa=cache[tuple(a)];q,fb=cache[tuple(b)]
        need(ex['complete_named_basis']==[[s&((1<<n)-1),s>>n] for s in B],'example basis')
        need(ex['F_T_vector']==[p.get(s,0) for s in B],'example full T vector')
        need(ex['F_U_vector']==[q.get(s,0) for s in B],'example full U vector')
        need(ex['difference_vector']==[p.get(s,0)-q.get(s,0) for s in B],'example exact difference')
        for x,c,a_count,b_count in ex['literal_queries']:
            need((fa[x|(c<<n)],fb[x|(c<<n)])==(a_count,b_count),'literal query battery')
        degrees=dict(Counter(s.bit_count() for s in B if p.get(s,0)!=q.get(s,0)))
        examples.append(dict(order=n,module_masks=ex['original_module_masks'],degree_counts=degrees,
                             independently_computed_full_truth_table_size=len(fa),published_queries_compared=len(ex['literal_queries'])))
        if ex['original_module_masks']==[3,12,48,63]:
            uniform=[sum(3<<(2*i) for i in range(3) if mask>>i&1) for mask in range(8)]
            need(all((fa[x|(c<<6)]-fb[x|(c<<6)])%4==0 for x,c in product(uniform,repeat=2)),'uniform mod4')
            need(sum(fa[x|(c<<6)]!=fb[x|(c<<6)] for x,c in product(uniform,repeat=2))==24,'uniform exact')
            need((fa[3|(12<<6)],fb[3|(12<<6)])==(509,545),'uniform separator')
            need((fa[1|(5<<6)]-fb[1|(5<<6)])%4==2,'nonuniform separator')
            need(degrees=={2:36,3:132,4:132},'six-vertex degree counts')
    return dict(status='PASS',tournaments=tournaments,families=families,families_by_type=dict(class_counts),
                exact_named_coordinates_compared=coordinates,pair_hulls=hulls,modules=module_count,
                census=census,examples=examples)


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
