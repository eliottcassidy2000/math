#!/usr/bin/env python3
"""Independent exact core-spine controls for two marked apices.
No producer imports. Uses 5**m <= H**3 instead of logarithmic arithmetic.
"""
import itertools,json

def need(c,msg):
    if not c:raise RuntimeError(msg)
def graph(word):
    m=len(word);out=[0]*(m+2)
    for i,c in enumerate(word):
        for j in range(i+1,m):out[i]|=1<<j
        for r in (0,1):
            if c>>r&1:out[m+r]|=1<<i
            else:out[i]|=1<<(m+r)
    out[m]|=1<<(m+1);return out

def paths(out):
    n=len(out);size=1<<n;d=[[0]*n for _ in range(size)];h=[1]+[0]*(size-1)
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,size):
        h[mask]=sum(d[mask])
        for i,v in enumerate(d[mask]):
            if not v:continue
            todo=out[i]&~mask
            while todo:
                b=todo&-todo;todo-=b;j=b.bit_length()-1;d[mask|b][j]+=v
    return h

def components(out):
    n=len(out);reach=[v|(1<<i) for i,v in enumerate(out)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    unseen=(1<<n)-1;comps=[]
    while unseen:
        b=unseen&-unseen;i=b.bit_length()-1
        C=sum(1<<j for j in range(n) if reach[i]>>j&1 and reach[j]>>i&1)
        comps.append(C);unseen&=~C
    return comps

def reduce_word(w):
    out=graph(w);m=len(w);cc=components(out);core=0
    for c in cc:
        if c.bit_count()>1:core|=c
    x,y=bool(core>>m&1),bool(core>>(m+1)&1)
    v=[]
    for i,c in enumerate(w):
        if core>>i&1:
            # Removed x is restored as universal source, removed y as sink.
            v.append((c&1 if x else 1)+(c&2 if y else 0))
    return tuple(v),core,cc

def main():
    cases=components_checked=restorations=0;largest=0
    for n in range(6):
        for w in itertools.product(range(4),repeat=n):
            h=paths(graph(w));v,core,cc=reduce_word(w);apex=(1<<n)|(1<<(n+1))
            for c in cc:
                if c.bit_count()>1:
                    need(c&apex,'strong SCC without apex')
                    need(5**(c.bit_count()-1)<=h[c]**3,'strong floor control')
                    components_checked+=1
            m=len(v)
            need(5**m<=h[-1]**3,'universal reduced-spine bound')
            need(paths(graph(v))[-1]==h[-1],'two-apex restoration changed H')
            preserved_apex_word=tuple(c for i,c in enumerate(w) if core>>i&1)
            need(paths(graph(preserved_apex_word))[-1]==h[-1],'retained-apex chain deletion changed H')
            need(len(preserved_apex_word)==m,'retained-apex spine mismatch')
            need(h[core]==h[-1],'singleton SCC deletion changed H')
            if core&apex!=apex:restorations+=1
            largest=max(largest,m);cases+=1
    # Scalar-core reduction must never be sold as continuation-profile reduction.
    empty=paths(graph(()))[-1];sinkword=(3,3,3,3)
    need(paths(graph(sinkword))[-1]==empty==1,'sink padding control')
    empty0=paths(graph((0,)))[-1];sink0=paths(graph(sinkword+(0,)))[-1]
    need(empty0!=sink0,'continuation loss hostile')
    print(json.dumps({'status':'PASS','all_two_apex_words_through_length':5,'word_cases':cases,'retained_original_apex_controls':cases,'nontrivial_SCC_floor_controls':components_checked,'missing_apex_restoration_cases':restorations,'largest_reduced_spine':largest,'exact_bound':'5**reduced_spine_length <= H**3','continuation_loss_hostile':{'empty_H':empty,'word':'3333','H':1,'empty_then_0_H':empty0,'3333_then_0_H':sink0},'scope':'Finite graph controls support the proved all-order SCC/Busch argument; code-size bound concerns existence of a scalar-H representative, not all profiles.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
