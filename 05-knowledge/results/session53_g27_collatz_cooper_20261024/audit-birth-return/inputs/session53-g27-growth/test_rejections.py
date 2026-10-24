#!/usr/bin/env python3
from check import *

def expect_reject(f):
    try:f()
    except (ValueError,TypeError,KeyError):return
    raise ValueError('hostile input accepted')

def audit():
    # Every labelled simple graph on1..5: compare endpoint-state birth to
    # independent degree-bounded acyclic spanning edge-subset enumeration.
    graphs=0;forests=0
    cache={}
    for n in range(1,6):
        E=list(combinations(range(1,n+1),2))
        for mask in range(1<<len(E)):
            G={v:[] for v in range(1,n+1)}
            for i,(a,b) in enumerate(E):
                if mask>>i&1:G[a].append(b);G[b].append(a)
            exact=Counter();poly=Counter();cover_list=[]
            for k in range(1,n+1):
                for c in edge_covers(G,k):
                    validate_cover(G,c);exact[endpoint_state(c)]+=1;poly[k]+=1;cover_list.append(c)
            if n==1:pred=state_transfer({():1},1,[])
            else:
                oldE=list(combinations(range(1,n),2));oldmask=sum((b in G[a])<<i for i,(a,b) in enumerate(oldE))
                olddist,oldcovers=cache[n-1,oldmask]
                pred=state_transfer(olddist,n,G[n]);polypred=Counter()
                m=min(len(c) for c in oldcovers);emin=jmin=jnext=0
                for c in oldcovers:
                    e,j=responses(c,G[n]);k=len(c)
                    polypred[k+1]+=1;polypred[k]+=e;polypred[k-1]+=j
                    if k==m:emin+=e;jmin+=j
                    elif k==m+1:jnext+=j
                polypred=+polypred
                require(poly==polypred,'polynomial birth mismatch')
                predicted_min=m-1 if jmin else m if emin or jnext else m+1
                require(min(poly)==predicted_min,'minimum-cover trichotomy mismatch')
            require(exact==pred,'endpoint transfer mismatch')
            cache[n,mask]=(exact,cover_list);graphs+=1;forests+=sum(exact.values())
    # Minimum endpoint response alone misses near-minimum insertion.
    P4={1:[2],2:[1,3],3:[2,4],4:[3]}
    require(sum(responses(c,{2,3})[0] for c in edge_covers(P4,1))==0,'P4 old endpoints')
    require(sum(responses(c,{2,3})[1] for c in edge_covers(P4,2))==1,'P4 near-minimum rescue')
    # Same exposed set, different pairing: future response differs.
    A=((1,2),(3,4));B=((1,3),(2,4))
    require(responses(A,{1,2})==(2,0),'same-path join hostile')
    require(responses(B,{1,2})==(2,1),'distinct-path join control')
    # A singleton has one usable vertex, not two oriented endpoint slots.
    require(sum(state_birth(((1,1),),2,{1}).values())==2,'singleton double count')
    # Degree-three birth can give two paths from the same old two-cover.
    require(responses(((1,2),(3,)),{1,2,3})==(3,2),'multiplicity hostile')
    expect_reject(lambda:graph(0));expect_reject(lambda:graph(True))
    expect_reject(lambda:validate_cover(graph(3),[(1,3),(3,2)]))
    expect_reject(lambda:validate_cover(graph(3),[(1,2,3)]))
    expect_reject(lambda:state_birth(((1,2),),2,{1}))
    print(json.dumps({'all_simple_graphs_orders_1_through_5':graphs,'linear_forests_checked':forests,'controls':'PASS'},sort_keys=True))
if __name__=='__main__':audit()
