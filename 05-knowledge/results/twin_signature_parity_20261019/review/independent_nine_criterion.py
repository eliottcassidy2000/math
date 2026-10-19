"""Independent complete order-nine criterion census using direct permutations."""
from independent_tree_capacity import trees_to, canon, automorphisms, adjacency
from independent_cherry_extensions import linear_extensions, arcs
from itertools import combinations, permutations
from math import factorial
import json,hashlib

def need(x,m):
    if not x:raise RuntimeError(m)

def subgraph(edges,vertices):
    vv=sorted(vertices);ix={x:i for i,x in enumerate(vv)}
    return tuple((ix[u],ix[v]) for u,v in edges if u in ix and v in ix)

def split(n,edges,edge):
    aa=adjacency(n,[e for e in edges if e!=edge]);seen={edge[0]};todo=[edge[0]]
    for u in todo:
        for v in aa[u]:
            if v not in seen:seen.add(v);todo.append(v)
    return seen,set(range(n))-seen

def main():
    tt=trees_to(9);auts={n:{c:len(automorphisms(n,e)) for c,e in row.items()} for n,row in tt.items()}
    need([len(tt[n]) for n in range(1,10)]==[1,1,1,2,3,6,11,23,47],'tree universe')
    good=[];bad=[];types={};padding={};hered=0
    for c,e in tt[8].items():
        aa=adjacency(8,e);same=True
        vals={linear_extensions(9,arcs(e,m))%2 for m in range(128)}
        need(len(vals)==1,'padding orientation independence control')
        if auts[8][c]>1:padding[c]=next(iter(vals))
    need(sum(padding.values())==2,'exact two odd symmetric padding types')
    for c,e in tt[9].items():
        aa=adjacency(9,e);fail=[]
        if auts[9][c]==1:fail.append('asymmetric tree')
        for v in range(9):
            if len(aa[v])==1:
                ee=subgraph(e,set(range(9))-{v})
                if auts[8][canon(8,ee)]==1:fail.append('asymmetric leaf deletion');break
        for cut in e:
            parts=split(9,e,cut)
            if sorted(map(len,parts))==[2,7]:
                part=next(p for p in parts if len(p)==7)
                if auts[7][canon(7,subgraph(e,part))]==1:fail.append('asymmetric seven-side cut');break
        H=True
        for k in (7,8,9):
            for vv in combinations(range(9),k):
                ee=subgraph(e,vv)
                if len(ee)==k-1 and auts[k][canon(k,ee)]==1:H=False;break
            if not H:break
        hered+=H
        if fail:bad.append({'canonical':c,'edges':e,'failed_conditions':fail})
        else:
            good.append(c);need(H or not H,'boolean')
            types[c]={'hereditary':H,'automorphisms':auts[9][c]}
    need(len(good)==36 and len(bad)==11 and hered==24,'nine census')
    # Independently verify the small all-orientation P4 input on all tournaments.
    p4=((0,1),(1,2),(2,3));pairs=list(combinations(range(4),2))
    for mask in range(64):
        hh=arcs(pairs,mask)
        for o in range(8):
            aa=arcs(p4,o)
            count=sum(all((p[u],p[v]) in hh for u,v in aa) for p in permutations(range(4)))
            need(count%2==1,'all P4 orientations odd')
    print(json.dumps({'status':'PASS','unlabelled_trees':47,'good_by_proved_criterion':36,'excluded':11,'hereditarily_symmetric':24,'nonhereditary_good':12,'symmetric_eight_padding_types':len(padding),'odd_symmetric_eight_padding_types':sum(padding.values()),'P4_orientation_host_pairs':512,'good_canonical_digest':hashlib.sha256('\n'.join(sorted(good)).encode()).hexdigest(),'excluded_classes':bad,'scope':'Exhaustive source-side criterion census. The theorem establishes universal host conclusions; no nine-tournament sampling is substituted.'},indent=2))

if __name__=='__main__':main()
