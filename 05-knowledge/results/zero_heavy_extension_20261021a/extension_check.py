#!/usr/bin/env python3
"""Two-unit-cell quotient and a literal three-presentation collision."""
from pathlib import Path
from importlib.util import module_from_spec, spec_from_file_location
from collections import Counter
from math import comb
import hashlib
import json

CANDIDATES = [Path(__file__).parent/'zero_heavy_check.py',
              Path(__file__).parent.parent/'session18-zero-heavy'/'zero_heavy_check.py']
SOURCE = next((p for p in CANDIDATES if p.is_file()), None)
if SOURCE is None:
    raise RuntimeError('Place the frozen zero_heavy_check.py beside this script')
if hashlib.sha256(SOURCE.read_bytes()).hexdigest() != '533fa6c595db3b3489825d1fb6ee48e1295f665c4ac059fcb441f4246e31b12c':
    raise RuntimeError('Frozen base checker changed')
spec = spec_from_file_location('Z', SOURCE); Z = module_from_spec(spec); spec.loader.exec_module(Z)


def choose(n, r):
    return 0 if r < 0 else comb(n, r)


def records(k, delta, strong):
    n=2*k; m=n-1 if strong else n+1
    for edges in Z.matchings(tuple(range(n))):
        if strong and not Z.indecomp(edges,n):
            continue
        for ds in Z.compositions(delta,m,k-1):
            yield 'A',edges,tuple(k-d for d in ds),n
        for hole in range(m):
            for ds in Z.compositions(delta-k,m-1,k-1):
                gaps=list(k-d for d in ds);gaps.insert(hole,0)
                yield 'B',edges,tuple(gaps),n
    if delta >= k+1:
        n=2*k-1;m=n-1 if strong else n+1
        for edges in Z.p3_graphs(k):
            if strong and not Z.indecomp(edges,n):
                continue
            for ds in Z.compositions(delta-k-1,m,k-1):
                yield 'C',edges,tuple(k-d for d in ds),n


def main():
    layers=[];relabels=0
    for k,delta in [(3,4),(3,5),(4,6)]:
        fs=1
        for a in range(1,k+1):fs*=2*a-1
        ms=list(Z.matchings(tuple(range(2*k))))
        ik=sum(Z.indecomp(e,2*k) for e in ms)
        jk=sum(sum(b==a+1 for a,b in e) for e in ms if Z.indecomp(e,2*k))
        pk=sum(Z.indecomp(e,2*k-1) for e in Z.p3_graphs(k))
        for strong in (False,True):
            m=2*k-1 if strong else 2*k+1
            E,J,P=(ik,jk,pk) if strong else (fs,fs,(k-1)*fs)
            expected=(E*(choose(m+delta-1,delta)-m*choose(m+delta-k-1,delta-k))
                     +(m*E-J)*choose(m+delta-k-2,delta-k)
                     +P*choose(m+delta-k-3,delta-k-1))
            seen=set();types=Counter();sizes=Counter();rejected=0
            for kind,edges,gaps,n in records(k,delta,strong):
                t,w=Z.model(n,edges,gaps,strong)
                clones=Z.clones(t)
                if kind=='B':
                    h=gaps.index(0)+(1 if strong else 0)
                    bad=0<h<n and (h-1,h) in edges
                    Z.need(bool(clones)==bad,'complete missing-gap clone filter with a unit gap')
                    if bad:rejected+=1;continue
                Z.need(not clones,'retained actual pair-free core')
                if strong:Z.need(Z.strongly_connected(t),'actual strongness')
                key,size=Z.canonical(t,w)
                Z.need(key not in seen,'two-unit-gap record injectivity')
                seen.add(key);types[kind]+=1;sizes[size]+=1
                if size>=4 and relabels<40:
                    order=tuple(reversed(range(len(t))))
                    tt=tuple(sum((t[i]>>j&1)<<p for p,j in enumerate(order)) for i in order)
                    ww=tuple(w[i] for i in order)
                    Z.need(Z.canonical(tt,ww)[0]==key,'exceptional-cell relabel covariance');relabels+=1
            Z.need(len(seen)==expected,'unified formula at enlarged boundary')
            layers.append(dict(k=k,delta=delta,strong=strong,classes=len(seen),types=dict(types),
                               unit_cells=dict(sizes),clone_rejections=rejected))
    collisions=[]
    for k in (2,3,4,5):
        # All three records have the same external marked C3, with an edge
        # from its root to the final endpoint in the backward graph.
        matching=((0,2),(1,2*k-1))+tuple((a,a+1) for a in range(3,2*k-2,2))
        d_gaps=[k]*(2*k+1);d_gaps[1]=d_gaps[2]=0
        p1=((0,1),(0,2*k-2))+tuple((a,a+1) for a in range(2,2*k-3,2))
        p2=((0,1),(1,2*k-2))+tuple((a,a+1) for a in range(2,2*k-3,2))
        c_gaps=[k]*(2*k);c_gaps[1]=1
        for strong in (False,True):
            dg=tuple(d_gaps[1:-1] if strong else d_gaps)
            cg=tuple(c_gaps[1:-1] if strong else c_gaps)
            examples=[Z.model(2*k,matching,dg,strong),Z.model(2*k-1,p1,cg,strong),
                      Z.model(2*k-1,p2,cg,strong)]
            keys=[]
            for t,w in examples:
                Z.need(not Z.clones(t),'collision uses valid cores, not clone merges')
                if strong:Z.need(Z.strongly_connected(t),'strong collision')
                Z.need(sum(w)==(2*k*k+k if strong else 2*k*k+3*k)-2*k,'first collision mass')
                keys.append(Z.canonical(t,w)[0])
            Z.need(keys[0]==keys[1]==keys[2],'literal three-record owner collision')
            collisions.append(dict(k=k,strong=strong,matching=matching,p3_first=p1,p3_second=p2,
                                   matching_gaps=dg,p3_gaps=cg,
                                   actual_tournaments=[{'rows':t,'weights':w} for t,w in examples]))
    print(json.dumps(dict(status='FINITE-EXACT enlarged corridor and first-collision controls',
                         layers=layers,relabel_controls=relabels,three_record_collisions=collisions),
                     indent=2,sort_keys=True))


if __name__=='__main__':main()
