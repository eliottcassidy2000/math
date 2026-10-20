"""Construct/replay all-order outer detours on complete small universes.

No assertions: checks remain active under python -O. This program's
finite scope is all fixed-path strong presentations through order seven;
the general quantifier is supplied by outer.md.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import json
from outer_probe import G, fixed_paths, paths_to, strong

HERE=Path(__file__).resolve().parent


def scc(a):
    n=len(a)
    reach=[]
    for v in range(n):
        r=1<<v
        while True:
            t=r
            for u in range(n):
                if r>>u&1:
                    t |= a[u]
            if t==r:
                break
            r=t
        reach.append(r)
    unused=set(range(n));components=[]
    while unused:
        u=min(unused)
        c=tuple(v for v in sorted(unused) if reach[u]>>v&1 and reach[v]>>u&1)
        components.append(c);unused.difference_update(c)
    # All inputs here have the displayed standard HP, hence SCC order is
    # the order of their first path position. Verify it rather than infer it.
    for i,c in enumerate(components):
        for d in components[i+1:]:
            G.need(all(a[u]>>v&1 for u in c for v in d),'SCC order failure')
    return components


def classify(b):
    r=len(b);z=r-1
    different=next((p for p in paths_to(b,z) if G.induced(b,p)!=b),None)
    if different is not None:
        return dict(kind='same_endpoint',path=different)
    G.need(b[z]&1,'endpoint-rigid HP does not close')
    G.need(G.induced(b,tuple(range(1,r))+(0,))!=b,'rotation guard failure')
    for i in range(1,r-2):
        if b[i+2]>>i&1:
            return dict(kind='triangle',triangle=i)
    # No consecutive cyclic triangle remains after the original first.
    start=0;depth=0
    while True:
        tail=tuple(range(start,r));d=G.induced(b,tail)
        G.need(strong(d),'recursive tail not strong')
        G.need(all(G.induced(d,p)==d for p in paths_to(d,len(d)-1)),
               'recursive tail not endpoint-rigid')
        if b[z]>>(start+1)&1:
            G.need(r-start>=4,'a C3 escaped the earlier triangle branch')
            return dict(kind='ear',start=start,depth=depth)
        # Prefix-rotation lemma: this first vertex dominates all but z.
        G.need(all(b[start]>>v&1 for v in range(start+1,z)),
               'proved near-source conclusion failed')
        remainder=G.induced(b,tuple(range(start+1,r)))
        parts=scc(remainder)
        G.need(len(parts)>=2,'deleting near-source left a strong remainder')
        # In this triangle-free branch all preceding rigid SCCs are singles.
        G.need(all(len(c)==1 for c in parts[:-1]),'non-singleton rigid prefix')
        last=parts[-1]
        if len(last)==1:
            k=r-start-3
            G.need(k>=1,'terminal triangle escaped initial branch')
            return dict(kind='five',start=start,depth=depth,k=k)
        start=start+1+last[0]
        depth+=1


def words(b,zmiddle,certificate):
    r=len(b);offset=r+len(zmiddle);n=offset+r
    A=tuple(range(offset));right=tuple(range(offset,n))
    typ=certificate['kind']
    if typ=='same_endpoint':
        alternate=A+tuple(offset+v for v in certificate['path'])
        original=tuple(range(n))
        return [alternate,alternate[1:]+alternate[:1],original[1:]+original[:1]]
    if typ=='triangle':
        i=certificate['triangle'];a=offset
        L=tuple(range(offset+1,offset+i))
        b0,c,d=range(offset+i,offset+i+3)
        R=tuple(range(offset+i+3,n))
        return [A+L+(b0,c,d)+R+(a,), (c,)+A+L+(b0,d)+R+(a,),
                (c,)+A+(a,)+L+(b0,d)+R]
    s=offset+certificate['start'];F=tuple(range(s));a=s;z=n-1
    if typ=='ear':
        b0=s+1;C=tuple(range(s+2,n-1))
        return [F+(b0,)+C+(z,a), C+F+(b0,z,a), C+F+(a,b0,z)]
    G.need(typ=='five','unknown certificate')
    V=tuple(range(s+1,n-2));c=n-2
    return [F+V+(c,z,a),(z,)+F+V+(c,a),(a,z)+F+V+(c,),
            (a,c,z)+F+V,(c,z)+F+(a,)+V]


def replay(b,zmiddle,certificate):
    source=G.join([G.transpose(b),zmiddle,b]);n=len(source)
    original=tuple(range(n))
    G.need(G.blue(source,original),'source not blue')
    state=source
    paths=words(b,zmiddle,certificate)
    for path in paths:
        G.need(G.is_path(state,path),'undirected path '+str((b,certificate,path)))
        G.need(not G.blue(state,path),'blue move '+str((b,certificate,path)))
        state=G.complement(state,path)
    G.need(state==G.complement(source,original),'literal target mismatch')
    return paths


def wrap_nearsource(b):
    # a -> singleton -> b, except b's terminal vertex beats a.
    out=list(G.join([(0,),(0,),b]));z=len(out)-1
    out[0]^=1<<z;out[z]^=1
    return tuple(out)


def family_controls():
    total=0;deep=[]
    # Arbitrary V is retained, not assumed rigid or transitive.
    for k in range(1,7):
        for V in fixed_paths(k):
            b=list(G.join([(0,),V,(0,),(0,)]));z=len(b)-1
            b[0]^=1<<z;b[z]^=1;b=tuple(b)
            for depth in (0,1,3):
                expanded=b
                for _ in range(depth):expanded=wrap_nearsource(expanded)
                cert=dict(kind='five',start=2*depth,depth=depth,k=k)
                for middle in ((),(0,)):
                    replay(expanded,middle,cert);total+=1
    for depth in range(5):
        b=(6,12,8,1)
        for _ in range(depth):b=wrap_nearsource(b)
        cert=classify(b)
        G.need(cert['kind']=='five' and cert['depth']==depth,'deep recursion control')
        replay(b,(),cert)
        deep.append(dict(order=len(b),depth=depth))
    return dict(arbitrary_V_controls=total,deep_recursion=deep)


def main():
    records=[];examples={};total=0
    middles=[(),(0,),(2,4,1)]
    for n in range(4,8):
        stats=Counter();depth=0
        for b in fixed_paths(n):
            if not strong(b):continue
            cert=classify(b);stats[cert['kind']]+=1
            depth=max(depth,cert.get('depth',0))
            for middle in middles:
                p=replay(b,middle,cert);total+=1
                if cert['kind'] not in examples:
                    examples[cert['kind']]=dict(B=b,middle=middle,certificate=cert,paths=p)
        records.append(dict(order=n,classes=dict(stats),maximum_recursion_depth=depth))
        print(records[-1],flush=True)
    family=family_controls()
    result=dict(universe='all fixed-path strong presentations of orders 4..7',
                middles=middles,total_literal_controls=total,records=records,examples=examples)
    result['family_controls']=family
    (HERE/'outer_verify.json').write_text(json.dumps(result,indent=2)+'\n')
    print('literal black detours passed',total,flush=True)


if __name__=='__main__':main()
