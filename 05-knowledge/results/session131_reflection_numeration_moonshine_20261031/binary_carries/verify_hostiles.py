#!/usr/bin/env python3
"""Independent literal-permutation verification of every hostile witness."""
import itertools as it
import json
from pathlib import Path

checks=0

def require(ok,msg):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(msg)


def recompute(d):
    n=d['n']; r=d['rows']; rho=d['rho']
    edge=lambda u,v: (r[u]>>v)&1
    require(all(edge(u,v)+edge(v,u)==1 for u,v in it.combinations(range(n),2)), 'tournament')
    h=sum(all(edge(u,v) for u,v in zip(p,p[1:])) for p in it.permutations(range(n)))
    require(h==d['H'], 'literal path count')
    cc=[]
    for size in range(3,n+1,2):
        for vertices in it.combinations(range(n),size):
            root=vertices[0]
            for tail in it.permutations(vertices[1:]):
                cyc=(root,)+tail
                if all(edge(cyc[k],cyc[(k+1)%size]) for k in range(size)):
                    cc.append(cyc)
    sup=[frozenset(c) for c in cc]
    m=[]
    for c in cc:
        raw=tuple(rho[u] for u in c[::-1]);k=raw.index(min(raw))
        m.append(cc.index(raw[k:]+raw[:k]))
    q=sum(v==m[v] for v in range(len(cc)))
    i=sum(v<m[v] and bool(sup[v]&sup[m[v]]) for v in range(len(cc)))
    j=sum(v<m[v] and not sup[v]&sup[m[v]] for v in range(len(cc)))
    pair=[s for s in it.combinations(range(len(cc)),2) if not sup[s[0]]&sup[s[1]]]
    nonfixed=[s for s in pair if tuple(sorted(m[v] for v in s))!=s]
    pair_orbits={min(s,tuple(sorted(m[v] for v in s))) for s in nonfixed}
    p=len(pair_orbits)
    triple=[s for s in it.combinations(range(len(cc)),3)
            if sum(len(sup[v]) for v in s)==len(set().union(*(sup[v] for v in s)))]
    t=sum(tuple(sorted(m[v] for v in s))==s for s in triple)
    for name,val in [('q',q),('i',i),('j',j),('p',p),('t',t)]:
        require(val==d[name],'independent '+name)
    alpha=[1,len(cc),len(pair),len(triple)]
    while not alpha[-1]:alpha.pop()
    require(alpha==d['alpha'],'independent OCF profile')
    require(h==sum(x*(1<<k) for k,x in enumerate(alpha)), 'OCF vs words')
    require(h%16==(1+2*q+4*i+8*(j+p+t))%16,'mod16')
    return {'n':n,'rows':r,'H':h,'q':q,'i':i,'j':j,'p':p,'t':t,'alpha':alpha}

x=json.loads((Path(__file__).parent/'normal.json').read_text())
seen={}; out={}
for label, data in x['hostiles'].items():
    batch=data if isinstance(data,list) else [data]
    vals=[]
    for d in batch:
        key=tuple(d['rows'])
        if key not in seen:seen[key]=recompute(d)
        vals.append(seen[key])
    out[label]=vals
    if label=='dropping_binary_carry':
        d=vals[0];require(d['H']%8!=(1+2*(d['q']%2)+4*(d['i']%2))%8,'mutation no carry')
    elif label=='residue7_is_attained':
        require(vals[0]['H']%8==7 and vals[0]['H']!=7,'exact vs residue')
    elif label.startswith('omitting_'):
        d=vals[0];term=label[-1]
        require((1+2*d['q']+4*d['i']+8*(d['j']+d['p']+d['t']-d[term]))%16!=d['H']%16,'mutation '+term)
    elif label=='same_q_i_different_mod16':
        a,b=vals;require((a['n'],a['q'],a['i'])==(b['n'],b['q'],b['i']) and a['H']%16!=b['H']%16,'q i insufficiency')
    elif label=='H_does_not_determine_OCF':
        a,b=vals;require(a['H']==b['H'] and a['alpha']!=b['alpha'],'OCF nonuniqueness')
print(json.dumps({'status':'PASS','checks':checks,'unique_tournaments':len(seen),'hostiles':out},sort_keys=True,indent=2))
