#!/usr/bin/env python3
"""Actual-source, named-port continuation certificates; no relabelling move.

Three temporary vertices implement a clone-compatible reversal of an exterior
arc. All path guards and every old named arc are checked during replay.
"""
from itertools import permutations,combinations
from pathlib import Path
import json
from scout import native,birth,clone,paths,hp,reverse,require,strong

HERE=Path(__file__).resolve().parent


def restrict(t,keep):
    return tuple(sum(1<<j for j,b in enumerate(keep) if t[a]>>b&1)for a in keep)


def flip(t,a,b):
    r=list(t);r[a]^=1<<b;r[b]^=1<<a;return tuple(r)


def exterior_recipe(t,v,a,b):
    """Construct base paths; expanding v in both words preserves its frame."""
    n=len(t)
    if not(t[a]>>b&1):a,b=b,a
    keep=[i for i in range(n)if i not in(a,b)]
    p0=next(paths(restrict(t,keep)))
    p0=tuple(keep[i]for i in p0)
    # z loses only to b; s is a sink; w beats only a.
    z,s,w=n,n+1,n+2
    masks=[1<<b,(1<<(n+1))-1,((1<<(n+2))-1)^(1<<a)]
    ext=t
    for mask in masks:ext=birth(ext,mask)
    p=(a,b,z)+p0+(s,w)
    q=(a,w,b,z)+p0+(s,)
    result=native(native(ext,p),q)
    require(restrict(result,range(n))==flip(t,a,b),'exterior flip endpoint')
    return {'a':a,'b':b,'masks':masks,'first_path':p,'second_path':q}


def expand_path(p,n,v,m):
    addresses=[(i,j)for i in range(n)for j in range(m if i==v else 1)]
    lookup={x:i for i,x in enumerate(addresses)}
    return tuple(lookup[(i,j)]for i in p for j in range(m if i==v else 1))


def replay(t,v,recipe,m):
    n=len(t);ext=t
    for mask in recipe['masks']:ext=birth(ext,mask)
    large=clone(ext,v,m)
    p=expand_path(recipe['first_path'],n+3,v,m)
    q=expand_path(recipe['second_path'],n+3,v,m)
    got=native(native(large,p),q)
    keep=range(n+m-1)
    require(restrict(got,keep)==clone(flip(t,recipe['a'],recipe['b']),v,m),'expanded native endpoint')


def rename(t,mapping):
    out=[0]*len(t)
    for i,row in enumerate(t):
        out[mapping[i]]=sum(1<<mapping[j]for j in range(len(t))if row>>j&1)
    return tuple(out)


def fit_target(source,v,target):
    """Choose a target's physical names at design time, fixing the old root star."""
    n,N=len(source),len(target);outold=[i for i in range(n)if source[v]>>i&1]
    inold=[i for i in range(n)if i!=v and not(source[v]>>i&1)]
    best=None
    for converse in (False,True):
        cand=reverse(target)if converse else target
        oo=[i for i in range(1,N)if cand[0]>>i&1];ii=[i for i in range(1,N)if not(cand[0]>>i&1)]
        if len(oo)<len(outold)or len(ii)<len(inold):continue
        newborn=list(range(n,N));k=len(oo)-len(outold)
        # New names are unconstrained; this fixed partition loses no old-arc minimum.
        outnames=outold+newborn[:k];innames=inold+newborn[k:]
        for po in permutations(outnames):
            for pi in permutations(innames):
                mapping=[0]*N;mapping[0]=v
                for i,j in zip(oo,po):mapping[i]=j
                for i,j in zip(ii,pi):mapping[i]=j
                named=rename(cand,mapping)
                diff=[(i,j)for i,j in combinations(range(n),2)if ((source[i]>>j)^(named[i]>>j))&1]
                require(all(v not in e for e in diff),'root star retained')
                key=(len(diff),named)
                if best is None or key<best[0]:best=(key,named,mapping,converse,diff)
    require(best is not None,'target star embeds source')
    _,named,mapping,converse,diff=best
    return named,mapping,converse,diff


def compile_return(source,v,target,offset):
    named,mapping,converse,diff=fit_target(source,v,target)
    t=source;births=[]
    for j in range(len(source),len(named)):
        mask=sum(1<<i for i in range(j)if named[i]>>j&1)
        births.append(mask);t=birth(t,mask)
    flips=[]
    for a,b in diff:
        recipe=exterior_recipe(t,v,a,b)
        for m in (1,2,3,5):replay(t,v,recipe,m)
        flips.append(recipe);t=flip(t,a,b)
    require(t==named,'actual target reached without relabelling or deleting old names')
    values=[]
    for m in range(1,6):
        h,k=hp(clone(source,v,m)),hp(clone(named,v,m))
        require(k==3*h+offset,'ternary scalar endpoint')
        values.append([m,h,k])
    d=len(flips);b=len(births)
    return {'source':source,'clone':v,'retained_ports':[5,1],'offset':offset,'target':named,
            'atlas_target':target,'design_mapping':mapping,'design_converse':converse,
            'permanent_birth_masks':births,'exterior_flips':flips,'values':values,
            'clock':{'permanent_births':b,'temporary_births':3*d,'native_complements':2*d,
                     'temporary_deletions':3*d,'peak_order_at_m1':len(named)+(3 if d else 0)}}


def peel(t):
    names=list(range(len(t)));history=[]
    while len(t)>3:
        n=len(t);choices=[]
        for v in range(n):
            incoming=[i for i in range(n)if t[i]>>v&1]
            if len(incoming)==1:
                keep=[i for i in range(n)if i!=v];sub=restrict(t,keep)
                if strong(sub):choices.append((v,incoming[0],keep,sub))
        if not choices:break
        require(len(choices)==1,'unique adaptive peel')
        v,r,keep,t=choices[0];history.append([names[v],names[r]]);names=[names[i]for i in keep]
    return {'core_names':names,'core':t,'core_H':hp(t),'reverse_history':history}


if __name__=='__main__':
    scout=json.loads((HERE/'scout31.json').read_text())['programmes']
    atlas=[json.loads(s)for f in ('catalog9_3.jsonl','catalog10_3.jsonl')for s in (HERE/f).read_text().splitlines()]
    targets={h:next(tuple(r['rows'])for r in atlas if r['H']==h)for h in (287,289,295)}
    results=[]
    for index,c,targetH in ((1,2,287),(1,4,289),(0,4,295)):
        r=compile_return(tuple(scout[index]['target']),3,targets[targetH],c)
        r['first_programme_index']=index;r['ancestry']=peel(tuple(r['target']))
        results.append(r)
    (HERE/'returns.json').write_text(json.dumps(results,indent=2)+'\n')
    for r in results:print(r['values'][0],r['clock'],r['ancestry'])
