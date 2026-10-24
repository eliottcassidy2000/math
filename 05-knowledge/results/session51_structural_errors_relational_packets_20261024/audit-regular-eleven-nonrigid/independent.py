#!/usr/bin/env python3
"""Independent bounded regular11 nonrigid generation, without a target count.
Matrices by literal row choices; tiny block types by degree constraints;
full isomorphisms by MRV adjacency backtracking, never target canonicalization.
"""
import argparse,collections,itertools,json,pathlib

def need(x,m):
    if not x:raise RuntimeError(m)
def cyclic(n):return tuple(sum(1<<((i+d)%n) for d in range(1,(n+1)//2)) for i in range(n))
def valid(r):
    n=len(r);need(all(x<(1<<n) and not (x>>i)&1 for i,x in enumerate(r)),'valid rows')
    need(all(((r[i]>>j)&1)+((r[j]>>i)&1)==1 for i in range(n) for j in range(i)),'one arc per pair')
def join(blocks):
    out=[];off=0
    for b in blocks:out.extend(x<<off for x in b);off+=len(b)
    return out
def arc(r,i,j):r[i]|=1<<j

def isomorphisms(a,b,first=False):
    n=len(a);da=[x.bit_count() for x in a];db=[x.bit_count() for x in b]
    if sorted(da)!=sorted(db):return []
    ca=[[(x&y).bit_count() for y in a] for x in a];cb=[[(x&y).bit_count() for y in b] for x in b]
    sa=[sorted((da[j],ca[i][j],(a[i]>>j)&1) for j in range(n) if j!=i) for i in range(n)]
    sb=[sorted((db[j],cb[i][j],(b[i]>>j)&1) for j in range(n) if j!=i) for i in range(n)]
    initial=[[j for j in range(n) if da[i]==db[j] and sa[i]==sb[j]] for i in range(n)]
    f=[-1]*n;out=[]
    def go(used):
        if len(used)==n:out.append(tuple(f));return first
        best=None;choices=None
        for u in range(n):
            if f[u]>=0:continue
            options=[v for v in initial[u] if v not in used and all(((a[u]>>x)&1)==((b[v]>>f[x])&1) and ca[u][x]==cb[v][f[x]] for x in range(n) if f[x]>=0)]
            if not options:return False
            if choices is None or len(options)<len(choices):best=u;choices=options
        for v in choices:
            f[best]=v
            if go(used|{v}):return True
            f[best]=-1
        return False
    go(set());return out

def orbits(group):
    remaining=set(range(len(group[0])));out=[]
    while remaining:
        v=min(remaining);o={p[v] for p in group};need(o<=remaining,'disjoint automorphism orbits');remaining-=o;out.append(sorted(o))
    return out

def matrices():
    rows=[sum(1<<j for j in s) for s in itertools.combinations(range(5),3)]
    return {x for x in itertools.product(rows,repeat=5) if all(sum((r>>j)&1 for r in x)==3 for j in range(5))}
def shift(x,u,v):
    out=[0]*5
    for i in range(5):
        for j in range(5):out[(i+u)%5]|=((x[i]>>j)&1)<<((j+v)%5)
    return tuple(out)

def q3():
    allm=matrices();remaining=set(allm);out=[];sector=collections.Counter()
    while remaining:
        m=min(remaining);orb={shift(m,u,v) for u in range(5) for v in range(5)};need(orb<=allm,'matrix shift stays in universe');remaining-=orb;stab=25//len(orb);sector[stab]+=1
        if stab==1:continue
        r=join([cyclic(5),cyclic(5),(0,)])
        for i in range(5):
            arc(r,10,i);arc(r,5+i,10)
            for j in range(5):arc(r,i,5+j) if (m[i]>>j)&1 else arc(r,5+j,i)
        out.append(tuple(r))
    return out,dict(labelled_degree_three_matrices=len(allm),translation_orbit_stabilizer_histogram=sorted(sector.items()))

def q7():
    out=[]
    for d in (0,1):
        remaining=set(itertools.combinations(range(5),4-d));subsets=[]
        while remaining:
            s=min(remaining);remaining-={tuple(sorted((x+z)%5 for x in s)) for z in range(5)};subsets.append(s)
        for s in subsets:
            for sign in ([1] if d==0 else [1,-1]):
                r=join([cyclic(3),cyclic(3),cyclic(5)])
                for i in range(3):
                    for j in range(3):arc(r,i,3+j) if d and j==sign*i%3 else arc(r,3+j,i)
                    for j in range(5):
                        if j in s:arc(r,i,6+j);arc(r,6+j,3+i)
                        else:arc(r,6+j,i);arc(r,3+i,6+j)
                out.append(tuple(r))
    return out

def strong(r):
    for v in range(len(r)):
        seen={v};todo=[v]
        while todo:
            u=todo.pop()
            for w in range(len(r)):
                if (r[u]>>w)&1 and w not in seen:seen.add(w);todo.append(w)
        if len(seen)!=len(r):return False
    return True

def q9():
    pairs=list(itertools.combinations(range(4),2));candidates=[]
    for mask in range(64):
        r=[0]*4
        for i,(u,v) in enumerate(pairs):arc(r,u,v) if (mask>>i)&1 else arc(r,v,u)
        if strong(r):candidates.append(tuple(r))
    small=candidates[0]
    need(all(isomorphisms(small,x,True) for x in candidates),'unique strong-four class by direct isomorphism')
    need(len(isomorphisms(small,small))==1,'strong-four rigid')
    low=[i for i in range(4) if small[i].bit_count()==1];high=[i for i in range(4) if small[i].bit_count()==2];out=[]
    for target in itertools.permutations(high):
        match=dict(zip(low,target));r=join([cyclic(3),small,small])
        for i in range(3):
            for j in range(4):arc(r,3+j,i);arc(r,i,7+j)
        for i in range(4):
            for j in range(4):arc(r,3+i,7+j) if match.get(i)==j else arc(r,7+j,3+i)
        out.append(tuple(r))
    return out,len(candidates)

def q1():
    out=[]
    for mask in range(1<<5):
        steps=[d if mask>>(d-1)&1 else 11-d for d in range(1,6)]
        out.append(tuple(sum(1<<((i+d)%11) for d in steps) for i in range(11)))
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--representatives',required=True);a=p.parse_args()
    c3,matrix=q3();c9,strong4=q9();families={'five-five-one':c3,'two-triples':q7(),'one-triple':c9,'circulants':q1()};reps=[];raw=[]
    for name,values in families.items():
        for r in values:
            valid(r);need(all(x.bit_count()==5 for x in r),'candidate regular eleven');g=isomorphisms(r,r);o=orbits(g);need(len(g)>1,'candidate nonrigid');q=len(o)
            match=next((j for j,item in enumerate(reps) if item['aut_order']==len(g) and item['orbit_sizes']==sorted(map(len,o)) and isomorphisms(r,tuple(item['rows']),True)),None)
            if match is None:
                match=len(reps);reps.append(dict(rows=list(r),aut_order=len(g),orbit_sizes=sorted(map(len,o)),pointed_types=q,error=11*(q-1),deficit=11-q))
            raw.append(dict(family=name,representative=match,aut_order=len(g),orbit_sizes=sorted(map(len,o))))
    error=collections.Counter(r['error'] for r in reps);aut=collections.Counter(r['aut_order'] for r in reps);delta=sum(r['deficit'] for r in reps)
    need(sum(r['pointed_types'] for r in reps)+delta==11*len(reps),'nonrigid pointing deficit')
    out=dict(status='PASS',universe='only nonrigid regular-eleven structural candidates; no full regular atlas',matrix_universe=matrix,labelled_strong_four=strong4,raw_family_counts={name:len(v) for name,v in families.items()},classes=len(reps),error_histogram=sorted(error.items()),aut_histogram=sorted(aut.items()),nonrigid_pointed_types=sum(r['pointed_types'] for r in reps),pointing_deficit=delta,representatives=reps,raw_presentations=raw)
    pathlib.Path(a.output).write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');pathlib.Path(a.representatives).write_text(str(len(reps))+'\n'+'\n'.join(' '.join(map(str,r['rows'])) for r in reps)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('representatives','raw_presentations')},sort_keys=True))
if __name__=='__main__':main()
